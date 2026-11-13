"""Audit a completed native targets store; write a plan without deleting files.

Run from the project root after exporting tar_manifest() to a CSV. Metadata may
contain superseded records and vector-slice child names without object files.
The last record for a name wins; only actual branch records need stored objects.
"""

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


def read_csv(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--source-root", help="Original project root for relocated absolute file paths")
    args = parser.parse_args()
    root = Path.cwd().resolve()
    store = root / "_targets"
    metadata = store / "meta/meta"
    audit = args.audit.resolve()
    audit.mkdir(parents=True, exist_ok=True)
    manifest = read_csv(args.manifest)
    records = read_csv(metadata, delimiter="|")
    by_name = {row["name"]: row for row in records}
    roots = {row["name"] for row in manifest}
    assert roots and roots.issubset(by_name), "Current target metadata is missing"
    keep = set(roots)
    pending = list(roots)
    while pending:
        record = by_name.get(pending.pop())
        if record is None:
            continue  # Native vector slices are indices, not persisted targets.
        for child in filter(None, record["children"].split("*")):
            if child not in keep:
                keep.add(child)
                pending.append(child)

    references = set()
    for name in keep:
        record = by_name.get(name)
        if record is None:
            continue
        if record["format"] == "file":
            for value in filter(None, record["path"].split("*")):
                value = value.replace("\\", "/")
                if args.source_root:
                    source_prefix = args.source_root.replace("\\", "/").rstrip("/") + "/"
                    if value.startswith(source_prefix):
                        value = value[len(source_prefix):]
                path = Path(value)
                if not path.is_absolute():
                    path = root / path
                path = path.resolve()
                assert path.exists(), f"Referenced file/directory missing: {path}"
                references.add(path)
        elif record["type"] in ("stem", "branch"):
            assert (store / "objects" / name).is_file(), f"Object missing: {name}"

    partition_roots = [
        by_name[name]
        for name in roots
        if name.startswith("partitions_") and not name.endswith("_keys")
    ]
    assert len(partition_roots) == 18
    assert all(row["type"] == "pattern" and row["format"] != "file"
               for row in partition_roots)
    pieces = [child for row in partition_roots
              for child in row["children"].split("*") if child]
    assert all(by_name[name]["type"] == "branch" for name in pieces)

    candidates = []
    inventory = []
    for category, directory in [
        ("object", store / "objects"),
        ("generated", store / "user/partitions"),
    ]:
        entries = directory.iterdir() if category == "object" else directory.rglob("*")
        for path in entries:
            if not path.is_file():
                continue
            path = path.resolve()
            size = path.stat().st_size
            inventory.append({"Kind": category, "Path": str(path), "Bytes": size})
            needed = category == "object" and path.name in keep
            referenced = any(path == ref or ref in path.parents for ref in references)
            if not needed and not referenced:
                candidates.append({"Kind": category, "Path": str(path), "Bytes": size})

    finals = sorted(name for name in roots
                    if name.startswith(("stats_", "clusters_stats_"))
                    and not name.endswith("_permuted")
                    and by_name[name]["type"] == "stem")
    hashes = [{"name": name, "SHA256": sha256(store / "objects" / name)}
              for name in finals]
    write_csv(audit / "deletion-plan.csv", candidates, ["Kind", "Path", "Bytes"])
    write_csv(audit / "inventory.csv", inventory, ["Kind", "Path", "Bytes"])
    write_csv(audit / "keep_names.csv", [{"name": name} for name in sorted(keep)], ["name"])
    write_csv(audit / "referenced_files.csv", [{"local": str(p)} for p in sorted(references)], ["local"])
    write_csv(audit / "results-hashes-before.csv", hashes, ["name", "SHA256"])
    shutil.copy2(metadata, audit / "meta.before-prune")
    tracked = [metadata, args.manifest.resolve(), Path(__file__).resolve(),
               root / "_targets.R", root / "R/partitioning.R",
               root / "R/tar_partitions.R", root / "R/partition_sources.R",
               audit / "deletion-plan.csv", audit / "keep_names.csv",
               audit / "referenced_files.csv"]
    snapshot = {
        "Created": datetime.now().isoformat(),
        "Files": [{"Path": str(p), "Hash": sha256(p)} for p in tracked],
        "Definitions": len(roots),
        "NativePartitionRoots": len(partition_roots),
        "NativePieces": len(pieces),
        "ResultHashes": len(hashes),
        "CandidateFiles": len(candidates),
        "CandidateGiB": sum(row["Bytes"] for row in candidates) / 2**30,
    }
    (audit / "audit-fingerprints.json").write_text(json.dumps(snapshot, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in snapshot.items() if k != "Files"}, indent=2))


if __name__ == "__main__":
    main()
