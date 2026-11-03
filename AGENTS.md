# Project conventions

- Use `air format` to format new or modified R code before finishing changes. Use `air format --check` to verify formatting.

<!-- CODEGRAPH_START -->
## CodeGraph

If `.codegraph/` exists at the repository root, use CodeGraph before searching for or reading code to understand or locate it. Prefer the `codegraph_explore` MCP tool when available; otherwise use `codegraph explore`. If neither is available, fall back to ordinary code search. Do not create an index automatically.
<!-- CODEGRAPH_END -->
