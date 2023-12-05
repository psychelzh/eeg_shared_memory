function cmbn = preallocate(val, opts)
%PREALLOCATE Preallocate results data
arguments (Repeating)
    val (:, 1)
end
arguments
    opts.VariableNames (1, :) {mustBeText}
    opts.AddFisherZ (1, 1) {mustBeNumericOrLogical} = true
    opts.ScalarFisherZ (1, 1) {mustBeNumericOrLogical} = true
end
cmbn = combinations(val{:});
if isfield(opts, "VariableNames")
    if width(val) ~= width(opts.VariableNames)
        error("Number of variable names does not equal to number of columns.")
    end
    cmbn.Properties.VariableNames = opts.VariableNames;
end
if opts.AddFisherZ
    if opts.ScalarFisherZ
        cmbn.fisher_z = nan(height(cmbn), 1);
    else
        cmbn.fisher_z = cell(height(cmbn), 1);
    end
end
end
