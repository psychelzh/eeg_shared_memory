load(fullfile("data", "cca_result_model.mat"), "Y")
num_times = 257;
num_trials = 150;
sz = size(Y);
Y = reshape(Y, num_times, num_trials, sz(2), sz(3));
num_comp = 3;
isc_comp = nan(sz(3), sz(3), num_comp, num_trials);
for i_comp = 1:num_comp % keep first 3 components only
    for i_trial = 1:num_trials
        isc_comp(:, :, i_comp, i_trial) = corr(squeeze(Y(:, i_trial, i_comp, :)), ...
            "Rows", "pairwise");
    end
end
save(fullfile("data", "cca_result_isc_trials"), "isc_comp")

[comp, trial, subj] = ndgrid(1:num_comp, 1:num_trials, 1:size(Y, 4));
Y_tbl = table(comp(:), trial(:), subj(:), VariableNames=["comp", "trial", "subj"]);
Y_tbl.ts = rowfun( ...
    @(comp, trial, subj) squeeze(Y(:, trial, comp, subj)), ...
    Y_tbl, "InputVariables", ["comp", "trial", "subj"], ...
    "OutputFormat", "cell");
parquetwrite(fullfile("data", "cca_y_first3comp.parquet"), Y_tbl)