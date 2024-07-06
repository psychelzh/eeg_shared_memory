data_source = load(fullfile("data", "grp_subjs206_nodemean_1000ms.mat"));
% raw data is channel x time x trial x subjects
sz = size(data_source.grp_data);
% we will collapse time and trial as time
% and transform as time x channel x subjects
X = permute(reshape(data_source.grp_data, sz(1), [], sz(end)), [2, 1, 3]);
clearvars data_source % raw data is very memory consuming

% shuffle to get p values for each component
fprintf("Begin corrca on shuffle data...\n")
num_surrogates = 100;
ISC_null = nan(1, num_surrogates);
% https://github.com/elgar328/matlab-code-examples/issues/2
PB = ProgressBar(num_surrogates, taskname = 'Shuffle Stats', ui = "cli");
NUM_WORKERS = 4;
parfor (i = 1:num_surrogates, NUM_WORKERS)
    [~, ISC] = cca.corrca(generate_surrogate(X));
    ISC_null(i) = ISC(1);
    count(PB)
end

fprintf("Begin corrca on real data...\n")
[W, ISC, Y, A] = cca.corrca(X);

fprintf("All done! Now saving results...\n")
save(fullfile("data", "cca_result_model"), "X", "W", "Y", "A", "-v7.3")
save(fullfile("data", "cca_result_stats"), "ISC", "ISC_null", "-v7.3")

function surrogate = generate_surrogate(X, opts)
arguments
    X {mustBeNumeric}
    opts.DimShift (1, 1) = 1;
    opts.DimLoop (1, 1) = 3;
end
surrogate = X;
dim_shift = opts.DimShift;
dim_loop  = opts.DimLoop;
shifts = randsample(size(X, dim_shift), size(X, dim_loop), true);
for i = 1:size(surrogate, dim_loop)
    surrogate(:, :, i) = circshift(surrogate(:, :, i), shifts(i), dim_shift);
end
end
