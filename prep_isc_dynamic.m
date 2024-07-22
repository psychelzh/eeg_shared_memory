% load(fullfile("data", "cca_result_model.mat"), "Y")
% num_times = 257; % means 1000 ms
% num_trials = 150;
% sz = size(Y);
% % time x trial x component x subject
% Y = reshape(Y, num_times, num_trials, sz(2), sz(3));
isc_time_win = calc_isc_time_win(Y);
save data\CorCAExtra\isc_time_win_subjs207 isc_time_win -v7.3

% save isc_time_win to parquet
sz = size(isc_time_win);
subj_pair_id = 1:sz(1);
comp_id = 1:sz(2);
trial_id = 1:sz(3);
window_id = 1:sz(4);
isc_tbl = combinations(window_id, trial_id, comp_id, subj_pair_id);
isc_tbl.isc = isc_time_win(:);
parquetwrite("data\CorCAExtra\cca_isc_time_win.parquet", isc_tbl)

function isc = calc_isc_time_win(Y, opts)
arguments
    Y
    opts.NComp = 3 % used components
    opts.SizeTimeWin = 51 % number of time points (default: 200 ms)
    opts.SizeTimeStep = 5 % number of time points (default: 20 ms)
end
[num_times, num_trials, num_comp, num_subj] = size(Y);
if opts.NComp > num_comp
    opts.NComp = num_comp;
end
num_comp = opts.NComp;
time_wins = arrayfun(@(start) start + (1:opts.SizeTimeWin), ...
    1:opts.SizeTimeStep:num_times, 'UniformOutput', false);
isc_cell = cell(size(time_wins));

pb1 = ProgressBar(length(time_wins), ...
    'UpdateRate', Inf, ...
    'Title', 'Time Windows');
pb1.setup([], [], []);

for i_time_win = 1:length(time_wins)
    time_win = time_wins{i_time_win};
    if any(time_win > num_times)
        continue
    end
    y_time_win = Y(time_win, :, :, :);
    isc_comp = nan(num_subj * (num_subj - 1) / 2, num_comp, num_trials);

    pb2 = ProgressBar(num_trials, ...
        'IsParallel', true, ...
        'Title', 'Trials');
    pb2.setup([], [], []);

    parfor i_trial = 1:num_trials
        for i_comp = 1:num_comp
            r = corr(squeeze(y_time_win(:, i_trial, i_comp, :)), "Rows", "pairwise");
            mask_lower = tril(true(size(r)), -1);
            isc_comp(:, i_comp, i_trial) = r(mask_lower);
        end
        updateParallel();
    end
    isc_cell{i_time_win} = isc_comp;

    pb2.release()
    pb1.step(1, [], []);
end

pb1.release()

isc = cat(4, isc_cell{:});
end
