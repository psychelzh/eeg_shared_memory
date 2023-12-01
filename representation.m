channel = readtable(fullfile("config", "eeg_channel_labels_64.csv"), TextType="string");
subjs_id = unique(readtable(fullfile("data", "group_task-wordencoding_events.csv")).subj_id, "stable");
load(fullfile("data", "grp_subjs206_nodemean_1000ms.mat"), "grp_data");
[~, len_time_point, len_trial, len_subj] = size(grp_data);
regions_id = 1:6;

% debug
% regions_id = 1:2;
% len_trial = 2;

% used to construct table pre-allocation
[regions_id_idx, trial_id_idx] = meshgrid(regions_id, 1:len_trial);
numel_region_trial = numel(regions_id_idx);

%% intersubject similarity

numel_cor_mat = len_subj * len_subj;
[subj_row, subj_col] = ind2sub([len_subj, len_subj], 1:(numel_cor_mat));

% acquire: trial-level
simi_inter_by_trial = array2table(...
    [repelem([regions_id_idx(:), trial_id_idx(:)], numel_cor_mat, 1), ...
    repmat([subjs_id(subj_row), subjs_id(subj_col), nan(numel_cor_mat, 1)], numel_region_trial, 1)], ...
    "VariableNames", ["region_id", "trial_id", "subj_id_row", "subj_id_col", "fisher_z"]);
% easier handle data store
store_start = 1;
store_end = numel_cor_mat;
fprintf("Processing trial level intersubject similarity...\n")
for i_region = regions_id
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    prog_bar = ProgressBar(len_trial, Title=['Region ', num2str(i_region)]);
    for i_trial = 1:len_trial
        % collapse channel and time (thus spatiotemporal pattern)
        cur_dat = reshape(grp_data(chan_in_reg, :, i_trial, :), ...
            [length(chan_in_reg) * len_time_point, len_subj]);
        cur_cor_mat = atanh(corr(cur_dat));
        simi_inter_by_trial.fisher_z(store_start:store_end) = cur_cor_mat(:);
        store_start = store_start + numel_cor_mat;
        store_end = store_end + numel_cor_mat;
        prog_bar([], [], [])
    end
    prog_bar.release()
end
parquetwrite(fullfile("data", "type-inter_acq-trial_rs.parquet"), ...
    simi_inter_by_trial(simi_inter_by_trial.subj_id_row < simi_inter_by_trial.subj_id_col, :))

% acquire: whole-time-series
simi_inter_by_whole = array2table(...
    [repelem(regions_id(:), numel_cor_mat, 1), ...
    repmat([subjs_id(subj_row), subjs_id(subj_col), nan(numel_cor_mat, 1)], length(regions_id), 1)], ...
    "VariableNames", ["region_id", "subj_id_row", "subj_id_col", "fisher_z"]);
store_start = 1;
store_end = numel_cor_mat;
fprintf("Processing whole time series intersubject similarity...\n")
for i_region = progress(regions_id)
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    cur_dat = reshape(grp_data(chan_in_reg, :, :, :), ...
        [length(chan_in_reg) * len_time_point * len_trial, len_subj]);
    cur_cor_mat = atanh(corr(cur_dat, rows="pairwise"));
    simi_inter_by_whole.fisher_z(store_start:store_end) = cur_cor_mat(:);
    store_start = store_start + numel_cor_mat;
    store_end = store_end + numel_cor_mat;
end
parquetwrite(fullfile("data", "type-inter_acq-whole_rs.parquet"), ...
    simi_inter_by_whole(simi_inter_by_whole.subj_id_row < simi_inter_by_whole.subj_id_col, :))

%% individual to group similarity

% acquire: trial-level
simi_grp_by_trial = array2table(...
    [repelem([regions_id_idx(:), trial_id_idx(:)], len_subj, 1), ...
    repmat([subjs_id, nan(len_subj, 1)], numel_region_trial, 1)], ...
    "VariableNames", ["region_id", "trial_id", "subj_id", "fisher_z"]);
% easier handle data store
store_start = 1;
store_end = len_subj;
fprintf("Processing trial level individual to group similarity...\n")
for i_region = regions_id
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    prog_bar = ProgressBar(len_trial, Title=['Region ', num2str(i_region)]);
    for i_trial = 1:len_trial
        % collapse channel and time (thus spatiotemporal pattern)
        cur_dat = reshape(grp_data(chan_in_reg, :, i_trial, :), ...
            [length(chan_in_reg) * len_time_point, len_subj]);
        simi_grp_by_trial.fisher_z(store_start:store_end) = utils.calc_simi_ind_to_grp(cur_dat);
        store_start = store_start + len_subj;
        store_end = store_end + len_subj;
        prog_bar([], [], [])
    end
    prog_bar.release()
end
parquetwrite(fullfile("data", "type-group_acq-trial_rs.parquet"), simi_grp_by_trial)

% acquire: whole-time-series
simi_grp_by_whole = array2table(...
    [repelem(regions_id(:), len_subj, 1), ...
    repmat([subjs_id, nan(len_subj, 1)], length(regions_id), 1)], ...
    "VariableNames", ["region_id", "subj_id", "fisher_z"]);
store_start = 1;
store_end = len_subj;
fprintf("Processing whole time series individual to group similarity...\n")
for i_region = progress(regions_id)
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    cur_dat = reshape(grp_data(chan_in_reg, :, :, :), ...
        [length(chan_in_reg) * len_time_point * len_trial, len_subj]);
    simi_grp_by_whole.fisher_z(store_start:store_end) = utils.calc_simi_ind_to_grp(cur_dat);
    store_start = store_start + len_subj;
    store_end = store_end + len_subj;
end
parquetwrite(fullfile("data", "type-group_acq-whole_rs.parquet"), simi_grp_by_whole)
