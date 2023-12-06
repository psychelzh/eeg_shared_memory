channel = readtable(fullfile("config", "eeg_channel_labels_64.csv"), TextType="string");
subjs_id = unique(readtable(fullfile("data", "group_task-wordencoding_events.csv")).subj_id, "stable");
load(fullfile("data", "grp_subjs206_nodemean_1000ms.mat"), "grp_data");
[~, len_time_point, len_trial, len_subj] = size(grp_data);
regions_id = 1:6;
% we keep the lower triangular correlation (no diagonal)
idx_keep_cors = tril(true(len_subj, len_subj), -1);

% debug
% regions_id = 1:2;
% len_trial = 2;

%% intersubject similarity

% acquire: trial-level
simi_inter_by_trial = utils.filter_triangular( ...
    utils.preallocate(regions_id, 1:len_trial, subjs_id, subjs_id, ...
    VariableNames=["region_id", "trial_id", "subj_id_col", "subj_id_row"]));
fprintf("Processing trial level intersubject similarity...\n")
for i_region = regions_id
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    fprintf("Region " + string(i_region) + "\n")
    for i_trial = progress(1:len_trial)
        % collapse channel and time (thus spatiotemporal pattern)
        cur_cor_mat = corr(reshape(grp_data(chan_in_reg, :, i_trial, :), ...
            length(chan_in_reg) * len_time_point, []));
        simi_inter_by_trial.fisher_z( ...
            simi_inter_by_trial.region_id == i_region & ...
            simi_inter_by_trial.trial_id == i_trial) = ...
            atanh(cur_cor_mat(idx_keep_cors));
    end
end
parquetwrite( ...
    fullfile("data", "type-inter_acq-trial_rs.parquet"), ...
    simi_inter_by_trial)
clearvars simi_inter_by_trial

% acquire: whole-time-series
simi_inter_by_whole = utils.filter_triangular( ...
    utils.preallocate(regions_id, subjs_id, subjs_id, ...
    VariableNames=["region_id", "subj_id_col", "subj_id_row"]));
fprintf("Processing whole time series intersubject similarity...\n")
for i_region = progress(regions_id)
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    cur_cor_mat = corr(reshape(grp_data(chan_in_reg, :, :, :), ...
        length(chan_in_reg) * len_time_point * len_trial, []), ...
        rows="pairwise");
    simi_inter_by_whole.fisher_z( ...
        simi_inter_by_whole.region_id == i_region) = ...
        atanh(cur_cor_mat(idx_keep_cors));
end
parquetwrite( ...
    fullfile("data", "type-inter_acq-whole_rs.parquet"), ...
    simi_inter_by_whole)
clearvars simi_inter_by_whole

%% individual to group similarity

% acquire: trial-level
simi_grp_by_trial = ...
    utils.preallocate(regions_id, 1:len_trial, subjs_id, ...
    VariableNames=["region_id", "trial_id", "subj_id"]);
fprintf("Processing trial level individual to group similarity...\n")
for i_region = regions_id
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    fprintf("Region " + string(i_region) + "\n")
    for i_trial = progress(1:len_trial)
        % collapse channel and time (thus spatiotemporal pattern)
        cur_dat = reshape(grp_data(chan_in_reg, :, i_trial, :), ...
            length(chan_in_reg) * len_time_point, []);
        simi_grp_by_trial.fisher_z(...
            simi_grp_by_trial.region_id == i_region & ...
            simi_grp_by_trial.trial_id == i_trial) = ...
            utils.calc_simi_ind_to_grp(cur_dat, FisherZ=true);
    end
end
parquetwrite( ...
    fullfile("data", "type-group_acq-trial_rs.parquet"), ...
    simi_grp_by_trial)
clearvars simi_grp_by_trial

% acquire: whole-time-series
simi_grp_by_whole = utils.preallocate(regions_id, subjs_id, ...
    VariableNames=["region_id", "subj_id"]);
fprintf("Processing whole time series individual to group similarity...\n")
for i_region = progress(regions_id)
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    cur_dat = reshape(grp_data(chan_in_reg, :, :, :), ...
        length(chan_in_reg) * len_time_point * len_trial, []);
    simi_grp_by_whole.fisher_z(...
        simi_grp_by_whole.region_id == i_region) = ...
        utils.calc_simi_ind_to_grp(cur_dat, FisherZ=true);
end
parquetwrite( ...
    fullfile("data", "type-group_acq-whole_rs.parquet"), ...
    simi_grp_by_whole)
clearvars simi_grp_by_whole

%% windowed results (separate regions)

% setup for windowed calculations
size_window = 26;
step = 5;
[window_start, window_end] = utils.setup_window(len_time_point, size_window, step);

% type: intersubject similarity
fprintf("Processing stepped window intersubject similarity...\n")
for i_region = regions_id
    cur_simi_inter_by_window = utils.filter_triangular( ...
        utils.preallocate(1:len_trial, 1:length(window_start), subjs_id, subjs_id, ...
        VariableNames=["trial_id", "window_id", "subj_id_col", "subj_id_row"]));
    cur_reg = "region" + string(i_region);
    chan_in_reg = channel.code(channel.(cur_reg) ~= 0);
    fprintf(cur_reg + "\n")
    for i_trial = progress(1:len_trial)
        cur_trial_dat = grp_data(chan_in_reg, :, i_trial, :);
        fisher_z = cell(length(window_start), 1);
        parfor i_win = 1:length(window_start)
            % collapse channel and time (thus spatiotemporal pattern)
            cur_cor_mat = corr(reshape( ...
                cur_trial_dat(:, window_start(i_win):window_end(i_win), :, :), ...
                length(chan_in_reg) * size_window, [])); %#ok<*PFBNS>
            fisher_z{i_win} = atanh(cur_cor_mat(idx_keep_cors));
        end
        cur_simi_inter_by_window.fisher_z( ...
            cur_simi_inter_by_window.trial_id == i_trial) = ...
            vertcat(fisher_z{:});
    end
    parquetwrite( ...
        fullfile("data", "type-inter_acq-window_region-" + cur_reg + "_rs.parquet"), ...
        cur_simi_inter_by_window)
end
clearvars cur_simi_inter_by_window

% type: individual to group similarity
simi_grp_by_window = nan(length(regions_id), len_trial, length(window_start), len_subj);
fprintf("Processing stepped window intersubject similarity...\n")
for i_region = regions_id
    cur_simi_grp_by_window = ...
        utils.preallocate(1:len_trial, 1:length(window_start), subjs_id, ...
        VariableNames=["trial_id", "window_id", "subj_id"]);
    cur_reg = "region" + string(i_region);
    chan_in_reg = channel.code(channel.(cur_reg) ~= 0);
    fprintf(cur_reg + "\n")
    for i_trial = progress(1:len_trial)
        cur_trial_dat = grp_data(chan_in_reg, :, i_trial, :);
        fisher_z = cell(length(window_start), 1);
        parfor i_win = 1:length(window_start)
            % collapse channel and time (thus spatiotemporal pattern)
            cur_dat = reshape( ...
                cur_trial_dat(:, window_start(i_win):window_end(i_win), :, :), ...
                length(chan_in_reg) * size_window, []);
            fisher_z{i_win} = utils.calc_simi_ind_to_grp(cur_dat, FisherZ=true);
        end
        cur_simi_grp_by_window.fisher_z(...
            cur_simi_grp_by_window.trial_id == i_trial) = ...
            vertcat(fisher_z{:});
    end
    parquetwrite( ...
        fullfile("data", "type-group_acq-window_region-" + cur_reg + "_rs.parquet"), ...
        cur_simi_grp_by_window)
end
clearvars cur_simi_grp_by_window

%% windowed results (separate trials)

% setup for windowed calculations
size_window = 26;
step = 5;
[window_start, window_end] = utils.setup_window(len_time_point, size_window, step);

% type: intersubject similarity
fprintf("Processing stepped window intersubject similarity...\n")
path = fullfile("data", "type-inter_acq-window");
if (~exist(path, "dir")), mkdir(path), end
for i_region = regions_id
    cur_reg = "region" + string(i_region);
    chan_in_reg = channel.code(channel.(cur_reg) ~= 0);
    fprintf(cur_reg + "\n")
    path_region = fullfile(path, "region-" + cur_reg);
    if (~exist(path_region, "dir")), mkdir(path_region), end
    for i_trial = progress(1:len_trial)
        cur_trial_dat = grp_data(chan_in_reg, :, i_trial, :);
        fisher_z = cell(length(window_start), 1);
        parfor i_win = 1:length(window_start)
            % collapse channel and time (thus spatiotemporal pattern)
            cur_cor_mat = corr(reshape( ...
                cur_trial_dat(:, window_start(i_win):window_end(i_win), :, :), ...
                length(chan_in_reg) * size_window, [])); %#ok<*PFBNS>
            fisher_z{i_win} = atanh(cur_cor_mat(idx_keep_cors));
        end
        cur_simi_inter_by_window = utils.filter_triangular( ...
            utils.preallocate(1:length(window_start), subjs_id, subjs_id, ...
            VariableNames=["winow_id", "subj_id_col", "subj_id_row"]));
        cur_simi_inter_by_window.fisher_z = vertcat(fisher_z{:});
        path_trial = fullfile(path_region, "trial-" + string(i_trial));
        if (~exist(path_trial, "dir")), mkdir(path_trial), end
        parquetwrite(fullfile(path_trial, "rs.parquet"), cur_simi_inter_by_window)
    end
end
clearvars cur_simi_inter_by_window

% type: individual to group similarity
fprintf("Processing stepped window intersubject similarity...\n")
path = fullfile("data", "type-group_acq-window");
for i_region = regions_id
    cur_reg = "region" + string(i_region);
    chan_in_reg = channel.code(channel.(cur_reg) ~= 0);
    fprintf(cur_reg + "\n")
    path_region = fullfile(path, "region-" + cur_reg);
    if (~exist(path_region, "dir")), mkdir(path_region), end
    for i_trial = progress(1:len_trial)
        cur_trial_dat = grp_data(chan_in_reg, :, i_trial, :);
        fisher_z = cell(length(window_start), 1);
        parfor i_win = 1:length(window_start)
            % collapse channel and time (thus spatiotemporal pattern)
            cur_dat = reshape( ...
                cur_trial_dat(:, window_start(i_win):window_end(i_win), :, :), ...
                length(chan_in_reg) * size_window, []);
            fisher_z{i_win} = utils.calc_simi_ind_to_grp(cur_dat, FisherZ=true);
        end
        cur_simi_grp_by_window = ...
            utils.preallocate(1:length(window_start), subjs_id, ...
            VariableNames=["window_id", "subj_id"]);
        cur_simi_grp_by_window.fisher_z = vertcat(fisher_z{:});
        path_trial = fullfile(path_region, "trial-" + string(i_trial));
        if (~exist(path_trial, "dir")), mkdir(path_trial), end
        parquetwrite(fullfile(path_trial, "rs.parquet"), cur_simi_grp_by_window)
    end
end
clearvars cur_simi_grp_by_window
