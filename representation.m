channel = readtable(fullfile("config", "eeg_channel_labels_64.csv"), TextType="string");
subjs_id = unique(readtable(fullfile("data", "group_task-wordencoding_events.csv")).subj_id, "stable");
load(fullfile("data", "grp_subjs206_nodemean_1000ms.mat"), "grp_data");
[~, len_time_point, len_trial, len_sub] = size(grp_data);
regions_id = 1:6;

% debug
% regions_id = 1:2;
% len_trial = 2;

%% intersubject similarity
[regions_id_idx, trial_id_idx] = meshgrid(regions_id, 1:len_trial);
numel_region_trial = numel(regions_id_idx);
numel_cor_mat = len_sub * len_sub;
[sub_row, sub_col] = ind2sub([len_sub, len_sub], 1:(numel_cor_mat));
simi_inter_subj = array2table(...
    [repelem([regions_id_idx(:), trial_id_idx(:)], numel_cor_mat, 1), ...
    repmat([subjs_id(sub_row), subjs_id(sub_col), nan(numel_cor_mat, 1)], numel_region_trial, 1)], ...
    "VariableNames", ["region_id", "trial_id", "subj_id_row", "subj_id_col", "r"]);
% easier handle data store
r_start = 1;
r_end = numel_cor_mat;
for i_region = regions_id
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    for i_trial = 1:len_trial
        % collapse channel and time (thus spatiotemporal pattern)
        cur_dat = reshape(grp_data(chan_in_reg, :, i_trial, :), ...
            [length(chan_in_reg) * len_time_point, len_sub]);
        cur_cor_mat = corr(cur_dat);
        simi_inter_subj.r(r_start:r_end) = cur_cor_mat(:);
        r_start = r_start + numel_cor_mat;
        r_end = r_end + numel_cor_mat;
    end
end
parquetwrite(fullfile("data", "task-rs_acq-inter_type-full.parquet"), simi_inter_subj)
parquetwrite(fullfile("data", "task-rs_acq-inter_type-triu.parquet"), ...
    simi_inter_subj(simi_inter_subj.subj_id_row < simi_inter_subj.subj_id_col, :))

%% individual to group similarity

