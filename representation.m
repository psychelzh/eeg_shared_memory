channel = readtable(fullfile("config", "eeg_channel_labels_64.csv"), TextType="string");
subjs_id = unique(readtable(fullfile("data", "group_task-wordencoding_events.csv")).subj_id, "stable");
load(fullfile("data", "grp_subjs206_nodemean_1000ms.mat"), "grp_data");
[~, len_time_point, len_trial, len_sub] = size(grp_data);
region_id = 1:6;

% debug
% region_id = 1:2;
% len_trial = 2;

%% intersubject similarity
numel_cor_mat = len_sub * len_sub;
[sub_row, sub_col] = ind2sub([len_sub, len_sub], 1:(numel_cor_mat));
simi_inter_subj = repmat(array2table( ...
    [subjs_id(sub_row), subjs_id(sub_col), nan(len_sub * len_sub, 1)], ...
    "VariableNames", ["subj_id_x", "subj_id_y", "r"]), ...
    length(region_id) * len_trial, 1);
% easier handle data store
r_start = 1;
r_end = numel_cor_mat;
for i_region = region_id
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
parquetwrite(fullfile("data", "type-inter_rsa.parquet"), simi_inter_subj)


%% individual to group similarity

