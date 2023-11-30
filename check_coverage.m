channel = readtable(fullfile("config", "eeg_channel_labels_64.csv"), TextType="string");
subjs_id = unique(readtable(fullfile("data", "group_task-wordencoding_events.csv")).subj_id, "stable");
load(fullfile("data", "grp_subjs206_nodemean_1000ms.mat"), "grp_data");
[~, len_time_point, len_trial, len_subj] = size(grp_data);
regions_id = 1:6;

[regions_id_idx, trial_id_idx, subj_id_idx] = ...
    meshgrid(regions_id, 1:len_trial, 1:len_subj);
coverage = array2table( ...
    [regions_id_idx(:), trial_id_idx(:), subjs_id(subj_id_idx(:)), ...
    nan(numel(regions_id_idx), 1)], ...
    VariableNames=["region_id", "trial_id", "subj_id", "prop"]);

start_idx = 1;
end_idx = len_trial * len_subj;
for i_region = regions_id
    chan_in_reg = channel.code(channel.("region" + string(i_region)) ~= 0);
    cur_dat = grp_data(chan_in_reg, :, :, :);
    prop_missing = mean(any(isnan(cur_dat), 2), 1);
    coverage.prop(start_idx:end_idx) = prop_missing(:);
    start_idx = start_idx + len_trial * len_subj;
    end_idx = end_idx + len_trial * len_subj;
end
