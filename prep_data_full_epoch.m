path_src = fullfile("data", "SM_RSA_rawEEG_preprocessed");
files_src = ls(fullfile(path_src, "*.mat"));
eeg = cell(1, size(files_src, 1));
subjs = strings(1, size(files_src, 1));
parfor i_file = 1:size(files_src, 1)
    file_src = fullfile(path_src, files_src(i_file, :));
    subjs(i_file) = extract(file_src, lookBehindBoundary("sub") + digitsPattern);
    eeg_src = load(file_src, "EEG_replace_badt")
    eeg{i_file} = eeg_src.EEG_replace_badt;
end
subjs = subjs(~cellfun(@isempty, eeg));
eeg = cat(4, eeg{:});
% save(fullfile("data", "EEG", "combined"), "eeg", "subjs", "-v7.3")

load(fullfile("data", "cca_result_model.mat"), "W")
X = reshape(permute(eeg, [2, 3, 1, 4]), [], size(eeg, 1), size(eeg, 4));
N = size(X, 3);
for l=N:-1:1, Y(:,:,l)=X(:,:,l)*W; end

% save(fullfile("data", "CorCA", "model_subjs207"), "X", "Y", "W", "-v7.3")

