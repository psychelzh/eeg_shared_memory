function [window_start, window_end] = setup_window(len, size, step)
%SETUP_WINDOW Setup window for further analysis
%   For each step, we need a duration (based on size) of timepoints to form
%   a window to calculate correlation.
arguments
    len (1, 1)
    size (1, 1) = 26
    step (1, 1) = 5
end
window_start = 1:step:len;
window_end = size:step:len;
window_start = window_start(1:length(window_end));
end

