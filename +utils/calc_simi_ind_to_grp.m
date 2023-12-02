function simi = calc_simi_ind_to_grp(dat, opts)
%CALC_SIMI_IND_TO_GRP Calculate individual to group similarity
%   The group representation is calcualted by averaging all the subjects
%   except the current calculating subject.

arguments
    dat (:, :) double
    opts.FisherZ (1, 1) logical = true
end

ncols = size(dat, 2);
simi = nan(ncols, 1);
for i_col = 1:ncols
    cur_rep = dat(:, i_col);
    grp_rep = mean(dat(:, setdiff(1:ncols, i_col)), 2, "omitmissing");
    r = corr(cur_rep, grp_rep, rows="pairwise");
    if opts.FisherZ
        simi(i_col) = atanh(r);
    else
        simi(i_col) = r;
    end
end
end
