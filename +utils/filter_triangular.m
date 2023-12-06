function T_filtered = filter_triangular(T, opts)
%FILTER_TRIANGULAR Filter out triangular part of table
arguments
    T {mustBeA(T, 'table')}
    opts.RowName {mustBeTextScalar} = "subj_id_row"
    opts.ColName {mustBeTextScalar} = "subj_id_col"
    opts.Part {mustBeMember(opts.Part, ["upper", "lower"])} = "lower"
end
switch opts.Part
    case "lower"
        T_filtered = T(T.(opts.RowName) > T.(opts.ColName), :);
    case "upper"
        T_filtered = T(T.(opts.RowName) < T.(opts.ColName), :);
end
end
