function [res, mats] = fn_FE_entry_point(mod, matls, steps, fe_options, varargin)
if numel(varargin) < 1
    solver = 'BristolFE';
else
    solver = varargin{1};
end

switch solver
    case 'Abaqus'
        fn_ABAQUS(mod, matls, steps, fe_options);
        res = [];
        mats = [];
    otherwise
        [res, mats] = fn_BristolFE_v2(mod, matls, steps, fe_options);
end


end



