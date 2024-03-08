function [res, mats] = fn_FE_entry_point(mod, matls, steps, fe_options)
default_options.solver = 'BristolFE';
fe_options = fn_set_default_fields(fe_options, default_options);


switch fe_options.solver
    case 'Abaqus'
        fn_ABAQUS(mod, matls, steps, fe_options);
        res = [];
        mats = [];
    case 'BristolFE'
        [res, mats] = fn_BristolFE_v2(mod, matls, steps, fe_options);
    otherwise
        error('Invalid FE solver');
end


end



