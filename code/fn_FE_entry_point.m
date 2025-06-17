function varargout = fn_FE_entry_point(mod, matls, steps, fe_options)
default_options.solver = 'BristolFE';
default_options.gl_mat_nds = [];
fe_options = fn_set_default_fields(fe_options, default_options);

%Set the solver
switch fe_options.solver
    case 'Abaqus'
        fn_solver = @fn_ABAQUS;
    case 'BristolFE'
        fn_solver = @fn_BristolFE_v2;
    case 'pogo'
        fn_solver = @fn_pogoFE;
    otherwise
        error('Invalid FE solver');
end

%Call selected solver with appropriate number of outputs
if nargout == 1
    varargout{1} = fn_solver(mod, matls, steps, fe_options);
elseif nargout == 2
    [varargout{1}, varargout{2}] = fn_solver(mod, matls, steps, fe_options);
end


end



