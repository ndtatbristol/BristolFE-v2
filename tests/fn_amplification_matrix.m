function [A, max_EV] = fn_amplification_matrix(M, C, K, dt, fe_options)
default_options.solver_mode = 'vel at last half time step';
fe_options = fn_set_default_fields(fe_options, default_options);

ndf = size(K, 1);
switch fe_options.solver_mode
    case 'vel at last half time step'
        if isdiag(M)
            diag_M = spdiags(sum(M).', 0, ndf, ndf);
            inv_M = spdiags(1 ./ sum(M).', 0, ndf, ndf);
            B1 = 2 * speye(ndf) - dt * (inv_M * C) - dt ^ 2 * (inv_M * K);
            B2 = dt * (inv_M * C) - speye(ndf);
        else
            B1 = 2 * speye(ndf) - dt * (M \ C) - dt ^ 2 * (M \ K);
            B2 = dt * (M \ C) - speye(ndf);
        end
        A = [B1,         B2
             speye(ndf), sparse(ndf, ndf)];

    case 'vel at curent time step'
        if isdiag(M)
            diag_M = spdiags(sum(M).', 0, ndf, ndf);
            inv_M = spdiags(1 ./ sum(M).', 0, ndf, ndf);
            B1 = 2 * speye(ndf) - dt ^ 2 * inv_M * K;
            B2 =    -speye(ndf) + dt / 2 * inv_M * C;
            B3 =     speye(ndf) + dt / 2 * inv_M * C;

        else
            B1 = 2 * speye(ndf) - dt ^ 2 * (M \ K);
            B2 =    -speye(ndf) + dt / 2 * (M \ C);
            B3 =     speye(ndf) + dt / 2 * (M \ C);
        end
        A = [B3 \ B1,   B3 \ B2
            speye(ndf), sparse(ndf, ndf)];

end
if prod(size(A)) < 500 ^ 2
    max_EV = max(abs(eig(full(A))), [], 'all');
else
    max_EV = max(abs(eigs(A, 1)), [], 'all');
end
end