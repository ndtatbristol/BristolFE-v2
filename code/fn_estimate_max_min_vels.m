function [mx, mn] = fn_estimate_max_min_vels(matls)
v = [0,0];
j = 1;
for i = 1:numel(matls)
    s = [];
    if isfield(matls(i), 'stiffness_matrix')
        s = matls(i).stiffness_matrix(:);
    end
    if isfield(matls(i), 'D')
        s = matls(i).D(:);
    end
    if isfield(matls(i), 'density') %deal with legacy naming
        rho = matls(i).density;
    else
        rho = matls(i).rho;
    end

    if ~isempty(s)
        s = s(abs(s) > 0);
        v(j,1) = sqrt(min(s) / rho);
        v(j,2) = sqrt(max(s) / rho);
        j = j + 1;
    end
end
mn = min(v(:,1));
mx = max(v(:,2));
end