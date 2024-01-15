function v = fn_get_plot_vals_v3(f_out, gl_lookup, els, M)

%In this version, one value per element is extracted (not per node). 

%Get KEs at each global DOF
KE_g = zeros(size(f_out));
diagM = diag(M);
for i = 1:size(f_out,2)
    KE_g(:,i) = abs(f_out(:,i) .^ 2 .* diagM);
end

%Need to convert them first to nodal values
KE_n = 0;
for i = 1:size(gl_lookup, 2) %loop over DoF
    j = gl_lookup(:, i);
    valid = j > 0;
    j(~valid) = 1; %to prevent indexing error
    KE_n = KE_n + KE_g(j, :) .* valid;
end

%Finally to element values
v = zeros(size(els, 1), size(f_out, 2));
for i = 1:size(els, 2)
    j = els(:, i);
    valid = j > 0;
    j(~valid) = 1;
    tmp = KE_n(j, :);
    tmp(~valid, :) = 0;
    v = v + tmp;
end

end