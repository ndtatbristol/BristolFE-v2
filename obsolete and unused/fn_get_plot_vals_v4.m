function v = fn_get_plot_vals_v4(f_out, gl_lookup, nds, els, M, field_output)

%In this version, one value per element is extracted (not per node) for
%compatibility with animation routine which animate element patch colours

switch field_output
    case {'v1', 'v2', 'v3', 'element curl', 'element div'}
        %vel components at element nodes
        nv = fn_convert_from_global_to_nodal(f_out, gl_lookup);

        %extract appropriate value for each elemeent
        switch field_output
            case {'v1', 'v2', 'v3'}
                k = str2num(field_output(2));
                v = fn_average_over_nodes_of_element(nv(:,:,k), els);
            case {'element curl', 'element div'}
                %velocity at element centre
                v_bar = fn_average_over_nodes_of_element(nv, els);
                el_ctrs = fn_calc_element_centres(nds, els);

                %create unit vectors pointing out at each element node
                unit_vecs = zeros([size(els), size(nds, 2)]);
                for i = 1:size(els, 2)
                    unit_vecs(:, i, :) = nds(els(:, i), :);
                end
                unit_vecs = unit_vecs - permute(el_ctrs, [1, 3, 2]);
                unit_vecs = unit_vecs ./ sqrt(sum(unit_vecs .^ 2, 3));
                
                switch field_output
                    case 'element div'
                        %do nothing, unit_vecs are already outward
                    case 'element curl'
                        %turn unit_vecs through 90 degs
                        unit_vecs = cat(3, unit_vecs(:,:,2), -unit_vecs(:,:,1));
                end
                
                %put sum of dot products around each element as the output
                v = zeros(size(els, 1), size(f_out,2));
                for i = 1:size(els, 2)
                    v = v + sum(nv(els(:, i), : , 1:size(unit_vecs, 3)) .* unit_vecs(:, i, :), 3);
                end
                %NOTE these are not strictly div or curl. They should be
                %proportional if elements are regular (probably need to
                %divide by area / vol of element) but will not be exactly
                %proportional for irregular elements
        end
        
    case 'element KE'


        %Get KEs at each global DOF
        KE_g = zeros(size(f_out));
        diagM = diag(M);
        for i = 1:size(f_out,2)
            KE_g(:,i) = abs(f_out(:,i) .^ 2 .* diagM);
        end

        %Need to convert them first to nodal values
        KE_n = fn_convert_from_global_to_nodal(KE_g, gl_lookup);
        KE_n = sum(KE_n, 3);

        %Finally to element values
        v = fn_average_over_nodes_of_element(KE_n, els);
end

end

%--------------------------------------------------------------------------

function nv = fn_convert_from_global_to_nodal(gv, gl_lookup)
nv = zeros(size(gl_lookup, 1), size(gv, 2), size(gl_lookup, 2));
for i = 1:size(gl_lookup, 2) %loop over DoF
    j = gl_lookup(:, i);
    valid = j > 0;
    j(~valid) = 1; %to prevent indexing error
    nv(:, :, i) = gv(j, :) .* valid;
end
end

%--------------------------------------------------------------------------

function ev = fn_average_over_nodes_of_element(nv, els)
sz = size(nv);
nv = reshape(nv, [sz(1), prod(sz(2:end))]); %flatten
ev = zeros(size(els, 1), size(nv, 2));
n = 0;
for i = 1:size(els, 2)
    j = els(:, i);
    valid = j > 0;
    j(~valid) = 1;
    ev = ev + nv(j, :) .* valid;
    n = n + double(els(:, i) > 0);
end
ev = ev ./ n;
ev = reshape(ev, [size(ev, 1), sz(2:end)]);
end