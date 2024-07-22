function v = fn_get_field_output(f_out, mod, mats, varargin)
%SUMMARY
%   Converts raw field output (displacement or velocity time-history at
%   every node) to single value per element to be used in animations etc.
%   depending on selected field_output_type
%INPUTS
%   f_out - raw field output from model
%   mod - model description with mod.nds and mod.els fields
%   mats - global matrices with mats.m and mats,gl_lookup fields
%   [field_output_type - type of field output:
%       'sqrt(KE)' - square root of kinetic energy of element (default)
%       'mean(ui)' - (i = 1,2,3) mean of displacement in i direction
%       'mean(vi)' - (i = 1,2,3) mean of velocity in i direction
%       'div(u)' or 'div(v)' - approximate divergence of displacement or 
%       velocity over element
%       'curl(u)' or 'curl(v)' - approximate curl of displacement or 
%       velocity over element (2D only)
%       'raw(u)' or 'raw(v)' - raw displacement of velocity field output per node
%OUTPUTS
%   v - field output of chosen quantity per element (unless 'raw' selected)
%--------------------------------------------------------------------------
if isempty(varargin)
    field_output_type = 'sqrt(KE)';
else
    field_output_type = varargin{1};
end

switch field_output_type
    case {'mean(u1)', 'mean(u2)', 'mean(u3)', 'curl(u)', 'div(u)', 'mean(v1)', 'mean(v2)', 'mean(v3)', 'curl(v)', 'div(v)'}
        %components at element nodes
        nodal_values = fn_convert_from_global_to_nodal(f_out, mats.gl_lookup);

        %extract appropriate value for each elemeent
        switch field_output_type
            case {'mean(u1)', 'mean(u2)', 'mean(u3)', 'mean(v1)', 'mean(v2)', 'mean(v3)'}
                k = str2num(field_output_type(7));
                v = fn_average_over_nodes_of_element(nodal_values(:,:,k), mod.els);
            case {'curl(u)', 'div(u)', 'curl(v)', 'div(v)'}
                %NOTE these are not strictly div or curl. They should be
                %proportional if elements are regular (probably need to
                %divide by area / vol of element) but will not be exactly
                %proportional for irregular elements
                switch field_output_type
                    case {'div(u)', 'div(v)'}
                        unit_vecs = fn_get_unit_vecs_div(mod);
                    case {'curl(u)', 'curl(v)'}
                        if size(mod.nds) > 2
                            error('Cannot output curl for 3D models')
                        end
                        unit_vecs = fn_get_unit_vecs_curl(mod);
                end
                %sum of dot products around each element is the output
                v = zeros(size(mod.els, 1), size(f_out,2));
                for i = 1:size(mod.els, 2)
                    j = mod.els(:, i);
                    j(~j) = 1; %avoid indexing error - should be OK because unit vecs are zero for these ones anyway
                    v = v + sum(nodal_values(j, : , 1:size(unit_vecs, 3)) .* unit_vecs(:, i, :), 3);
                end
        end
        
    case {'raw(u)', 'raw(v)'}
        v = field_output_type;
    
    otherwise %KE is the default

        %Get KEs at each global DOF
        KE_g = zeros(size(f_out));
        diagM = diag(mats.M);
        for i = 1:size(f_out,2)
            KE_g(:,i) = abs(f_out(:,i) .^ 2 .* diagM) / 2;
        end

        %Need to convert them first to nodal values
        KE_n = fn_convert_from_global_to_nodal(KE_g, mats.gl_lookup);
        KE_n = sum(KE_n, 3);

        %Finally to sqrt(KE) for element
        v = sqrt(fn_average_over_nodes_of_element(KE_n, mod.els));
end

end

%--------------------------------------------------------------------------

function unit_vecs = fn_get_unit_vecs_div(mod)
%create unit vectors pointing out at each element node
unit_vecs = zeros([size(mod.els), size(mod.nds, 2)]);
valid_els = ones(size(mod.els, 1), 1);
for i = 1:size(mod.els, 2)
    j = mod.els(:, i);
    valid_els = valid_els & j;
    j(~j) = 1;
    unit_vecs(:, i, :) = mod.nds(j, :);
end
el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
unit_vecs = unit_vecs - permute(el_ctrs, [1, 3, 2]);
unit_vecs = unit_vecs ./ sqrt(sum(unit_vecs .^ 2, 3));
unit_vecs(~valid_els, :) = 0;
end

function unit_vecs = fn_get_unit_vecs_curl(mod)
%rotate outward unit vectors by 90degs to get curl ones
unit_vecs = fn_get_unit_vecs_div(mod);
unit_vecs = cat(3, unit_vecs(:,:,2), -unit_vecs(:,:,1));
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