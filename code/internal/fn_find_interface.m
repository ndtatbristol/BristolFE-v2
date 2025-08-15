function [interface_facets, e1, e2] = fn_find_interface(mod, el_i1, el_i2)
%Returns matrix of 2-node line segments (2D) or n-node polygons (3D) at
%interface between mod.els(el_i1,:) and mod.els(el_i2,:). Node numbering of
%facets is NOT guaranteed to be consistent

%Used by fn_add_fluid_solid_interface to generate nodes of interface
%elements and fn_add_crack_{2/3}d to indentify facets of crack surface

if any(el_i1 & el_i2)
    error('Element sets must not contain any common elements')
end
if ~any(el_i1) || ~any(el_i2)
    warning('One or both element sets is empty')
    interface_facets = [];
    e1 = [];
    e2 = [];
    return
end

%Get lists of faces for both element sets

fcs(1).data = fn_faces_from_els(mod.els(el_i1,:), mod.el_typ_i(el_i1), mod.el_types);
fcs(2).data = fn_faces_from_els(mod.els(el_i2,:), mod.el_typ_i(el_i2), mod.el_types);

%Loop over every face in set 1 and flag any that are in set 2
% max_nds_per_facet = 0; %also get maximum number of nodes per facet
% no_facets = 0;
% for i1 = 1:numel(fcs(1).data)
%     max_nds_per_facet = max(max_nds_per_facet, size(fcs(1).data{i1}.fcs, 2));
%     fcs(1).data{i1}.el_in_set2 = zeros(size(fcs(1).data{i1}.fcs, 1), 1);
%     for j1 = 1:size(fcs(1).data{i1}.fcs, 1)
%         fc1 = sort(fcs(1).data{i1}.fcs(j1, :), 2);
%         for i2 = 1:numel(fcs(2).data)
%             if size(fcs(2).data{i2}.fcs, 2) ~= numel(fc1)
%                 continue
%             end
%             fc2 = sort(fcs(2).data{i1}.fcs, 2);
%             k1 = all(fc2 == fc1, 2);
%             if nnz(k1) > 0
%                 fcs(1).data{i1}.el_in_set2(j1) = fcs(2).data{i2}.el_i(k1);
%             end
%         end
%     end
%     no_facets = no_facets + nnz(fcs(1).data{i1}.el_in_set2);
% end

%new version
max_nds_per_facet = 0; %also get maximum number of nodes per facet
no_facets = 0;
for i1 = 1:numel(fcs(1).data)
    max_nds_per_facet = max(max_nds_per_facet, size(fcs(1).data{i1}.fcs, 2));
    fcs(1).data{i1}.el_in_set2 = zeros(size(fcs(1).data{i1}.fcs, 1), 1);
    for i2 = 1:numel(fcs(2).data)
        if size(fcs(1).data{i1}.fcs, 2) ~= size(fcs(2).data{i2}.fcs, 2)
            continue
        end
        [k1, k2] = ismember(sort(fcs(1).data{i1}.fcs, 2), sort(fcs(2).data{i2}.fcs, 2), 'rows');
        fcs(1).data{i1}.el_in_set2(k1) = fcs(2).data{i2}.el_i(k2(k1));
    end
    no_facets = no_facets + nnz(fcs(1).data{i1}.el_in_set2);
end

%Pull out the facets
interface_facets = zeros(no_facets, max_nds_per_facet);
e1 = zeros(no_facets, 1);
e2 = zeros(no_facets, 1);
k1 = 1;
for i1 = 1:numel(fcs(1).data)
    j1 = fcs(1).data{i1}.el_in_set2 > 0;
    tmp = fcs(1).data{i1}.fcs(j1, :);
    k2 = k1 + size(tmp, 1) - 1;
    interface_facets(k1:k2, 1:size(tmp, 2)) = tmp;
    e1(k1:k2) = fcs(1).data{i1}.el_i(j1);
    e2(k1:k2) = fcs(1).data{i1}.el_in_set2(j1);
    k1 = k2 + 1;
end
end