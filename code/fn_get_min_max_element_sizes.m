function [min_el_size, max_el_size] = fn_get_min_max_element_sizes(mod)
%USAGE
%   [min_el_size, max_el_size] = fn_get_min_max_element_sizes(mod)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns minimum and maximum side lengths of elements in a model
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, element nodes, mod.els, element types,
%   mod.el_types, and element type indices, mod.el_typ_i
%OUTPUT
%   min_el_size - minimum element side length found in model
%   max_el_size - maximum element side length found in model
%--------------------------------------------------------------------------

el_edges = fn_edges_from_els(mod.els, 1:size(mod.els, 1), mod.el_typ_i, mod.el_types);
min_el_size = inf;
max_el_size = 0;
for i = 1:numel(el_edges)
    r = sqrt(sum((mod.nds(el_edges{i}.eds(:,1), :) - mod.nds(el_edges{i}.eds(:,2), :)) .^ 2, 2));
    min_el_size = min(min_el_size, min(r));
    max_el_size = max(max_el_size, max(r));
end

end