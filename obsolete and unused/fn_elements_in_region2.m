function [in, out] = fn_elements_in_region2(nds, els, region)
%SUMMARY
%   Returns logical n_els x 1 vectors indicating whether elements in model
%   are inside or outside the specified region.

el_centres = fn_calc_element_centres(nds, els);
in = inpolygon(el_centres(:,1), el_centres(:,2), region(:,1), region(:,2));
out = ~in;
end