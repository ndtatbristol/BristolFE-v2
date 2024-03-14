function [in, out] = fn_elements_in_region(mod, region)
%SUMMARY
%   Returns logical n_els x 1 vectors indicating whether elements in model
%   are inside or outside the specified region.

in = inpolygon(mod.el_centres(:,1), mod.el_centres(:,2), region(:,1), region(:,2));
out = ~in;
end