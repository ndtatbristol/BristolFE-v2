function [in, out] = fn_elements_in_region(nds, els, region)
el_centres = fn_calc_element_centres(nds, els);
in = inpolygon(el_centres(:,1), el_centres(:,2), region(:,1), region(:,2));
out = ~in;
end