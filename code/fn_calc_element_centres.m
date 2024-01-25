function el_ctrs = fn_calc_element_centres(nds, els)
%SUMMARY
%   Returns n_els x 2 matrix of the coordinates of all element centres in a
%   model.

el_ctrs = 0;
n = 0;
for i = 1:size(els, 2)
    j = els(:, i);
    v = j > 0;
    n = n + v;
    j(~v) = 1;
    el_ctrs = el_ctrs + nds(j, :) .* v;
end

el_ctrs = el_ctrs ./ n;

end