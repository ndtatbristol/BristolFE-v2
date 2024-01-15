function el_centres = fn_calc_element_centres(nds, els)

% el_centres = [...
%     mean([nds(els(:,1), 1), nds(els(:,2), 1), nds(els(:,3), 1)], 2) ,...
%     mean([nds(els(:,1), 2), nds(els(:,2), 2), nds(els(:,3), 2)], 2)];

el_centres = 0;
n = 0;
for i = 1:size(els, 2)
    j = els(:, i);
    v = j > 0;
    n = n + v;
    j(~v) = 1;
    el_centres = el_centres + nds(j, :) .* v;
end

el_centres = el_centres ./ n;

end