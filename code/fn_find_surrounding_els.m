function i = fn_find_surrounding_els(els, nd)
%Returns elements surrounding node
[i, ~] = find(els == nd);

%Remove any interface elements - this is a hack only for 2D
j = sum(els(i, :) > 0, 2) < 3;
i(j) = [];
end
