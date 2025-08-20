function nds = fn_find_surrounding_nds(nd, els)
%Returns nodes surrounding nd
i = fn_find_surrounding_els(els, nd);
nds = els(i, :);
nds = nds(:);
nds(nds == nd) = [];
nds = unique(nds);
nds(nds == 0) = [];
end