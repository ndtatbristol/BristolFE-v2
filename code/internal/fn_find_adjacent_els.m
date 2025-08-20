function adj_els = fn_find_adjacent_els(els, nd1, nd2)
%returns elements from els with common edges to edge defined ny nd1 and nd2
%(which obviously must be adjacent nodes)
adj_els = find(sum(els == nd1, 2) & sum(els == nd2, 2));
end