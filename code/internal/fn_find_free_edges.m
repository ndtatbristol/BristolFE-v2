function free_ed = fn_find_free_edges(els)
% %add dummy edges at start and end so subsequent logic is general
% ed = [[0, 0]; ed; [0, 0]];
% ind = find((ed(2:end-1, 1) ~= ed(1:end-2, 1) | ed(2:end-1, 2) ~= ed(1:end-2, 2)) & (ed(2:end-1, 1) ~= ed(3:end, 1) | ed(2:end-1, 2) ~= ed(3:end, 2)));
% free_ed = ed(ind + 1, :);


[free_ed, els_with_ed] = fn_get_edges(els);
j = (sum(els_with_ed > 0, 2) == 1) & (sum(free_ed > 0, 2) == 2);
free_ed = free_ed(j, :);
end