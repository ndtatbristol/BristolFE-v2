function free_ed_nds = fn_find_free_edge_nds(free_ed)
%works though list of free edges and gets ordered list of nds
free_ed_nds = zeros(size(free_ed, 1), 1);
free_ed_no = 1;
end_no = 1;
other_end_no = [2, 1];
for i = 1:size(free_ed, 1)
    %add node at one end of edge to the list
    free_ed_nds(i) = free_ed(free_ed_no, end_no);
    %get node at other end of edge
    k = free_ed(free_ed_no, other_end_no(end_no));
    %zero the edge being processed so not picked up later
    free_ed(free_ed_no, :) = 0;
    %find the other free edge containing node j
    [free_ed_no, end_no] = find(free_ed == k);
end
    
end