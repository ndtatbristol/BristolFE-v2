function [ed, els_with_ed] = fn_get_edges(els)
%Add extra return arguments with element indices associated with each edge
%- TOFINISH!
ed = zeros(numel(els), 2);
el_i = zeros(numel(els), 1);

%First just extract a list of all edges for all elements (which will
%contain duplicates of any common edges). Also record the elements to which
%the edges belong in el_i.
kk = 1;
for ii = 1:size(els, 1)
    tmp = els(ii, :);
    tmp = tmp(find(tmp));
    tmp = [tmp, tmp(1)];
    tmp_ed = [tmp(1:end-1); tmp(2:end)]';
    ed(kk: kk + size(tmp_ed,1) - 1, :) = tmp_ed;
    el_i(kk: kk + size(tmp_ed,1) - 1) = ii;
    kk = kk + size(tmp_ed,1);
end

%Trim the lists
ed = ed(1: kk - 1, :);
el_i = el_i(1: kk - 1);
%Sort the edges so that nodes are increasing order for each edge (sort(ed, 2)) and then
%by order of lowest node (sortrows(...))
[ed, j] = sortrows(sort(ed, 2));
el_i = el_i(j);

if nargout > 1
    [ed, ~, ci] = unique(ed, 'rows');
    els_with_ed = sparse(ci, el_i, el_i);
    els_with_ed = sort(els_with_ed, 2, 'descend');
    els_with_ed = full(els_with_ed(:, 1:max(find(sum(els_with_ed) > 0))));
end

end