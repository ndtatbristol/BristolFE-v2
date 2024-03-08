function [nds2, els2, old_nds, new_nds] = fn_remove_unused_nodes(nds, els)
%list of original node numbers
new_nds = [1:size(nds, 1)]';
old_nds = [1:size(nds, 1)]';
in_use = zeros(size(old_nds));

u = unique(els(:));
u(u == 0 | isnan(u)) = [];
in_use(u) = 1;
in_use = find(in_use);
nds2 = nds(in_use, :);
old_nds = old_nds(in_use);

nd_nos_new = [1:numel(old_nds)]';
els2 = interp1(old_nds, nd_nos_new, els, 'nearest');

new_nds = zeros(size(nds, 1), 1);
new_nds(old_nds) = nd_nos_new;

end