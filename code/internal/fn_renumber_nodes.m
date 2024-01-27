function [nds2, els2, old_nd_i] = fn_renumber_nodes(nds, els, new_nd_i)
%Reorder node coordinate list and original node index list
nds2 = nds(new_nd_i, :);

%List of original node numbers
old_nd_i = [1:size(nds, 1)]';

%Update elements
z = find(els == 0);
els2 = interp1(new_nd_i, old_nd_i, els, 'nearest'); %NB counter-intuitive way around!
els2(z) = 0;
end