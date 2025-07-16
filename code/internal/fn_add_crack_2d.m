function [mod, el_cents, ep, en, crack_nds] = fn_add_crack_2d(mod, crack_vtcs, tol)

%This works but
%   tol needs to be at least el_size, but currently crack length is
%   effectively increasd by this. Really need two separate things:
%       tol that can be set conservatively that just determines number of
%       elements in crack vicinity to be considered, and
%       a separate calc to remove crack nodes that are beyond project of
%       desired crack line. Modify fn_signed_dist_to_bdry to return some
%       extra parameters, e.g. nearest point on bdry - then drop all crack
%       nds where nearest point is one of crack end vertices?

%2d crack assumed to be line of coordinates in order!
crack_fcs = [1:size(crack_vtcs, 1) - 1; 2:size(crack_vtcs, 1)]';

%Get signed distance of each element from crack
el_cents = fn_calc_element_centres(mod.nds, mod.els);
d = fn_signed_dist_to_bdry(el_cents, crack_vtcs, crack_fcs);

%identify elements on either side of crack and within tol of crack
e0 = abs(d) < tol;
ep = (sign(d)) > 0 & e0;
en = (sign(d)) < 0 & e0;

%List all edges of all the pos and neg side elements
ep_edges = fn_edges(mod.els(ep, :));
en_edges = fn_edges(mod.els(en, :));

%find elements in ep that share edge with an element in en as these will
%make up mesh edges that define crack
k = 1;
crack_nds = zeros(size(ep_edges, 1), 2);
for i = 1:size(ep_edges, 1)
    j = all(ep_edges(i, :) == en_edges, 2);
    if any(j)
        crack_nds(k, :) = ep_edges(i, :);
        k = k + 1;
    end
end
crack_nds = crack_nds(1:k - 1, :);
crack_nds = unique(crack_nds(:)); %these are the nodes that need to be duplicated

%duplicate crack nodes
new_node_indices = [1:numel(crack_nds)]' + size(mod.nds,1);
mod.nds = [mod.nds; mod.nds(crack_nds, :)];

%finally loop through crack_nds and for any occurences inelements on -ve 
%side, change to equivalent new nd
tmp = mod.els(en, :);
for i = 1:numel(crack_nds)
    tmp(tmp == crack_nds(i)) = new_node_indices(i);
end
mod.els(en, :) = tmp;
end

function edges = fn_edges(els)
edges = sort([els(:,[1,2]); els(:,[2,3]); els(:,[3,1])], 2);
end
