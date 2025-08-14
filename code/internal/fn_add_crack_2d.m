function [mod, el_cents, ep, en, crack_nd_i] = fn_add_crack_2d(mod, crack_vtcs, varargin)
if isempty(varargin)
    cod = 0;
else
    cod = varargin{1}
end

%2d crack assumed to be line of coordinates in order!
crack_fcs = [1:size(crack_vtcs, 1) - 1; 2:size(crack_vtcs, 1)]';

%Get signed distance of each element from crack
el_cents = fn_calc_element_centres(mod.nds, mod.els);
d = fn_signed_dist_to_bdry(el_cents, crack_vtcs, crack_fcs);

%identify elements on either side of crack and within tol of crack
[min_el_size, max_el_size] = fn_get_min_max_element_sizes(mod);

e0 = abs(d) < max_el_size * 2; %this is purely for computational efficiency to minimuise number of elements considered
ep = (sign(d)) >= 0 & e0;
en = (sign(d)) <  0 & e0;

%List all edges of all the pos and neg side elements
ep_edges = fn_edges(mod.els(ep, :));
en_edges = fn_edges(mod.els(en, :));

%find elements in ep that share edge with an element in en as these will
%make up mesh edges that define crack
k = 1;
crack_nd_i = zeros(size(ep_edges, 1), 2);
for i = 1:size(ep_edges, 1)
    j = all(ep_edges(i, :) == en_edges, 2);
    if any(j)
        crack_nd_i(k, :) = ep_edges(i, :);
        k = k + 1;
    end
end
crack_nd_i = crack_nd_i(1:k - 1, :);
crack_nd_i = unique(crack_nd_i(:)); %these are the nodes that need to be duplicated
crack_nds = mod.nds(crack_nd_i, :);
% crack_nds = crack_nds + randn(size(crack_nds)) * cod / 10; %this adds a small random perturbation to crack nodes if a COD is required to minimise chances of any falling exactly on crack line as this means cod-direction

%Identify any crack nodes where nearest point on desired crack is one of
%end points
[d, nearest_pts, norm_vecs] = fn_signed_dist_to_bdry(crack_nds, crack_vtcs, crack_fcs);
d1 = (sqrt(sum((nearest_pts - crack_vtcs(1,:)) .^ 2, 2)) < eps) | ...
     (sqrt(sum((nearest_pts - crack_vtcs(end,:)) .^ 2, 2)) < eps);
%Find any crack nodes where nearest point is within half an element of end
%nodes (because crack will effectively extend to next node beyond last
%disconnected one)
d2 = (sqrt(sum((nearest_pts - crack_vtcs(1,  :)) .^ 2, 2)) < min_el_size / 2) | ...
     (sqrt(sum((nearest_pts - crack_vtcs(end,:)) .^ 2, 2)) < min_el_size / 2);


%Remove these points
crack_nd_i(d1 | d2) = [];
crack_nds(d1 | d2, :) = [];
norm_vecs(d1 | d2, :) = [];

mod.nds(crack_nd_i, :) = mod.nds(crack_nd_i, :) + norm_vecs .* cod;
crack_nds = crack_nds - norm_vecs .* cod;

%duplicate crack nodes
new_node_indices = [1:numel(crack_nd_i)]' + size(mod.nds,1);
mod.nds = [mod.nds; crack_nds];

%finally loop through crack_nds and for any occurences inelements on -ve 
%side, change to equivalent new nd
tmp = mod.els(en, :);
for i = 1:numel(crack_nd_i)
    tmp(tmp == crack_nd_i(i)) = new_node_indices(i);
end
mod.els(en, :) = tmp;
end

function edges = fn_edges(els)
edges = sort([els(:,[1,2]); els(:,[2,3]); els(:,[3,1])], 2);
end
