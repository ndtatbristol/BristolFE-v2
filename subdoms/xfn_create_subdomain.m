function dom = xfn_create_subdomain(mod, matls, inner_bdry, abs_layer_thick, varargin)

%Last optional argument allows the start of the absorbing region to be
%specified. If not specified it will be penultimate boundary layer.

interface_el_name = 'ASI2D2';

dom.mod.nds = mod.nds;
dom.mod.els = mod.els;
dom.mod.el_mat_i = mod.el_mat_i;
dom.mod.el_abs_i = mod.el_abs_i;
dom.mod.el_typ_i = mod.el_typ_i;
dom.mod.bdry_lyrs = zeros(size(mod.nds, 1), 1);
dom.mod.el_abs_i = zeros(size(mod.els, 1), 1);

%add 4 layers of elements to get the bdry set
max_abs = 0;
edge_break_detected = 0;
warning_given = 0;
for i = 1:4
    if i == 1
        %identify nodes on inner bdry first
        [bdry_nds, free_ed_nds] = fn_nds_on_bdry(inner_bdry, dom.mod.nds, dom.mod.els);
    else
        %work out from there to get next 3
        [bdry_nds, layer_els, edge_break_detected] = fn_get_next_layer(dom.mod.nds, dom.mod.els, bdry_nds, free_ed_nds);
    end
        dom.mod.bdry_lyrs(bdry_nds) = i;
        if edge_break_detected
            error('Boundary layers break edge of model');
        end
    if i == 3 %for 3rd one, use boundary as start of absorbing region
        if isempty(varargin)
            abs_layer_start_bdry = [dom.mod.nds(bdry_nds, 1), dom.mod.nds(bdry_nds, 2)];
        else
            abs_layer_start_bdry = varargin{1};
        end
    end
end

%Delete original interface elements and then regenerate later in this 
%function - this is to avoid potential instability at edge of domain where 
%original interface elements copied from main mesh may now be on free edges
els_in_use = ~strcmp(dom.mod.el_typ_i, interface_el_name);
[dom.mod.main_el_i, ~, dom.mod.els, dom.mod.el_mat_i, dom.mod.el_abs_i, dom.mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dom.mod.els, dom.mod.el_mat_i, dom.mod.el_abs_i, dom.mod.el_typ_i);


%Add the absorbing layers by working out from centre of region
[~, cand_els] = fn_elements_in_region2(dom.mod.nds, dom.mod.els, abs_layer_start_bdry);
dom.mod.el_abs_i(cand_els) = fn_dist_point_to_bdry_2D(fn_calc_element_centres(dom.mod.nds, dom.mod.els(cand_els, :)), abs_layer_start_bdry) / abs_layer_thick;
els_in_use = ones(size(dom.mod.els, 1), 1);
els_in_use(dom.mod.el_abs_i > 1) = 0;

[dom.mod.main_el_i, ~, dom.mod.els, dom.mod.el_mat_i, dom.mod.el_abs_i, dom.mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dom.mod.els, dom.mod.el_mat_i, dom.mod.el_abs_i, dom.mod.el_typ_i);
[dom.mod.nds, dom.mod.els, old_nds, ~] = fn_remove_unused_nodes(dom.mod.nds, dom.mod.els);
dom.mod.main_nd_i = old_nds;
dom.mod.bdry_lyrs = dom.mod.bdry_lyrs(old_nds);

dom.mod = fn_add_fluid_solid_interface_els(dom.mod, matls);

free_ed = fn_find_free_edges(dom.mod.els);

dom.mod.outer_bndry_pts = [dom.mod.nds(free_ed, 1), dom.mod.nds(free_ed, 2)];
dom.mod.inner_bndry_pts = inner_bdry;

dom.mod = rmfield(dom.mod, 'main_el_i');%not needed anywhere?
end

%--------------------------------------------------------------------------

function [layer_nds, layer_els, edge_break_detected] = fn_get_next_layer(nds, els, bdry_nds, free_ed_nds)
%Returns ordered list of nodes and elements on next layer out from bdry_nds

%Find first node on next layer - either adjoining edge node if bdry started
%on edge of model on any (external) adjoinng node if it was internal bdry
nd_count = 1;
ly_count = 1;
surr_nds = [];
b = 1;
while isempty(surr_nds)
    surr_nds = fn_find_surrounding_nds(bdry_nds(b), els);
    out = ~inpolygon(nds(surr_nds, 1), nds(surr_nds, 2), nds(bdry_nds, 1), nds(bdry_nds, 2));
    surr_nds = surr_nds(out);
    b = b + 1;
end
j = find(ismember(surr_nds, free_ed_nds));
if any(j) %Case boundary starts at edge
    %work out which of the two is outside and make that the first
    %node in the new layer. First element is the one with those
    %nodes in common
    first_nd = surr_nds(j(1));
else %Case boundary is internal
    %pick any of the surrounding nodes
    first_nd = surr_nds(1);
end

%Find all elements adjoining bndry
[layer_els, ~] = find(ismember(els, bdry_nds));
%Find elements outside bdry
[~, out] = fn_elements_in_region2(nds, els(layer_els,:), nds(bdry_nds, :));
layer_els = layer_els(out);

%Find associated nodes
layer_nds = els(layer_els, :);
layer_nds = layer_nds(:);
layer_nds = setdiff(layer_nds, bdry_nds); 

%OK to here - now need to sort nds (and elements although that is less
%critical

tmp = zeros(size(layer_nds));
tmp(1) = first_nd;
layer_nds(layer_nds == first_nd) = 0;
edge_break_detected = 0;
for i = 2:numel(layer_nds)
    last_nd = tmp(i - 1);
    %Next node has to be (a) one of the surrounding nodes, (b) in the list
    %of remaining nodes
    surr_nds = fn_find_surrounding_nds(last_nd, els);
    [next_nd, j] = intersect(layer_nds, surr_nds);
    if isempty(next_nd)
        %This can occur either on an internal boundary where it is the
        %correct stopping criterion because the next node is the first node
        %OR if the boundary intersects the model edge in a second place.
        if any(layer_nds)
            q = layer_nds(layer_nds > 0);
            r2 = sum((nds(q, :) - nds(last_nd, :)) .^ 2, 2);
            [~, j] = min(r2);
            next_nd = q(j); %restart at nearest remaining one to last one - probably not very reliable
            edge_break_detected = 1;
        else
            break %case for internal bdries
        end
    end
    tmp(i) = next_nd(1);
    layer_nds(j(1)) = 0;
end
layer_nds = tmp;
layer_nds(layer_nds == 0) = [];
end

function [bdry_nds, free_ed_nds] = fn_nds_on_bdry(inner_bdry, nds, els)
%Returns ordered list of nodes tracking boundary

%find which bdry nds are within overall model and pick one of these as the
%first node
free_ed = fn_find_free_edges(els);
free_ed_nds = fn_find_free_edge_nds(free_ed);
bdry_start_pt = inpolygon(inner_bdry(:,1), inner_bdry(:,2), nds(free_ed_nds, 1), nds(free_ed_nds, 2));
bdry_start_pt = min(find(bdry_start_pt));
start_node = fn_find_node_at_point(nds, inner_bdry(bdry_start_pt,:), inf);

bdry_nds = zeros(size(nds, 1), 1);
[bdry_nds, exit_cond] = fn_track_bdry(nds, els, inner_bdry, bdry_nds, bdry_start_pt, 1, start_node, free_ed_nds);
if strcmp(exit_cond, 'edge_node')
    %trim the list, reverse it and start again from what was the first node
    %(now the last node) in the opposite direction
    bdry_nds = bdry_nds(1:min(find(bdry_nds == 0)) - 1);
    bdry_nds = [flipud(bdry_nds); zeros(size(nds, 1), 1)];
    [bdry_nds, exit_cond] = fn_track_bdry(nds, els, inner_bdry, bdry_nds, bdry_start_pt, -1, start_node, free_ed_nds);
    if ~strcmp(exit_cond, 'edge_node')
        error('Failed to find other end of boundary')
    end
end
% figure; plot(nds(:,1), nds(:,2), 'y.'); hold on; plot([inner_bdry(:,1); inner_bdry(1,1)], [inner_bdry(:,2); inner_bdry(1,2)], 'r:'); xlim([min(inner_bdry(:,1)) - 1e-3, max(inner_bdry(:,1)) + 1e-3]); ylim([min(inner_bdry(:,2))- 1e-3, max(inner_bdry(:,2))+ 1e-3]);
bdry_nds = bdry_nds(1:min(find(bdry_nds == 0)) - 1);



end

function [bdry_nds, exit_cond] = fn_track_bdry(nds, els, inner_bdry, bdry_nds, bdry_start_pt, bdry_dir, start_node, free_ed_nds)
bdry_pt = bdry_start_pt;
bdry_nd = start_node;
b = min(find(bdry_nds == 0)); %index of next one
%track boundary in direction stated
while 1
    bdry_nds(b) = bdry_nd;
    % plot(nds(bdry_nds(1:b),1), nds(bdry_nds(1:b), 2), 'r');
    b = b + 1;
    surr_nds = fn_find_surrounding_nds(bdry_nd, els);
    surr_els = fn_find_surrounding_els(els, bdry_nd);
    %increment bdry_pt until it goes outside polygon of surrounding nodes
    while any(fn_point_in_element(nds, els(surr_els, :), inner_bdry(bdry_pt, :)))
        last_bdry_pt = bdry_pt;
        bdry_pt = bdry_pt + bdry_dir;
        if bdry_pt > size(inner_bdry,1)
            bdry_pt = 1;
        end
        if bdry_pt < 1
            bdry_pt = size(inner_bdry,1);
        end
    end
    if b > 2
        surr_nds(surr_nds == bdry_nds(b - 2)) = [];
    end
    % plot(inner_bdry([last_bdry_pt, bdry_pt], 1), inner_bdry([last_bdry_pt, bdry_pt], 2), 'k');
    % plot(nds(surr_nds,1), nds(surr_nds, 2), 'co')
    %find nearest node in s to line between bndry pts
    [node_list, s, r] = fn_find_nodes_on_line(nds(surr_nds,:), inner_bdry(last_bdry_pt, :), inner_bdry(bdry_pt, :), inf);
    node_list = surr_nds(node_list);
    k = s>0;
    node_list = node_list(k);
    if isempty(node_list)
        exit_cond = 'edge_node';
        break
    end
    r = r(k);
    [~, i] = min(abs(r));
    bdry_nd = node_list(i);
    if bdry_nd == bdry_nds(1)
        exit_cond = 'return_to_start';
        break
    end
    if ismember(bdry_nd, free_ed_nds)
        exit_cond = 'edge_node';
        bdry_nds(b) = bdry_nd;
        break
    end
end
end
