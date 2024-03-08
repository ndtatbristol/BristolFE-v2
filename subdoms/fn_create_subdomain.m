function dm_mod = fn_create_subdomain(mn_mod, matls, inner_bdry, abs_layer_thick, varargin)
%New version - based on tidier way of getting layers
%Needs to return something like this
    %             nds: [11424×2 double]
    %             els: [22516×3 double]
    %        el_mat_i: [22516×1 double]
    %        el_abs_i: [22516×1 double]
    %        el_typ_i: {22516×1 cell}
    %       bdry_lyrs: [11424×1 double]
    %       main_nd_i: [11424×1 double]
    % outer_bndry_pts: [882×2 double]
    % inner_bndry_pts: [361×2 double]

%Last optional argument allows the start of the absorbing region to be
%specified. If not specified it will be penultimate boundary layer.

interface_el_name = 'ASI2D2';

if ~isfield(mn_mod, 'el_typ_i')
    mn_mod.el_typ_i = {matls(mn_mod.el_mat_i).el_typ};
    mn_mod.el_typ_i = mn_mod.el_typ_i(:);
end

dm_mod.nds = mn_mod.nds;
dm_mod.els = mn_mod.els;
dm_mod.el_mat_i = mn_mod.el_mat_i;
if isfield(mn_mod, 'el_abs_i')
    dm_mod.el_abs_i = mn_mod.el_abs_i;
else
    dm_mod.el_abs_i = zeros(size(mn_mod.el_mat_i));
end
dm_mod.el_typ_i = mn_mod.el_typ_i;
dm_mod.bdry_lyrs = zeros(size(mn_mod.nds, 1), 1);
dm_mod.el_abs_i = zeros(size(mn_mod.els, 1), 1);
dm_mod.inner_bndry_pts = inner_bdry;

%       main_nd_i: [11424×1 double]
% outer_bndry_pts: [882×2 double]

%Get elements in region

el_used = fn_elements_in_region(dm_mod, inner_bdry);
dm_mod.main_int_el_i = find(el_used); %the indices of the internal els in the main model for this sub-domain (need to be identifiable when doing validation models)

%Work out and assign bdry nodes to layers
for i = 1:4
    [el_i, bdry_nds] = fn_find_adjacent_els_to_els(dm_mod.els, el_used, ~el_used);
    dm_mod.bdry_lyrs(bdry_nds)= i;
    el_used(el_i) = 1;

    if i == 3 %for 3rd one, use boundary as start of absorbing region
        if isempty(varargin)
            abs_layer_start_bdry = dm_mod.nds(bdry_nds, :);
        else
            abs_layer_start_bdry = varargin{1};
        end
    end
end

% figure;c = 'rgbm'; h = fn_show_geometry(dm_mod, matls, []);for i = 1:4; fn_plot_line(dm_mod.nds(dm_mod.bdry_lyrs == i,:), [c(i), '.']); end


%Delete original interface elements and then regenerate later in this 
%function - this is to avoid potential instability at edge of domain where 
%original interface elements copied from main mesh may now be on free edges
els_in_use = ~strcmp(dm_mod.el_typ_i, interface_el_name);
[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i, el_used] = fn_remove_unused_elements(els_in_use, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i, el_used);


%Add the absorbing layers by working out from centre of region
% [~, cand_els] = fn_elements_in_region2(dm_mod.nds, dm_mod.els, abs_layer_start_bdry);
cand_els = ~el_used;
dm_mod.el_abs_i(cand_els) = fn_dist_point_to_bdry_2D(fn_calc_element_centres(dm_mod.nds, dm_mod.els(cand_els, :)), abs_layer_start_bdry) / abs_layer_thick;
els_in_use = ones(size(dm_mod.els, 1), 1);
els_in_use(dm_mod.el_abs_i > 1) = 0;

[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(els_in_use, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);
[dm_mod.nds, dm_mod.els, old_nds, ~] = fn_remove_unused_nodes(dm_mod.nds, dm_mod.els);
dm_mod.main_nd_i = old_nds;
dm_mod.bdry_lyrs = dm_mod.bdry_lyrs(old_nds);

dm_mod = fn_add_fluid_solid_interface_els(dm_mod, matls);

free_ed = fn_find_free_edges(dm_mod.els);

dm_mod.outer_bndry_pts = [dm_mod.nds(free_ed, 1), dm_mod.nds(free_ed, 2)];
dm_mod.int_el_i = fn_elements_in_region(dm_mod, dm_mod.inner_bndry_pts);

end

%--------------------------------------------------------------------------

function [adj_els, common_nds] = fn_find_adjacent_els_to_els(els, els_to_consider, els_to_choose_from)
%returns logical array [size(els,1)x1] of elements in
%els(els_to_choose_from,:) that share common nodes with
%els(els_to_consider,:)

%unique nodes in els
tmp = els(els_to_consider,:);
un_nds = unique(tmp(:));
un_nds(un_nds == 0) = [];

%find rows in els_to_choose_from that contain un_nds
% nds_to_choose_from = els(els_to_choose_from, :);
tmp = ismember(els, un_nds);
tmp(~els_to_choose_from, :) = 0;

common_nds = els(tmp);
adj_els = sum(tmp, 2) > 0;

end
