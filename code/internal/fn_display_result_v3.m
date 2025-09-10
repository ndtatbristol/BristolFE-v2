function h_patch = fn_display_result_v3(mod, display_options)
%SUMMARY
%   Displays mesh from 2D or 3D model, returning handle to patches for later
%   animations
%USAGE
%   fn_display_result(mod, display_options) to display mesh OR
%INPUTS
%   nodes - n x {2|3} matrix of nodal coordinates. The row number is the node
%   number; columns are the x and y (and z) coordinates of the node.
%   elements - m x n matrix of element nodes. The row number is the element
%   number; columns 1, 2 and 3 are the node numbers of the nodes for each
%   element
%   display_options - structured variable allowing optional plotting properties to
%   be set. See below for defaults. In particular:
%   default_options.node_sets_to_plot - allows specific nodes to be plotted
%   in a particular color. It is a vector of structured variables with
%   fields nd and col. nd is a vector of node indices and col is the color
%   (e.g. 'r') in which nodes in that set will be plotted.
%OUTPUTS
%   none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default_options.draw_elements = 0;
default_options.element_edge_color = [1,1,1] * 0.95;
default_options.element_face_colour = [1,1,1] * 0.75;
default_options.mesh_edge_color = 'k';
default_options.draw_mesh_edges = 1;
default_options.node_sets_to_plot = [];
default_options.scale_factor = [];
default_options.show_abs = 1;
default_options.matl_cols = [];
default_options.interface_el_col = [0,0,1];
default_options.transparency = 0.5;
default_options.scale = 1;
default_options.offset = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_options = fn_set_default_fields(display_options, default_options);

% if ~isfield(mod, 'el_abs_i')
%     mod.el_abs_i = zeros(size(mod.els, 1), 1);
% end

if isempty(display_options.matl_cols)
    no_matls = numel(unique(mod.el_mat_i));
    display_options.matl_cols = linspace(0.5, 0.75, no_matls)' * [1, 1, 1];
end

mod.el_mat_i(mod.el_mat_i == 0) = 1; %needed to prevent errors for elements with no material
base_cdata = permute(display_options.matl_cols(mod.el_mat_i, :), [1, 3, 2]);

if size(mod.nds, 2) == 2
    %2D CASE
    
    if numel(display_options.offset) ~= 2
        display_options.offset = [0,0];
    end

    %Draw elements always to create the patch object. Then turn edges on or off
    %etc to get desired effects
    hold on;
    mod.els(mod.els == 0) = NaN;
    h_patch = patch('Faces', mod.els, 'Vertices', display_options.scale * mod.nds + display_options.offset, 'CData', base_cdata, 'FaceColor', 'flat');

    if ~isempty(display_options.interface_el_col)
        [i, ~] = find(mod.els == 0 | isnan(mod.els));
        x = [mod.nds(mod.els(i, 1), 1), mod.nds(mod.els(i, 2), 1)]';
        y = [mod.nds(mod.els(i, 1), 2), mod.nds(mod.els(i, 2), 2)]';
        plot(x * display_options.scale + display_options.offset(1), ...
             y * display_options.scale + display_options.offset(2), ...
             'Color', display_options.interface_el_col);
    end

    if display_options.show_abs
        set(h_patch, 'CData', base_cdata .* (1 - mod.el_abs_i / 2))
    end

    if display_options.draw_elements
        set(h_patch, 'EdgeColor', display_options.element_edge_color);
    else
        set(h_patch, 'EdgeColor', 'none');
    end

    if display_options.draw_mesh_edges
        %find edges that only occur once (i.e. they are the free edges)
        free_ed = fn_find_free_edges(mod.els);
        %plot them
        hold on;
        plot(reshape(mod.nds(free_ed, 1), size(free_ed))' * display_options.scale + display_options.offset(1) , ...
             reshape(mod.nds(free_ed, 2), size(free_ed))' * display_options.scale + display_options.offset(2) , ...
             display_options.mesh_edge_color);
    end

    if ~isempty(display_options.node_sets_to_plot)
        hold on;
        for ii = 1:length(display_options.node_sets_to_plot)
            plot(mod.nds(display_options.node_sets_to_plot(ii).nd, 1) * display_options.scale + display_options.offset(1) , ...
                 mod.nds(display_options.node_sets_to_plot(ii).nd, 2) * display_options.scale + display_options.offset(2) , ...
                 display_options.node_sets_to_plot(ii).col);
        end
    end

else
    %3D CASE
    if numel(display_options.offset) ~= 3
        display_options.offset = [0,0,0];
    end
    el_faces = fn_faces_from_els(mod.els, 1:size(mod.els,1), mod.el_typ_i, mod.el_types);
    for i = 1:numel(el_faces)
        cdata = base_cdata(el_faces{i}.el_i, :, :);
        j = fn_exterior_faces(el_faces{i}.fcs);
        fcs = el_faces{i}.fcs(j,:);
        cdata = cdata(j,:,:);
        el_i = el_faces{i}.el_i(j);
        h_patch = patch('Faces', fcs, 'Vertices', ...
            mod.nds * display_options.scale + display_options.offset, ...
            'FaceColor', 'flat', ...
            'FaceAlpha', display_options.transparency,  'CData', cdata, 'EdgeColor', 'none');

        if display_options.show_abs
            set(h_patch, 'CData', cdata .* (1 - mod.el_abs_i(el_i) / 2))
        end

        if display_options.draw_elements
            set(h_patch, 'EdgeColor', display_options.element_edge_color);
        else
            set(h_patch, 'EdgeColor', 'none');
        end
        %Mesh edges
        ext_e = fn_external_edges(fcs, el_i);
        hold on;
        plot3([mod.nds(ext_e(:, 1), 1), mod.nds(ext_e(:, 2), 1)]' * display_options.scale + display_options.offset(1), ...
              [mod.nds(ext_e(:, 1), 2), mod.nds(ext_e(:, 2), 2)]' * display_options.scale + display_options.offset(2), ...
              [mod.nds(ext_e(:, 1), 3), mod.nds(ext_e(:, 2), 3)]' * display_options.scale + display_options.offset(3), ...
              display_options.mesh_edge_color);
    end
    view(3)

    if ~isempty(display_options.node_sets_to_plot)
        for ii = 1:length(display_options.node_sets_to_plot)
            plot3(mod.nds(display_options.node_sets_to_plot(ii).nd, 1) * display_options.scale + display_options.offset(1), ...
                mod.nds(display_options.node_sets_to_plot(ii).nd, 2) * display_options.scale + display_options.offset(2), ...
                mod.nds(display_options.node_sets_to_plot(ii).nd, 3) * display_options.scale + display_options.offset(3), ...
                display_options.node_sets_to_plot(ii).col);
        end
    end

end

axis equal;
axis off;
end


%Following used in 3D case - could move these out to separate functions for
%consistency with fn_find_free_edges in 2D case?

function ext_e = fn_external_edges(fcs, el_i)
%Criterion: any face edge that is shared by >1 faces of same element
el_i2 = reshape(el_i * ones(1, size(fcs, 2)), [], 1);
tmp2 = [fcs(:,2:end), fcs(:, 1)];
all_e = [fcs(:), tmp2(:)];
all_e = sort(all_e, 2);
all_e = [el_i2, all_e];
[~, i, j] = unique(all_e, 'rows');
k = accumarray(j, 1);
ext_e = all_e(i(k > 1), :);
ext_e = ext_e(:, 2:end);
end

function i = fn_exterior_faces(fcs)
%Identify exterior faces. Criterion: any faces that are not shared by different elements
tmp = sort(fcs,2); %so each row has nodes numbered in ascending order
[~, i, j] = unique(tmp, 'rows');
k = accumarray(j, 1);
i = i(k == 1);
end
