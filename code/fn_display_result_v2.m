function h_patch = fn_display_result_v2(nds, els, display_options)
%SUMMARY
%   Displays mesh from 2D or 3D model, returning handle to patches for later
%   animations
%USAGE
%   fn_display_result(nodes, elements, display_options) to display mesh OR
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
default_options.el_mat_i = ones(size(els, 1), 1);
default_options.el_abs_i = zeros(size(els, 1), 1);
default_options.show_abs = 1;
default_options.matl_cols = [];
default_options.interface_el_col = [0,0,1];
default_options.transparency = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display_options = fn_set_default_fields(display_options, default_options);

if isempty(display_options.matl_cols)
    no_matls = numel(unique(display_options.el_mat_i));
    display_options.matl_cols = linspace(0.5, 0.75, no_matls)' * [1, 1, 1];
end

display_options.el_mat_i(display_options.el_mat_i == 0) = 1; %needed to prevent errors for elements with no material
base_cdata = permute(display_options.matl_cols(display_options.el_mat_i, :), [1, 3, 2]);

if size(nds, 2) == 2
    %2D CASE

    %Draw elements always to create the patch object. Then turn edges on or off
    %etc to get desired effects
    hold on;
    els(els == 0) = NaN;
    h_patch = patch('Faces', els, 'Vertices', nds, 'CData', base_cdata, 'FaceColor', 'flat');

    if ~isempty(display_options.interface_el_col)
        [i, ~] = find(els == 0 | isnan(els));
        x = [nds(els(i, 1), 1), nds(els(i, 2), 1)]';
        y = [nds(els(i, 1), 2), nds(els(i, 2), 2)]';
        plot(x, y, 'Color', display_options.interface_el_col);
    end

    if display_options.show_abs
        set(h_patch, 'CData', base_cdata .* (1 - display_options.el_abs_i / 2))
    end

    if display_options.draw_elements
        set(h_patch, 'EdgeColor', display_options.element_edge_color);
    else
        set(h_patch, 'EdgeColor', 'none');
    end

    if display_options.draw_mesh_edges
        %find edges that only occur once (i.e. they are the free edges)
        free_ed = fn_find_free_edges(els);
        %plot them
        hold on;
        plot(reshape(nds(free_ed, 1), size(free_ed))', reshape(nds(free_ed, 2), size(free_ed))', display_options.mesh_edge_color);
    end

    if ~isempty(display_options.node_sets_to_plot)
        hold on;
        for ii = 1:length(display_options.node_sets_to_plot)
            plot(nds(display_options.node_sets_to_plot(ii).nd, 1), nds(display_options.node_sets_to_plot(ii).nd, 2), display_options.node_sets_to_plot(ii).col);
        end
    end

else
    %3D CASE
    switch size(els, 2)
        case 4 %tetrahedron
            %TODO
        case 8 % cube
            fc_i = [
                1,2,3,4
                1,2,6,5
                2,3,7,6
                3,4,8,7
                4,1,5,8
                5,6,7,8
                ];

    end
    base_cdata = reshape(repmat(base_cdata, [size(fc_i,1), 1]), [], 1, 3);
    fcs = reshape(els(:, fc_i), [], size(fc_i, 2));
    el_i = reshape((1:size(els, 1))' * ones(1, size(fc_i, 1)), [], 1);
    i = fn_exterior_faces(fcs);
    fcs = fcs(i,:);
    base_cdata = base_cdata(i,:,:);
    el_i = el_i(i);


    h_patch = patch('Faces', fcs, 'Vertices', nds, 'FaceColor', 'flat', ...
        'FaceAlpha', display_options.transparency,  'CData', base_cdata, 'EdgeColor', 'none');

    if display_options.show_abs
        set(h_patch, 'CData', base_cdata .* (1 - display_options.el_abs_i(el_i) / 2))
    end

    if display_options.draw_elements
        set(h_patch, 'EdgeColor', display_options.element_edge_color);
    else
        set(h_patch, 'EdgeColor', 'none');
    end

    %Mesh edges
    ext_e = fn_external_edges(fcs, el_i);
    hold on;
    plot3([nds(ext_e(:, 1), 1), nds(ext_e(:, 2), 1)]' , ...
        [nds(ext_e(:, 1), 2), nds(ext_e(:, 2), 2)]' , ...
        [nds(ext_e(:, 1), 3), nds(ext_e(:, 2), 3)]', display_options.mesh_edge_color);
    view(3)

    if ~isempty(display_options.node_sets_to_plot)
        for ii = 1:length(display_options.node_sets_to_plot)
            plot3(nds(display_options.node_sets_to_plot(ii).nd, 1), ...
                nds(display_options.node_sets_to_plot(ii).nd, 2), ...
                nds(display_options.node_sets_to_plot(ii).nd, 3), ...
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
