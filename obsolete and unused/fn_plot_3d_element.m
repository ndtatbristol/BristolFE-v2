function fcs = fn_plot_3d_element(mod, cols, alpha)
%Assumes Abaqus ordering of nodes!
switch size(mod.els, 2)
    case 4 %tetrahedron
    case 8 % cube
        cols = reshape(repmat(cols, [6,1]), [], 1, 3);
        fc_i = [
            1,2,3,4
            1,2,6,5
            2,3,7,6
            3,4,8,7
            4,1,5,8
            5,6,7,8
            ];
        fcs = reshape(mod.els(:, fc_i), [], size(fc_i, 2));
        el_i = reshape((1:size(mod.els, 1))' * ones(1, size(fc_i, 1)), [], 1);

end
i = fn_exterior_faces(fcs);
fcs = fcs(i,:); 
cols = cols(i,:,:);
el_i = el_i(i);
patch('Faces', fcs, 'Vertices', mod.nds, 'FaceColor', 'flat', 'facealpha', alpha,  'CData', cols, 'EdgeColor', 'None');

%Mark external edges. Criterion: any face edge that is shared by >1 faces of same element
ext_e = fn_external_edges(fcs, el_i);
hold on;
plot3([mod.nds(ext_e(:, 1), 1), mod.nds(ext_e(:, 2), 1)]' , ...
    [mod.nds(ext_e(:, 1), 2), mod.nds(ext_e(:, 2), 2)]' , ...
    [mod.nds(ext_e(:, 1), 3), mod.nds(ext_e(:, 2), 3)]', 'g');

view(3);


end

function ext_e = fn_external_edges(fcs, el_i)
%Criterion: any face edge that is shared by >1 faces of same element
el_i2 = reshape(el_i * ones(1, size(fcs, 2)), [], 1);
tmp2 = [fcs(:,2:end), fcs(:, 1)];
all_e = [fcs(:), tmp2(:)];
all_e = sort(all_e, 2);
all_e = [el_i2, all_e];
[tmp2, i, j] = unique(all_e, 'rows');
k = accumarray(j, 1);
ext_e = all_e(i(k > 1), :);
ext_e = ext_e(:, 2:end);
end

function i = fn_exterior_faces(fcs)
%Identify exterior faces. Criterion: any faces that are not shared by different elements
tmp = sort(fcs,2); %so each row has nodes numbered in ascending order
[tmp2, i, j] = unique(tmp, 'rows');
k = accumarray(j, 1);
i = i(k == 1);
end