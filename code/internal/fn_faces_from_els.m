function el_faces = fn_faces_from_els(els, el_typ_i, el_types)
%this returns a cell array (one cell per element type), each cell
%containing a matrix of element face nodes for elements of that type in the
%model.

un_el_typ_i = unique(el_typ_i);
el_i = (1:size(els, 1))';

for i = 1:numel(un_el_typ_i)
    ei = el_types{un_el_typ_i(i)};

    switch ei
        case {'CPE3', 'AC2D3'} %2D triangles
            fc_i = [
                1,2
                2,3
                3,1];
        case 'ASI2D2' %2D interface element
            fc_i = [
                1,2];
        case 4 %3D Tetrahedron
            fc_i = [
                1,2,3
                1,2,4
                2,3,4
                1,3,4];
        case 'C3D8R' %3D cubes
            fc_i = [
                1,2,3,4
                1,2,6,5
                2,3,7,6
                3,4,8,7
                4,1,5,8
                5,6,7,8
                ];

    end
    j = el_typ_i == un_el_typ_i(i);

    el_faces{i}.el_typ_i = ei;
    el_faces{i}.fcs = reshape(els(j, fc_i), [], size(fc_i, 2));
    el_faces{i}.el_i = reshape(el_i(j) * ones(1, size(fc_i, 1)), [], 1);
end

end
