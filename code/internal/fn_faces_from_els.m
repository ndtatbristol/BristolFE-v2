function el_faces = fn_faces_from_els(els, el_i, el_typ_i, el_types)
%this returns a cell array (one cell per element type), each cell
%containing a matrix of element face nodes for elements of that type in the
%model. Structure will look like:
%   el_faces{q}.el_typ = name of element given by el_types{q}
%   el_faces{q}.fcs = matrix of nodes of faces of each element of el_types{q} with 
%       size [n_els_of_typ_q * n_fcs_per_el_of_type_q, n_dim - 1]
%   el_faces{q}.el_i = vector of element numbers associated with each row
%       of el_faces{q}.fcs of size [n_els_of_typ_q * n_fcs_per_el_of_type_q, 1]

un_el_typ_i = unique(el_typ_i);
el_i = el_i(:);

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
    el_faces{i}.fcs = reshape(els(j, fc_i')', size(fc_i, 2), [])';
    el_faces{i}.el_i = reshape(reshape(el_i(j), 1, []) .* ones(size(fc_i, 1), 1), [], 1);

end

end
