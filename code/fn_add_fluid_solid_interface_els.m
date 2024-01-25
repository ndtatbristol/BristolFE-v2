function mod = fn_add_fluid_solid_interface_els(mod, matls)
%SUMMARY
%   Adds the necessary interface elements between all solid and fluid
%   elements in a model. Without these there is no coupling between the
%   solid and fluid domains.


fluid_el_name = 'AC2D3';
solid_el_name = 'CPE3';
interface_el_name = 'ASI2D2';

[ed, els_with_ed] = fn_get_edges(mod.els);

if ~isfield(mod, 'el_typ_i')
    mod.el_typ_i = {matls(mod.el_mat_i).el_typ};
    mod.el_typ_i = mod.el_typ_i(:);
end

invalid = ~els_with_ed;
els_with_ed(invalid) = 1;
el_typ = mod.el_typ_i(els_with_ed);
el_typ(invalid) = {''};

j = (strcmp(el_typ(:,1), fluid_el_name) & strcmp(el_typ(:,2), solid_el_name)) | ...
    (strcmp(el_typ(:,2), fluid_el_name) & strcmp(el_typ(:,1), solid_el_name));

el_typ = el_typ(j, :);
els_with_ed = els_with_ed(j, :);
ed = ed(j, :);

%New_els needs ordering so solid and fluid are on correct sides for all
%elements - this is why mod.nds data is necessary
no_int_els = size(ed,1);
for i = 1:no_int_els
    %work out centre of fluid element (typ == 2) adjoiing this edge
    % e = els_with_ed(i, typ(i, :) == 2);
    e = els_with_ed(i, strcmp(el_typ(i, :), fluid_el_name));
    ec = fn_calc_element_centres(mod.nds, mod.els(e,:));
    %line between nodes
    a = mod.nds(ed(i, 2), :) - mod.nds(ed(i, 1), :);
    %line at right angle to line between nodes
    b = [a(2), -a(1)];
    %line from first node to ec
    c = ec - mod.nds(ed(i, 1), :);
    %check sign of dot product
    if dot(c, b) < 0
        ed(i, :) = fliplr(ed(i, : ));
    end
end

%Add the new interface elements to the model
mod.els = [mod.els; [ed, zeros(no_int_els, size(mod.els, 2)- size(ed,2))]];
mod.el_typ_i = [mod.el_typ_i; repmat({interface_el_name}, [no_int_els, 1])];
mod.el_mat_i = [mod.el_mat_i; zeros(no_int_els, 1)];

if isfield(mod, 'el_abs_i')
    mod.el_abs_i = [mod.el_abs_i; zeros(no_int_els, 1)];
end



end