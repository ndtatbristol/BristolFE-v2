function mod = fn_add_fluid_solid_interface_els(mod, matls, varargin)
%SUMMARY
%   Adds the necessary interface elements between all solid and fluid
%   elements in a model. Without these there is no coupling between the
%   solid and fluid domains.
if numel(varargin) < 1
    options = [];
else
    options = varargin{1};
end
default_options.interface_el_name = 'ASI2D2';
default_options.fluid_el_names = {'AC2D3'};
default_options.solid_el_names = {'CPE3'};
options = fn_set_default_fields(options, default_options);

[ed, els_with_ed] = fn_get_edges(mod.els);

if ~isfield(mod, 'el_typ_i')
    mod.el_typ_i = {matls(mod.el_mat_i).el_typ};
    mod.el_typ_i = mod.el_typ_i(:);
end

invalid = ~els_with_ed;
els_with_ed(invalid) = 1;
el_typ = mod.el_typ_i(els_with_ed);
el_typ(invalid) = {''};

j = (ismember(el_typ(:,1), options.fluid_el_names) & ismember(el_typ(:,2), options.solid_el_names)) | ...
    (ismember(el_typ(:,2), options.fluid_el_names) & ismember(el_typ(:,1), options.solid_el_names));

el_typ = el_typ(j, :);
els_with_ed = els_with_ed(j, :);
ed = ed(j, :);

%New_els needs ordering so solid and fluid are on correct sides for all
%elements - this is why mod.nds data is necessary
no_int_els = size(ed,1);
for i = 1:no_int_els
    %work out centre of fluid element (typ == 2) adjoiing this edge
    % e = els_with_ed(i, typ(i, :) == 2);
    e = els_with_ed(i, ismember(el_typ(i, :), options.fluid_el_names));
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
mod.el_typ_i = [mod.el_typ_i; repmat({options.interface_el_name}, [no_int_els, 1])];
mod.el_mat_i = [mod.el_mat_i; zeros(no_int_els, 1)];

if isfield(mod, 'el_abs_i')
    mod.el_abs_i = [mod.el_abs_i; zeros(no_int_els, 1)];
end

end