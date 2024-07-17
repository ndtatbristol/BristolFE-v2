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
default_options.rayleigh_damping_level = 0;
options = fn_set_default_fields(options, default_options);

%Add el_type_i field if it's not there already
if ~isfield(mod, 'el_typ_i')
    mod.el_typ_i = {matls(mod.el_mat_i).el_typ};
    mod.el_typ_i = mod.el_typ_i(:);
end

%Get list of all (internal and external) unique element edges in model and 
%the elements adjoining each edge (zero in second column means that there 
%is only element on one side of that edge, so it is external edge)
[ed, els_with_ed] = fn_get_edges(mod.els);

%Build equivalent cell array to els_with_ed which contains the element types
invalid = ~els_with_ed;
els_with_ed(invalid) = 1;
el_typ = mod.el_typ_i(els_with_ed);
el_typ(invalid) = {''};

%Look for where the elements on either side of an edge are fluid and solid
%or vice versa (i.e. where the interface elements need to go)
j = (ismember(el_typ(:,1), options.fluid_el_names) & ismember(el_typ(:,2), options.solid_el_names)) | ...
    (ismember(el_typ(:,2), options.fluid_el_names) & ismember(el_typ(:,1), options.solid_el_names));

%Restrict the list of edges to just these ones
el_typ = el_typ(j, :);
els_with_ed = els_with_ed(j, :);
ed = ed(j, :);

%New_els needs ordering so solid and fluid are on correct sides for all
%elements - this is why mod.nds data is necessary
no_int_els = size(ed,1);
%Loop through each interface edge in turn and flip node order if necessary
%to they are all same way around
for i = 1:no_int_els
    %work out centre of fluid element adjoining this edge
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
mod.els = [mod.els; [ed, zeros(no_int_els, size(mod.els, 2) - size(ed, 2))]];
mod.el_typ_i = [mod.el_typ_i; repmat({options.interface_el_name}, [no_int_els, 1])];
mod.el_mat_i = [mod.el_mat_i; zeros(no_int_els, 1)]; %interface elements have no material

if isfield(mod, 'el_abs_i')
    mod.el_abs_i = [mod.el_abs_i; zeros(no_int_els, 1)];
end

end