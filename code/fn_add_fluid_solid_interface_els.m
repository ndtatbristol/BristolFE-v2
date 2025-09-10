function [mod, el_types] = fn_add_fluid_solid_interface_els(mod, el_types, varargin)
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

%Get lists of indices of fluid and solid element types
solid_el_i = find(ismember(el_types, options.solid_el_names));
fluid_el_i = find(ismember(el_types, options.fluid_el_names));

if isempty(solid_el_i) || isempty(fluid_el_i)
    %model has no solid or no fluid element types 
    return
end

%Add interface element to list of element types if not already there
if ~any(strcmp(el_types, options.interface_el_name))
    el_types{end + 1} = options.interface_el_name;
end

%Get interface element index
interface_el_i = find(strcmp(el_types, options.interface_el_name));

%New method using find interface function (should work for 2D and 3D up
%to and including this function)
[interface_facets, interface_el_solid, interface_el_fluid] = fn_find_interface(mod, ismember(mod.el_typ_i, solid_el_i), ismember(mod.el_typ_i, fluid_el_i), el_types);

if isempty(interface_facets)
    %No interface found, so do nothing
    return
end

%Nodes in interface_facets need ordering so solid and fluid are 
%on correct sides for all elements - this is why mod.nds data is necessary
%for this function to work. Currently this part is specific to 2D.
no_int_els = size(interface_facets,1);
%Loop through each interface edge in turn and flip node order if necessary
%to they are all same way around
for i = 1:no_int_els
    %work out centre of fluid element adjoining this edge
    % e = els_adjoining_fluid_solid_interface_edges(i, fluid_or_solid(i, :) == 2);
    e = interface_el_fluid(i);
    ec = fn_calc_element_centres(mod.nds, mod.els(e,:));
    %line between nodes
    a = mod.nds(interface_facets(i, 2), :) - mod.nds(interface_facets(i, 1), :);
    %line at right angle to line between nodes
    b = [a(2), -a(1)];
    %line from first node to ec
    c = ec - mod.nds(interface_facets(i, 1), :);
    %check sign of dot product
    if dot(c, b) < 0
        interface_facets(i, :) = fliplr(interface_facets(i, : ));
    end
end

%Add the new interface elements to the model
mod.els = [mod.els; [interface_facets, zeros(no_int_els, size(mod.els, 2) - size(interface_facets, 2))]];
mod.el_typ_i = [mod.el_typ_i; repmat(interface_el_i, [no_int_els, 1])];

%Extend material and absorbing indices to include new elements
mod.el_mat_i = [mod.el_mat_i; zeros(no_int_els, 1)];
mod.el_abs_i = [mod.el_abs_i; zeros(no_int_els, 1)];

end