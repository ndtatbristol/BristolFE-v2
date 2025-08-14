function mod = fn_add_fluid_solid_interface_els(mod, varargin)
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

%Add interface element to list of element types if not already there
if ~any(strcmp(mod.el_types, options.interface_el_name))
    mod.el_types{end + 1} = options.interface_el_name;
end
%Get interface element index
interface_el_i = find(strcmp(mod.el_types, options.interface_el_name));
%Get lists of indices of fluid and solid element types
solid_el_i = find(ismember(mod.el_types, options.solid_el_names));
fluid_el_i = find(ismember(mod.el_types, options.fluid_el_names));

%First get list of all (internal and external) unique element edges in model and 
%the elements adjoining each edge (zero in second column means that there 
%is only element on one side of that edge, so it is external edge)
[fluid_solid_interface_edges, els_adjoining_fluid_solid_interface_edges] = fn_get_edges(mod.els);

%Eliminate rows associated with edge elements
j = ~any(els_adjoining_fluid_solid_interface_edges == 0, 2);
fluid_solid_interface_edges = fluid_solid_interface_edges(j, :);
els_adjoining_fluid_solid_interface_edges = els_adjoining_fluid_solid_interface_edges(j, :);

%Create matrix describing whether elements on either side of each edge are
%fluid (2) or solid (1)
el_typ = mod.el_typ_i(els_adjoining_fluid_solid_interface_edges);
fluid_or_solid = zeros(size(el_typ));
fluid_or_solid(ismember(el_typ, solid_el_i)) = 1; %solid = 1
fluid_or_solid(ismember(el_typ, fluid_el_i)) = 2; %fluid = 2

%Identify edges where material is fluid on one side an solid on the other
j = fluid_or_solid(:, 1) ~= fluid_or_solid(:, 2);

%Restrict the list of edges to just these ones
% el_typ = el_typ(j, :);
els_adjoining_fluid_solid_interface_edges = els_adjoining_fluid_solid_interface_edges(j, :);
fluid_solid_interface_edges = fluid_solid_interface_edges(j, :);
fluid_or_solid = fluid_or_solid(j, :);

%Nodes in fluid_solid_interface_edges need ordering so solid and fluid are 
%on correct sides for all elements - this is why mod.nds data is necessary
%for this function to work
no_int_els = size(fluid_solid_interface_edges,1);
%Loop through each interface edge in turn and flip node order if necessary
%to they are all same way around
for i = 1:no_int_els
    %work out centre of fluid element adjoining this edge
    e = els_adjoining_fluid_solid_interface_edges(i, fluid_or_solid(i, :) == 2);
    ec = fn_calc_element_centres(mod.nds, mod.els(e,:));
    %line between nodes
    a = mod.nds(fluid_solid_interface_edges(i, 2), :) - mod.nds(fluid_solid_interface_edges(i, 1), :);
    %line at right angle to line between nodes
    b = [a(2), -a(1)];
    %line from first node to ec
    c = ec - mod.nds(fluid_solid_interface_edges(i, 1), :);
    %check sign of dot product
    if dot(c, b) < 0
        fluid_solid_interface_edges(i, :) = fliplr(fluid_solid_interface_edges(i, : ));
    end
end

%Add the new interface elements to the model
mod.els = [mod.els; [fluid_solid_interface_edges, zeros(no_int_els, size(mod.els, 2) - size(fluid_solid_interface_edges, 2))]];
mod.el_typ_i = [mod.el_typ_i; repmat(interface_el_i, [no_int_els, 1])];
mod.el_mat_i = [mod.el_mat_i; zeros(no_int_els, 1)]; %interface elements have no material

%Extend absorbing indices
mod.el_abs_i = [mod.el_abs_i; zeros(no_int_els, 1)];

end