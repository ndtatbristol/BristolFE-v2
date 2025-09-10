function mod = fn_remove_fluid_solid_interface_els(mod, el_types, varargin)
%SUMMARY
%   Removes interface elements from model. Typically used when creating
%   subdomains from another model and introducing scatterers.
if numel(varargin) < 1
    options = [];
else
    options = varargin{1};
end
default_options.interface_el_name = 'ASI2D2';
options = fn_set_default_fields(options, default_options);

interface_el_i = find(strcmp(el_types, options.interface_el_name));
if isempty(interface_el_i)
    %No interface elements defined, so none to remove!
    return
end

els_in_use = ~(mod.el_typ_i == interface_el_i);
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);
end
