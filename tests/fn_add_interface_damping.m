function C = fn_add_interface_damping(mod, M, C, gl_lookup, level)
interface_el_name = 'ASI2D2';
pressure_dof = 4;
%Adds Rayleigh damping to interface nodes in model

%Find the interface elements and associated nodes
nds = mod.els(strcmp(mod.el_typ_i, interface_el_name),:);
nds = nds(nds > 0);

% keyboard
%Find the relevant entries in global matrix
gi = gl_lookup(nds, pressure_dof);

C(gi,gi) = C(gi,gi) + M(gi,gi) * level;

end