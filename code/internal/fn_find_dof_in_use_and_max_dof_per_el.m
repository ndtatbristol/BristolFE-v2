function [dof_in_use, max_dof_per_el] = fn_find_dof_in_use_and_max_dof_per_el(unique_typs, varargin)
%This will aways return dof_in_use = dof_to_use if latter is specified 
%(even if some of those are not actually used in model - this is intended
%behaviour for cases with subdomains where extra DoF may be required in
%subdomain models built from same main one); if dof_to_use not specified
%than it will return all the unique DoF associated with elements in the model

if numel(varargin) < 1
    dof_to_use = [];
else
    dof_to_use = varargin{1};
end

%work out max_dof needed in global matrix and max dof per element
dof_in_use = dof_to_use;
max_dof_per_el = 0;
for t = 1:numel(unique_typs)
    fn_el_mats = str2func(['fn_el_', unique_typs{t}]);
    [~, ~, ~, ~, loc_df] = fn_el_mats([], [], [], [], dof_to_use);
    dof_in_use = [dof_in_use, loc_df];
    max_dof_per_el = max(max_dof_per_el, numel(loc_df));
end
dof_in_use = unique(dof_in_use);

end