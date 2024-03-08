function [dof_in_use, max_dof_per_el] = fn_find_dof_in_use_and_max_dof_per_el(unique_typs, varargin)
if numel(varargin) < 1
    dof_to_use = [];
else
    dof_to_use = varargin{1};
end

%find unique types
% un_typs = unique(el_typ_i);

%work out max_dof needed in global matrix and max dof per element
dof_in_use = [];
max_dof_per_el = 0;
for t = 1:numel(unique_typs)
    fn_el_mats = str2func(['fn_el_', unique_typs{t}]);
    [~, ~, ~, ~, loc_df] = fn_el_mats([], [], [], [], dof_to_use);
    dof_in_use = [dof_in_use, loc_df];
    max_dof_per_el = max(max_dof_per_el, numel(loc_df));
end
dof_in_use = unique(dof_in_use);

end