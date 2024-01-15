function [loc_nd, loc_df, varargout] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use, varargin)

j = ismember(loc_df, dofs_to_use);

% if el_dim_last
% 
%     for i = 1:length(varargin)
%         varargout{i} = varargin{i}(j, j, :);
%     end
% else

    for i = 1:length(varargin)
        varargout{i} = varargin{i}(:, j, j);
    end
    
% end
    loc_nd = loc_nd(j);
    loc_df = loc_df(j);

end