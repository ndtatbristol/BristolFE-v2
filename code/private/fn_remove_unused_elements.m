function [old_els, new_els, varargout] = fn_remove_unused_elements(in_use, varargin)
in_use = find(in_use);
old_els = [1:size(varargin{1}, 1)]';
old_els = old_els(in_use);
new_els = zeros(size(varargin{1}, 1), 1);
new_els(old_els) = 1:numel(old_els);
for i = 1:length(varargin)
    varargout{i} = varargin{i}(in_use, :);
end
end