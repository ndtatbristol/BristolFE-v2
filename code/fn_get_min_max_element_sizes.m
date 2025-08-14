function [min_el_size, max_el_size] = fn_get_min_max_element_sizes(mod)

%Finds min and max edge length of elements in model

%This will only work for linear elements where distance between adjacent
%nodes is edge length.

n = [1: size(mod.els, 2), 1];

all_eds = unique([reshape(mod.els(:,1:end-1), [], 1), reshape(mod.els(:,2:end), [], 1)], 'rows');

%will need to remove any with zeros in (e.g. for 2-noded interface elements
%in 3-noded tri mesh

p1 = mod.nds(all_eds(:, 1), :);
p2 = mod.nds(all_eds(:, 2), :);
d_sq = sum((p1 - p2) .^ 2, 2);
min_el_size = sqrt(min(d_sq));
max_el_size = sqrt(max(d_sq));
end