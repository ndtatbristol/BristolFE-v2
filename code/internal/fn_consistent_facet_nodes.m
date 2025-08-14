function [facet_nds, eds, ed_fcs] = fn_consistent_facet_nodes(facet_nds)

%gets matrix of facet nodes into consistent ordering (say clockwise) around
%each facet, so surface normals of adjacent facets are in same direction

if size(facet_nds, 2) ~= 3
    error('Facet nodes must have 3 columns (triangular facets')
end

eds = [facet_nds(:,1), facet_nds(:,2)
    facet_nds(:,2), facet_nds(:,3)
    facet_nds(:,3), facet_nds(:,1)];
ed_fcs = repmat((1:size(facet_nds, 1))', [3, 1]);

[~,ia,ic] = unique(sort(eds, 2), 'rows');
occs = accumarray(ic,1);
if any(occs) > 2
    error('More than two facets at same edge')
end
occs = find(occs > 1);
for i = 1:numel(occs)
    j = find(ic == occs(i));
    if all(eds(j(1), :) == eds(j(2), :))
        %Flip order of nodes in associated face
        flip_face = ed_fcs(j(1));
        facet_nds(flip_face,:) = fliplr(facet_nds(flip_face,:));
        eds(j(1), :) = fliplr(eds(j(1), :));
    end
end
end