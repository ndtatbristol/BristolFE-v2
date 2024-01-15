function isit = fn_point_in_element(nds, els, pts)
%els or pt can be multiple but not both
if size(pts, 1) >= 1 && size(els,1) == 1
    isit = inpolygon(pts(:,1), pts(:,2), nds(els, 1), nds(els, 2));
end

if size(pts, 1) == 1 && size(els,1) >= 1
    isit = zeros(size(els, 1), 1);
    for i = 1:size(els, 1)
        j = els(i, :) > 0;
        isit(i) = inpolygon(pts(:,1), pts(:,2), nds(els(i, j), 1), nds(els(i, j), 2));
    end
end

if size(pts, 1) > 1 && size(els,1) > 1
    error('Either pts or els can be more than one but not both')
end

end