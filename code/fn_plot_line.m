function fn_plot_line(pts, col, varargin)
if isempty(varargin)
    closed = 0;
else
    closed = varargin{1};
end
if closed
    pts = [pts; pts(1,:)];
end

hold on;
plot(pts(:, 1), pts(:, 2), col);
end