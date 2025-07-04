function pts = fn_create_smooth_random_blob(min_rad_frac, complexity, no_pts)
%Use complexity = 0 to get a circle
a = linspace(0, 2 * pi, no_pts + 1)';
a = a(1:end - 1);
if complexity > 0
    n = [0: complexity];
    phi = rand(size(n)) * 2 * pi;
    r = sum(sin(a * n + phi), 2);
    r = (r - min(r)) / (max(r) - min(r)) * (1 - min_rad_frac) + min_rad_frac;
else
    r = 1;
end
pts = r .* [cos(a), sin(a)];
end
