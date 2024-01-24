function time_step = fn_get_suitable_time_step(matls, el_size, varargin)
if isempty(varargin)
    safety_factor = 3;
else
    safety_factor = varargin{1};
end
mx_vel = fn_estimate_max_min_vels(matls);

time_step = el_size / mx_vel / safety_factor;
end
