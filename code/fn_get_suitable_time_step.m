function time_step = fn_get_suitable_time_step(matls, el_size, varargin)
%SUMMARY
%   Returns what should be a stable time step by calculating the fastest
%   possible wavespeed in the materials, workout out how fast such a wave
%   traverses the specified element size and dividing that by a safety
%   factor (default = 3, specify a alternative value as 3rd optional
%   argument if desired

if isempty(varargin)
    safety_factor = 3;
else
    safety_factor = varargin{1};
end
mx_vel = fn_estimate_max_min_vels(matls);

time_step = el_size / mx_vel / safety_factor;
end
