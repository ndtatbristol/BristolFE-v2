function [force_inp, force_set, force_in_set, force_out_set] = fn_convert_disps_to_forces_v2(K_sub, C_sub, M_sub, time_step, disp_inp, lyrs, in_or_out, solver_mode)
force_in_set = lyrs == 2;
force_out_set = lyrs == 3;
force_set = force_in_set | force_out_set;
switch in_or_out
    case 'in'
        freeze_set = [lyrs == 3 | lyrs == 4];
    case 'out'
        freeze_set = [lyrs == 1 | lyrs == 2];
end

disp_inp(freeze_set, :) = 0;
M_sub = spdiags(sum(M_sub).', 0, size(M_sub,1), size(M_sub,2));

tmp = [zeros(size(disp_inp, 1), 2), disp_inp];
accn = (tmp(:, 3:end) - 2 * tmp(:, 2:end-1) + tmp(:, 1:end-2)) / time_step ^ 2;
switch solver_mode
    case 'vel at last half time step'
        vel = (tmp(:, 2:end - 1) - tmp(:, 1:end - 2)) / time_step;
    case 'vel at curent time step'
        vel = (tmp(:, 3:end) - tmp(:, 1:end - 2)) / (2 * time_step);
end
disp = [zeros(size(disp_inp, 1), 1), disp_inp(:, 1:end - 1)];
% disp = [zeros(size(disp_inp, 1), 2), disp_inp(:, 1:end - 2)];
% disp = disp_inp;
force_inp  = ...
        M_sub(force_set, :) * accn + ...
        C_sub(force_set, :) * vel + ...
        K_sub(force_set, :) * disp;


end