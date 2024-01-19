function [force_inp, force_set, force_in_set, force_out_set] = fn_convert_disps_to_forces_v2(K_sub, C_sub, M_sub, time_step, disp_inp, lyrs, in_or_out)
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

% Original version (that works!) is in loop here
% tmp = disp_inp(force_set, :);
% accn = zeros(size(tmp));
% vel = zeros(size(tmp));
% accn(:, 3:end) = (tmp(:, 3:end) - 2 * tmp(:, 2:end-1) + tmp(:, 1:end-2)) / time_step ^ 2;
% vel(:, 2:end) = (tmp(:, 2:end) - tmp(:, 1:end-1)) / time_step;
% force_inp = zeros(size(accn));
% tic
% for i = 2:size(accn, 2)
%     force_inp(:, i)  = ...
%         M_sub(force_set, force_set) * accn(:, i) + ...
%         C_sub(force_set, force_set) * vel(:, i) + ...
%         K_sub(force_set, :) * disp_inp(:, i - 1);
% end
% % [vel(150,300:305) / max(vel, [], 'all');
% %     accn(150,300:305) / max(accn, [], 'all');
% %     disp_inp(150,300:305) / max(disp_inp, [], 'all');
% %     force_inp(150,300:305) / max(force_inp, [], 'all')]
% toc
% 
% dumdum = force_inp;

%Faster matrix version here - seems to be identical to original
tmp = [zeros(sum(force_set), 2), disp_inp(force_set, :)];
accn = (tmp(:, 3:end) - 2 * tmp(:, 2:end-1) + tmp(:, 1:end-2)) / time_step ^ 2;
vel = (tmp(:, 3:end) - tmp(:, 2:end - 1)) / time_step;
disp = [zeros(size(disp_inp, 1), 1), disp_inp(:, 1:end -1)];
force_inp  = ...
        M_sub(force_set, force_set) * accn + ...
        C_sub(force_set, force_set) * vel + ...
        K_sub(force_set, :) * disp;

end