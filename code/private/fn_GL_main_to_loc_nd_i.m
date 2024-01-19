function loc_nd_i = fn_GL_main_to_loc_nd_i(main_nd_i_to_convert, equiv_main_nd_i_for_loc)

[~, loc_nd_i] = ismember(main_nd_i_to_convert, equiv_main_nd_i_for_loc);

end

% [~, loc_nd_i] = ismember(main.doms{d}.res{tx}.frc_gl_nds, main.doms{d}.scats{s}.mod.main_nd_i);
% main.doms{d}.scats{s}.steps{tx}.load.frc_nds = loc_nd_i;
% main.doms{d}.scats{s}.steps{tx}.load.frc_dfs = main.doms{d}.res{tx}.frc_dfs;
% main.doms{d}.scats{s}.steps{tx}.load.frcs = main.doms{d}.res{tx}.frcs;
% main.doms{d}.scats{s}.steps{tx}.load.time = main.mod.time;
% 
% 
% [~, loc_nd_i] = ismember(main.doms{d}.res{tx}.dsp_gl_nds, main.doms{d}.scats{s}.mod.main_nd_i);
% main.doms{d}.scats{s}.steps{tx}.mon.nds = loc_nd_i;
% main.doms{d}.scats{s}.steps{tx}.mon.dfs = main.doms{d}.res{tx}.dsp_dfs;
% main.doms{d}.scats{s}.steps{tx}.mon.field_output_every_n_frames = options.field_output_every_n_frames;