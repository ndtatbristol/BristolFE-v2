function main = fn_run_GL_whole_model(main, time_pts, options)
default_options.centre_freq = main.mod.design_centre_freq;
default_options.number_of_cycles = 5;
default_options.time_step = main.mod.max_safe_time_step;
% default_options.tx_to_run = 1:numel(main.steps);
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [1:4];
default_options.movie_mode = 0; %when this is one, the input signal is set to be the toneburst rather than an impulse, so movies look nice. Other results will be wrong!!!
options = fn_set_default_fields(options, default_options);

%Runs everything necessary
time_step = options.time_step;
if isempty(time_pts)
    %choose  number of points for slowest wave to traverse largest
    %dimension of model
    time_pts = round(sqrt(sum((max(main.mod.nds) - min(main.mod.nds)) .^ 2)) / main.mod.design_min_vel / time_step);
end
main.inp.time = [0:time_pts - 1] * time_step;

%Force re-run if time axis has changed
% if ~isfield(main.mod, 'time') || time_pts ~= numel(main.mod.time)  || time_step ~= (main.mod.time(2) - main.mod.time(1))
%     fprintf('All models will be run');
%     if isfield(main.mod, 'time') && (time_pts ~= numel(main.mod.time)  || time_step ~= (main.mod.time(2) - main.mod.time(1)))
%         fprintf(' (time axis has changed)');
%     end
%     fprintf('\n');
%     main = fn_clear_fields_in_GL_model(main, 'res', 1, 1);
%     % main = fn_clear_fields_in_GL_model(main, 'mats', 1, 1);
%     %Create new time axis
%     if time_step > main.mod.max_safe_time_step
%         warning('Specified time step is greater than maximum recommended time step')
%     end
%     main.mod.time = [0:time_pts - 1] * time_step;
% end

%The desired output only affects how final A-scans are filtered; it does
%not change the FE runs, which are performed with impulse excitation.
main.inp.sig = fn_gaussian_pulse(main.inp.time, options.centre_freq, options.number_of_cycles);

%Main model is always run with impulse excitation unless in movie mode
if options.movie_mode
    inp = main.inp.sig;
else
    inp = zeros(size(main.inp.time));
    inp(1) = 1;
end

%For each domain, we need to know which nds/dfs we need to monitor at in
%main model
mon_nds = zeros(size(main.mod.nds, 1), 1);
for d = 1:numel(main.doms)
    mon_nds(main.doms{d}.mod.main_nd_i) = 1;
end

%Also need to include the transducer nodes so we can obtain the pristine
%response
%NB - no facility at the moment to do transmit on some and receive on
%others, it is always all on both

for t = 1:numel(main.tx)
    mon_nds(main.tx{t}.nds) = 1;
end

if ~isfield(main, 'rx')
    % main.rx = main.tx;
    separate_rx = 0;
else
    separate_rx = 1;
    for r = 1:numel(main.rx)
        mon_nds(main.rx{r}.nds) = 1;
    end
end

mon_nds = find(mon_nds);
[mon_nds, mon_dfs] = meshgrid(mon_nds, options.dof_to_use);

%Copy force input and monitoring info to each step (bit of duplication, but 
%allows general function to be used where this information is potentially 
%different for each step
 
for t = 1:numel(main.tx)
    steps{t}.load.time = main.inp.time;
    steps{t}.load.frcs = inp;
    steps{t}.load.frc_nds = main.tx{1}.nds;
    steps{t}.load.frc_dfs = main.tx{1}.dfs;
    steps{t}.mon.nds = mon_nds;
    steps{t}.mon.dfs = mon_dfs;
    steps{t}.mon.field_output_every_n_frames = options.field_output_every_n_frames;
end

%Case of separate rx needs some thought as to how to do - should be all
%part of same steps but need to sort out results at end
% if separate_rx
%     for t = 1:numel(main.rx)
%         steps{r}.load.time = main.inp.time;
%         steps{t}.load.frcs = inp;
%         steps{t}.load.frc_nds = main.tx{1}.nds;
%         steps{t}.load.frc_dfs = main.tx{1}.dfs;
%         steps{t}.mon.nds = mon_nds;
%         steps{t}.mon.dfs = mon_dfs;
%         steps{t}.mon.field_output_every_n_frames = options.field_output_every_n_frames;
%     end
% end

%Actually run the model for each transmitter
[res, mats] = fn_BristolFE_v2(main.mod, main.matls, steps, options);

%Parse pristine results to main.res{tx}.rx{rx}
[main.res.fmc.tx, main.res.fmc.rx] =meshgrid(1:numel(main.tx), 1:numel(main.rx));
main.res.fmc.tx = main.res.fmc.tx(:)';
main.res.fmc.rx = main.res.fmc.rx(:)';
main.res.fmc.time_data = zeros(numel(main.inp.time), numel(main.res.fmc.rx(:)));
for t = 1:numel(res)
    %copy field results
    main.res.tx{t}.fld = res{t}.fld;
    main.res.tx{t}.fld_time = res{t}.fld_time;

    for r = 1:numel(main.rx)
        k = find(t == main.res.fmc.tx & r == main.res.fmc.rx);
        i = ismember(res{t}.dsp_nds, mon_nds(main.rx{r}.nds));
        gl_nds = res{t}.dsp_nds(i);
        tmp = res{t}.dsps(i, :);
        if isfield(main.rx{r}, 'wt')
            tmp = main.rx{r}.wt(:)' * tmp;
        else
            tmp = sum(tmp);
        end
        main.res.fmc.time_data(:, k) = tmp(:); 
        % main.res{t}.rx{r}.dsps = main.res{t}.dsps(i, :);
        % main.res{t}.rx{r}.dsp_nds = main.res{t}.dsp_nds(i);
        % main.res{t}.rx{r}.dsp_dfs = main.res{t}.dsp_dfs(i);
    end
end

%At this point all the disp data is in res{tx}

%Parse domain results back to each domain. Need:
%   1. The incident displacements for all active DoFs on layers 2 and 3, 
%   which will be used in subsequent reciprocity calcs to project the 
%   scattered fields back to transducer
%   1. The forces, local nodes and DFs to apply in subsequent scatter
%   models (these will be all active DoFs on layers 2 and 3)
main.res.mats = mats;
for t = 1:numel(main.tx)
    for d = 1:numel(main.doms)
        %Work out the local nodes etc and how they map to global ones
        i = main.doms{d}.mod.bdry_lyrs > 0;
        loc_ndi = find(i);
        loc_lyrs = main.doms{d}.mod.bdry_lyrs(i);
        gl_nds = main.doms{d}.mod.main_nd_i(i);
        [i, j] = ismember(res{t}.dsp_nds, gl_nds);
        main.doms{d}.tx{t}.mn_nd_i = res{t}.dsp_nds(i);
        main.doms{d}.tx{t}.dsps = res{t}.dsps(i, :);
        %TO HERE
        main.doms{d}.tx{t}.dfs = res{t}.dsp_dfs(i);
        main.doms{d}.tx{t}.bdry_lyrs = loc_lyrs(j(i));
        % main.doms{d}.tx{t}.dsp_nds = loc_ndi(j(i));%because indexes into .main_nd_i are the local node numbers
        % main.doms{d}.tx{t}.dsp_lys = loc_lyrs(j(i));
        %At this point, we have the incident displacement field at the boundary nodes

        %Now convert it into a forcing field, firstly in terms of the nodes
        %in the domain model (which are not nesc same numbering as those in
        %scatterer models ... hence also keep track of global node numbers
        %for when those models are run with these forces)
        % gl_gi = fn_gl_ind_for_nd_and_dof(...
        %     main.res.mats.gl_lookup, main.doms{d}.tx{t}.mn_nd_i, ...
        %     main.doms{d}.tx{t}.dfs );
        % 
        % [main.doms{d}.tx{t}.frcs, force_set] = fn_convert_disps_to_forces_v2(...
        %     main.res.mats.K(gl_gi,gl_gi), main.res.mats.C(gl_gi,gl_gi), ...
        %     main.res.mats.M(gl_gi,gl_gi), time_step, main.doms{d}.tx{t}.dsps, ...
        %     main.doms{d}.tx{t}.bdry_lyrs, 'in');
        % main.doms{d}.res{t}.frc_nds = main.doms{d}.res{t}.dsp_nds(force_set);
        % main.doms{d}.res{t}.frc_gl_nds = main.doms{d}.res{t}.dsp_gl_nds(force_set);
        % main.doms{d}.res{t}.frc_dfs = main.doms{d}.res{t}.dsp_dfs(force_set);        %

    end
end


%New version end


%Generate main model global matrices if not already there
% if ~isfield(main, 'mats')
%     %build model stiffness matrices etc
%     % [main.mats.K, main.mats.C, main.mats.M, main.mats.Q, main.mats.gl_lookup] = fn_build_global_matrices_v3(...
%     %     main.mod.nds, main.mod.els, main.mod.el_mat_i, main.mod.el_abs_i, main.mod.matls);
%     [main.mats.K, main.mats.C, main.mats.M, main.mats.gl_lookup] = fn_build_global_matrices_v4(...
%         main.mod.nds, main.mod.els, main.mod.el_mat_i, main.mod.el_abs_i, main.mod.el_typ_i, main.mod.matls, options);
% 
%     %Work out global indices and weights for transducer positions
%     %this has to be done for all tx_rx because we still need the rx transfer functions even if only doing 1 tx
%     for r = 1:numel(main.mod.tx_rx) 
%         [main.mod.tx_rx{r}.gl_ind, gl_nds, gl_dofs, nd_i] = fn_global_indices_for_all_dof_at_nodes(main.mats.gl_lookup, main.mod.tx_rx{r}.nds);
%         main.mod.tx_rx{r}.wt = zeros(size(main.mod.tx_rx{r}.gl_ind));
%         for i = 1:numel(main.mod.tx_rx{r}.wt)
%             main.mod.tx_rx{r}.wt(i) = main.mod.tx_rx{r}.wt_by_dof(nd_i(i), gl_dofs(i));
%         end
%     end
% 
%     %force main model to be re-run
%     main = fn_clear_fields(main, 'res', 1, 1);
% end
% 
% %Run main model if necessary
% if ~isfield(main, 'res')
%     [main.res.hist_gi, main.doms] = fn_find_gl_inds_for_domains(main.mats.gl_lookup, main.doms);
% 
%     h_ref = zeros(size(main.res.hist_gi)); %this is just used for decoding the history responses immediately after running model
%     %work out which GL indices need to be monitored in full model for the tx rx
%     %positions (we need the time signal at each eventually)
%     for t = 1:numel(main.mod.tx_rx)
%         tmp = main.mod.tx_rx{t}.gl_ind;
%         main.res.hist_gi = [main.res.hist_gi; tmp];
%         h_ref  = [h_ref; ones(size(tmp)) * t];
%     end
%     i = 1; %iteration - always 1 for pristine model run
%     for t = 1:numel(tx_to_run)
%         fprintf('Full run of main model for excitation %i ', tx_to_run(t));
%         [h_out, main.res.tx_rx{tx_to_run(t)}.f_out, ~] = fn_explicit_dynamic_solver_v5(...
%             main.mats.K, main.mats.C, main.mats.M, main.mod.time, ...
%             main.mod.tx_rx{tx_to_run(t)}.gl_ind, main.mod.tx_rx{tx_to_run(t)}.wt * main.mod.input, ... %forcing input functions
%             [], [], ... %disp input functions
%             main.res.hist_gi, options.field_output_every_n_frames, options.use_gpu_if_available);
        %Extract displacement histories at monitoring nodes - thes need to
        %be deconvolved with input signal
    %     main.res.tx_rx{tx_to_run(t)}.hist = h_out(h_ref == 0, :); %h_ref == 0 flags the history outputs on boundary nodes of domains - this is what is extracted for running scatterer models
    %     % main.res.tx_rx{tx_to_run(t)}.deconv_hist = fn_deconvolve(h_out(h_ref == 0, :), main.mod.input, 2, 1e-3);
    % 
    %     %Extract A-scans for all receive elements
    %     main.res.tx_rx{tx_to_run(t)}.rx = zeros(numel(main.mod.tx_rx), numel(main.mod.time));
    %     for r = 1:numel(main.mod.tx_rx)
    %         main.res.tx_rx{tx_to_run(t)}.rx(r, :) = main.mod.tx_rx{r}.wt.' * h_out(h_ref == r, :);
    %     end
    % end
end

