function main = fn_run_main_G_to_L_model(main, time_pts, options)
default_options.centre_freq = main.mod.design_centre_freq;
default_options.number_of_cycles = 5;
default_options.time_step = main.mod.max_safe_time_step;
default_options.tx_to_run = 1:numel(main.mod.tx_rx);
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
default_options.movie_mode = 0; %when this is one, the input signal is set to be the toneburst rather than an impulse, so movies look nice. Other results will be wrong!!!
options = fn_set_default_fields(options, default_options);

%Runs everything necessary
time_step = options.time_step;
if isempty(time_pts)
    %choose  number of points for slowest wave to traverse largest
    %dimension of model
    time_pts = round(sqrt(sum((max(main.mod.nds) - min(main.mod.nds)) .^ 2)) / main.mod.design_min_vel / time_step);
end

%Force re-run if time axis has changed
if ~isfield(main.mod, 'time') || time_pts ~= numel(main.mod.time)  || time_step ~= (main.mod.time(2) - main.mod.time(1))
    warning('Time axis has changed - all models will be re-run');
    main = fn_clear_fields(main, 'res', 1, 1);
    main = fn_clear_fields(main, 'mats', 1, 1);
    %Create new time axis
    if time_step > main.mod.max_safe_time_step
        warning('Specified time step is greater than maximum recommended time step')
    end
    main.mod.time = [0:time_pts - 1] * time_step;
end

%Main model is always run with impulse excitation
main.mod.input = zeros(size(main.mod.time));
main.mod.input(1) = 1;

%The desired output only affects how final A-scans are filtered; it does
%not change the FE runs, which are performed with impulse excitation.
main.mod.desired_input = fn_gaussian_pulse(main.mod.time, options.centre_freq, options.number_of_cycles);

if options.movie_mode
    main.mod.input = main.mod.desired_input;
end

tx_to_run = options.tx_to_run;

%Generate main model global matrices if not already there
if ~isfield(main, 'mats')
    %build model stiffness matrices etc
    % [main.mats.K, main.mats.C, main.mats.M, main.mats.Q, main.mats.gl_lookup] = fn_build_global_matrices_v3(...
    %     main.mod.nds, main.mod.els, main.mod.el_mat_i, main.mod.el_abs_i, main.mod.matls);
    [main.mats.K, main.mats.C, main.mats.M, main.mats.gl_lookup] = fn_build_global_matrices_v4(...
        main.mod.nds, main.mod.els, main.mod.el_mat_i, main.mod.el_abs_i, main.mod.el_typ_i, main.mod.matls, options);
    
    %Work out global indices and weights for transducer positions
    %this has to be done for all tx_rx because we still need the rx transfer functions even if only doing 1 tx
    for r = 1:numel(main.mod.tx_rx) 
        [main.mod.tx_rx{r}.gl_ind, gl_nds, gl_dofs, nd_i] = fn_global_indices_for_all_dof_at_nodes(main.mats.gl_lookup, main.mod.tx_rx{r}.nds);
        main.mod.tx_rx{r}.wt = zeros(size(main.mod.tx_rx{r}.gl_ind));
        for i = 1:numel(main.mod.tx_rx{r}.wt)
            main.mod.tx_rx{r}.wt(i) = main.mod.tx_rx{r}.wt_by_dof(nd_i(i), gl_dofs(i));
        end
    end

    %force main model to be re-run
    main = fn_clear_fields(main, 'res', 1, 1);
end

%Run main model if necessary
if ~isfield(main, 'res')
    [main.res.hist_gi, main.doms] = fn_find_gl_inds_for_domains(main.mats.gl_lookup, main.doms);

    h_ref = zeros(size(main.res.hist_gi)); %this is just used for decoding the history responses immediately after running model
    %work out which GL indices need to be monitored in full model for the tx rx
    %positions (we need the time signal at each eventually)
    for t = 1:numel(main.mod.tx_rx)
        tmp = main.mod.tx_rx{t}.gl_ind;
        main.res.hist_gi = [main.res.hist_gi; tmp];
        h_ref  = [h_ref; ones(size(tmp)) * t];
    end
    i = 1; %iteration - always 1 for pristine model run
    for t = 1:numel(tx_to_run)
        fprintf('Full run of main model for excitation %i ', tx_to_run(t));
        [h_out, main.res.tx_rx{tx_to_run(t)}.f_out, ~] = fn_explicit_dynamic_solver_v5(...
            main.mats.K, main.mats.C, main.mats.M, main.mod.time, ...
            main.mod.tx_rx{tx_to_run(t)}.gl_ind, main.mod.tx_rx{tx_to_run(t)}.wt * main.mod.input, ... %forcing input functions
            [], [], ... %disp input functions
            main.res.hist_gi, options.field_output_every_n_frames, options.use_gpu_if_available);
        %Extract displacement histories at monitoring nodes - thes need to
        %be deconvolved with input signal
        main.res.tx_rx{tx_to_run(t)}.hist = h_out(h_ref == 0, :); %h_ref == 0 flags the history outputs on boundary nodes of domains - this is what is extracted for running scatterer models
        % main.res.tx_rx{tx_to_run(t)}.deconv_hist = fn_deconvolve(h_out(h_ref == 0, :), main.mod.input, 2, 1e-3);

        %Extract A-scans for all receive elements
        main.res.tx_rx{tx_to_run(t)}.rx = zeros(numel(main.mod.tx_rx), numel(main.mod.time));
        for r = 1:numel(main.mod.tx_rx)
            main.res.tx_rx{tx_to_run(t)}.rx(r, :) = main.mod.tx_rx{r}.wt.' * h_out(h_ref == r, :);
        end
    end
end

