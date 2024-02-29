function main = fn_run_main_model(main, time_pts, fe_options)
default_options.doms_to_run = []; %only relevant in validation mode
default_options.centre_freq = main.mod.design_centre_freq;
default_options.number_of_cycles = 5;
default_options.time_step = main.mod.max_safe_time_step;
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
default_options.tx_trans = 1:numel(main.trans); %by default all transducers are transmitters
default_options.rx_trans = 1:numel(main.trans); %by default all transducers are also receivers
default_options.validation_mode = 0;
fe_options = fn_set_default_fields(fe_options, default_options);
if isempty(fe_options.doms_to_run)
    fe_options.doms_to_run = 1:numel(main.doms);
end

fe_options.dof_to_use = fn_find_dof_in_use_and_max_dof_per_el(unique(main.mod.el_typ_i), fe_options.dof_to_use);

%Runs everything necessary
time_step = fe_options.time_step;
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
main.inp.sig = fn_gaussian_pulse(main.inp.time, fe_options.centre_freq, fe_options.number_of_cycles);

%Main model is always run with impulse excitation unless in movie mode
if ~isinf(fe_options.field_output_every_n_frames)
    inp = main.inp.sig;
    warning('Movie mode - toneburst input will be used to generate nice animations but time-domain outputs will be suppressed')
else
    inp = zeros(size(main.inp.time));
    inp(1) = 1;
end

if ~fe_options.validation_mode
    %For each domain, we need to know which nds/dfs we need to monitor at in
    %main model. Not needed in validation mode
    mon_nds = zeros(size(main.mod.nds, 1), 1);
    for d = 1:numel(main.doms)
        mon_nds(main.doms{d}.mod.main_nd_i(main.doms{d}.mod.bdry_lyrs > 0)) = 1;
    end
end

%Also need to include the transducer nodes so we can obtain the pristine
%response
%NB - no facility at the moment to do transmit on some and receive on
%others, it is always all on both

for e = 1:numel(main.trans)
    mon_nds(main.trans{e}.nds) = 1;
end

mon_nds = find(mon_nds);
[mon_nds, mon_dfs] = meshgrid(mon_nds, fe_options.dof_to_use);
mon_nds = mon_nds(:);
mon_dfs = mon_dfs(:);

%Copy force input and monitoring info to each step (bit of duplication, but
%allows general function to be used where this information is potentially
%different for each step

for e = 1:numel(main.trans)
    steps{e} = fn_convert_to_step_data(...
        main.inp.time, inp, ...
        main.trans{e}.nds, main.trans{e}.dfs, ...
        mon_nds, mon_dfs, fe_options.field_output_every_n_frames);
end

if ~fe_options.validation_mode
    %NORMAL MODE

    %Run the model for each transducer (need boundary data whether
    %transmitter or receiver anyway
    [res, main.res.mats] = fn_BristolFE_v2(main.mod, main.matls, steps, fe_options);

    %Parse the essential data (the displacements at the subdomain boundaries)
    %which is needed regardless of whether looking at field animations or FMC
    valid_mon_dsps = res{1}.valid_mon_dsps; %same for all steps
    main.res.dsp_gi = res{1}.dsp_gi;
    main.res.mon_nds = mon_nds(valid_mon_dsps);
    main.res.mon_dfs = mon_dfs(valid_mon_dsps);
    for e = 1:numel(main.trans)
        %copy field results
        % main.res.trans{e} = res{e};
        main.res.trans{e}.dsps = res{e}.dsps(valid_mon_dsps,:);
    end

    if ~isinf(fe_options.field_output_every_n_frames)
        %Parse pristine field results to main.res{tx}.trans{t}
        for e = 1:numel(main.trans)
            %copy field results
            main.res.trans{e}.fld = res{e}.fld;
            main.res.trans{e}.fld_time = res{e}.fld_time;
        end
        main.res.fmc = [];
    else
        %Parse the pristine FMC results
        [main.res.fmc.tx, main.res.fmc.rx] =meshgrid(fe_options.tx_trans, fe_options.rx_trans);
        main.res.fmc.tx = main.res.fmc.tx(:)';
        main.res.fmc.rx = main.res.fmc.rx(:)';
        main.res.fmc.time_data = zeros(numel(main.inp.time), numel(main.res.fmc.rx(:)));

        %TODO Need to add array details in usual format here

        for t = fe_options.tx_trans
            for r = fe_options.rx_trans
                k = find(t == main.res.fmc.tx & r == main.res.fmc.rx);
                % i = ismember(res{t}.dsp_nds, mon_nds(main.trans{r}.nds));
                i = ismember(steps{t}.mon.nds, main.trans{r}.nds);

                % gl_nds = res{t}.dsp_nds(i);
                tmp = res{t}.dsps(i, :);
                if isfield(main.trans{r}, 'wt')
                    tmp = main.trans{r}.wt(:)' * tmp;
                else
                    tmp = sum(tmp);
                end
                main.res.fmc.time_data(:, k) = tmp(:);
            end
        end
    end
else
    %VALIDATION MODE
    for d = fe_options.doms_to_run
        main.doms{d}.val_mod = fn_insert_subdomain_model_into_main(main.mod, main.doms{d}.mod);

        %Run the model for each transducer (need boundary data whether
        %transmitter or receiver anyway
        [res, main.res.mats] = fn_BristolFE_v2(main.doms{d}.val_mod, main.matls, steps, fe_options);
    
        if ~isinf(fe_options.field_output_every_n_frames)
            %Parse pristine field results to main.val.trans{t}
            for e = 1:numel(main.trans)
                %copy field results
                main.doms{d}.val.trans{e}.fld = res{e}.fld;
                main.doms{d}.val.trans{e}.fld_time = res{e}.fld_time;
            end
            main.doms{d}.val.fmc = [];
        else
            %Parse the pristine FMC results
            [main.doms{d}.val.fmc.tx, main.doms{d}.val.fmc.rx] =meshgrid(fe_options.tx_trans, fe_options.rx_trans);
            main.doms{d}.val.fmc.tx = main.doms{d}.val.fmc.tx(:)';
            main.doms{d}.val.fmc.rx = main.doms{d}.val.fmc.rx(:)';
            main.doms{d}.val.fmc.time_data = zeros(numel(main.inp.time), numel(main.doms{d}.val.fmc.rx(:)));
    
            %TODO Need to add array details in usual format here
    
            for t = fe_options.tx_trans
                for r = fe_options.rx_trans
                    k = find(t == main.doms{d}.val.fmc.tx & r == main.doms{d}.val.fmc.rx);
                    % i = ismember(res{t}.dsp_nds, mon_nds(main.trans{r}.nds));
                    i = ismember(steps{t}.mon.nds, main.trans{r}.nds);
    
                    % gl_nds = res{t}.dsp_nds(i);
                    tmp = res{t}.dsps(i, :);
                    if isfield(main.trans{r}, 'wt')
                        tmp = main.trans{r}.wt(:)' * tmp;
                    else
                        tmp = sum(tmp);
                    end
                    main.doms{d}.val.fmc.time_data(:, k) = tmp(:);
                end
            end
        end
    end
end

end

function step = fn_convert_to_step_data(t, inp, frc_nds, frc_dfs, mon_nds, mon_dfs, f_every)
step.load.time = t;
step.load.frcs = inp;
step.load.frc_nds = frc_nds;
step.load.frc_dfs = frc_dfs;
step.mon.nds = mon_nds;
step.mon.dfs = mon_dfs;
step.mon.field_output_every_n_frames = f_every;
end

