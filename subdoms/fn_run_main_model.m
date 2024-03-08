function main = fn_run_main_model(main, fe_options)

%New version to sort out deconvolution issue

default_options.doms_to_run = []; %only relevant in validation mode
default_options.centre_freq = main.mod.design_centre_freq;
default_options.number_of_cycles = 5;
default_options.time_step = main.mod.max_safe_time_step;
default_options.time_pts = 1000;
default_options.max_time = [];
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
default_options.tx_trans = 1:numel(main.trans); %by default all transducers are transmitters
default_options.rx_trans = 1:numel(main.trans); %by default all transducers are also receivers
default_options.validation_mode = 0;
default_options.max_damping = [];
default_options.damping_power_law = 3;
default_options.max_stiffness_reduction = 0.01;
fe_options = fn_set_default_fields(fe_options, default_options);
if isempty(fe_options.doms_to_run)
    fe_options.doms_to_run = 1:numel(main.doms);
end
if isempty(fe_options.max_damping)
    fe_options.max_damping = fe_options.centre_freq * 2 * pi;
end
if ~isfield(main.mod, 'el_typ_i')
    main.mod.el_typ_i = {main.matls(main.mod.el_mat_i).el_typ};
    main.mod.el_typ_i = main.mod.el_typ_i(:);
end



fe_options.dof_to_use = fn_find_dof_in_use_and_max_dof_per_el(unique(main.mod.el_typ_i), fe_options.dof_to_use);

%Input signal and time-axis used for all simulations
if ~isempty(fe_options.max_time)
    fe_options.time_pts = ceil(fe_options.max_time / fe_options.time_step);
end
main.inp.time = [0:fe_options.time_pts - 1] * fe_options.time_step;
main.inp.sig = fn_gaussian_pulse(main.inp.time, fe_options.centre_freq, fe_options.number_of_cycles);

%Specify which full-domain models need to be run depending on options. Note that
%requesting field output causes the number of runs to be doubled as one is
%run with impulse excitation to obtain transfer functions and second is run
%with toneburst excitation to generate nice field outputs for movies
if fe_options.validation_mode
    main_modes = {'validation'};
    main = fn_clear_fields(main, 'val', 1, 1);
else
    if ~isinf(fe_options.field_output_every_n_frames)
        main_modes = {'impulse response', 'field output'};
    else
        main_modes = {'impulse response'};
        main = fn_clear_fields(main, 'res', 1, 1);
    end
end

%Run the full-domain model for each case required
for m = 1:numel(main_modes)
    %Flag monitoring nodes:
    mon_nds = zeros(size(main.mod.nds, 1), 1);

    %First, for sub-domains bdries (only needed in impulse mode)
    if strcmp(main_modes{m}, 'impulse response')
        for d = 1:numel(main.doms)
            mon_nds(main.doms{d}.mod.main_nd_i(main.doms{d}.mod.bdry_lyrs > 0)) = 1;
        end
    end

    %Second, for the transducer nodes
    for e = 1:numel(main.trans)
        mon_nds(main.trans{e}.nds) = 1;
    end

    %Create list of monitoring nodes and DoFs
    mon_nds = find(mon_nds);
    [mon_nds, mon_dfs] = meshgrid(mon_nds, fe_options.dof_to_use);
    mon_nds = mon_nds(:);
    mon_dfs = mon_dfs(:);

    %Choose appropriate input signal
    if strcmp(main_modes{m}, 'impulse response')
        inp = zeros(size(main.inp.time));
        inp(1) = 1;
    else
        %For validation and field output, the actual input signal is used
        inp = main.inp.sig;
    end

    %Prepare the basic info for FMC output
    if any(strcmp(main_modes{m}, {'impulse response', 'validation'}))
        if ~isfield(main, 'res')
            main.res.fmc_template = fn_create_fmc_template(main, fe_options);
        end
    end

    %Prepare the steps
    for e = 1:numel(main.trans)
        steps{e} = fn_convert_to_step_data(...
            main.inp.time, inp, ...
            main.trans{e}.nds, main.trans{e}.dfs, ...
            mon_nds, mon_dfs, fe_options.field_output_every_n_frames);
    end

    switch main_modes{m}
         case {'impulse response', 'field output'}
            %Run the model for each transducer (need boundary data whether
            %transmitter or receiver anyway
            % [fe_res, main.res.mats] = fn_BristolFE_v2(main.mod, main.matls, steps, fe_options);
            [fe_res, main.res.mats] = fn_FE_entry_point(main.mod, main.matls, steps, fe_options);
            

            %Parse the impulse data
            if strcmp(main_modes{m}, 'impulse response')
                %Sub-domain boundary displacements - these ALWAYS hold
                %impulse responses
                main.res = fn_parse_to_bdry_nds(main.res, fe_res, mon_nds, mon_dfs);
                %Pristine FMC results
                % main.res.fmc = fn_parse_to_fmc(steps, fe_res, main.trans, fe_options, main.inp.sig);
                main.res.fmc = fn_parse_to_fmc(main.res.fmc_template, steps, fe_res, main.trans, main.inp.sig);
            end

            %Parse the pristine toneburst field data for movies
            if strcmp(main_modes{m}, 'field output')
                main.res = fn_parse_to_movie(main.res, fe_res);
            end

        case 'validation'
            for d = fe_options.doms_to_run
                %Create a new model of whole domain containing the scatterer
                main.doms{d}.val_mod = fn_insert_subdomain_model_into_main(main.mod, main.doms{d}.mod, main.matls);

                %Run the model for each transducer
                % [fe_res, main.res.mats] = fn_BristolFE_v2(main.doms{d}.val_mod, main.matls, steps, fe_options);
                [fe_res, main.res.mats] = fn_FE_entry_point(main.doms{d}.val_mod, main.matls, steps, fe_options);

                %Parse the field data for movies if requested
                if ~isinf(fe_options.field_output_every_n_frames)
                    main.doms{d}.val = fn_parse_to_movie([], fe_res);
                end

                %Parse the validation FMC results
                % main.doms{d}.val.fmc = fn_parse_to_fmc(steps, fe_res, main.trans, fe_options, []);
                main.doms{d}.val.fmc = fn_parse_to_fmc(main.res.fmc_template, steps, fe_res, main.trans, []);
            end
       
    end
end
end


%--------------------------------------------------------------------------
function res = fn_parse_to_movie(res, fe_res)
for e = 1:numel(fe_res)
    res.trans{e}.fld = fe_res{e}.fld;
    res.trans{e}.fld_time = fe_res{e}.fld_time;
end
end

function res = fn_parse_to_bdry_nds(res, fe_res, mon_nds, mon_dfs)
valid_mon_dsps = fe_res{1}.valid_mon_dsps; %same for all steps
res.dsp_gi = fe_res{1}.dsp_gi;
res.mon_nds = mon_nds(valid_mon_dsps);
res.mon_dfs = mon_dfs(valid_mon_dsps);
for e = 1:numel(fe_res)
    res.trans{e}.dsps = fe_res{e}.dsps(valid_mon_dsps,:);
end

end

function fmc_template = fn_create_fmc_template(main, fe_options)
[fmc_template.tx, fmc_template.rx] = meshgrid(fe_options.tx_trans, fe_options.rx_trans);
fmc_template.tx = fmc_template.tx(:)';
fmc_template.rx = fmc_template.rx(:)';
fmc_template.time_data = zeros(fe_options.time_pts, numel(fmc_template.rx(:)));
fmc_template.time = main.inp.time(:);
ne = numel(main.trans);

%TODO Need to add other array geom details in usual format here
fmc_template.array.centre_freq = fe_options.centre_freq;
fmc_template.array.el_xc = zeros(1, ne);
fmc_template.array.el_yc = zeros(1, ne);
fmc_template.array.el_zc = zeros(1, ne);
fmc_template.array.el_x1 = zeros(1, ne);
fmc_template.array.el_y1 = zeros(1, ne);
fmc_template.array.el_z1 = zeros(1, ne);
fmc_template.array.el_x2 = zeros(1, ne);
fmc_template.array.el_y2 = zeros(1, ne);
fmc_template.array.el_z2 = zeros(1, ne);

fmc_template.array.el_type = 'rectangular';
for t = 1:ne
    c = mean(main.mod.nds(main.trans{t}.nds, :));
    m = max(main.mod.nds(main.trans{t}.nds, :));
    if numel(c) == 2
        %2D models
        fmc_template.array.el_xc(t) = c(1);
        fmc_template.array.el_zc(t) = c(2);

        %in 2d model, point 1 is in xz plane at extremum of element
        fmc_template.array.el_x1(t) = m(1);
        fmc_template.array.el_z1(t) = m(2);

        %in 2d model point in centre of element in xy plane but offset in
        %y, here by unit distance as there is no physical dimension of
        %elements in this direction in FE model
        fmc_template.array.el_x2(t) = c(1);
        fmc_template.array.el_y2(t) = 1;
        fmc_template.array.el_z2(t) = c(2);
    else
        %3D models
        fmc_template.array.el_xc(t) = c(1);
        fmc_template.array.el_yc(t) = c(2);
        fmc_template.array.el_zc(t) = c(3);

        %following should be done properly - here they assume elements lie
        %in x-y plane, but really should do PCA of nodes to determine
        %pinciple axes of element and go from there.
        fmc_template.array.el_x1(t) = m(1);
        fmc_template.array.el_y1(t) = c(2);
        fmc_template.array.el_z1(t) = c(3);

        fmc_template.array.el_x2(t) = c(1);
        fmc_template.array.el_y2(t) = m(2);
        fmc_template.array.el_z2(t) = c(3);
    end
end
end


function fmc = fn_parse_to_fmc(fmc_template, steps, res, trans, in_sig)
fmc = fmc_template;
for k = 1:numel(fmc.tx)
    i = ismember(steps{fmc.tx(k)}.mon.nds, trans{fmc.rx(k)}.nds);
    tmp = res{fmc.tx(k)}.dsps(i, :);
    if isfield(trans{fmc.rx(k)}, 'wt')
        tmp = trans{fmc.rx(k)}.wt(:)' * tmp;
    else
        tmp = sum(tmp);
    end
    %Convolved with input if required (which is case for pristine
    %results which are generated for impulse response)
    fmc.time_data(:, k) = tmp(:);
    if ~isempty(in_sig)
        fmc.time_data(:, k) = fn_convolve(fmc.time_data(:, k), in_sig(:), 1);
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

function main = fn_clear_fields(main, fieldname, clear_main, clear_doms)
if clear_main
    main = fn_safe_clear(main, fieldname);
end
if clear_doms
    for d = 1:numel(main.doms)
            main.doms{d} = fn_safe_clear(main.doms{d}, fieldname);
    end
end
end

function x = fn_safe_clear(x, fieldname)
if isfield(x, fieldname)
    x = rmfield(x, fieldname);
end
end