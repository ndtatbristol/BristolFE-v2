function main = fn_run_subdomain_model(main, fe_options)
default_options.doms_to_run = [];
default_options.tx_trans = []; %by default all transducers are transmitters
default_options.rx_trans = []; %by default all transducers are also receivers
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
fe_options = fn_set_default_fields(fe_options, default_options);
fe_options.dof_to_use = fn_find_dof_in_use_and_max_dof_per_el(unique(main.mod.el_typ_i), fe_options.dof_to_use);

if isempty(fe_options.tx_trans)
    fe_options.tx_trans = 1:numel(main.trans); %by default all transducers are transmitters
end
if isempty(fe_options.rx_trans)
    fe_options.rx_trans = 1:numel(main.trans); %by default all transducers are transmitters
end
if isempty(fe_options.doms_to_run)
    fe_options.doms_to_run = 1:numel(main.doms);
end

time_step = main.inp.time(2) - main.inp.time(1);

%Run the scatterer models
for d = fe_options.doms_to_run
    main.doms{d} = fn_clear_fields(main.doms{d}, 'res', 1, 0);
    for si = 1:numel(fe_options.tx_trans) %si is step counter for FE, one step per transmitting transducer
        t = fe_options.tx_trans(si);

        %Get the mappings between nodes etc to main domain
        [mn_res_i, gl_i, bdry_nds, bdry_dfs, bdry_lyrs] = fn_get_main_subdomain_mappings(main.doms{d}.mod, main.res, main.res.mats.gl_lookup);

        %Get relevant sub matrices
        K_sub = main.res.mats.K(gl_i, gl_i);
        C_sub = main.res.mats.C(gl_i, gl_i);
        M_sub = main.res.mats.M(gl_i, gl_i);

        %Get relevant incident displacements
        bdry_dsps = main.res.trans{t}.dsps(mn_res_i, :);

        %Sub-domain models are run with true input signal rather than
        %impulse response, which is what is recorded
        bdry_dsps = fn_convolve(bdry_dsps, main.inp.sig, 2);

        %Convert to forces
        [frcs, frce_set] = fn_convert_disps_to_forces_v2(...
            K_sub, C_sub, M_sub, time_step, bdry_dsps, bdry_lyrs, 'in');

        %Create the load step - forcing
        steps{si}.load.frc_nds = bdry_nds(frce_set);
        steps{si}.load.frc_dfs = bdry_dfs(frce_set);
        steps{si}.load.frcs = frcs;
        steps{si}.load.time = main.inp.time;

        %Create the load step - monitoring
        steps{si}.mon.nds = bdry_nds;
        steps{si}.mon.dfs = bdry_dfs;
        steps{si}.mon.field_output_every_n_frames = fe_options.field_output_every_n_frames;

    end

    %Run the model for all incident fields
    % res = fn_BristolFE_v2(main.doms{d}.mod, main.matls, steps, fe_options);
    res = fn_FE_entry_point(main.doms{d}.mod, main.matls, steps, fe_options);

    %Copy the pristine FMC data into this domain's results - scatterered
    %results will be added on to this
    main.doms{d}.res.fmc = main.res.fmc;


    %Parse the field data for movies if requested
    if ~isinf(fe_options.field_output_every_n_frames)
        for si = 1:numel(fe_options.tx_trans)
            t = fe_options.tx_trans(si);
            main.doms{d}.res.trans{t}.fld = res{si}.fld;
            main.doms{d}.res.trans{t}.fld_time = res{si}.fld_time;
        end
    end

    %Parse the history results (these are not saved, just used to calculate
    %the received signals in the main model)
    mn_all_nds_dfs = [main.res.mon_nds, main.res.mon_dfs];

    for si = 1:numel(fe_options.tx_trans)
        t = fe_options.tx_trans(si);

        %Convert to forces
        [frcs, frce_set] = fn_convert_disps_to_forces_v2(...
            K_sub, C_sub, M_sub, time_step, res{si}.dsps, bdry_lyrs, 'out');

        %Loop over receivers
        for r = fe_options.rx_trans
            %Main nodes and DoFs associated with forcing points
            mn_nds_i = main.doms{d}.mod.main_nd_i(steps{t}.mon.nds(frce_set));
            mn_bdry_nds_dfs = [mn_nds_i, steps{t}.mon.dfs(frce_set)];

            %USE THE PARSE TO FMC DATA unction here instead!
            %Convolve with relevant receiver transfer function
            i = ismember(mn_all_nds_dfs, mn_bdry_nds_dfs, 'rows');
            tmp = sum(fn_convolve(main.res.trans{r}.dsps(i, :), frcs, 2));

            %Add it onto the pristine FMC data already copied into this
            %domain's results
            k = find(t == main.res.fmc.tx & r == main.res.fmc.rx);
            main.doms{d}.res.fmc.time_data(:, k) = main.doms{d}.res.fmc.time_data(:, k) + tmp(:);
        end
    end
end

end


function [mn_res_i, gl_i, bdry_nds, bdry_dfs, bdry_lyrs] = fn_get_main_subdomain_mappings(dom_mod, main_res, main_gl_lookup)

%Find boundary nodes in subdomain
bdry_nds = find(dom_mod.bdry_lyrs > 0);

%Equivalent ones in main
nb_m = dom_mod.main_nd_i(bdry_nds);

%Global matrix indices associated with these nodes
G_mon = main_gl_lookup(sub2ind(size(main_gl_lookup),main_res.mon_nds,main_res.mon_dfs));

mn_res_i = find(ismember(main_res.mon_nds, nb_m));
gl_i = main_res.dsp_gi(mn_res_i);

Nb_m = main_res.mon_nds(mn_res_i);
Db_m = main_res.mon_dfs(mn_res_i);

%Now need to work back to find equvalent sub-domain nodes
bdry_nds = interp1(nb_m, bdry_nds, Nb_m, 'nearest');

%And then the layers
bdry_lyrs = dom_mod.bdry_lyrs(bdry_nds);

%Finally DoFs which are same in main and sub-domain
bdry_dfs = Db_m;

end

