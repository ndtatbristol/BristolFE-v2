function main = fn_run_subdomain_model(main, options)
default_options.doms_to_run = 1:numel(main.doms);
default_options.tx_trans = 1:numel(main.trans); %by default all transducers are transmitters
default_options.rx_trans = 1:numel(main.trans); %by default all transducers are also receivers
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
options = fn_set_default_fields(options, default_options);

time_step = main.inp.time(2) - main.inp.time(1);

%Run main model if not already run
if ~isfield(main, 'res')
    main = fn_run_GL_whole_model(main, options);
end

%Run the scatterer models
for d = options.doms_to_run
    main.doms{d}.res.fmc.time_data = zeros(size(main.res.fmc.time_data));
    for si = 1:numel(options.tx_trans) %si is step counter for FE, one step per transmitting transducer
        t = options.tx_trans(si);

        %Figure out which nodes in domain we need displacements for
        dm_bdry_nds_i = find(main.doms{d}.mod.bdry_lyrs > 0);
        
        %And the equivalent main nodes
        mn_bdry_nds_i = main.doms{d}.mod.main_nd_i(dm_bdry_nds_i);

        %Get the corresponding indices into the main results data
        %(mn_res_i) as well as the corresponding domain nodes / DoFs
        [mn_res_i, dm_res_nds_i] = fn_dm_to_mn(mn_bdry_nds_i, dm_bdry_nds_i, main.res.mon_nds);
        %Get the displacements at the boudary from the main model as well
        %as DoF and layer indices
        bdry_dsps = main.res.trans{t}.dsps(mn_res_i, :);
        bdry_dfs = main.res.mon_dfs(mn_res_i);
        bdry_lyrs = main.doms{d}.mod.bdry_lyrs(dm_res_nds_i);

        gl_i = fn_gl_ind_for_nd_and_dof(...
            main.res.mats.gl_lookup, ...
            main.res.mon_nds(mn_res_i), ...
            main.res.mon_dfs(mn_res_i));

        %Get relevant sub matrices
        K_sub = main.res.mats.K(gl_i, gl_i);
        C_sub = main.res.mats.C(gl_i, gl_i);
        M_sub = main.res.mats.M(gl_i, gl_i);

        %Convert to forces
        [frcs, frce_set] = fn_convert_disps_to_forces_v2(...
            K_sub, C_sub, M_sub, time_step, bdry_dsps, bdry_lyrs, 'in');

        %Create the load step - forcing
        steps{si}.load.frc_nds = dm_res_nds_i(frce_set);
        steps{si}.load.frc_dfs = bdry_dfs(frce_set);
        steps{si}.load.frcs = frcs;
        steps{si}.load.time = main.inp.time;

        %Create the load step - monitoring
        steps{si}.mon.nds = dm_res_nds_i;
        steps{si}.mon.dfs = bdry_dfs;
        steps{si}.mon.field_output_every_n_frames = options.field_output_every_n_frames;

    end

    %Run the model for all incident fields
    res = fn_BristolFE_v2(main.doms{d}.mod, main.matls, steps, options);

    %Parse the field results
    for si = 1:numel(options.tx_trans)
        t = options.tx_trans(si);
        main.doms{d}.res.trans{t}.fld = res{si}.fld;
        main.doms{d}.res.trans{t}.fld_time = res{si}.fld_time;
    end

    %Parse the history results (these are not saved, just used to calcualte
    %the received signals in the main model
    mn_all_nds_dfs = [main.res.mon_nds, main.res.mon_dfs];

    for si = 1:numel(options.tx_trans)
        t = options.tx_trans(si);

        %Convert to forces
        [frcs, frce_set] = fn_convert_disps_to_forces_v2(...
             K_sub, C_sub, M_sub, time_step, res{si}.dsps, bdry_lyrs, 'out');

        %Loop over receivers
        for r = options.rx_trans
            %Main nodes and DoFs associated with forcing points
            mn_nds_i = main.doms{d}.mod.main_nd_i(steps{t}.mon.nds(frce_set));
            mn_bdry_nds_dfs = [mn_nds_i, steps{t}.mon.dfs(frce_set)];
            
            %Convolve with relevant receiver transfer function
            i = ismember(mn_all_nds_dfs, mn_bdry_nds_dfs, 'rows');
            tmp = sum(fn_convolve(main.res.trans{r}.dsps(i, :), frcs, 2));
            
            %Stick it in the FMC data for this domain
            k = find(t == main.res.fmc.tx & r == main.res.fmc.rx);
            main.doms{d}.res.fmc.time_data(:, k) = tmp(:);
        end
    end
    
end

end

function [mn_res_i, dm_res_nds_i] = fn_dm_to_mn(mn_bdry_nds_i, dm_bdry_nds_i, mn_res_nds_i)
mn_res_i = ismember(mn_res_nds_i, mn_bdry_nds_i);
[~, j] = ismember(mn_res_nds_i(mn_res_i), mn_bdry_nds_i);
dm_res_nds_i = dm_bdry_nds_i(j);
end