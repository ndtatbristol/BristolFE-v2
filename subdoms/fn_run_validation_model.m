function main = fn_run_validation_model(main, fe_options)
default_options.doms_to_run = []; %by default run for all
default_options.tx_trans = []; %by default all transducers are transmitters
default_options.rx_trans = []; %by default all transducers are also receivers
default_options.field_output_every_n_frames = inf;
default_options.use_gpu_if_available = 1;
default_options.dof_to_use = [];
fe_options = fn_set_default_fields(fe_options, default_options);

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
    %Create whole domain mesh that includes the scatterer
    

end

end