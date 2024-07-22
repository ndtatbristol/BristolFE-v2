%Run all examples - used prior to build to check no errors
clear;
close all
restoredefaultpath;

scripts_to_run = {
    'solid_example.m',
    'fluid_example.m',
    'fluid_example_with_abs_layer.m',
    'solid_example_angled_excitation.m',
    'coupled_solid_fluid_example.m',
    'absorbing_layer_example.m',
    'subdomain_example.m',
    'subdomain_array_example.m'
    };

for i = 1:numel(scripts_to_run)
    fprintf([scripts_to_run{i}, '\n']);
    run(scripts_to_run{i});
end