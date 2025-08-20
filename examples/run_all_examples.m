%Run all examples - used prior to build to check no errors
clear all;
close all
restoredefaultpath;

scripts_to_run = {
    'solid_example.m',
    'fluid_example.m',
    'fluid_solid_with_absorbing_layer_example.m',
    'solid_example_angled_excitation.m',
    'subdomain_example.m',
    'subdomain_array_example.m'
    };

for i = 1:numel(scripts_to_run)
    fprintf('\n\n--------------------------------------------------------------\n')
    fprintf([scripts_to_run{i}, '\n']);
    fprintf('--------------------------------------------------------------\n\n')
    run(scripts_to_run{i});
end

clear all;
