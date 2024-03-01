clear
close all

%try and get deconvolution working properly

load('validation_result.mat');
validation_result = tmp;

load('toneburst_result.mat');
toneburst_result = tmp;

load('impulse_result.mat')
impulse_result = tmp;

load('inp.mat')

convolved_result = fn_convolve(impulse_result, inp, 2);
deconvolved_result = fn_deconvolve(toneburst_result, inp, 2);

% figure;
% plot(inp);

figure;
plot(validation_result);
hold on;
plot(convolved_result);
ylim([-1,1]*20000)
