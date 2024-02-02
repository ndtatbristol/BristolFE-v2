function fmc = fn_extract_FMC_from_main(main)
fmc = main.res.fmc;
fmc.time = main.inp.time;
fmc.time_data = fn_convolve(fmc.time_data, main.inp.sig(:), 1);
end