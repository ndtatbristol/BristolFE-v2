function [fmc, val_fmc] = fn_extract_FMC_from_subdomain(main, d)
fmc = main.res.fmc;
fmc.time = main.inp.time;
fmc.time_data = fmc.time_data + main.doms{d}.res.fmc.time_data;
fmc.time_data = fn_convolve(fmc.time_data, main.inp.sig(:), 1);
if isfield(main.doms{d}, 'val') && isfield(main.doms{d}.val, 'fmc')
    val_fmc = main.res.fmc;
    val_fmc.time = main.inp.time;
    val_fmc.time_data = main.doms{d}.val.fmc.time_data;
    val_fmc.time_data = fn_convolve(val_fmc.time_data, main.inp.sig(:), 1);
end

end