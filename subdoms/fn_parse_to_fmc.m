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