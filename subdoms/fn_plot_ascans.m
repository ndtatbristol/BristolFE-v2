function doms_and_scats_to_show = fn_plot_ascans(main, options)
default_options.doms_and_scats_to_show = [];
default_options.show_what = {'val', 'di', 'res'};
default_options.convolve_with_desired_input = 1;
options = fn_set_default_fields(options, default_options);
doms_and_scats_to_show = options.doms_and_scats_to_show;
show_what = options.show_what;
filter_on = options.convolve_with_desired_input;
%--------------------------------------------------------------------------
if isempty(doms_and_scats_to_show)
    doms_and_scats_to_show = [0, 0];
    i = 1;
    for d = 1:numel(main.doms)
        for s = 1:numel(main.doms{1}.scats)
            doms_and_scats_to_show(i, :) = [d, s];
            i = i + 1;
        end
    end
end
unique_doms = unique(doms_and_scats_to_show(:,1));
unique_scats = unique(doms_and_scats_to_show(:,2));
clf;
t = 1; r = 1; 
for i = 1:size(doms_and_scats_to_show, 1)
    rw = find(unique_doms == doms_and_scats_to_show(i, 1));
    cl = find(unique_scats == doms_and_scats_to_show(i, 2));
    subplot(numel(unique_doms), numel(unique_scats), (rw - 1) *numel(unique_scats) + cl);
    d = doms_and_scats_to_show(i,1);
    s = doms_and_scats_to_show(i,2);
    leg = {};
    for j = 1:numel(show_what)
        switch show_what{j}
            case 'res'
                if isfield(main.doms{d}.scats{s}, show_what{j})
                    data = main.doms{d}.scats{s}.res.tx_rx{t}.rx(r, :);
                    col = 'r';
                    leg{end + 1} = 'New method';
                else 
                    data = [];
                end
            case 'val'
                if isfield(main.doms{d}.scats{s}, show_what{j})
                    data = main.doms{d}.scats{s}.val.res.tx_rx{1}.rx(r, :);
                    col = [1,1,1] * 0.5;
                    leg{end + 1} = 'Validation';
                else 
                    data = [];
                end
            case 'di'
                if isfield(main.doms{d}.scats{s}, show_what{j})
                    data = main.doms{d}.scats{s}.di.res.tx_rx{t}.rx(r, :);
                    col = 'b';
                    leg{end + 1} = 'Direct inj.';
                else 
                    data = [];
                end
        end
        if ~isempty(data)
            hold on;
            if filter_on
                data = fn_convolve(data, main.mod.desired_input, 2);
            end
            plot(main.mod.time, data, 'Color', col);
        end
    end
    legend(leg);
end

end