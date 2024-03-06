function main = fn_clear_fields(main, fieldname, clear_main, clear_scats)
if clear_main
    main = fn_safe_clear(main, fieldname);
end
if clear_scats
    for d = 1:numel(main.doms)
        for s = 1:numel(main.doms{d}.scats)
            main.doms{d}.scats{s} = fn_safe_clear(main.doms{d}.scats{s}, fieldname);
        end
    end
end
end

function x = fn_safe_clear(x, fieldname)
if isfield(x, fieldname)
    x = rmfield(x, fieldname);
end
end
