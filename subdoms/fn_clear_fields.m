function main = fn_clear_fields(main, fieldname, clear_main, clear_doms)
if clear_main
    main = fn_safe_clear(main, fieldname);
end
if clear_doms
    for d = 1:numel(main.doms)
        main.doms{d} = fn_safe_clear(main.doms{d}, fieldname);
    end
end
end

function x = fn_safe_clear(x, fieldname)
if isfield(x, fieldname)
    x = rmfield(x, fieldname);
end
end
