function fn_run_animation(h_patch, fld, anim_options)
%SUMMARY
%   Animates field output on a previously displayed model geometry. The
%   wave intensity is presented as the local kinetic energy on a dB scale
%   (default range is 40dB, normalised to peak amplitude over all locations
%   and all times). The range and normalisation value can be over-ridden by
%   setting the appropriate options (see below).

default_options.pause_value = 0.001;
default_options.max_color = [1,1,1];
default_options.min_color = [0,0,0];
default_options.mp4_out = [];
default_options.frame_rate = 48;
default_options.db_range = [-40, 0];
default_options.min_ti = 1;
default_options.max_ti = inf;
default_options.ti_step = 1;
default_options.norm_val = [];
default_options.repeat_n_times = 1;
default_options.db_or_linear = 'db';
default_options.clear_at_end = 1;
anim_options = fn_set_default_fields(anim_options, default_options);

if ~iscell(h_patch)
    h_patch = {h_patch};
end

if ~iscell(fld)
    fld = {fld};
end

if isempty(anim_options.norm_val)
    anim_options.norm_val = 0;
    for f = 1:numel(fld)
        if ~isempty(fld{f})
            anim_options.norm_val = max(max(abs(fld{f}(:)), [], 'all'), anim_options.norm_val);
        end
    end
end


if ~isempty(anim_options.mp4_out)
    if isinf(anim_options.repeat_n_times)
        %Avoid making infinitely long videos!
        anim_options.repeat_n_times = 1;
    end
    video_profile = 'MPEG-4';
    vidObj = VideoWriter(anim_options.mp4_out, video_profile);
    vidObj.FrameRate = anim_options.frame_rate;
    vidObj.Quality = 100;
    open(vidObj);

    switch video_profile
        case 'MPEG-4'
            opengl software;
            set(gcf, 'Color', [1,1,1]);
            set(gcf,'Renderer','OpenGL');
        otherwise
            set(gcf,'Renderer','painters');
    end
    set(gcf, 'Color', [1,1,1]);
end


for a = 1:anim_options.repeat_n_times
    for ti = anim_options.min_ti:anim_options.ti_step:min(anim_options.max_ti, size(fld{1}, 2))
        for i = 1:numel(h_patch)
            if ~isempty(fld{i})
                if ti == anim_options.min_ti && a == 1
                    base_cdata{i} = get(h_patch{i}, 'CData');
                    max_cdata = permute(anim_options.max_color, [3,1,2]);
                    min_cdata = permute(anim_options.min_color, [3,1,2]);
                end
                %Note animation data is already in energy, hence dB is 10 * log10(.),
                %not 20 * log10(.)
                switch anim_options.db_or_linear
                    case 'db'
                        v = (20 * log10(abs(fld{i}(:, ti)) / anim_options.norm_val) - anim_options.db_range(1)) / (anim_options.db_range(2) - anim_options.db_range(1));
                        % v(isnan(v)) = 0;
                        v(v > 1) = 1;
                        v(v < 0) = 0;
                        v = v .* sign(fld{i}(:, ti));
                    case 'linear'
                        v = fld{i}(:, ti) / anim_options.norm_val;
                        v(v < -1) = -1;
                        v(v > 1) = 1;
                end
                % set(h_patch{i}, 'CData', base_cdata{i} .* (1 - abs(v)) + abs(v) .* (v > 0) .* max_cdata + abs(v) .* (v < 0) .* min_cdata);
                set(h_patch{i}, 'CData', base_cdata{i} .* (1 - v) + v .* max_cdata);
            end
        end
        title(sprintf('Frame %i', ti))
        pause(anim_options.pause_value);
        if ~isempty(anim_options.mp4_out)
            writeVideo(vidObj,getframe(gcf));
        end
    end
    %Reset figure
    if anim_options.clear_at_end
        for i = 1:numel(h_patch)
            set(h_patch{i}, 'CData', base_cdata{i});
        end
    end
end

if ~isempty(anim_options.mp4_out)
    close(vidObj);
end

end

