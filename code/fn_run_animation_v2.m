function fn_run_animation_v2(h_patch, fld, options)
default_options.pause_value = 0.001;
default_options.wave_color = [1,1,1];
default_options.mp4_out = [];
default_options.frame_rate = 48;
default_options.db_range = [-40, 0];
default_options.min_ti = 1;
default_options.max_ti = inf;
default_options.ti_step = 1;
options = fn_set_default_fields(options, default_options);

if ~iscell(h_patch)
    % tmp = animation_data;
    % clear('animation_data');
    % animation_data{1} = tmp;
    h_patch = {h_patch};
    fld = {fld};
end

% max_ti = 0;
% min_ti = inf;
% for i = 1:numel(h_patch)
%     max_ti = max(max_ti, animation_data{i}.max_ti);
%     min_ti = min(min_ti, animation_data{i}.min_ti);
% end

if ~isempty(options.mp4_out)
    video_profile = 'MPEG-4';
    vidObj = VideoWriter(options.mp4_out, video_profile);
    vidObj.FrameRate = options.frame_rate;
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

for ti = options.min_ti:options.ti_step:min(options.max_ti, size(fld{1}, 2))
    for i = 1:numel(h_patch)
        if ti == options.min_ti
            base_cdata{i} = get(h_patch{i}, 'CData');
            max_cdata = permute(options.wave_color, [3,1,2]);
        end
        %Note animation data is already in energy, hence dB is 10 * log10(.),
        %not 20 * log10(.)
        v = (10 * log10(fld{i}(:, ti)) - options.db_range(1)) / (options.db_range(2) - options.db_range(1));
        v(v < 0) = 0;
        v(v > 1) = 1;
        set(h_patch{i}, 'CData', base_cdata{i} .* (1 - v) + v .* max_cdata);
    end
    pause(options.pause_value);
    if ~isempty(options.mp4_out)
        writeVideo(vidObj,getframe(gcf));
    end
end

if ~isempty(options.mp4_out)
    close(vidObj);
end

%Reset figure
for i = 1:numel(h_patch)
    set(h_patch{i}, 'CData', base_cdata{i});
end
end

