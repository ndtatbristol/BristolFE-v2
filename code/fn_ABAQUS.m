function [res, mats] = fn_ABAQUS(mod, matls, steps, fe_options)
default_options.abaqus_jobname = 'test';
default_options.abaqus_exe = 'C:\SIMULIA\Commands\abaqus';
default_options.abaqus_working_folder = 'C:\temp\abaqus';
fe_options = fn_set_default_fields(fe_options, default_options);

if exist(fe_options.abaqus_working_folder) ~= 7
    mkdir(fe_options.abaqus_working_folder);
end

fn_create_abaqus_input_file(mod, matls, steps, [fe_options.abaqus_working_folder, filesep, fe_options.abaqus_jobname, '.inp'], fe_options);
cdir = pwd;
cd(fe_options.abaqus_working_folder)
cmd = [fe_options.abaqus_exe, ' job=', fe_options.abaqus_jobname, ' interactive'];
% cmd = [fe_options.abaqus_exe, ' job=', fe_options.abaqus_jobname];
% fprintf(cmd)
system(cmd);
end

function fn_create_abaqus_input_file(mod, matls, steps, fname, fe_options)
fid = fopen(fname, 'wt');
fprintf(fid, '*HEADING\nFile created in Matlab\n');

%Materials
for i = 1:numel(matls)
    fprintf(fid, ['*MATERIAL, NAME=MATL%i\n'], i);
    fprintf(fid, '*ELASTIC, TYPE=ANISOTROPIC\n');
    utD = fn_expand_upper_triangle_of_D(matls(i).D);
    fprintf(fid, fn_abaqus_fmt_str('%e', 8), utD(1:8));
    fprintf(fid, fn_abaqus_fmt_str('%e', 8), utD(9:16));
    fprintf(fid, fn_abaqus_fmt_str('%e', 8), [utD(17:21), 0, 0, 0]);
    fprintf(fid, '*DENSITY\n');
    fprintf(fid, '%e\n', matls(i).rho);
end

%Nodes
fprintf(fid, '*NODE\n');
fprintf(fid, ['%i, ', fn_abaqus_fmt_str('%e', size(mod.nds, 2))], [[1:size(mod.nds, 1)]', mod.nds]');

%Elements
if ~isfield(mod, 'el_typ_i')
    mod.el_typ_i = repmat({'CPE3'}, size(mod.els,1), 1);
end
un_els = unique(mod.el_typ_i);
for j = 1:numel(un_els)
    fprintf(fid, ['*ELEMENT, TYPE=', un_els{j}, '\n']);
    i = strcmp(mod.el_typ_i, un_els{j});
    fprintf(fid, ['%i, ', fn_abaqus_fmt_str('%i', size(mod.els, 2))], [find(i), mod.els(i,:)]');
end

%Assign elements of each material to a different elset
un_matls = unique(mod.el_mat_i);
for i = 1:numel(un_matls)
    j = mod.el_mat_i == un_matls(i);
    fprintf(fid, '*ELSET, ELSET=ELSET%i\n', i);
    fn_write_abaqus_number_list(fid, find(j), '%i', 16);
end

%Sections
for i = 1:numel(un_matls)
    fprintf(fid, '*SOLID SECTION, ELSET=ELSET%i, MATERIAL=MATL%i, ORIENTATION=ORIENT_GLOBAL\n', i, i);
end

fprintf(fid, '*ORIENTATION, NAME=ORIENT_GLOBAL, SYSTEM=RECTANGULAR\n');
fprintf(fid, '1.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.00000E+00, 0.00000E+00\n');
fprintf(fid, '1, 0.00000E+00\n');

for s = 1:numel(steps)
    for n = 1:size(steps{s}.load.frcs, 1)
        fprintf(fid, '*AMPLITUDE, NAME=S%iN%i, DEFINITION=EQUALLY SPACED, FIXED INTERVAL=%e\n', s, n, steps{s}.load.time(2) - steps{s}.load.time(1));
        fn_write_abaqus_number_list(fid, steps{s}.load.frcs(n,:), '%i', 8);
    end
    fprintf(fid, '*STEP, NLGEOM=NO\n');
    fprintf(fid, '*DYNAMIC, EXPLICIT, DIRECT USER CONTROL\n');
    fprintf(fid, '%e, %e\n', steps{s}.load.time(2) - steps{s}.load.time(1), steps{s}.load.time(end));
    for n = 1:size(steps{s}.load.frc_nds, 1)
        if size(steps{s}.load.frcs, 1) == 1
            fprintf(fid, '*CLOAD, AMPLITUDE=S%iN%i\n', s, 1);
        else
            fprintf(fid, '*CLOAD, AMPLITUDE=S%iN%i\n', s, n);
        end
        fprintf(fid, '%i, %i, %e\n', steps{s}.load.frc_nds(n), steps{s}.load.frc_dfs(n), 1.0);
    end
    fprintf(fid, '*NSET, NSET=NSET_S%i\n', s);
    %TODO - list monitoring nodes here
    %This works to here - file actually executres in ABAQUS and results in
    %OBD file look plausible. Need to add outputs to text file here so they
    %can be read back into Matlab. Acoustic = DoF 8 in Abaqus
    % fprintf(fid, '*OUTPUT, FIELD, NUMBER INTERVAL=%i\n', fe_options.field_output_every_n_frames);
    fprintf(fid, '*OUTPUT, HISTORY\n', fe_options.field_output_every_n_frames);
    fprintf(fid, '*NODE OUTPUT, NSET=NSET_S%i\n', s);
    fprintf(fid, 'U\n');
% *NODE OUTPUT, NSET=NSET_PDW_ALL
% U
    fprintf(fid, '*END STEP\n');
end

% *BOUNDARY
% NSET_PDW_EDGE_XMIN, XSYMM
% NSET_PDW_EDGE_XMAX, XSYMM
% *CLOAD, AMPLITUDE=5_CYCLE_1_MHZ_HANNING
% 10511, 2, 1.00000E+00
% *OUTPUT, FIELD, NUMBER INTERVAL=10
% *NODE OUTPUT, NSET=NSET_PDW_ALL
% U
% *END STEP

fclose(fid);
end

function str = fn_abaqus_fmt_str(typ, no)
str = repmat([typ,', '], 1, no - 1);
str = [str, typ, '\n'];
end

function utD = fn_expand_upper_triangle_of_D(D)
if all(size(D) == [3,3])
    D = [D(:, 1:2), zeros(3, 3), D(:,3)];
    D = [D(1:2, :); zeros(3, 6); D(3, :)];
end
utD = zeros(1, 21);
k = 1;
for i = 1:6
    for j = 1:i
        utD(k) = D(j,i);
        k = k + 1;
    end
end

end

function fn_write_abaqus_number_list(fid, data, type_str, no_per_line)
left_overs = rem(numel(data), no_per_line);
k = numel(data) - left_overs;
fprintf(fid, fn_abaqus_fmt_str(type_str, no_per_line), data(1:k));
if left_overs
    fprintf(fid, fn_abaqus_fmt_str(type_str, left_overs), data(k+1:end));
end
end