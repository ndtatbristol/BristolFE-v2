function varargout = fn_pogoFE(mod, matls, steps, fe_options)
default_options.pogo_path = 'C:\Program Files\Pogo\windows';
default_options.pogo_matlab_path = 'C:\Program Files\Pogo\matlab';
default_options.pogo_verbosity = -1;
default_options.pogo_compression = 0;
default_options.max_diff_absorbing_levels = 50;
default_options.solver_precision = 'single';
fe_options = fn_set_default_fields(fe_options, default_options);
addpath(genpath(fe_options.pogo_matlab_path));

if fe_options.pogo_verbosity < 0
    verb_flag = ' -o >NUL';
else
    verb_flag = sprintf(' -o --setVerbosity %i', fe_options.pogo_verbosity);
end
if fe_options.pogo_compression == 0
    solver_flags = ' --setCompressOff';
else
    solver_flags = '';
end
if isinf(fe_options.field_output_every_n_frames)
    solver_flags = [solver_flags, ' --setFieldSaveOff'];
end

pogo_model = fn_convert_to_pogo_model(mod, matls, steps, fe_options);

%Save model
fname = 'pogoPW';
dummy_fname = [fname, '.txt'];
warning('off', 'all');
savePogoInp(fname, pogo_model);
warning('on', 'all');

if pogo_model.nDims == 3
    pogo_solver = 'pogoSolve3D';
else
    switch fe_options.solver_precision
        case 'single'
            pogo_solver = 'pogoSolve';
        case 'double'
            pogo_solver = 'pogoSolve64';
    end
end
%Blocking
system(['"', fe_options.pogo_path, filesep, 'pogoBlock" ',fname, '.pogo-inp', verb_flag]);
%Solving
system(['"', fe_options.pogo_path, filesep, pogo_solver, '" ',fname, solver_flags, verb_flag]);

for s = 1:numel(steps)
    if numel(steps) > 1
        fn = sprintf([fname,'-%i'], s);
    else
        fn = fname;
    end
    h = loadPogoHist(fn);
    res{s}.dsps = h.sets(1).MeasureSet1.histTraces';
    delete([fn, '.*']);
end

varargout{1} = res;

if nargout > 1
    mats = [];
    if isempty(fe_options.gl_mat_nds)
        %default is all nodes
        fe_options.gl_mat_nds = 1:size(mod.nds, 1);
    end
    no_nds = size(pogo_model.nodePos, 2);
    total_df = no_nds * pogo_model.nDofPerNode;
    mats.M = spalloc(total_df, total_df, no_nds * pogo_model.nDofPerNode);
    mats.C = spalloc(total_df, total_df, no_nds * pogo_model.nDofPerNode);
    mats.K = spalloc(total_df, total_df, no_nds * pogo_model.nDofPerNode);
    df = 1:pogo_model.nDofPerNode;
    gl_nd_i = ones(pogo_model.nDofPerNode, 1) * (1:no_nds);
    gl_nd_i = gl_nd_i(:);
    df = 1:pogo_model.nDofPerNode;
    gl_df_i = df' * ones(1, no_nds);
    gl_df_i = gl_df_i(:);
    mats.gl_lookup = fn_create_fast_lookup(gl_nd_i, gl_df_i, 0, 0);
    %get stiffness matrix
    steps = steps(1);
    steps{1}.load.time = [0, 1];
    steps{1}.load.frcs = [0, 0];
    pogo_model = fn_convert_to_pogo_model(mod, matls, steps, fe_options);
    savePogoInp(fname, pogo_model);
    % % keyboard
    system(['"', fe_options.pogo_path, filesep, 'pogoBlock" ',fname, '.pogo-inp', verb_flag]);
    % system(['"', fe_options.pogo_path, filesep, pogo_solver, '" ',fname, solver_flags, verb_flag]);
    for i = 1:numel(fe_options.gl_mat_nds)
        system(['"', fe_options.pogo_path, filesep, pogo_solver, '" ',fname, sprintf(' --printStiff %i ', fe_options.gl_mat_nds(i)), '-o --setFieldSaveOff --setVerbosity 0 > ', dummy_fname]);
        [K, C, M, gl_i] = fn_scrape_global_matrix_values(dummy_fname, mats.gl_lookup);
        %these need to be put (actually added) into actual global matrices
        %mats.M etc
        mats.K(gl_i, gl_i) = mats.K(gl_i, gl_i) + K;
        mats.M(gl_i, gl_i) = mats.M(gl_i, gl_i) + M;
        mats.C(gl_i, gl_i) = mats.C(gl_i, gl_i) + C;
        disp(i)
    end
    varargout{2} = mats;
end
end


function model = fn_convert_to_pogo_model(mod, matls, steps, fe_options)
model.nDims = size(mod.nds, 2);
model.nDofPerNode = model.nDims;                    %This is not general enough - but how do you get DoF in my code?
%model.runName = "";                                 %not important (unused at present)
for n = 1:numel(steps)
    dt(n) = steps{n}.load.time(2) - steps{n}.load.time(1);
    nt(n) = numel(steps{n}.load.time);
end
if numel(unique(dt)) > 1 || numel(unique(nt)) > 1
    warning('Not sure how Pogo interprets different time steps per shot')
end
model.dt = dt(1);                                   %time step
model.nt = nt(1);                                   %number of increments
model.nodePos = mod.nds';                            %nodal locations; size nDims x nNodes
model.elNodes = mod.els';                           %nodes for each element; size nNodesPerElMax x nEls
%In bristol FE elements are identified with material and element type is
%specified in material definition. There is no concept of material orientation in
%BristolFE (different orientations require different materials)
model.elTypeRefs = mod.el_mat_i';                   %- which of the element types each element refers to, length nEls
model.matTypeRefs = mod.el_mat_i';                  % - which of the material types each element refers to, length nEls
model.orientRefs = zeros(size(model.matTypeRefs));                   %- which of the orientations each element refers to, length nEls

for n = 1:numel(matls)
    model.elTypes{n}.name = matls(n).el_typ;       % - element name (matching Abaqus library typically)
    model.elTypes{n}.paramsType = 0;                %This will need changing for damping - parameters associated with the element type - usually just 0
    if isfield(matls(n), 'alpha')
        model.matTypes{n} = fn_pogo_matl(matls(n).D, matls(n).rho, matls(n).alpha, 0);
    else
        model.matTypes{n} = fn_pogo_matl(matls(n).D, matls(n).rho, 0, 0);
    end
    % model.matTypes{n}.parent = 0;                   %- what is the parent material (0 if no parent) - used in absorbing boundaries
    % model.matTypes{n}.paramsType = 2;               %2 is anisoptropic with 21 params
    % model.matTypes{n}.paramValues = [matls(n).D(find(triu(ones(size(matls(n).D)))))', matls(n).rho];%- the parameters mentioned above
end

%Deal with absorbing materials
%First work out levels of discrete absorption levels to use in Pogo (each
%one requires a separate material definition)
if isfield(mod, 'el_abs_i')
    unique_abs_levels = unique(mod.el_abs_i);
    if fe_options.max_diff_absorbing_levels < numel(unique_abs_levels)
        %pogo_abs_levels = interp1(unique_abs_levels, linspace(0, 1, numel(unique_abs_levels)), linspace(0, 1, fe_options.max_diff_absorbing_levels));
        pogo_abs_levels = linspace(0, 1, fe_options.max_diff_absorbing_levels + 1); %plus one because first one is zero
        % pogo_abs_levels = pogo_abs_levels(2:end);
    else
        pogo_abs_levels = unique_abs_levels;
    end
    mod.el_abs_i_pogo = interp1(pogo_abs_levels, 1:numel(pogo_abs_levels) ,mod.el_abs_i, 'nearest');
    %convert Bristol FE's continuous abs levels to Pogo indices for discrete levels
    new_mat_ind = numel(matls) + 1; %indices of new materials to be added with damping
    no_matls_initially = numel(matls);
    for n = 1:no_matls_initially
        %max_diff_absorbing_levels
        if all(mod.el_abs_i(mod.el_mat_i == n) == 0)
            %no absorbing elements for this material, move onto next one
            continue
        end
        for m = 2:numel(pogo_abs_levels) %start at 2 because 1 is zero absorption case
            %Figure out the absorbing material parameters
            if isfield(matls(n), 'alpha')
                alpha = matls(n).alpha;
            else
                alpha = 0.0;
                alpha = alpha + pogo_abs_levels(m) ^ fe_options.damping_power_law *  fe_options.max_damping;
                D = matls(n).D * exp(log(fe_options.max_stiffness_reduction) * pogo_abs_levels(m) ^ (fe_options.damping_power_law + 1));
                %Add new material to Pogo model
                model.matTypes{new_mat_ind} = fn_pogo_matl(D, matls(n).rho, alpha, n);
                %Assign all relevant elements in pogo model to this material
                model.matTypeRefs((mod.el_mat_i == n) & (mod.el_abs_i_pogo == m)) = new_mat_ind;
                new_mat_ind = new_mat_ind + 1;
            end
        end
    end
end
%Following probably not needed - will define forcing directly

%Following will need to be pulled from displacement control?
% model.fixNodes - nodes with DOFs to be fixed
% model.fixDof - DOF corresponding to the nodes above. Nodes can appear in
% fixNodes multiple times to specify different DOF.

%will need to switch for disp excitation as well at some point
for n = 1:numel(steps)

    model.shots{n}.ntSig = numel(steps{n}.load.time); %- number of time points for all the signals (commonly just set to model.nt)
    model.shots{n}.dtSig = steps{n}.load.time(2) - steps{n}.load.time(1); %- time step for the signals (commonly just set to model.dt)
    ones_for_time_pts = ones(1, model.shots{n}.ntSig);
    if size(steps{n}.load.frcs, 1) == 1
        model.shots{n}.sigs{1}.isDofGroup = 0;                      % - 0 if not, 1 if dofGroups are referred to in dofSpec (NB only from v1.08)
        model.shots{n}.sigs{1}.nodeSpec = steps{n}.load.frc_nds(:)';    % - nodes the signal is applied to.
        model.shots{n}.sigs{1}.dofSpec = steps{n}.load.frc_dfs(:)';     %- DOF to apply the signal to (matches nodeSpec) in (range: 1 to model.nDofPerNode)
        model.shots{n}.sigs{1}.sigType = 0;                         %- 0 force, 1 displacement or (unused at present) 2 velocity
        model.shots{n}.sigs{1}.sigAmps = ones(size(model.shots{n}.sigs{1}.dofSpec));                         %- amplitudes the signals are multiplied by for each dof specified
        model.shots{n}.sigs{1}.sig = steps{n}.load.frcs(1, :);
    else
        for m = 1:size(steps{n}.load.frcs, 1)
            model.shots{n}.sigs{m}.isDofGroup = 0;                      % - 0 if not, 1 if dofGroups are referred to in dofSpec (NB only from v1.08)
            model.shots{n}.sigs{m}.nodeSpec = steps{n}.load.frc_nds(m);    % - nodes the signal is applied to.
            model.shots{n}.sigs{m}.dofSpec = steps{n}.load.frc_dfs(m);     %- DOF to apply the signal to (matches nodeSpec) in (range: 1 to model.nDofPerNode)
            model.shots{n}.sigs{m}.sigType = 0;                         %- 0 force, 1 displacement or (unused at present) 2 velocity
            model.shots{n}.sigs{m}.sigAmps = 1;                         %- amplitudes the signals are multiplied by for each dof specified
            model.shots{n}.sigs{m}.sig = steps{n}.load.frcs(m, :);
        end
    end

    if numel(steps) > 1
        for n = 2:numel(steps)
            if ~isequal(steps{1}.mon, steps{n}.mon)
                warning('Using monitoring information for step 1 only')
            end
        end
    end
    model.measFreq = 1;%- number of time increments between when history measurements are taken
    model.measStart = 1;%- starting increment (1 indexed)
    % model.measSets{n} - structure for each of the measured node sets
    model.measSets{1}.name = 'MeasureSet1'; %- string name for the set
    model.measSets{1}.isDofGroup = 0; %- do we refer to DoF groups or do we use nodal values (0 or 1)
    %if 0:
    model.measSets{1}.measNodes =steps{1}.mon.nds;  %- nodes to take history measurements from
    model.measSets{1}.measDof = steps{1}.mon.dfs; %- degrees of freedom to take history measurements from
    %if 1:
    % model.measSets{n}.dofGroup - which Dof groups
    %
    if fe_options.field_output_every_n_frames < inf
        warning('Pogo field output not supported')
        % model.fieldStoreIncs = fe_options.field_output_every_n_frames; %- which increments to output the field at. Omit or empty array if none.
    end
    model.fieldStoreIncs = [];
    %
    % model.metadata - struct containing metadata, stored as
    %   model.metadata.temperature = 20
    % or
    %   model.metadata.owner = 'Boris Johnson'
    % Note that both text or numerical values should work OK - both will be stored
    % as text then interpreted as numbers when loaded back in.
    %
    % model.mpiDivNum - what is the MPI number associated with this model (i.e. how
    %   will the other models refer to this one)
    % model.mpi{X}.outNodes and model.mpi{Y}.inNodes - which nodal values should I pass
    %   out to process X, and where should I place (overwrite) the nodal values
    %   from process Y
end

end

function pogo_matl = fn_pogo_matl(D, rho, alpha, parent_matl)
pogo_matl.parent = parent_matl;                   %- what is the parent material (0 if no parent) - used in absorbing boundaries
pogo_matl.paramsType = 2;               %2 is anisoptropic with 21 params
pogo_matl.paramValues = [D(find(triu(ones(size(D)))))', rho];%- the parameters mentioned above
if alpha > 0
    pogo_matl.paramValues = [pogo_matl.paramValues, alpha];
end
end
