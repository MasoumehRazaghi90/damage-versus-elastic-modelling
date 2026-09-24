
clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8; %transparency
markerSize=40; %For plotted points
markerSize2=10; %For nodes on patches
lineWidth1=1; %For meshes
lineWidth2=2; %For boundary edges
cMap=spectral(250); %colormap

%%

% Geometry parameters
length_cut = 29; % Uneven and larger than 2
length_spacing = 3;
num_repeat_x = 15;
num_repeat_y = 2;

num_nodes_x = 1+(length_spacing+1)*(num_repeat_x+1);
num_nodes_y = num_repeat_y*length_cut + length_spacing*(num_repeat_y+1);

plateEl = [num_nodes_x-1,num_nodes_y-1];
plateDim = plateEl;

layerThickness = 3;
elementType = 'hex8';%'quad4';
numRefine = 6;


appliedStretch = 1.40;       % 40% engineering strain


%Material parameter set

%% Material parameters from the second code


k_factor = 1000;

% 1-term Ogden matrix
c1_ini = 2.7271875051643106e+00;
m1_ini = 6.9427891844363305e+00;
k_ini  = c1_ini*k_factor;

% GOH fibers
k1_ini    = 2.1009948723657132e+01;
k2_ini    = 6.9399954857954005e-01;
kappa_ini = 1.4455130336209238e-01;
gamma     = 41;

% Very small GOH ground-matrix contribution
cGOH_ini = 1e-4;
kGOH_ini = cGOH_ini*k_factor;

% Ground-matrix damage
mu_matrix_ini    = 1.0302688108433917e+00;
sigma_matrix_ini = 3.7822643808779288e-01;
Dmax_matrix_ini  = 1;

% GOH damage
mu_GOH_ini       = 5.5063703735016223e-01;
sigma_GOH_ini    = 4.7249662329555610e-02;
Dmax_GOH_ini     = 1;



LangerAngle_perp = 0/180*pi; % 0 deg: e1 along X, perpendicular to Y loading
LangerAngle_para = 90/180*pi; % 90 deg: e1 along Y, parallel to Y loading
angleSelection=1;%1 perpendicular to LL and 2 parallel
LangerAngles = [LangerAngle_perp, LangerAngle_para];


% FEA control settings
numTimeSteps=80; %Number of time steps desired
max_refs=800; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=15; %Optimum number of iterations
max_retries=50; %Maximum number of retires
dtmin=(1/numTimeSteps)/500; %Minimum time step size
dtmax=(1/numTimeSteps); %Maximum time step size
runMode='internal';
runFEBio = true; % true: run FEBio and directly post-process log outputs
% false: post-process existing FEBio log outputs

% Path names
% This MATLAB code is stored in:
defaultFolder = fileparts(mfilename('fullpath'));

% Zimmer Damage2 results are stored in:
savePath = 'C:\Users\Sahar\Downloads\chapter5\damage\damage2\data\zimmermesher';

if ~exist(savePath,'dir')
    mkdir(savePath);
end

fprintf('Zimmer Damage2 results will be saved in:\n%s\n',savePath);


% Defining file names
febioFebFileNamePart='tempModel';
%%

[F,V] = quadPlate(plateDim,plateEl);

%%
ind = [];
indSplit = [];

s = ceil((length_cut+length_spacing)/2);
i_split1 = [];
i_split2 = [];

p = length_spacing+2:length_cut-2+length_spacing+1;
i1 = [];
for i_s = 1:1:num_repeat_y+1
    ii = p + ((i_s-1)*(length_spacing+length_cut));
    i_add = ii(ii<num_nodes_y);
    i1 = [i1 i_add];

    if ~isempty(i_add)
        i_split_add = [i_add(1)-1 i_add(end)+1];
        i_split1 = [i_split1 i_split_add(i_split_add>0 & i_split_add<num_nodes_y)];
    end
end

p = length_spacing+2-s:length_cut-2+length_spacing+1-s;
i2 = [];
for i_s = 1:1:num_repeat_y+1
    ii = p + ((i_s-1)*(length_spacing+length_cut));
    i_add = ii(ii<=num_nodes_y & ii>0);
    if ~isempty(i_add)
        i_split_add = [i_add(1)-1 i_add(end)+1];
        i_split2 = [i_split2 i_split_add(i_split_add>0 & i_split_add<num_nodes_y)];
    end

    i2 = [i2 i_add];
end

% ------------------------------------------------------------
% Save the Y coordinates of the slit-tip nodes that define
% ------------------------------------------------------------

yBCLevels = sort(unique(V(i_split2,2)));

if numel(yBCLevels) ~= 4
    error('Expected four BC-defining slit-tip Y levels, but found %d.', ...
        numel(yBCLevels));
end

fprintf('\nBC levels obtained automatically from slit-tip nodes:\n');
disp(yBCLevels(:)');

c = 1;
for j = (length_spacing+2):length_spacing+1:num_nodes_x-1

    if iseven(c)
        indAdd = sub2ind([num_nodes_y,num_nodes_x],i2,j*ones(size(i2)));
        indSplitAdd = sub2ind([num_nodes_y,num_nodes_x],i_split2,j*ones(size(i_split2)));
    else
        indAdd = sub2ind([num_nodes_y,num_nodes_x],i1,j*ones(size(i1)));
        indSplitAdd = sub2ind([num_nodes_y,num_nodes_x],i_split1,j*ones(size(i_split1)));
    end
    ind = [ind; indAdd(:)];
    indSplit = [indSplit; indSplitAdd(:)];
    c = c+1;
end

%%

cFigure; hold on;
title('grid');
hp=gpatch(F,V,'bw','k',0.9);
plotV(V(ind,:),'r.','MarkerSize',markerSize)
plotV(V(indSplit,:),'g.','MarkerSize',markerSize)
axisGeom; camlight headlight;
gdrawnow;

%%

for q = 1:1:numRefine
    indFacesRemove = [];
    C = patchConnectivity(F,V,{'vf'});

    for i_vs = 1:length(indSplit)
        i_v = indSplit(i_vs);
        V = [V; V(i_v,:)]; % Append copy of point
        indFacesNow = C.vertex.face(i_v,:);
        indFacesNow = indFacesNow(indFacesNow>0);
        indFacesRemove = [indFacesRemove indFacesNow];
        for i_f = indFacesNow
            f = F(i_f,:);
            % Reorder so split point is first
            iStart = find(f == i_v); % Index of split point
            if iStart>1
                f = [f(iStart:end) f(1:iStart-1)]; % reorder
            end
            e1 = [f(1),f(2)];
            e2 = [f(2),f(3)];
            e3 = [f(3),f(4)];
            e4 = [f(4),f(1)];
            if any(ismember(e1,ind))
                ind=[ind; length(V)+1];
            elseif any(ismember(e4,ind))
                ind=[ind; length(V)+3];
            end
            F_add = [f(1) length(V)+1 length(V)+2 length(V)+3;...
                length(V)+1 f(2) f(3) length(V)+2;...
                length(V)+3 length(V)+2 f(3) f(4)];
            F = [F; F_add];
            V = [V; mean(V(e1,:),1); mean(V(f,:),1); mean(V(e4,:),1)];
        end
    end
    L = true(size(F,1),1);
    L(indFacesRemove) = 0;
    F = F(L,:);

    [F,V,~,indFix] = mergeVertices(F,V);
    ind = indFix(ind);
    indSplit = indFix(indSplit);

    [F,V,indFix] = patchCleanUnused(F,V);
    ind = indFix(ind);
    indSplit = indFix(indSplit);

    ind = unique(ind);
end

%% Smooth refined mesh like the GIBBON expanding-lattice demo
Eb = patchBoundary(F,V);

smoothOpt.Method = 'HC';
smoothOpt.n = 250;
smoothOpt.RigidConstraints = unique(Eb(:));

V = patchSmooth(F,V,[],smoothOpt);



cFigure; hold on;
title('grid');
hp=gpatch(F,V,'bw','k',0.9);
plotV(V(ind,:),'r.','MarkerSize',markerSize)
plotV(V(indSplit,:),'g.','MarkerSize',markerSize)
axisGeom; camlight headlight;
gdrawnow;


%%

C = patchConnectivity(F,V,{'vf'});
VF = patchCentre(F,V);

indNew = [];
for i_v = ind(:)'
    V = [V; V(i_v,:)]; % Append copy of point
    indFacesNow = C.vertex.face(i_v,:);
    indFacesNow = indFacesNow(indFacesNow>0);
    for i_f = indFacesNow(:)'
        if VF(i_f,1)>V(i_v,1) % On the right side
            f = F(i_f,:);
            f(f==i_v)=length(V);
            F(i_f,:)=f;
            indNew = [indNew; length(V)];
        end
    end
end


%% Rotate complete Zimmer geometry by +90 degrees about Z
%
% old X -> new Y
% old Y -> new -X
%
% This is a rigid rotation: mesh connectivity F is unchanged.

V_old = V;

V(:,1) = -V_old(:,2);
V(:,2) =  V_old(:,1);

% Actual reference length in the Y loading direction
sampleWidth = max(V(:,2)) - min(V(:,2));

% Prescribed displacement required for exactly 40% engineering strain
displacementMagnitude = ...
    sampleWidth*(appliedStretch-1);

fprintf('Reference Y loading length = %.6g mm\n',sampleWidth);
fprintf('Prescribed Y displacement = %.6g mm\n',displacementMagnitude);
fprintf('Target stretch = %.4f\n',appliedStretch);

% The BC interval coordinates were originally measured in old Y.
% After rotation they are coordinates in new X.
xBCLevels = sort(-yBCLevels);

fprintf('\nZimmer geometry rotated +90 degrees.\n');
fprintf('Rotated BC X levels:\n');
disp(xBCLevels(:)');
%%

cFigure; hold on;
title('grid');
hp=gpatch(F,V,'bw','k',0.9);
plotV(V(ind,:),'r.','MarkerSize',markerSize)
plotV(V(indSplit,:),'g.','MarkerSize',markerSize)
plotV(V(indNew,:),'y.','MarkerSize',markerSize/2)
axisGeom; camlight headlight;
gdrawnow;

%%

%% Define element type and boundary conditions
searchtol=1e-3;

switch elementType
    case 'quad4'
        E = F;
        bcFixList = find(V(:,1)<(min(V(:,1))+searchtol));
        bcPrescribeList = find(V(:,1)>(max(V(:,1))-searchtol));

    case 'hex8'

        [E,V] = patchThick(F,V,1,layerThickness,ceil(layerThickness));
        F = element2patch(E);

        % Check initial hexahedral elements
        [elementVolume0,elementValid0] = hexVol(E,V);

        fprintf('Invalid initial elements: %d\n',sum(~elementValid0));
        fprintf('Minimum initial element volume: %.6g\n',min(elementVolume0));

        if any(~elementValid0)
            error('The initial mesh contains inverted or invalid hex elements.');
        end


        % No BC dimensions are hard-coded here.
        % ------------------------------------------------------------

        % ------------------------------------------------------------
        % LOCAL ZIMMER BC PATCHES AFTER 90-DEGREE ROTATION
        %
        % old RIGHT side -> new TOP
        % old LEFT side  -> new BOTTOM
        %
        % xBCLevels are obtained automatically from the original
        % slit-tip geometry.
        % ------------------------------------------------------------

        yMin = min(V(:,2));
        yMax = max(V(:,2));

        % Two X intervals corresponding to the original BC Y intervals
        bcXMask = ...
            (V(:,1) >= xBCLevels(1)-searchtol & ...
            V(:,1) <= xBCLevels(2)+searchtol) | ...
            ...
            (V(:,1) >= xBCLevels(3)-searchtol & ...
            V(:,1) <= xBCLevels(4)+searchtol);

        % TOP: prescribed/loading patches
        bcPrescribeList = find( ...
            abs(V(:,2)-yMax) <= searchtol & ...
            bcXMask);

        % BOTTOM: support patches
        bcSupportList = find( ...
            abs(V(:,2)-yMin) <= searchtol & ...
            bcXMask);

        % ------------------------------------------------------------
        % One reference node on the bottom-left support patch
        % ------------------------------------------------------------

        bcSupportLeftList = bcSupportList( ...
            V(bcSupportList,1) >= xBCLevels(1)-searchtol & ...
            V(bcSupportList,1) <= xBCLevels(2)+searchtol);

        if isempty(bcSupportLeftList)
            error('No nodes found in the bottom-left BC patch.');
        end

        zMid = 0.5*(min(V(:,3))+max(V(:,3)));
        xLeftMid = mean(xBCLevels(1:2));

        [~,indPin] = min( ...
            (V(bcSupportLeftList,1)-xLeftMid).^2 + ...
            (V(bcSupportLeftList,3)-zMid).^2);

        bcPinList = bcSupportLeftList(indPin);

        % ------------------------------------------------------------
        % Verification
        % ------------------------------------------------------------

        fprintf('\n--- Rotated Zimmer BC verification ---\n');

        fprintf('Specimen Y range       : %.4f to %.4f mm\n', ...
            yMin,yMax);

        fprintf('Left BC X interval     : %.4f to %.4f mm\n', ...
            xBCLevels(1),xBCLevels(2));

        fprintf('Right BC X interval    : %.4f to %.4f mm\n', ...
            xBCLevels(3),xBCLevels(4));

        fprintf('Top prescribed nodes   : %d\n', ...
            numel(bcPrescribeList));

        fprintf('Bottom support nodes   : %d\n', ...
            numel(bcSupportList));

        fprintf('Reference nodes        : %d\n', ...
            numel(bcPinList));

        fprintf('Reference node xyz     : %.4f %.4f %.4f mm\n', ...
            V(bcPinList,1), ...
            V(bcPinList,2), ...
            V(bcPinList,3));

        if isempty(bcPrescribeList)
            error('No top prescribed nodes were selected.');
        end

        if isempty(bcSupportList)
            error('No bottom support nodes were selected.');
        end

        if numel(bcPinList) ~= 1
            error('Exactly one reference node is required.');
        end

        if ~isempty(intersect(bcSupportList,bcPrescribeList))
            error('Top and bottom BC node sets overlap.');
        end


end

fprintf('\nZimmer Damage2 boundary check:\n');
fprintf('  Top prescribed nodes = %d\n',numel(bcPrescribeList));
fprintf('  Bottom support nodes = %d\n',numel(bcSupportList));
fprintf('  Reference nodes      = %d\n\n',numel(bcPinList));
%%

% V(:,3)=V(:,3)+1e-2*randn(size(V(:,3)));

%%
% Visualizing boundary conditions.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(F,V,'w','none',0.5);


hl(1)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize);
hl(2)=plotV(V(bcSupportList,:),'g.','MarkerSize',markerSize);
hl(3)=plotV(V(bcPinList,:),'b.','MarkerSize',markerSize*1.5);

legend(hl,{ ...
    'Top patches: prescribed u_y, u_z=0', ...
    'Bottom patches: u_y=0, u_z=0', ...
    'Reference node: u_x=0'}, ...
    'Location','best');

axisGeom(gca,fontSize);
view(2);
camlight headlight;
drawnow;

%% DEFINE FIBRE DIRECTIONS
nElem=length(E);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%%
setenv('GIBBON_cFigure_close','0');

results = struct( ...
    'lambda',[], ...
    'force',[], ...
    'energy_mean',[], ...
    'energy_sum',[], ...
    'energy_damage_mean',[], ...
    'energy_intact_mean',[], ...
    'damage',[], ...
    'damage_max',[], ...
    'damage_fraction',[], ...   % volume fraction with aggregate mixture D >= 0.99
    'p90_vm',[], ...            % volume-weighted p90 von Mises stress [MPa]
    'nu_true',[], ...           % apparent/global true Poisson ratio
    'area_ratio',[]);           % projected bounding-box area ratio A/A0

results(2) = results;

angles_deg = [0 90];  % Loop from 0 to 90 degrees

for ia = 1:numel(angles_deg)
    ang = angles_deg(ia);

    if ang == 0
        angLabel = 'Perpendicular';
    elseif ang == 90
        angLabel = 'Parallel';
    else
        error('Only 0 and 90 degrees are supported in this setup.');
    end

    % Material directions: SAME convention as AUXETIC
    R = euler2DCM([0,0,-ang*pi/180]);

    e1 = (R*[1 0 0]')';  % Rotated fibre/Langer direction

    e1_dir = repmat(e1,nElem,1);
    e2_dir = repmat([0 0 1],nElem,1);  % Same as auxetic
    e3_dir = cross(e1_dir,e2_dir,2);    % Same as auxetic
    % Plot material directions
    VE = patchCentre(E, V);  % Make sure VE is defined

    cFigure; hold on;
    title(sprintf('Fiber direction = %d°', ang), 'FontSize', 20);
    gpatch(F, V, 'kw', 'none', 0.25);
    quiverVec(VE, e1_dir, 1, 'r');
    quiverVec(VE, e2_dir, 0.5, 'g');
    quiverVec(VE, e3_dir, 0.5, 'b');
    axisGeom(gca, fontSize);
    camlight headlight;
    drawnow;

    %%
    % Visualizing material directions

    [VE]=patchCentre(E,V);

    hf=cFigure; hold on;
    title(sprintf('Material directions  |  Fiber Direction = %s', angLabel), 'FontSize', fontSize);

    gpatch(F,V,'kw','none',0.25);
    hf(1)=quiverVec(VE,e1_dir,1,'r');
    hf(2)=quiverVec(VE,e2_dir,0.5,'g');
    hf(3)=quiverVec(VE,e3_dir,0.5,'b');
    legend(hf,{'e1-direction','e2-direction','e3-direction'});
    axisGeom(gca,fontSize);
    camlight headlight;


    legend(hf,{'e1-direction (fiber)','e2-direction','e3-direction'});
    axisGeom(gca,fontSize);
    camlight headlight;
    drawnow;


    %% Defining the FEBio input structure
    % See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
    % manual.

    %Get a template with default settings
    [febio_spec]=febioStructTemplate;

    %febio_spec version
    febio_spec.ATTR.version='4.0';

    %Module section
    febio_spec.Module.ATTR.type='solid';

    %% Control settings for loading and unloading steps

    stepStruct.Control.analysis='STATIC';
    stepStruct.Control.time_steps=numTimeSteps;
    stepStruct.Control.step_size=1/numTimeSteps;

    stepStruct.Control.solver.max_refs=max_refs;
    stepStruct.Control.solver.qn_method.max_ups=max_ups;

    stepStruct.Control.time_stepper.dtmin=dtmin;
    stepStruct.Control.time_stepper.dtmax=dtmax;
    stepStruct.Control.time_stepper.max_retries=max_retries;
    stepStruct.Control.time_stepper.opt_iter=opt_iter;

    % Complete missing settings using the FEBio template
    stepStruct.Control = structComplete( stepStruct.Control,febio_spec.Control,1);

    % Remove the original global Control section
    febio_spec = rmfield(febio_spec,'Control');

    % Step 1: loading
    febio_spec.Step.step{1}.ATTR.id=1;
    febio_spec.Step.step{1}.ATTR.name='Loading';
    febio_spec.Step.step{1}.Control=stepStruct.Control;

    % Step 2: unloading
    febio_spec.Step.step{2}.ATTR.id=2;
    febio_spec.Step.step{2}.ATTR.name='Unloading';
    febio_spec.Step.step{2}.Control=stepStruct.Control;

    %% Material section: damaged Ogden matrix + two damaged HGO families

    materialName1 = 'Material1';

    febio_spec.Material.material{1}.ATTR.id   = 1;
    febio_spec.Material.material{1}.ATTR.name = materialName1;
    febio_spec.Material.material{1}.ATTR.type = 'solid mixture';

    %% Constituent 1: damaged Ogden ground matrix

    febio_spec.Material.material{1}.solid{1}.ATTR.type = 'uncoupled elastic damage';
    febio_spec.Material.material{1}.solid{1}.k = k_ini;

    % Elastic matrix
    febio_spec.Material.material{1}.solid{1}.elastic.ATTR.type = 'Ogden';
    febio_spec.Material.material{1}.solid{1}.elastic.density = 1;
    febio_spec.Material.material{1}.solid{1}.elastic.c1 = c1_ini;
    febio_spec.Material.material{1}.solid{1}.elastic.m1 = m1_ini;
    febio_spec.Material.material{1}.solid{1}.elastic.k  = k_ini;

    % Matrix damage
    febio_spec.Material.material{1}.solid{1}.damage.ATTR.type = 'CDF log-normal';
    febio_spec.Material.material{1}.solid{1}.damage.mu = mu_matrix_ini;
    febio_spec.Material.material{1}.solid{1}.damage.sigma = sigma_matrix_ini;
    febio_spec.Material.material{1}.solid{1}.damage.Dmax = Dmax_matrix_ini;
    febio_spec.Material.material{1}.solid{1}.criterion.ATTR.type = 'DC max normal Lagrange strain';


    %% Constituent 2: damaged HGO fibre family at +gamma

    febio_spec.Material.material{1}.solid{2}.ATTR.type = 'elastic damage';
    febio_spec.Material.material{1}.solid{2}.elastic.ATTR.type = 'HGO unconstrained';
    febio_spec.Material.material{1}.solid{2}.elastic.c = cGOH_ini/2;
    febio_spec.Material.material{1}.solid{2}.elastic.k1 = k1_ini;
    febio_spec.Material.material{1}.solid{2}.elastic.k2 = k2_ini;
    febio_spec.Material.material{1}.solid{2}.elastic.kappa = kappa_ini;
    febio_spec.Material.material{1}.solid{2}.elastic.gamma = +gamma;
    febio_spec.Material.material{1}.solid{2}.elastic.k = kGOH_ini/2;

    % Fibre damage
    febio_spec.Material.material{1}.solid{2}.damage.ATTR.type = 'CDF log-normal';
    febio_spec.Material.material{1}.solid{2}.damage.mu = mu_GOH_ini;
    febio_spec.Material.material{1}.solid{2}.damage.sigma = sigma_GOH_ini;
    febio_spec.Material.material{1}.solid{2}.damage.Dmax = Dmax_GOH_ini;
    febio_spec.Material.material{1}.solid{2}.criterion.ATTR.type = 'DC max normal Lagrange strain';


    %% Constituent 3: damaged HGO fibre family at -gamma

    febio_spec.Material.material{1}.solid{3}.ATTR.type = 'elastic damage';
    febio_spec.Material.material{1}.solid{3}.elastic.ATTR.type = 'HGO unconstrained';
    febio_spec.Material.material{1}.solid{3}.elastic.c = cGOH_ini/2;
    febio_spec.Material.material{1}.solid{3}.elastic.k1 =k1_ini;
    febio_spec.Material.material{1}.solid{3}.elastic.k2 = k2_ini;
    febio_spec.Material.material{1}.solid{3}.elastic.kappa = kappa_ini;
    febio_spec.Material.material{1}.solid{3}.elastic.gamma =-gamma;
    febio_spec.Material.material{1}.solid{3}.elastic.k = kGOH_ini/2;

    % Same fibre-damage parameters as the +gamma family
    febio_spec.Material.material{1}.solid{3}.damage.ATTR.type = 'CDF log-normal';
    febio_spec.Material.material{1}.solid{3}.damage.mu = mu_GOH_ini;
    febio_spec.Material.material{1}.solid{3}.damage.sigma = sigma_GOH_ini;
    febio_spec.Material.material{1}.solid{3}.damage.Dmax = Dmax_GOH_ini;
    febio_spec.Material.material{1}.solid{3}.criterion.ATTR.type = 'DC max normal Lagrange strain';


    % Mesh section
    % -> Nodes

    %%Area of interest
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

    % -> Elements
    partName1='Part1';
    febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type=elementType; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

    % -> NodeSets
    nodeSetName1='bcPrescribeList1';
    nodeSetName2='bcSupportList2';
    nodeSetName3='bcPinList3';

    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcPrescribeList);

    febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
    febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList);

    febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
    febio_spec.Mesh.NodeSet{3}.VAL=mrow(bcPinList);

    %MeshDomains section
    switch elementType
        case 'quad4'
            febio_spec.MeshDomains.ShellDomain{1}.ATTR.name=partName1;
            febio_spec.MeshDomains.ShellDomain{1}.ATTR.mat=materialName1;
            febio_spec.MeshDomains.ShellDomain{1}.shell_thickness=layerThickness;
        case 'hex8'
            febio_spec.MeshDomains.SolidDomain{1}.ATTR.name=partName1;
            febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat=materialName1;
    end

    febio_spec.MeshData.ElementData{1}.ATTR.elem_set = partName1;
    febio_spec.MeshData.ElementData{1}.ATTR.type = 'mat_axis';

    for q = 1:size(E,1)
        febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid = q;
        febio_spec.MeshData.ElementData{1}.elem{q}.a = e1_dir(q,:);
        febio_spec.MeshData.ElementData{1}.elem{q}.d = e2_dir(q,:);
    end
    %% Boundary conditions -- Y-loading, consistent with AUXETIC

    % TOP prescribed patches:
    % ux free, uy prescribed in Step 1, uz fixed
    febio_spec.Boundary.bc{1}.ATTR.name = ...
        'zero_displacement_z_top';
    febio_spec.Boundary.bc{1}.ATTR.type = ...
        'zero displacement';
    febio_spec.Boundary.bc{1}.ATTR.node_set = nodeSetName1;

    febio_spec.Boundary.bc{1}.x_dof = 0;
    febio_spec.Boundary.bc{1}.y_dof = 0;
    febio_spec.Boundary.bc{1}.z_dof = 1;


    % BOTTOM support patches:
    % ux free, uy fixed, uz fixed
    febio_spec.Boundary.bc{2}.ATTR.name = ...
        'zero_displacement_yz_bottom';
    febio_spec.Boundary.bc{2}.ATTR.type = ...
        'zero displacement';
    febio_spec.Boundary.bc{2}.ATTR.node_set = nodeSetName2;

    febio_spec.Boundary.bc{2}.x_dof = 0;
    febio_spec.Boundary.bc{2}.y_dof = 1;
    febio_spec.Boundary.bc{2}.z_dof = 1;


    % One reference node:
    % ux fixed only to prevent rigid-body translation in X
    febio_spec.Boundary.bc{3}.ATTR.name = ...
        'zero_displacement_x_pin';
    febio_spec.Boundary.bc{3}.ATTR.type = ...
        'zero displacement';
    febio_spec.Boundary.bc{3}.ATTR.node_set = nodeSetName3;

    febio_spec.Boundary.bc{3}.x_dof = 1;
    febio_spec.Boundary.bc{3}.y_dof = 0;
    febio_spec.Boundary.bc{3}.z_dof = 0;


    %% STEP 1: displacement-controlled loading in Y

    febio_spec.Step.step{1}.Boundary.bc{1}.ATTR.name = ...
        'prescribed_displacement_y';

    febio_spec.Step.step{1}.Boundary.bc{1}.ATTR.type = ...
        'prescribed displacement';

    febio_spec.Step.step{1}.Boundary.bc{1}.ATTR.node_set = ...
        nodeSetName1;

    febio_spec.Step.step{1}.Boundary.bc{1}.dof = 'y';

    febio_spec.Step.step{1}.Boundary.bc{1}.value.ATTR.lc = 1;
    febio_spec.Step.step{1}.Boundary.bc{1}.value.VAL = ...
        displacementMagnitude;

    febio_spec.Step.step{1}.Boundary.bc{1}.relative = 0;


    %% STEP 2: unloading only

    febio_spec.Step.step{2}.Loads.nodal_load{1}.ATTR.name = 'zero_force_unloading';
    febio_spec.Step.step{2}.Loads.nodal_load{1}.ATTR.type = 'nodal_target_force';
    febio_spec.Step.step{2}.Loads.nodal_load{1}.ATTR.node_set = nodeSetName1;
    febio_spec.Step.step{2}.Loads.nodal_load{1}.force = [0 0 0];
    febio_spec.Step.step{2}.Loads.nodal_load{1}.scale.ATTR.lc = 2;
    febio_spec.Step.step{2}.Loads.nodal_load{1}.scale.VAL = 1;


    %% Load controllers

    % Step 1 loading
    febio_spec.LoadData.load_controller{1}.ATTR.id = 1;
    febio_spec.LoadData.load_controller{1}.ATTR.name = 'LC_Loading';
    febio_spec.LoadData.load_controller{1}.ATTR.type = 'loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate = 'LINEAR';
    febio_spec.LoadData.load_controller{1}.extend = 'CONSTANT';
    febio_spec.LoadData.load_controller{1}.points.pt.VAL = ...
        [0 0; 1 1];

    % Step 2 unloading
    febio_spec.LoadData.load_controller{2}.ATTR.id = 2;
    febio_spec.LoadData.load_controller{2}.ATTR.name = 'LC_Unloading';
    febio_spec.LoadData.load_controller{2}.ATTR.type = 'loadcurve';
    febio_spec.LoadData.load_controller{2}.interpolate = 'LINEAR';
    febio_spec.LoadData.load_controller{2}.extend = 'CONSTANT';
    febio_spec.LoadData.load_controller{2}.points.pt.VAL = ...
        [1 0; 2 1];



    %Output section
    % -> log file
    logNamePrefix = sprintf('ang_%03d', ang);

    febioLogFileName      = [logNamePrefix, '.txt'];
    febioLogFileName_disp = [logNamePrefix, '_disp_out.txt'];
    febioLogFileName_stress_prin = [logNamePrefix, '_stress_prin_out.txt'];
    febioLogFileName_force = [logNamePrefix, '_force_out.txt'];
    % febioLogFileName_energy = [logNamePrefix, '_energy_out.txt'];
    % febioLogFileName_damage = [logNamePrefix, '_damage_out.txt'];
    febioLogFileName_stress_y = ...
        [logNamePrefix, '_stress_y_out.txt'];

   

    febio_spec.Output.logfile.ATTR.file=febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

    febio_spec.Output.logfile.node_data{2}.ATTR.file = ...
        febioLogFileName_force;

    febio_spec.Output.logfile.node_data{2}.ATTR.data = ...
        'Rx;Ry;Rz';

    febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

    % Export nodal reactions from the top prescribed/loading patches.
    febio_spec.Output.logfile.node_data{2}.VAL = ...
        mrow(bcPrescribeList);


    febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress_prin;
    febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

    % strain energy density
    febioLogFileName_energy = [logNamePrefix, '_energy_out.txt'];
    febio_spec.Output.logfile.element_data{2}.ATTR.file = febioLogFileName_energy;
    febio_spec.Output.logfile.element_data{2}.ATTR.data = 'sed';
    febio_spec.Output.logfile.element_data{2}.ATTR.delim = ',';

    % damage
    febioLogFileName_damage = [logNamePrefix, '_damage_out.txt'];
    febio_spec.Output.logfile.element_data{3}.ATTR.file = febioLogFileName_damage;
    febio_spec.Output.logfile.element_data{3}.ATTR.data = 'D';
    febio_spec.Output.logfile.element_data{3}.ATTR.delim = ',';

    % Axial Cauchy stress in the Y-direction
    febio_spec.Output.logfile.element_data{4}.ATTR.file = febioLogFileName_stress_y;
    febio_spec.Output.logfile.element_data{4}.ATTR.data = 'sy';
    febio_spec.Output.logfile.element_data{4}.ATTR.delim = ',';



    % Plotfile section
    febio_spec.Output.plotfile.compression=0;


    %% Running the FEBio analysis

    febFileName = fullfile( ...
        savePath, ...
        sprintf('%s_ang_%03d.feb',febioFebFileNamePart,ang));

    if runFEBio

        % Delete only files that will be regenerated during this simulation pass.
        if exist(fullfile(savePath,febioLogFileName),'file')
            delete(fullfile(savePath,febioLogFileName));
        end

        if exist(fullfile(savePath,febioLogFileName_disp),'file')
            delete(fullfile(savePath,febioLogFileName_disp));
        end

        if exist(fullfile(savePath,febioLogFileName_stress_prin),'file')
            delete(fullfile(savePath,febioLogFileName_stress_prin));
        end

        if exist(fullfile(savePath,febioLogFileName_force),'file')
            delete(fullfile(savePath,febioLogFileName_force));
        end

        if exist(fullfile(savePath,febioLogFileName_energy),'file')
            delete(fullfile(savePath,febioLogFileName_energy));
        end

        if exist(fullfile(savePath,febioLogFileName_damage),'file')
            delete(fullfile(savePath,febioLogFileName_damage));
        end

        if exist(fullfile(savePath,febioLogFileName_stress_y),'file')
            delete(fullfile(savePath,febioLogFileName_stress_y));
        end

       

        % Export and run the full loading-unloading Damage2 model.
        febioStruct2xml(febio_spec,febFileName);

        febioAnalysis.run_filename = febFileName;
        febioAnalysis.run_logname  = fullfile(savePath,febioLogFileName);
        febioAnalysis.disp_on = 1;
        febioAnalysis.runMode = runMode;
        febioAnalysis.maxLogCheckTime = 300;

        runFlag = runMonitorFEBio(febioAnalysis);

        if runFlag ~= 1
            error('FEBio failed for Zimmer Damage2 angle %d degrees.',ang);
        end

        fprintf('\nFEBio completed for Zimmer Damage2 angle %d degrees.\n',ang);

    else

        % Post-processing pass: use the existing FEBio logfile outputs.
        % Do not create, overwrite, or rerun FEBio files.
        runFlag = 1;

    end

    %% Import FEBio results


    if runFlag==1 %i.e. a succesful run

        %%

        % Importing nodal displacements from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

        %Access data
        N_disp_mat=dataStruct.data; %Displacement
        timeVec=dataStruct.time; %Time

        %Create deformed coordinate set
        V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);

        %% Import principal stresses

        dataStruct = importFEBio_logfile( ...
            fullfile(savePath,febioLogFileName_stress_prin),0,1);

        E_stress_mat = dataStruct.data;

        %%
        % Import axial Cauchy stress in Y
        dataStruct = importFEBio_logfile( ...
            fullfile(savePath,febioLogFileName_stress_y),0,1);

        E_stress_y_mat = dataStruct.data;

      
        % Import strain-energy density
        dataStruct = importFEBio_logfile( ...
            fullfile(savePath,febioLogFileName_energy),0,1);

        E_energy_mat = dataStruct.data;

        % Import damage
        dataStruct = importFEBio_logfile( ...
            fullfile(savePath,febioLogFileName_damage),0,1);

        E_damage_mat = dataStruct.data;


       %% Import TOP-PATCH reaction forces directly from FEBio logfile

        dataStruct_force = importFEBio_logfile( ...
            fullfile(savePath,febioLogFileName_force),0,1);

        N_force_mat = dataStruct_force.data;
        timeForceVec = dataStruct_force.time(:);

        % Verify force-log dimensions
        if size(N_force_mat,2) ~= 3
            error('Expected FEBio reaction-force output [Rx Ry Rz].');
        end

        if size(N_force_mat,1) ~= numel(bcPrescribeList)
            error('Number of reaction-force nodes does not match the top prescribed node set.');
        end

        %% Y-DIRECTION POST-PROCESSING
        % Same coordinate convention as AUXETIC:
        % Y = axial/loading direction
        % X = in-plane transverse direction
        % Z = through-thickness direction

        % Mean axial Cauchy stress in Y
        stress_y_sim = squeeze(mean(E_stress_y_mat,1,'omitnan'));
        stress_y_sim = stress_y_sim(:);

        % Hydrostatic pressure calculated from principal Cauchy stresses
        pressure_sim = squeeze( ...
            -1/3 .* mean(sum(E_stress_mat,2),1,'omitnan'));
        pressure_sim = pressure_sim(:);


        %% Sum top prescribed-patch Y reactions

        % N_force_mat columns are [Rx Ry Rz]; Y loading therefore uses Ry.
        Fy_signed = squeeze(sum(N_force_mat(:,2,:),1));
        Fy_signed = Fy_signed(:);

        if numel(Fy_signed) ~= numel(timeForceVec)
            error('Reaction-force and time vectors have inconsistent lengths.');
        end

        % Loading step ends at global analysis time approximately 1.
        [~,iPeakForce] = min(abs(timeForceVec(:)-1));

        % Report tensile force as positive irrespective of FEBio reaction sign.
        if Fy_signed(iPeakForce) < 0
            Fy_total = -Fy_signed;
        else
            Fy_total = Fy_signed;
        end

        fprintf('Peak total force: %.6g N\n',max(Fy_total));
        fprintf('Final unloaded force: %.6g N\n',Fy_total(end));


        %% Global specimen stretches -- same convention as AUXETIC

        w0 = max(V(:,1)) - min(V(:,1));   % transverse X dimension
        h0 = max(V(:,2)) - min(V(:,2));   % axial/loading Y dimension
        A0 = w0*h0;

        nSteps = size(V_DEF,3);

        lambda_x_vec = zeros(nSteps,1);   % transverse stretch
        lambda_y_vec = zeros(nSteps,1);   % axial stretch
        area_ratio_vec = zeros(nSteps,1);

        for i = 1:nSteps

            V_temp = V_DEF(:,:,i);

            w = max(V_temp(:,1)) - min(V_temp(:,1));
            h = max(V_temp(:,2)) - min(V_temp(:,2));
            A = w*h;

            lambda_x_vec(i) = w/w0;
            lambda_y_vec(i) = h/h0;
            area_ratio_vec(i) = A/A0;

        end


        %% Match Y reaction force to displacement-output times

        if numel(timeForceVec) < 2 || any(diff(timeForceVec(:)) <= 0)
            error('The reaction-force logfile time vector is not strictly increasing.');
        end

        Fy_vec = interp1( ...
            timeForceVec(:),Fy_total,timeVec(:), ...
            'linear',NaN);

        if any(~isfinite(Fy_vec))
            error(['Reaction-force logfile does not cover the complete FEBio ' ...
                'displacement time history.']);
        end


        %% Axial Cauchy stress versus global axial stretch

        nStress = min(numel(lambda_y_vec),numel(stress_y_sim));

        cFigure; hold on;

        plot(lambda_y_vec(1:nStress), ...
            stress_y_sim(1:nStress), ...
            'r-','LineWidth',3);

        title(sprintf( ...
            'Axial stress-stretch | %s',angLabel), ...
            'FontSize',fontSize);

        xlabel('\lambda_y','FontSize',fontSize);
        ylabel('\sigma_{yy} [MPa]','FontSize',fontSize);

        grid on;
        box on;
        axis square;
        set(gca,'FontSize',fontSize);


        %% Pressure versus global axial stretch

        nPressure = min(numel(lambda_y_vec),numel(pressure_sim));

        cFigure; hold on;

        plot(lambda_y_vec(1:nPressure), ...
            pressure_sim(1:nPressure), ...
            'r-','LineWidth',3);

        title(sprintf( ...
            'Pressure-stretch | %s',angLabel), ...
            'FontSize',fontSize);

        xlabel('\lambda_y','FontSize',fontSize);
        ylabel('Pressure [MPa]','FontSize',fontSize);

        grid on;
        box on;
        axis square;
        set(gca,'FontSize',fontSize);


        %% Time-reaction-force curve

        cFigure; hold on;

        plot(timeForceVec(:),Fy_total, ...
            'b-','LineWidth',3);

        title(sprintf( ...
            'Reaction force history | %s',angLabel), ...
            'FontSize',fontSize);

        xlabel('Time','FontSize',fontSize);
        ylabel('F_y [N]','FontSize',fontSize);

        grid on;
        box on;
        axis square;
        set(gca,'FontSize',fontSize);


        %% Volume-weighted energy and damage histories
        nHistorySteps = min(size(E_energy_mat,3),size(E_damage_mat,3));

        energy_total_vec = zeros(nHistorySteps,1);
        energy_total_sum = zeros(nHistorySteps,1);

        % Add these two
        energy_damage_vec = nan(nHistorySteps,1);
        energy_intact_vec = nan(nHistorySteps,1);

        damage_vec        = zeros(nHistorySteps,1);
        damage_max_vec    = zeros(nHistorySteps,1);
        damage_fraction_vec = zeros(nHistorySteps,1);

        if ~strcmpi(elementType,'hex8')
            error(['This post-processing uses reference element volumes and is ' ...
                'implemented for the current hex8 Zimmer model.']);
        end

        elementWeights = elementVolume0(:); % reference element volumes [mm^3]

        if numel(elementWeights) ~= size(E_energy_mat,1) || ...
                numel(elementWeights) ~= size(E_damage_mat,1)
            error(['Number of reference element volumes does not match the ' ...
                'number of FEBio element output values.']);
        end

        totalVolume = sum(elementWeights,'omitnan');
        if ~isfinite(totalVolume) || totalVolume <= 0
            error('Invalid total reference volume for post-processing.');
        end

        % FEBio reports D for a solid mixture as the aggregate/sum of constituent
        % damage variables.  Therefore D>=0.99 is treated here as a HIGH-DAMAGE
        % threshold, not as proof that every constituent is fully damaged.
        % Same small damage threshold as auxetic code
        damageOnsetThreshold = 0.01;

        % High-damage threshold
        highDamageThreshold = 0.99;

        for i = 1:nHistorySteps
            E_step = E_energy_mat(:,1,i);
            D_step = E_damage_mat(:,1,i);

            goodE = isfinite(E_step) & isfinite(elementWeights) & elementWeights>0;
            goodD = isfinite(D_step) & isfinite(elementWeights) & elementWeights>0;

            % Integrated strain energy: sed [N/mm^2] * reference volume [mm^3]
            energy_total_sum(i) = sum(E_step(goodE).*elementWeights(goodE),'omitnan');
            energy_total_vec(i) = energy_total_sum(i) / sum(elementWeights(goodE),'omitnan');

            % Same damaged/intact threshold as auxetic code
            damaged = goodE & goodD & ...
                (D_step > damageOnsetThreshold);

            intact = goodE & goodD & ...
                (D_step <= damageOnsetThreshold);

            % Mean strain-energy density in damaged region
            if any(damaged)

                energy_damage_vec(i) = ...
                    sum(E_step(damaged).*elementWeights(damaged),'omitnan') / ...
                    sum(elementWeights(damaged),'omitnan');

            end

            % Mean strain-energy density in intact region
            if any(intact)

                energy_intact_vec(i) = ...
                    sum(E_step(intact).*elementWeights(intact),'omitnan') / ...
                    sum(elementWeights(intact),'omitnan');

            end

            % Volume-weighted damage statistics.
            damage_vec(i) = ...
                sum(D_step(goodD).*elementWeights(goodD),'omitnan') / ...
                sum(elementWeights(goodD),'omitnan');
            damage_max_vec(i) = max(D_step,[],'omitnan');

            isHighDamage = goodD & (D_step >= highDamageThreshold);
            damage_fraction_vec(i) = ...
                sum(elementWeights(isHighDamage),'omitnan') / totalVolume;
        end

        %% Damage versus actual global axial stretch

        nDamageHist = min(numel(lambda_y_vec),numel(damage_vec));

        cFigure; hold on;

        plot(lambda_y_vec(1:nDamageHist), ...
            damage_vec(1:nDamageHist), ...
            'm-','LineWidth',3);

        plot(lambda_y_vec(1:nDamageHist), ...
            damage_max_vec(1:nDamageHist), ...
            'k--','LineWidth',3);

        legend({'Mean damage','Maximum element damage'}, ...
            'Location','best');

        xlabel('\lambda_y');
        ylabel('Mean damage');

        title(sprintf( ...
            'Damage versus stretch | %s',angLabel));

        grid on;
        box on;
        axis square;


        %% High aggregate-damage volume fraction versus axial stretch

        nDamage = min(numel(lambda_y_vec),numel(damage_fraction_vec));

        cFigure; hold on;

        plot(lambda_y_vec(1:nDamage), ...
            damage_fraction_vec(1:nDamage), ...
            'LineWidth',3);

        xlabel('\lambda_y');
        ylabel(sprintf('Volume fraction with D \\geq %.2f',highDamageThreshold));

        title(sprintf( ...
            'Spatial spread of high aggregate damage | %s',angLabel));

        grid on;
        box on;
        axis square;


        %% Reaction force versus actual global axial stretch

        cFigure; hold on;

        [~,iPeakStretch] = max(lambda_y_vec);

        hLoad = plot( ...
            lambda_y_vec(1:iPeakStretch), ...
            Fy_vec(1:iPeakStretch), ...
            'b-', ...
            'LineWidth',3);

        hUnload = plot( ...
            lambda_y_vec(iPeakStretch:end), ...
            Fy_vec(iPeakStretch:end), ...
            'b--', ...
            'LineWidth',3);

        legend([hLoad,hUnload], ...
            {'Loading','Unloading'}, ...
            'Location','best');

        xlabel('\lambda_y');
        ylabel('Reaction force F_y [N]');

        title(sprintf( ...
            'Loading-unloading response | %s',angLabel));

        grid on;
        box on;
        axis square;


        %% Von Mises stress from principal Cauchy stresses

        E_stress_mat_VM = sqrt(( ...
            (E_stress_mat(:,1,:)-E_stress_mat(:,2,:)).^2 + ...
            (E_stress_mat(:,2,:)-E_stress_mat(:,3,:)).^2 + ...
            (E_stress_mat(:,1,:)-E_stress_mat(:,3,:)).^2 )/2);


        %% Volume-weighted 90th percentile of von Mises stress

        p90_vec = zeros(size(E_stress_mat_VM,3),1);

        for qt = 1:size(E_stress_mat_VM,3)
            sigmaVM = E_stress_mat_VM(:,1,qt);
            p90_vec(qt) = weightedPercentile( ...
                sigmaVM(:),elementWeights,0.90);
        end

        nP90 = min(numel(lambda_y_vec),numel(p90_vec));

        cFigure; hold on;

        plot(lambda_y_vec(1:nP90), ...
            p90_vec(1:nP90), ...
            'r-o','LineWidth',3,'MarkerFaceColor','r');

        xlabel('\lambda_y');
        ylabel('p90 stress');
        title(sprintf('p90 stress vs stretch | %s',angLabel));
        axis square;
        grid on;
        box on;


        %% Apparent/global true Poisson ratio
        % AUXETIC convention: axial = Y, transverse = X.

        nu_true = nan(size(lambda_y_vec));

        validNu = ...
            lambda_y_vec > 0 & ...
            lambda_x_vec > 0 & ...
            abs(log(lambda_y_vec)) > 1e-10;

        nu_true(validNu) = ...
            -log(lambda_x_vec(validNu)) ./ ...
            log(lambda_y_vec(validNu));

        cFigure; hold on;

        plot(lambda_y_vec,nu_true,'b-','LineWidth',5);
        xlabel('\lambda_y');
        ylabel('\nu_{true}');
        title(sprintf('Poisson ratio vs stretch | %s',angLabel));
        axis square;
        grid on;
        box on;


        %% Projected area ratio

        cFigure; hold on;

        plot(lambda_y_vec,area_ratio_vec,'k-','LineWidth',5);
        xlabel('\lambda_y');
        ylabel('A/A_0');
        title(sprintf('Projected area ratio vs stretch | %s',angLabel));
        axis square;
        grid on;
        box on;

        % Damage contour
        CVd = faceToVertexMeasure(E,V,E_damage_mat(:,:,end));

        hf = cFigure;
        title(sprintf('Damage contour | %s',angLabel));
        hp = gpatch(F,V_DEF(:,:,end),CVd,'none',1,lineWidth1);
        hp.FaceColor = 'interp';
        axisGeom(gca,fontSize);
        colormap(parula); colorbar;
        cMax = max(E_damage_mat(:),[],'omitnan');
        if ~isfinite(cMax) || cMax <= 0
            cMax = 1;
        end
        caxis([0 cMax]);
        camlight headlight;

        % Damage animation
        hf = cFigure;
        title(sprintf('Damage animation | %s',angLabel));

        CVd = faceToVertexMeasure(E,V,E_damage_mat(:,:,end));
        hp = gpatch(F,V_DEF(:,:,end),CVd,'none',1,lineWidth1);
        hp.FaceColor = 'interp';

        axisGeom(gca,fontSize);
        colormap(parula); colorbar;
        cMax = max(E_damage_mat(:),[],'omitnan');

        if isempty(cMax) || ~isfinite(cMax) || cMax <= 0
            cMax = 1;
        end

        caxis([0 cMax]);
        axis(axisLim(V_DEF));
        camlight headlight;

        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)
            CVd = faceToVertexMeasure(E,V,E_damage_mat(:,:,qt));

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt} = {'Vertices','CData'};
            animStruct.Set{qt} = {V_DEF(:,:,qt),CVd};
        end

        anim8(hf,animStruct);

        hf = cFigure;
        title(sprintf('p90 clipped stress animation | %s',angLabel));

        CV = faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,end));
        CV = min(CV,p90_vec(end));

        hp = gpatch(F,V_DEF(:,:,end),CV,'none',1,lineWidth1);
        hp.FaceColor = 'interp';

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        cMax = max(p90_vec,[],'omitnan');
        if ~isfinite(cMax) || cMax <= 0
            cMax = 1;
        end
        caxis([0 cMax]);
        axis(axisLim(V_DEF));
        camlight headlight;

        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)
            sigmaVM = E_stress_mat_VM(:,:,qt);
            CV = faceToVertexMeasure(E,V,sigmaVM);
            CV = min(CV,p90_vec(qt));

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt} = {'Vertices','CData'};
            animStruct.Set{qt} = {V_DEF(:,:,qt),CV};
        end

        anim8(hf,animStruct);

        % x-displacement animation
        hf = cFigure;
        title(sprintf('Transverse x-displacement animation | %s',angLabel));

        CV = abs(N_disp_mat(:,1,end));
        hp = gpatch(F,V_DEF(:,:,end),CV,'none',1,lineWidth1);
        hp.FaceColor = 'interp';

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([0 max(abs(N_disp_mat(:,1,:)),[],'all')]);
        axis(axisLim(V_DEF));
        camlight headlight;

        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)
            CV = abs(N_disp_mat(:,1,qt));

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt} = {'Vertices','CData'};
            animStruct.Set{qt} = {V_DEF(:,:,qt),CV};
        end

        anim8(hf,animStruct);

        %% Total displacement-magnitude animation

        DN_magnitude = sqrt( ...
            sum(N_disp_mat(:,:,end).^2,2));

        hf = cFigure;

        title(sprintf( ...
            'Displacement magnitude | %s',angLabel), ...
            'FontSize',fontSize);

        hp = gpatch( ...
            F,V_DEF(:,:,end),DN_magnitude, ...
            'none',1,lineWidth1);

        hp.FaceColor='interp';

        axisGeom(gca,fontSize);
        colormap(cMap);
        colorbar;

        DN_all = sqrt(sum(N_disp_mat.^2,2));
        cMax = max(DN_all,[],'all','omitnan');

        if ~isfinite(cMax) || cMax <= 0
            cMax = 1;
        end

        caxis([0 cMax]);
        axis(axisLim(V_DEF));
        camlight headlight;

        animStruct=[];
        animStruct.Time=timeVec;

        for qt=1:size(N_disp_mat,3)

            DN_magnitude = sqrt( ...
                sum(N_disp_mat(:,:,qt).^2,2));

            animStruct.Handles{qt}=hp;
            animStruct.Props{qt}={'Vertices','CData'};
            animStruct.Set{qt} = ...
                {V_DEF(:,:,qt),DN_magnitude};
        end

        anim8(hf,animStruct);



        %% Axial Cauchy-stress animation

        CV = faceToVertexMeasure( ...
            E,V,E_stress_y_mat(:,:,end));

        hf = cFigure;

        title(sprintf( ...
            '\\sigma_{yy} animation | %s',angLabel), ...
            'FontSize',fontSize);

        hp = gpatch( ...
            F,V_DEF(:,:,end),CV, ...
            'none',1,lineWidth1);

        hp.FaceColor = 'interp';

        axisGeom(gca,fontSize);
        colormap(cMap);
        colorbar;

        cMin = min(E_stress_y_mat(:),[],'omitnan');
        cMax = max(E_stress_y_mat(:),[],'omitnan');

        if ~isfinite(cMin)
            cMin = 0;
        end

        if ~isfinite(cMax) || cMax <= cMin
            cMax = cMin + 1;
        end

        caxis([cMin cMax]);
        axis(axisLim(V_DEF));
        camlight headlight;

        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)

            CV = faceToVertexMeasure( ...
                E,V,E_stress_y_mat(:,:,qt));

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt} = {'Vertices','CData'};
            animStruct.Set{qt} = ...
                {V_DEF(:,:,qt),CV};
        end

        anim8(hf,animStruct);



        %%
        % Plotting the simulated results using |anim8| to visualize and animate

        [CV]=faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,end));

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure;

        gtitle([febioFebFileNamePart,': Press play to animate']);
        title('$\sigma_{vm}$ [MPa]','Interpreter','Latex')
        hp1=gpatch(F,V_DEF(:,:,end),CV,'k',1,0.5); %Add graphics object to animate
        hp1.FaceColor='interp';

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(E_stress_mat_VM(:))]);
        axis(axisLim(V_DEF)); %Set axis limits statically
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments

            [CV]=faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,qt));

            %Set entries in animation structure
            animStruct.Handles{qt}=hp1;
            animStruct.Props{qt}={'Vertices','CData'};
            animStruct.Set{qt}={V_DEF(:,:,qt),CV};
        end
        anim8(hf,animStruct); %Initiate animation feature

        results(ia).lambda          = lambda_y_vec;
        results(ia).force           = Fy_vec;
        results(ia).energy_mean = energy_total_vec;
        results(ia).energy_sum  = energy_total_sum;

        % Add these
        results(ia).energy_damage_mean = energy_damage_vec;
        results(ia).energy_intact_mean = energy_intact_vec;

        results(ia).damage = damage_vec;
        results(ia).damage_max      = damage_max_vec;
        results(ia).damage_fraction = damage_fraction_vec;
        results(ia).p90_vm          = p90_vec;
        results(ia).nu_true         = nu_true;
        results(ia).area_ratio      = area_ratio_vec;

        % Explicit per-orientation processed result file (not a whole-workspace dump).
        resultOrientation = results(ia); %#ok<NASGU>
        save(fullfile(savePath, ...
            sprintf('ZimmerDamage_%s_processed.mat',angLabel)), ...
            'resultOrientation','ang','angLabel','-v7.3');
        gdrawnow;

    end
end

%% Final comparison plots and combined save: only after post-processing pass

if all(arrayfun(@(s) ~isempty(s.lambda), results))

    plotComparison(results, 'force', ...
        'Reaction Force F_y [N]', ...
        'Force vs Axial Stretch (Perpendicular vs Parallel)');

    plotComparison(results, 'energy_sum', ...
        'Integrated Strain Energy [N mm]', ...
        'Strain Energy vs Axial Stretch (Perpendicular vs Parallel)');

    plotComparison(results, 'damage', ...
        'Volume-weighted mean aggregate damage D', ...
        'Mean Aggregate Damage vs Axial Stretch (Perpendicular vs Parallel)');

    plotComparison(results, 'damage_max', ...
        'Maximum aggregate damage D', ...
        'Maximum Aggregate Damage vs Axial Stretch (Perpendicular vs Parallel)');

    plotComparison(results, 'damage_fraction', ...
        sprintf('Volume fraction with aggregate D \geq %.2f',highDamageThreshold), ...
        'Spatial Spread of High Aggregate Damage (Perpendicular vs Parallel)');

    plotComparison(results, 'p90_vm', ...
        'Volume-weighted p90 von Mises stress [MPa]', ...
        'p90 von Mises Stress vs Axial Stretch (Perpendicular vs Parallel)');

    % Save both orientations plus essential model metadata.
    modelInfo = struct();
    modelInfo.model = 'Zimmer damage';
    modelInfo.loadingDirection = 'y';
    modelInfo.orientationConvention = '0 deg = Perpendicular; 90 deg = Parallel';
    modelInfo.appliedStretch = appliedStretch;
    modelInfo.displacementMagnitude = displacementMagnitude;
    modelInfo.layerThickness = layerThickness;
    modelInfo.elementType = elementType;
    modelInfo.numRefine = numRefine;
    modelInfo.k_factor = k_factor;
    modelInfo.c1 = c1_ini;
    modelInfo.m1 = m1_ini;
    modelInfo.k1 = k1_ini;
    modelInfo.k2 = k2_ini;
    modelInfo.kappa = kappa_ini;
    modelInfo.gamma = gamma;
    modelInfo.mu_matrix = mu_matrix_ini;
    modelInfo.sigma_matrix = sigma_matrix_ini;
    modelInfo.Dmax_matrix = Dmax_matrix_ini;
    modelInfo.mu_GOH = mu_GOH_ini;
    modelInfo.sigma_GOH = sigma_GOH_ini;
    modelInfo.Dmax_GOH = Dmax_GOH_ini;
    modelInfo.damageOnsetThreshold = damageOnsetThreshold;
    modelInfo.highDamageThreshold = highDamageThreshold;

    resultsFile = fullfile(savePath,'ZimmerDamage_processed_results.mat');
    save(resultsFile,'results','modelInfo','angles_deg','-v7.3');

    fprintf('\nProcessed Zimmer damage results saved as:\n%s\n',resultsFile);

else

    warning(['Final comparison plots were skipped because at least one ' ...
        'orientation did not produce processed results.']);

end

function plotComparison(results, fieldName, yLabelText, titleText)
% plotComparison  Compare 0-deg Perpendicular and 90-deg Parallel responses.

cFigure; hold on;

% 0 degrees = Perpendicular
lambdaPerp = results(1).lambda(:);
dataPerp   = results(1).(fieldName);
dataPerp   = dataPerp(:);

% 90 degrees = Parallel
lambdaPara = results(2).lambda(:);
dataPara   = results(2).(fieldName);
dataPara   = dataPara(:);

nPerp = min(numel(lambdaPerp),numel(dataPerp));
nPara = min(numel(lambdaPara),numel(dataPara));

lambdaPerp = lambdaPerp(1:nPerp);
dataPerp   = dataPerp(1:nPerp);
lambdaPara = lambdaPara(1:nPara);
dataPara   = dataPara(1:nPara);

if nPerp < 2 || nPara < 2
    warning('Insufficient data for comparison field %s.',fieldName);
    return;
end

[~,iPeakPerp] = max(lambdaPerp);
[~,iPeakPara] = max(lambdaPara);

% Perpendicular: blue
h1 = plot(lambdaPerp(1:iPeakPerp),dataPerp(1:iPeakPerp), ...
    '-','Color',[0.15 0.35 0.85],'LineWidth',5);

h2 = plot(lambdaPerp(iPeakPerp:end),dataPerp(iPeakPerp:end), ...
    '--','Color',[0.15 0.35 0.85],'LineWidth',5);

% Parallel: red
h3 = plot(lambdaPara(1:iPeakPara),dataPara(1:iPeakPara), ...
    '-','Color',[0.85 0.15 0.15],'LineWidth',5);

h4 = plot(lambdaPara(iPeakPara:end),dataPara(iPeakPara:end), ...
    '--','Color',[0.85 0.15 0.15],'LineWidth',5);

xlabel('\lambda_y (Axial Stretch)', ...
    'Interpreter','tex','FontSize',20,'FontWeight','bold');
ylabel(yLabelText, ...
    'Interpreter','tex','FontSize',20,'FontWeight','bold');
title(titleText, ...
    'Interpreter','tex','FontSize',20,'FontWeight','bold');

legend([h1 h2 h3 h4], ...
    {'Perpendicular loading','Perpendicular unloading', ...
    'Parallel loading','Parallel unloading'}, ...
    'FontSize',14,'Location','best');

set(gca,'FontSize',16,'FontWeight','bold');
lambdaAll = [lambdaPerp;lambdaPara];
lambdaAll = lambdaAll(isfinite(lambdaAll));

if ~isempty(lambdaAll)
    xlim([min(lambdaAll) max(lambdaAll)]);
end

axis square;
box on;
grid on;
end

function q = weightedPercentile(x,w,p)
% weightedPercentile  Weighted empirical percentile.
%   q = weightedPercentile(x,w,p) returns the value q for which the
%   cumulative normalized positive weight first reaches p (0 <= p <= 1).

x = x(:);
w = w(:);

if numel(x) ~= numel(w)
    error('weightedPercentile:SizeMismatch', ...
        'Data and weights must have the same number of entries.');
end

if ~isscalar(p) || ~isfinite(p) || p < 0 || p > 1
    error('weightedPercentile:InvalidProbability', ...
        'p must be a finite scalar between 0 and 1.');
end

good = isfinite(x) & isfinite(w) & w > 0;
x = x(good);
w = w(good);

if isempty(x)
    q = NaN;
    return;
end

[x,idx] = sort(x);
w = w(idx);
cw = cumsum(w) ./ sum(w);

iq = find(cw >= p,1,'first');
if isempty(iq)
    iq = numel(x);
end
q = x(iq);
end

