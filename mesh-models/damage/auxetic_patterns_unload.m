
clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8; %transparency
markerSize=40; %For plotted points
markerSize2=10; %For nodes on patches
lineWidth1=1; %For meshes
cMap=spectral(250); %colormap


%% Control parameters

meshSearchTolerance=1e-3; %Tolerance for logic on mesh searching thresholds


%% Material parameters from the second code


k_factor = 100;

% 1-term Ogden matrix
c1_ini = 2.7271875051643106e+00;
m1_ini = 6.9427891844363305e+00;
k_matrix_ini = c1_ini*k_factor;

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

LangerAngle_perp = 0/180*pi; % Angle of Langer line with loading direction
LangerAngle_para = 90/180*pi; % Angle of Langer line with loading direction
angleSelection=1;%1 perpendicular to LL and 2 parallel
LangerAngles = [LangerAngle_perp, LangerAngle_para];

% FEA control settings
numTimeSteps=20; %Number of time steps desired
max_refs=500; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=20; %Maximum number of retires
dtmin=(1/numTimeSteps)/50; %Minimum time step size
dtmax=(1/numTimeSteps); %Maximum time step size
runMode='internal';
runFEBio = true; %true manual extraction of CSV file and then false for running the full simulation

% Path names
defaultFolder = fileparts(mfilename('fullpath'));
savePath=fullfile(defaultFolder,'data','temp');

% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting principal stress
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%% geometries

% pointSpacing =2;
% distRefine = [3];

% pointSpacing =2;
% distRefine = [4 3];

% pointSpacing =1;
% distRefine = [1];

% pointSpacing =1;
% distRefine = [3 2];


pointSpacing =2;
distRefine = [2];


geometryID = 5;   % Use 1, 3, 4, or 5

switch geometryID
    case 1
        geometryLabel = 'AS1';
    case 3
        geometryLabel = 'G3';
    case 4
        geometryLabel = 'G4';
    case 5
        geometryLabel = 'HS';
    otherwise
        error('Only geometries 1, 3, 4, and 5 are selected.');
end

geometryTag = sprintf('G%02d_%s',geometryID,geometryLabel);

geo = allGeometries(geometryID,pointSpacing,distRefine);

savePath = fullfile(defaultFolder,'data',geometryTag);

if ~exist(savePath,'dir')
    mkdir(savePath);
end
if ~exist(savePath,'dir'); mkdir(savePath); end

F = geo.F;
V = geo.V;
indRefine = geo.indRefine;
V_refine = geo.V_refine;
V1 = geo.V1;
p = geo.p;
d = geo.d;
t = 3;

%FEBio Parameters
h0 = max(V(:,2)) - min(V(:,2));  % Initial mesh height
displacementMagnitude = 61.59*0.4;  % Slightly boost to avoid lambda_y = 1

%% Creating triangle Mesh

for i = 1:1:length(distRefine)
    % VF = patchCentre(F,V);
    % [distCorners,indMin] = minDist(VF,V_refine);

    %Compute distances on mesh description
    distCorners_V=meshDistMarch(F,V,indRefine);
    distCorners=vertexToFaceMeasure(F,distCorners_V);

    % cFigure;
    % hold on;
    % % % title('Colormapped face colors');
    % Hp=gpatch(F,V, distCorners_V,'k');
    % Hp.FaceColor='interp';
    % plotV(V(indRefine,:),'r.','Markersize',25);
    % plotV(V_refine,'b.','Markersize',25);
    % axisGeom;
    % drawnow;


    logicRefine = distCorners<=distRefine(i);
    logicRefine = triSurfLogicSharpFix(F, logicRefine, 3);

    % [F,V,C_type,indIni,C]=subTriDual(F,V,logicRefine);
    inputStruct.F=F;
    inputStruct.V=V;
    inputStruct.indFaces=find(logicRefine);
    % inputStruct.f=f;
    [outputStruct]=subTriLocal(inputStruct);
    F=outputStruct.F;
    V=outputStruct.V;

    %Smooth with constraints
    Eb = patchBoundary(F);
    indRigid = unique(Eb(:));
    smoothParameters.n=25; %Number of iterations
    smoothParameters.Method='LAP'; %Smoothing method
    smoothParameters.RigidConstraints=indRigid; %Vertices to not keep constant
    V=patchSmooth(F,V,[],smoothParameters);
end
V = V(:,[1,2]);

% cFigure;
% hold on;
% % title('Colormapped face colors');
% gpatch(F,V,'gw','k');
% gpatch(F,V,'w','k');
% pointAnnotate(V1,1:1:length(V1),'fontsize',15)
% axisGeom;
% drawnow;


%%
F=[F;fliplr(F)+length(V)];
Vc=V;
Vc(:,1)=-Vc(:,1);
V=[V;Vc];


F=[F;fliplr(F)+length(V)];
Vc=V;
Vc(:,2)=-Vc(:,2);
V=[V;Vc];


%% replicate mesh

[FT,VT,CT]=replicatemesh(F,V,p,d);

%%

%Smooth with constraints
Eb = patchBoundary(FT);
indRigid = unique(Eb(:));
smoothParameters2.n=25; %Number of iterations
smoothParameters2.Method='LAP'; %Smoothing method
smoothParameters2.RigidConstraints=indRigid; %Vertices to not keep constant
VT=patchSmooth(FT,VT,[],smoothParameters2);

%% Visualizing triangle Mesh
% cFigure;
% hold on;
% % title('Colormapped face colors');
% gpatch(FT,VT,CT,'k');
% plotV(V1,'k.-','linewidth',3,'markersize', 15)
% % pointAnnotate(V1,1:1:length(V1),'fontsize',15)
% axisGeom;
% colormap(gjet(prod(p))); icolorbar;
% drawnow;


%% Thickened the solidmesh
numSteps=ceil(t/pointSpacing);
[E,V,Fp1,Fp2]=patchThick(FT,VT,1,t,numSteps);

%Use element2patch to get patch data
F=element2patch(E,[],'penta6');



%%

%Get boundary faces (two sets due to pentahedra)
indb = tesBoundary(F); %Cell containing boundary face indices
Fb = {F{1}(indb{1},:),F{2}(indb{2},:)}; %Cell containing boundary faces

%%
% Plotting meshed model
% cFigure; hold on;
% title('The meshed model','FontSize',fontSize);
%
% gpatch(F,V,'g','k',0.5); %All faces
% gpatch(Fb,V,'gw','k',faceAlpha1); %Boundary faces
%
% patchNormPlot(F,V); %Visualise normal directions
%
% axisGeom(gca,fontSize);
% camlight headlight;
% drawnow;


%% Defining the boundary conditions

%Prescribed displacement nodes
bcPrescribeList=find(V(:,2)> max(V(:,2))-meshSearchTolerance);
bcSupportList=find(V(:,2)< min(V(:,2))+meshSearchTolerance);

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

% hf=cFigure;
% title('Boundary conditions','FontSize',fontSize);
% xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
% hold on;
%
% gpatch(Fb,V,'kw','k',0.5);
%
% hl(1)=plotV(V(bcPrescribeList,:),'r.','MarkerSize',markerSize);
% hl(2)=plotV(V(bcSupportList,:),'g.','MarkerSize',markerSize);
%
% legend(hl,{'Tensile Positive', 'Tensile Negative'});
%
% axisGeom(gca,fontSize);
% camlight headlight;
% drawnow;


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
    'damage',[], ...
    'damage_max',[], ...
    'damage_fraction',[]);
results(2) = results;
angles_deg = [0 90];  % Loop from 0 to 90 degrees
% angles_deg = [0];  % Loop from 0 to 90 degrees

% for ang = angles_deg
% %     % Define fiber direction by rotating vector [1 0 0] around z-axis
% %     R = euler2DCM([0,0,-ang*pi/180]);
for ia = 1:numel(angles_deg)
    ang = angles_deg(ia);

    if ang == 0
        %         mu_ini_local    = mu_perp_ini;
        %         sigma_ini_local = sigma_perp_ini;
        %         Dmax_ini_local  = Dmax_perp_ini;
        angLabel        = 'Perpendicular';
    elseif ang == 90
        %         mu_ini_local    = mu_para_ini;
        %         sigma_ini_local = sigma_para_ini;
        %         Dmax_ini_local  = Dmax_para_ini;
        angLabel        = 'Parallel';
    else
        error('Only 0 and 90 degrees are supported in this setup.');
    end
    %
    %     % Define fiber direction by rotating vector [1 0 0] around z-axis
    R = euler2DCM([0,0,-ang*pi/180]);
    e1 = (R*[1 0 0]')';  % Rotated fiber direction

    e1_dir = repmat(e1, nElem, 1);
    e2_dir = repmat([0 0 1], nElem, 1);  % Keep this fixed
    e3_dir = cross(e1_dir, e2_dir, 2);   % Orthogonal vector

    % Plot material directions
    VE = patchCentre(E, V);  % Make sure VE is defined

    % cFigure; hold on;
    % title(sprintf('Fiber direction = %d°', ang), 'FontSize', 20);
    % gpatch(Fb, V, 'kw', 'none', 0.25);
    % quiverVec(VE, e1_dir, mean(pointSpacing), 'r');
    % quiverVec(VE, e2_dir, mean(pointSpacing)/2, 'g');
    % quiverVec(VE, e3_dir, mean(pointSpacing)/2, 'b');
    % axisGeom(gca, fontSize);
    % camlight headlight;
    % drawnow;



    %%
    % Visualizing material directions

    % [VE]=patchCentre(E,V);
    %
    % hf=cFigure; hold on;
    % title(sprintf('Material directions  |  Fiber Angle = %d°', ang), 'FontSize', fontSize);
    %
    % gpatch(Fb,V,'kw','none',0.25);
    % hf(1)=quiverVec(VE,e1_dir,mean(pointSpacing),'r');
    % hf(2)=quiverVec(VE,e2_dir,mean(pointSpacing)/2,'g');
    % hf(3)=quiverVec(VE,e3_dir,mean(pointSpacing)/2,'b');

    % legend(hf,{'e1-direction','e2-direction','e3-direction'});
    % axisGeom(gca,fontSize);
    % camlight headlight;
    %
    %
    % legend(hf,{'e1-direction (fiber)','e2-direction','e3-direction'});
    % axisGeom(gca,fontSize);
    % camlight headlight;
    % drawnow;


    %% Defining the FEBio input structure
    % See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
    % manual.

    %Get a template with default settings
    [febio_spec]=febioStructTemplate;

    %febio_spec version
    febio_spec.ATTR.version='4.0';

    %Module section
    febio_spec.Module.ATTR.type='solid';

    %% Control settings for loading and unloading

    stepStruct.Control.analysis = 'STATIC';
    stepStruct.Control.time_steps = numTimeSteps;
    stepStruct.Control.step_size = 1/numTimeSteps;

    stepStruct.Control.solver.max_refs = max_refs;
    stepStruct.Control.solver.qn_method.max_ups = max_ups;

    stepStruct.Control.time_stepper.dtmin = dtmin;
    stepStruct.Control.time_stepper.dtmax = dtmax;
    stepStruct.Control.time_stepper.max_retries = max_retries;
    stepStruct.Control.time_stepper.opt_iter = opt_iter;

    % Complete missing FEBio control settings
    stepStruct.Control = structComplete( ...
        stepStruct.Control,febio_spec.Control,1);

    % Remove the original global control section
    febio_spec = rmfield(febio_spec,'Control');

    % Step 1: loading
    febio_spec.Step.step{1}.ATTR.id = 1;
    febio_spec.Step.step{1}.ATTR.name = 'Loading';
    febio_spec.Step.step{1}.Control = stepStruct.Control;

    % Step 2: unloading
    febio_spec.Step.step{2}.ATTR.id = 2;
    febio_spec.Step.step{2}.ATTR.name = 'Unloading';
    febio_spec.Step.step{2}.Control = stepStruct.Control;

    %% Material definition corresponding to inverse analysis

    materialName1 = 'Material1';

    febio_spec.Material.material{1}.ATTR.id   = 1;
    febio_spec.Material.material{1}.ATTR.name = materialName1;
    febio_spec.Material.material{1}.ATTR.type = 'solid mixture';


    %% Component 1: damaged Ogden ground matrix

    febio_spec.Material.material{1}.solid{1}.ATTR.type = ...
        'uncoupled elastic damage';

    % Bulk modulus of uncoupled damage wrapper
    febio_spec.Material.material{1}.solid{1}.k = ...
        k_matrix_ini;

    % Elastic Ogden matrix
    febio_spec.Material.material{1}.solid{1}.elastic.ATTR.type = ...
        'Ogden';

    febio_spec.Material.material{1}.solid{1}.elastic.density = 1;
    febio_spec.Material.material{1}.solid{1}.elastic.c1 = c1_ini;
    febio_spec.Material.material{1}.solid{1}.elastic.m1 = m1_ini;
    febio_spec.Material.material{1}.solid{1}.elastic.k  = ...
        k_matrix_ini;

    % Matrix damage law
    febio_spec.Material.material{1}.solid{1}.damage.ATTR.type = ...
        'CDF log-normal';

    febio_spec.Material.material{1}.solid{1}.damage.mu = ...
        mu_matrix_ini;

    febio_spec.Material.material{1}.solid{1}.damage.sigma = ...
        sigma_matrix_ini;

    febio_spec.Material.material{1}.solid{1}.damage.Dmax = ...
        Dmax_matrix_ini;

    febio_spec.Material.material{1}.solid{1}.criterion.ATTR.type = ...
        'DC max normal Lagrange strain';


    %% Component 2: damaged HGO fibre family at +gamma

    febio_spec.Material.material{1}.solid{2}.ATTR.type = ...
        'elastic damage';

    % Elastic HGO material
    febio_spec.Material.material{1}.solid{2}.elastic.ATTR.type = ...
        'HGO unconstrained';

    febio_spec.Material.material{1}.solid{2}.elastic.c = ...
        cGOH_ini/2;

    febio_spec.Material.material{1}.solid{2}.elastic.k1 = ...
        k1_ini;

    febio_spec.Material.material{1}.solid{2}.elastic.k2 = ...
        k2_ini;

    febio_spec.Material.material{1}.solid{2}.elastic.kappa = ...
        kappa_ini;

    febio_spec.Material.material{1}.solid{2}.elastic.gamma = ...
        +gamma;

    febio_spec.Material.material{1}.solid{2}.elastic.k = ...
        kGOH_ini/2;

    % HGO damage law
    febio_spec.Material.material{1}.solid{2}.damage.ATTR.type = ...
        'CDF log-normal';

    febio_spec.Material.material{1}.solid{2}.damage.mu = ...
        mu_GOH_ini;

    febio_spec.Material.material{1}.solid{2}.damage.sigma = ...
        sigma_GOH_ini;

    febio_spec.Material.material{1}.solid{2}.damage.Dmax = ...
        Dmax_GOH_ini;

    febio_spec.Material.material{1}.solid{2}.criterion.ATTR.type = ...
        'DC max normal Lagrange strain';


    %% Component 3: damaged HGO fibre family at -gamma

    febio_spec.Material.material{1}.solid{3}.ATTR.type = ...
        'elastic damage';

    % Elastic HGO material
    febio_spec.Material.material{1}.solid{3}.elastic.ATTR.type = ...
        'HGO unconstrained';

    febio_spec.Material.material{1}.solid{3}.elastic.c = ...
        cGOH_ini/2;

    febio_spec.Material.material{1}.solid{3}.elastic.k1 = ...
        k1_ini;

    febio_spec.Material.material{1}.solid{3}.elastic.k2 = ...
        k2_ini;

    febio_spec.Material.material{1}.solid{3}.elastic.kappa = ...
        kappa_ini;

    febio_spec.Material.material{1}.solid{3}.elastic.gamma = ...
        -gamma;

    febio_spec.Material.material{1}.solid{3}.elastic.k = ...
        kGOH_ini/2;

    % Same HGO damage law as the +gamma family
    febio_spec.Material.material{1}.solid{3}.damage.ATTR.type = ...
        'CDF log-normal';

    febio_spec.Material.material{1}.solid{3}.damage.mu = ...
        mu_GOH_ini;

    febio_spec.Material.material{1}.solid{3}.damage.sigma = ...
        sigma_GOH_ini;

    febio_spec.Material.material{1}.solid{3}.damage.Dmax = ...
        Dmax_GOH_ini;

    febio_spec.Material.material{1}.solid{3}.criterion.ATTR.type = ...
        'DC max normal Lagrange strain';
    % Mesh section
    % -> Nodes

    %%Area of interest
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

    % -> Elements
    partName1='KirigamiGripper';
    febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type='penta6'; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

    % -> NodeSets
    nodeSetName1='bcPrescribeList1';
    nodeSetName2='bcSupportList2';
    nodeSetName3='bcSupportList3';

    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcPrescribeList);

    febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
    febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList);

    febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
    febio_spec.Mesh.NodeSet{3}.VAL=bcSupportList(1);

    %MeshDomains section
    febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
    febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

    %MeshData section
    % -> ElementData
    febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
    febio_spec.MeshData.ElementData{1}.ATTR.type='mat_axis';

    for q=1:1:size(E,1)
        febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
        febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir(q,:);
        febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir(q,:);
    end
    %Boundary condition section
    % -> Fix boundary conditions
    febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_xz';
    febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
    febio_spec.Boundary.bc{1}.x_dof=0;
    febio_spec.Boundary.bc{1}.y_dof=0;
    febio_spec.Boundary.bc{1}.z_dof=1;

    febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement';
    febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
    febio_spec.Boundary.bc{2}.x_dof=0;
    febio_spec.Boundary.bc{2}.y_dof=1;
    febio_spec.Boundary.bc{2}.z_dof=1;

    %% STEP 1: displacement-controlled loading

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

    %% STEP 2: unload the upper boundary toward zero force

    febio_spec.Step.step{2}.Loads.nodal_load{1}.ATTR.name = ...
        'zero_force_unloading';

    febio_spec.Step.step{2}.Loads.nodal_load{1}.ATTR.type = ...
        'nodal_target_force';

    febio_spec.Step.step{2}.Loads.nodal_load{1}.ATTR.node_set = ...
        nodeSetName1;

    febio_spec.Step.step{2}.Loads.nodal_load{1}.force = [0 0 0];

    febio_spec.Step.step{2}.Loads.nodal_load{1}.scale.ATTR.lc = 2;
    febio_spec.Step.step{2}.Loads.nodal_load{1}.scale.VAL = 1;

    febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement';
    febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
    febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
    febio_spec.Boundary.bc{3}.x_dof=1;
    febio_spec.Boundary.bc{3}.y_dof=0;
    febio_spec.Boundary.bc{3}.z_dof=0;


    %% Load controllers

    % Step 1: loading
    febio_spec.LoadData.load_controller{1}.ATTR.id = 1;
    febio_spec.LoadData.load_controller{1}.ATTR.name = 'LC_Loading';
    febio_spec.LoadData.load_controller{1}.ATTR.type = 'loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate = 'LINEAR';
    febio_spec.LoadData.load_controller{1}.extend = 'CONSTANT';
    febio_spec.LoadData.load_controller{1}.points.pt.VAL = ...
        [0 0; 1 1];

    % Step 2: unloading
    febio_spec.LoadData.load_controller{2}.ATTR.id = 2;
    febio_spec.LoadData.load_controller{2}.ATTR.name = 'LC_Unloading';
    febio_spec.LoadData.load_controller{2}.ATTR.type = 'loadcurve';
    febio_spec.LoadData.load_controller{2}.interpolate = 'LINEAR';
    febio_spec.LoadData.load_controller{2}.extend = 'CONSTANT';
    febio_spec.LoadData.load_controller{2}.points.pt.VAL = ...
        [1 0; 2 1];

    %Output section
    % -> log file
    logNamePrefix = sprintf('%s_ang_%03d',geometryTag,ang);
    febioLogFileName        = [logNamePrefix, '.txt'];
    febioLogFileName_disp   = [logNamePrefix, '_disp_out.txt'];
    febioLogFileName_stress_prin = [logNamePrefix, '_stress_prin_out.txt'];
    febioLogFileName_force  = [logNamePrefix, '_force_out.txt'];

    forceExportFile = fullfile( ...
        savePath, ...
        sprintf('%s_ang_%03d_Ry_xplt.csv',geometryTag,ang));

    febio_spec.Output.logfile.ATTR.file = febioLogFileName;

    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

    febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
    febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
    febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

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


    % Plotfile section
    febio_spec.Output.plotfile.compression=0;

    %% Quick viewing of the FEBio input file structure
    % The |febView| function can be used to view the xml structure in a MATLAB
    % figure window.

    %%
    %%|febView(febio_spec); %Viewing the febio file|

    %% Exporting the FEBio input file
    % Exporting the febio_spec structure to an FEBio input file is done using
    % the |febioStruct2xml| function.

    febFileName = fullfile(savePath, ...
        sprintf('%s_ang_%03d.feb',geometryTag,ang));

    if runFEBio

        % febio_spec.Step.step = febio_spec.Step.step(1);
        % Create the FEB file only during the simulation pass
        febioStruct2xml(febio_spec,febFileName);

        febioAnalysis.run_filename = febFileName;
        febioAnalysis.run_logname = ...
            fullfile(savePath,febioLogFileName);

        febioAnalysis.disp_on = 0;
        febioAnalysis.runMode = runMode;
        febioAnalysis.maxLogCheckTime = 300;

        runFlag = runMonitorFEBio(febioAnalysis);

        if runFlag ~= 1
            error('FEBio failed for %s at %d degrees.', ...
                geometryTag,ang);
        end

        fprintf(['FEBio completed for %s at %d degrees.\n' ...
            'Manually export Ry from the XPLT file as:\n%s\n'], ...
            geometryTag,ang,forceExportFile);

        % Continue to the next angle without reading a CSV
        continue

    else

        % Post-processing only: do not generate or overwrite FEBio files
        runFlag = 1;

    end

    %% Import FEBio results

    if runFlag==1 %i.e. a succesful run
        %%
        %% Import reaction forces manually exported from the XPLT file

        if ~isfile(forceExportFile)
            error(['Cannot find the manually exported force file:' newline ...
                '%s' newline ...
                'Export Ry from FEBio Studio using this exact filename.'], ...
                forceExportFile);
        end

        % FEBio Studio CSV export contains two header lines
        forceExport = readmatrix( ...
            forceExportFile, ...
            'FileType','text', ...
            'NumHeaderLines',2);

        % Remove completely empty rows and columns
        forceExport = forceExport( ...
            ~all(isnan(forceExport),2),:);

        forceExport = forceExport(:, ...
            ~all(isnan(forceExport),1));

        if size(forceExport,2) < 2
            error('The exported reaction-force CSV was not imported correctly.');
        end

        % Column 1 is time
        timeForceVec = forceExport(:,1);

        % Remaining columns are nodal y-reaction forces
        Ry_selected_nodes = forceExport(:,2:end);

        fprintf('Imported force time points: %d\n',numel(timeForceVec));
        fprintf('Imported nodal Ry curves: %d\n', ...
            size(Ry_selected_nodes,2));

        %%
        % Importing nodal displacements from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

        %Access data
        N_disp_mat=dataStruct.data; %Displacement

        %Create deformed coordinate set
        V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
        timeVec = dataStruct.time;

        %%
        %% True Poisson's Ratio vs. Axial Stretch

        % Reference dimensions
        w0 = max(V(:,1)) - min(V(:,1));  % Initial width (x-direction)
        h0 = max(V(:,2)) - min(V(:,2));  % Initial height (y-direction)
        A0 = w0 * h0;                    % Initial area

        nSteps = size(V_DEF,3);
        lambda_x_vec = zeros(nSteps,1);
        lambda_y_vec = zeros(nSteps,1);
        nu_true      = zeros(nSteps,1);
        J_vec        = zeros(nSteps,1);

        for i = 1:nSteps
            V_temp = V_DEF(:,:,i);

            w = max(V_temp(:,1)) - min(V_temp(:,1)); % Current width
            h = max(V_temp(:,2)) - min(V_temp(:,2)); % Current height
            A = w * h;

            lambda_x_vec(i) = w / w0;
            lambda_y_vec(i) = h / h0;
            J_vec(i)        = A / A0;%area ratio

            % True Poisson’s ratio (log-based definition)
            if lambda_y_vec(i) > 0 && lambda_x_vec(i) > 0
                nu_true(i) = -log(lambda_x_vec(i)) / log(lambda_y_vec(i));
            else
                nu_true(i) = NaN;
            end
        end

        %% Sum the XPLT nodal reactions over the selected boundary

        Fy_signed = sum(Ry_selected_nodes,2,'omitnan');

        % End of loading occurs approximately at time = 1
        [~,iPeakForce] = min(abs(timeForceVec(:)-1));

        % Display tensile reaction force as positive
        if Fy_signed(iPeakForce) < 0
            Fy_total = -Fy_signed;
        else
            Fy_total = Fy_signed;
        end

        fprintf('Peak total force: %.6g N\n',max(Fy_total));
        % fprintf('Final unloaded force: %.6g N\n',Fy_total(end));
        % %Unloading
        % fprintf('Final loading force: %.6g N\n',Fy_total(end));
        fprintf('Final unloaded force: %.6g N\n',Fy_total(end));
        % Match reaction force to the displacement-output time points
        Fy_vec = interp1( ...
            timeForceVec(:), ...
            Fy_total(:), ...
            timeVec(:), ...
            'linear','extrap');

        % Plot true Poisson’s ratio
        cFigure; hold on;
        plot(lambda_y_vec, nu_true, 'b-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('\nu_{\rm true} (Poisson''s Ratio)', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Poisson''s Ratio vs. Stretch  |  Fiber Direction = %s', angLabel), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;
        %% Area Ratio J vs. Axial Stretch

        cFigure; hold on;
        plot(lambda_y_vec, J_vec, 'k-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('J (Area Ratio)', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Area Ratio J vs. Axial Stretch  |  Fiber Direction = %s', angLabel), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;

        %% Loading-unloading force versus axial stretch

        cFigure; hold on;

        % Maximum stretch separates loading from unloading
        [~,iPeakStretch] = max(lambda_y_vec);

        % Loading: solid
        hLoad = plot( ...
            lambda_y_vec(1:iPeakStretch), ...
            Fy_vec(1:iPeakStretch), ...
            'k-', ...
            'LineWidth',5);

        % Unloading: dashed
        hUnload = plot( ...
            lambda_y_vec(iPeakStretch:end), ...
            Fy_vec(iPeakStretch:end), ...
            'k--', ...
            'LineWidth',5);

        xlabel('\lambda_y (Axial Stretch)', ...
            'Interpreter','tex','FontSize',20,'FontWeight','bold');

        ylabel('Reaction Force F_y [N]', ...
            'Interpreter','tex','FontSize',20,'FontWeight','bold');

        title(sprintf( ...
            'Loading-unloading response | %s | %s', ...
            geometryTag,angLabel), ...
            'FontSize',20,'FontWeight','bold');

        legend([hLoad,hUnload], ...
            {'Loading','Unloading'}, ...
            'Location','best');

        set(gca,'FontSize',16,'FontWeight','bold');
        axis square;
        box on;
        grid on;

        %%
        % Importing element stress from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress_prin),0,1);

        %Access data
        E_stress_mat=dataStruct.data;

        % Import strain energy density
        dataStruct = importFEBio_logfile(fullfile(savePath,febioLogFileName_energy),0,1);
        E_energy_mat = dataStruct.data;

        % Import damage
        dataStruct = importFEBio_logfile(fullfile(savePath,febioLogFileName_damage),0,1);
        E_damage_mat = dataStruct.data;

        % Convert element damage to vertex values
        [CVd] = faceToVertexMeasure(E, V, E_damage_mat(:,:,end));

        hf=cFigure;
        title(sprintf('Damage contour | Fiber Direction = %s', angLabel), 'FontSize', 20)

        hp=gpatch(Fb, V_DEF(:,:,end), CVd, 'none', 1, lineWidth1);

        for qp=1:numel(hp)
            hp(qp).FaceColor='interp';
        end

        axisGeom(gca,fontSize);
        colormap(parula);
        colorbar;

        Dmax_plot = max(E_damage_mat(:), [], 'omitnan');

        if isempty(Dmax_plot) || ~isfinite(Dmax_plot) || Dmax_plot <= 0
            Dmax_plot = 1e-6;
        end

        caxis([0 Dmax_plot]);

        view(140,30);
        camlight headlight;

        % Convert element results to one value per step (mean over elements)

        nSteps = size(E_energy_mat,3);

        energy_total_vec   = zeros(nSteps,1);
        energy_damage_vec  = zeros(nSteps,1);
        energy_intact_vec  = zeros(nSteps,1);
        damage_vec         = zeros(nSteps,1);
        energy_total_vec   = zeros(nSteps,1);
        energy_damage_vec  = zeros(nSteps,1);
        energy_intact_vec  = zeros(nSteps,1);
        damage_vec         = zeros(nSteps,1);
        damage_max_vec     = zeros(nSteps,1);
        damage_fraction_vec = zeros(nSteps,1);

        % Damage values near 1 are considered fully damaged
        fullDamageThreshold = 0.99;

        % Reference volume of each penta6 element
        % Decompose each straight-sided wedge into three tetrahedra

        V1e = V(E(:,1),:);
        V2e = V(E(:,2),:);
        V3e = V(E(:,3),:);
        V4e = V(E(:,4),:);
        V5e = V(E(:,5),:);
        V6e = V(E(:,6),:);

        tetVolume = @(A,B,C,D) ...
            abs(dot(B-A,cross(C-A,D-A,2),2))/6;

        elementVolume0 = ...
            tetVolume(V1e,V2e,V3e,V4e) + ...
            tetVolume(V2e,V3e,V4e,V5e) + ...
            tetVolume(V3e,V4e,V5e,V6e);

        % Force a column vector for logical indexing
        elementWeights = elementVolume0(:);

        % Mesh-volume validity checks
        if numel(elementWeights) ~= size(E,1)
            error(['Element-volume count does not match the number of ' ...
                'penta6 elements.']);
        end

        if any(~isfinite(elementWeights))
            error('Non-finite penta6 reference volumes were detected.');
        end

        if any(elementWeights <= 0)
            error('Zero or negative penta6 reference volumes were detected.');
        end

        totalVolume = sum(elementWeights);

        fprintf('Number of penta6 elements: %d\n',size(E,1));
        fprintf('Reference mesh volume: %.6g\n',totalVolume);
        fprintf('Minimum element volume: %.6g\n',min(elementWeights));
        fprintf('Maximum element volume: %.6g\n',max(elementWeights));

        % for i = 1:nSteps
        %
        %     E_step = E_energy_mat(:,1,i);
        %     D_step = E_damage_mat(:,1,i);
        %
        %     damaged = D_step > 0.01;
        %     intact  = ~damaged;
        %
        %     % Total strain energy
        %     energy_total_vec(i) = mean(E_step);
        %
        %     % Energy in damaged region
        %     if any(damaged)
        %         energy_damage_vec(i) = mean(E_step(damaged));
        %     else
        %         energy_damage_vec(i) = NaN;
        %     end
        %
        %     % Energy in intact region
        %     if any(intact)
        %         energy_intact_vec(i) = mean(E_step(intact));
        %     else
        %         energy_intact_vec(i) = NaN;
        %     end
        %
        %     damage_vec(i) = mean(D_step);
        %
        %     % --- total energy as SUM, like the first script ---
        %     energy_total_sum   = zeros(nSteps,1);
        %     energy_damage_sum  = zeros(nSteps,1);
        %     energy_intact_sum  = zeros(nSteps,1);
        %
        %     for k = 1:nSteps
        %         E_step = E_energy_mat(:,1,k);
        %         D_step = E_damage_mat(:,1,k);
        %
        %         damaged = D_step > 0.01;
        %         intact  = ~damaged;
        %
        %         energy_damage_sum(k) = sum(E_step(damaged));
        %         energy_intact_sum(k) = sum(E_step(intact));
        %         energy_total_sum(k)  = energy_damage_sum(k) + energy_intact_sum(k);
        %     end
        %
        %     % --- save results for final comparison plots ---
        %     idx = find(angles_deg == ang);
        %
        %     results(idx).lambda      = lambda_y_vec;
        %     results(idx).force       = Fy_vec;
        %     results(idx).energy_mean = energy_total_vec;
        %     results(idx).energy_sum  = energy_total_sum;
        %     results(idx).damage      = damage_vec;
        %
        % end
        for i = 1:nSteps

            E_step = E_energy_mat(:,1,i);
            D_step = E_damage_mat(:,1,i);

            E_step = E_step(:);
            D_step = D_step(:);

            if numel(D_step) ~= numel(elementWeights)
                error(['Damage-result count does not match the number of ' ...
                    'element-volume weights.']);
            end

            damaged = D_step > 0.01;
            intact  = ~damaged;

            % Total strain energy
            energy_total_vec(i) = mean(E_step);

            % Energy in damaged region
            if any(damaged)
                energy_damage_vec(i) = mean(E_step(damaged));
            else
                energy_damage_vec(i) = NaN;
            end

            % Energy in intact region
            if any(intact)
                energy_intact_vec(i) = mean(E_step(intact));
            else
                energy_intact_vec(i) = NaN;
            end

            damage_vec(i) = mean(D_step);
            damage_max_vec(i) = max(D_step,[],'omitnan');
            isFullyDamaged = D_step >= fullDamageThreshold;

            damage_fraction_vec(i) = ...
                sum(elementWeights(isFullyDamaged),'omitnan') / totalVolume;


        end

        % --- total energy as SUM, like the first script ---
        energy_total_sum   = zeros(nSteps,1);
        energy_damage_sum  = zeros(nSteps,1);
        energy_intact_sum  = zeros(nSteps,1);

        for k = 1:nSteps
            E_step = E_energy_mat(:,1,k);
            D_step = E_damage_mat(:,1,k);

            damaged = D_step > 0.01;
            intact  = ~damaged;

            energy_damage_sum(k) = sum(E_step(damaged));
            energy_intact_sum(k) = sum(E_step(intact));
            energy_total_sum(k)  = energy_damage_sum(k) + energy_intact_sum(k);
        end

        % --- save results for final comparison plots ---

        results(ia).lambda      = lambda_y_vec;
        results(ia).force       = Fy_vec;
        results(ia).energy_mean = energy_total_vec;
        results(ia).energy_sum  = energy_total_sum;
        results(ia).damage      = damage_vec;
        results(ia).damage_max = damage_max_vec;
        results(ia).damage_fraction = damage_fraction_vec;

        if length(lambda_y_vec) ~= length(energy_total_vec)
            error('Size mismatch between stretch and energy results');
        end

        if length(lambda_y_vec) ~= length(damage_vec)
            error('Size mismatch between stretch and damage results');
        end

        E_stress_mat_VM=sqrt(( (E_stress_mat(:,1,:)-E_stress_mat(:,2,:)).^2 + ...
            (E_stress_mat(:,2,:)-E_stress_mat(:,3,:)).^2 + ...
            (E_stress_mat(:,1,:)-E_stress_mat(:,3,:)).^2  )/2); %Von Mises stress

        cFigure; hold on;
        plot(lambda_y_vec, energy_total_vec, 'k-', 'LineWidth', 5);
        plot(lambda_y_vec, energy_damage_vec, 'r-', 'LineWidth', 5);
        plot(lambda_y_vec, energy_intact_vec, 'b--', 'LineWidth', 5);
        legend('Total energy','Energy in damaged region','Energy in intact region', ...
            'Location','best')

        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('Strain Energy Density', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

        title(sprintf('Strain Energy Components | Fiber Direction = %s', angLabel), ...
            'FontSize',20,'FontWeight','bold')

        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;


        cFigure; hold on;
        plot(lambda_y_vec, damage_vec, 'r-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('Damage', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Damage vs Stretch  |  Fiber Direction = %s', angLabel), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;

        % ===== DAMAGE CONTOUR ANIMATION =====

        hf=cFigure;
        title(sprintf('Damage evolution | Fiber Direction = %s', angLabel), 'FontSize', 20)

        [CVd]=faceToVertexMeasure(E,V,E_damage_mat(:,:,end));

        hp=gpatch(Fb,V_DEF(:,:,end),CVd,'none',1,lineWidth1);

        for qp=1:numel(hp)
            hp(qp).FaceColor='interp';
        end

        axisGeom(gca,fontSize);
        colormap(parula);
        colorbar;
        Dmax_plot = max(E_damage_mat(:), [], 'omitnan');

        if isempty(Dmax_plot) || ~isfinite(Dmax_plot) || Dmax_plot <= 0
            Dmax_plot = 1e-6;
        end

        caxis([0 Dmax_plot]);

        view(140,30);
        camlight headlight;

        animStruct.Time=timeVec;

        for qt=1:size(N_disp_mat,3)

            [CVd]=faceToVertexMeasure(E,V,E_damage_mat(:,:,qt));

            animStruct.Handles{qt}=[hp(1) hp(1) hp(2) hp(2)];
            animStruct.Props{qt}={'Vertices','CData','Vertices','CData'};
            animStruct.Set{qt}={V_DEF(:,:,qt),CVd,V_DEF(:,:,qt),CVd};

        end

        anim8(hf,animStruct);

        %%
        % Compute element area (only once)
        A = patchArea(E, V_DEF(:,:,1));
        A_tot = sum(A);

        % Storage for p90 values for each time step
        p90_vec = zeros(size(N_disp_mat,3),1);

        %%
        % Plotting the simulated results using |anim8| to visualize and animate
        % deformations

        [CV]=faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,end));

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy

        gtitle([febioFebFileNamePart,': Press play to animate']);
        title('$\sigma_{vm}$ [MPa]','Interpreter','Latex')
        hp=gpatch(Fb,V_DEF(:,:,end),CV,'none',1,lineWidth1); %Add graphics object to animate

        for qp=1:1:numel(hp) %For all graphics objects e.g. triangles/quads
            %         hp(qp).Marker=".";
            %         hp(qp).MarkerSize=markerSize2;
            hp(qp).FaceColor='interp';
        end

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(E_stress_mat_VM(:))]);
        axis(axisLim(V_DEF)); %Set axis limits statically
        view(140,30);
        camlight headlight;

        %     % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:size(N_disp_mat,3) %Loop over time increments

            % Compute 90% stress (area-weighted)
            sigmaVM = E_stress_mat_VM(:,:,qt);
            Q = (sigmaVM(:) .* A) ./ A_tot;
            p90_vec(qt) = prctile(Q,90);

            [CV]=faceToVertexMeasure(E,V,E_stress_mat_VM(:,:,qt));

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp(1) hp(1) hp(2) hp(2)]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        end
        anim8(hf,animStruct); %Initiate animation feature
        drawnow;

        figure; hold on;
        plot(lambda_y_vec, p90_vec, 'r-o', 'LineWidth', 3, 'MarkerFaceColor', 'r');
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('90% Stress (Normalized)', 'FontSize', 18, 'FontWeight', 'bold');
        title(sprintf('p_{90} Stress vs Axial Stretch  |  Fiber Direction = %s', angLabel), 'FontSize', 20, 'FontWeight', 'bold');
        grid on; axis square;
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');

        hf = cFigure;
        tiledlayout(1,2);

        % === Left: Max Stress Animation ===
        ax1 = nexttile;
        title(ax1, ...
            sprintf('Max Stress $\\sigma_{vm}$ [MPa] | Fiber Direction = %s | Elements = %d', ...
            angLabel, size(E,1)), ...
            'Interpreter','latex','FontSize',18);
        hp1 = gpatch(Fb, V_DEF(:,:,end), CV, 'none', 1, lineWidth1);
        axisGeom; colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(E_stress_mat_VM(:))]);
        view(140,30); camlight headlight;

        % === Right: 90% Stress Animation ===
        ax2 = nexttile;
        title(ax2, ...
            sprintf('90%% Stress $\\sigma_{vm}$ [MPa] | Fiber Direction = %s | Elements = %d', ...
            angLabel, size(E,1)), ...
            'Interpreter','latex','FontSize',18);
        hp2 = gpatch(Fb, V_DEF(:,:,end), CV, 'none', 1, lineWidth1);
        axisGeom; colormap(cMap); colorbar;
        caxis([min(E_stress_mat_VM(:)) max(p90_vec)]);
        view(140,30); camlight headlight;




        %%
        % Plotting the simulated results using |anim8| to visualize and animate
        % deformations

        CV=abs(N_disp_mat(:,1,end)); %Current displacement magnitude

        % Create basic view and store graphics handle to initiate animation
        % hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy

        % gtitle([febioFebFileNamePart,': Press play to animate']);
        % title(sprintf('$u_x$ [mm]  |  Fiber Angle = %d°', ang), 'Interpreter', 'latex');
        % hp=gpatch(Fb,V_DEF(:,:,end),CV,'none',1,lineWidth1); %Add graphics object to animate
        %
        % for qp=1:1:numel(hp) %For all graphics objects e.g. triangles/quads
        %     %         hp(qp).Marker=".";
        %     %         hp(qp).MarkerSize=markerSize2;
        %     hp(qp).FaceColor='interp';
        % end
        %
        % axisGeom(gca,fontSize);
        % colormap(cMap); colorbar;
        % caxis([0     2.0668]);
        % % caxis([0 max(CV)]);
        % % axis(axisLim(V_DEF)); %Set axis limits statically
        % axis([-5 30 -5 90 -10 10])
        % view(140,30);
        % camlight headlight;
        %
        % % Set up animation features
        % animStruct.Time=timeVec; %The time vector
        % for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        %
        %     CV=abs(N_disp_mat(:,1,qt)); %Current displacement magnitude
        %
        %     %Set entries in animation structure
        %     animStruct.Handles{qt}=[hp(1) hp(1) hp(2) hp(2)]; %Handles of objects to animate
        %     animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
        %     animStruct.Set{qt}={V_DEF(:,:,qt),CV,V_DEF(:,:,qt),CV}; %Property values for to set in order to animate
        % end
        % anim8(hf,animStruct); %Initiate animation feature
        % drawnow;
    end

end

%% Final comparison plots: 0° vs 90°

if all(arrayfun(@(s) ~isempty(s.lambda), results))

    plotComparison(results, 'force', ...
        'Reaction Force F_y', ...
        'Force vs Axial Stretch (Perpendicular vs Parallel)');

    plotComparison(results, 'energy_sum', ...
        'Total Strain Energy', ...
        'Energy vs Axial Stretch (Perpendicular vs Parallel)');

    plotComparison(results, 'damage', ...
        'Damage', ...
        'Damage vs Axial Stretch (Perpendicular vs Parallel)');
    plotComparison(results, 'damage_max', ...
        'Maximum damage D', ...
        'Maximum Damage vs Axial Stretch (Perpendicular vs Parallel)');
    plotComparison(results, 'damage_fraction', ...
        'Fully damaged volume fraction', ...
        'Spatial Spread of Full Damage (Perpendicular vs Parallel)');

    % Save both orientation results for this geometry
    resultsFile = fullfile( ...
        savePath, ...
        sprintf('%s_processed_results.mat',geometryTag));

    save(resultsFile, ...
        'results', ...
        'geometryID', ...
        'geometryLabel', ...
        'geometryTag', ...
        'angles_deg', ...
        '-v7.3');

    fprintf('Processed results saved as:\n%s\n',resultsFile);

else
    warning('Skipping comparison plots because one angle failed.');
end
%%
function [FT,VT,CT]=replicatemesh(F,V,p,d)
numcopies=prod(p);
numfaces=length(F);
numvertices=length(V);
C=ones(numfaces,1);
FT=repmat(F,numcopies,1);
VT=repmat(V,numcopies,1);
CT=repmat(C,numcopies,1);
c=1;
cf=1;
cv=1;
indexoffset=0;
for i=1:1:p(1)
    for j=1:1:p(2)
        FT(cf:cf+numfaces-1,:)=FT(cf:cf+numfaces-1,:)+indexoffset;
        VT(cv:cv+numvertices-1,:)=VT(cv:cv+numvertices-1,:)+[(i-1)*d(1),(j-1)*d(2)];
        CT(cf:cf+numfaces-1,:)=c;
        indexoffset=indexoffset+numvertices;
        c=c+1;
        cf=cf+numfaces;
        cv=cv+numvertices;
    end
end

[FT,VT]=mergeVertices(FT,VT);

end

function plotComparison(results, fieldName, yLabelText, titleText)

cFigure; hold on;

% 0 degrees: Perpendicular
lambdaPerp = results(1).lambda(:);
dataPerp   = results(1).(fieldName);
dataPerp   = dataPerp(:);

% 90 degrees: Parallel
lambdaPara = results(2).lambda(:);
dataPara   = results(2).(fieldName);
dataPara   = dataPara(:);

% Separate loading and unloading at maximum stretch
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
    'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

ylabel(yLabelText, ...
    'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

title(titleText, ...
    'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

legend([h1 h2 h3 h4], ...
    {'Perpendicular loading', ...
    'Perpendicular unloading', ...
    'Parallel loading', ...
    'Parallel unloading'}, ...
    'FontSize',14, ...
    'Location','best');

set(gca, 'FontSize', 16, 'FontWeight', 'bold');

lambdaAll = [results(1).lambda(:); results(2).lambda(:)];
xlim([1 max(lambdaAll)]);

axis square;
box on;
grid on;

end

%%
% %% Select geometry: use 1, 3, 4, or 5
% geometryID = 5;
% 
% switch geometryID
%     case 1
%         geometryLabel = 'AS1';
%     case 3
%         geometryLabel = 'G3';
%     case 4
%         geometryLabel = 'G4';
%     case 5
%         geometryLabel = 'HS';
%     otherwise
%         error('geometryID must be 1, 3, 4, or 5.');
% end
% 
% geometryTag = sprintf('G%02d_%s',geometryID,geometryLabel);
% 
% Locate the folder containing auxetic_patterns_unload.m
% codeFolder = fileparts(which('auxetic_patterns_unload'));
% 
% if isempty(codeFolder)
%     error(['MATLAB cannot locate auxetic_patterns_unload.m. ' ...
%         'Open its folder or add that folder to the MATLAB path.']);
% end
% 
% savePath = fullfile(codeFolder,'data',geometryTag);
% 
% filePerp = fullfile(savePath, ...
%     sprintf('%s_ang_000_Ry_xplt.csv',geometryTag));
% 
% filePara = fullfile(savePath, ...
%     sprintf('%s_ang_090_Ry_xplt.csv',geometryTag));
% 
% if ~isfile(filePerp)
%     error('Cannot find perpendicular CSV:\n%s',filePerp);
% end
% 
% if ~isfile(filePara)
%     error('Cannot find parallel CSV:\n%s',filePara);
% end
% 
% % Read perpendicular result: angle 0 degrees
% 
% A0 = readmatrix( ...
%     filePerp, ...
%     'FileType','text', ...
%     'NumHeaderLines',2);
% 
% A0 = A0(~all(isnan(A0),2),:);
% A0 = A0(:,~all(isnan(A0),1));
% 
% t0 = A0(:,1);
% F0 = sum(A0(:,2:end),2,'omitnan');
% 
% End of loading is at approximately time = 1
% [~,i0] = min(abs(t0-1));
% 
% Display tensile force as positive
% if F0(i0) < 0
%     F0 = -F0;
% end
% 
% % Read parallel result: angle 90 degrees
% 
% A90 = readmatrix( ...
%     filePara, ...
%     'FileType','text', ...
%     'NumHeaderLines',2);
% 
% A90 = A90(~all(isnan(A90),2),:);
% A90 = A90(:,~all(isnan(A90),1));
% 
% t90 = A90(:,1);
% F90 = sum(A90(:,2:end),2,'omitnan');
% 
% End of loading is at approximately time = 1
% [~,i90] = min(abs(t90-1));
% 
% Display tensile force as positive
% if F90(i90) < 0
%     F90 = -F90;
% end
% 
% % Plot the two orientations
% 
% figure;
% hold on;
% 
% Perpendicular, angle 0 degrees: blue
% h1 = plot( ...
%     t0(1:i0),F0(1:i0), ...
%     '-','Color',[0.15 0.35 0.85], ...
%     'LineWidth',4);
% 
% h2 = plot( ...
%     t0(i0:end),F0(i0:end), ...
%     '--','Color',[0.15 0.35 0.85], ...
%     'LineWidth',4);
% 
% Parallel, angle 90 degrees: red
% h3 = plot( ...
%     t90(1:i90),F90(1:i90), ...
%     '-','Color',[0.85 0.15 0.15], ...
%     'LineWidth',4);
% 
% h4 = plot( ...
%     t90(i90:end),F90(i90:end), ...
%     '--','Color',[0.85 0.15 0.15], ...
%     'LineWidth',4);
% 
% xlabel('Time');
% ylabel('Reaction force F_y [N]');
% 
% title(sprintf( ...
%     'Reaction Force History: %s',geometryTag), ...
%     'Interpreter','none');
% 
% legend([h1 h2 h3 h4], ...
%     {'Perpendicular loading', ...
%      'Perpendicular unloading', ...
%      'Parallel loading', ...
%      'Parallel unloading'}, ...
%     'Location','best');
% 
% xline(1,'k:','End of loading', ...
%     'LabelVerticalAlignment','bottom');
% 
% grid on;
% box on;
% axis square;
% 
% set(gca, ...
%     'FontSize',16, ...
%     'FontWeight','bold');
% 
% fprintf('\nGeometry plotted: %s\n',geometryTag);
% fprintf('Perpendicular peak force: %.6g N\n',max(F0));
% fprintf('Parallel peak force: %.6g N\n',max(F90));
% fprintf('Perpendicular final force: %.6g N\n',F0(end));
% fprintf('Parallel final force: %.6g N\n',F90(end));