%% Human-back_damage

clear; close all; clc;

%% GIBBON path + check FEBio import functions
gibbonFolder = 'C:\Users\G00425456@atu.ie\Downloads\GIBBON-master (1)\GIBBON-master\docs\1damage-2goh';
addpath(genpath(gibbonFolder));

which importFEBio_xplt -all
which importFEBio_plotfile -all
which importFEBio_logfile -all

%%
% Plot settings
fontSize=20;
faceAlpha1=0.8;
faceAlpha2=1;
edgeColor=0.25*ones(1,3);
edgeWidth=1.5;
markerSize1=15;
markerSize2=40;
markerSize3=50;
lineWidth1=3;
lineWidth2=3;
cMap=viridis(20);
numPlotPointsGraphs = 100;

%% Control parameters

% Path names
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
savePath = fullfile(pwd,'tempFEBio');
if ~exist(savePath,'dir')
    mkdir(savePath);
end


% Defining file names
febioFebFileNamePart='tempModel_twoDamage_03';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; %Log file name for exporting stress sigma_z
febioLogFileName_stretch=[febioFebFileNamePart,'_stretch_out.txt']; %Log file name for exporting stretch U_z

% Define data paths
loadPath_experimental = fullfile(fileparts(defaultFolder),'data','Ni_Annaidh_2012');
dataName_1 = fullfile(loadPath_experimental,'HumanBack_stress_stretch_perp.mat');
dataName_2 = fullfile(loadPath_experimental,'HumanBack_stress_stretch_para.mat');

%Specifying dimensions and number of elements
sampleWidth=10;
sampleThickness=10;
sampleHeight=10;
pointSpacings=10*ones(1,3);
initialArea=sampleWidth*sampleThickness;

numElementsWidth=round(sampleWidth/pointSpacings(1));
numElementsThickness=round(sampleThickness/pointSpacings(2));
numElementsHeight=round(sampleHeight/pointSpacings(3));


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
mu_matrix_ini    = 1.0222386595083708e+00;
sigma_matrix_ini = 1.4169480827757489e-01;
Dmax_matrix_ini  = 1;

% GOH damage
mu_GOH_ini       = 5.5063703735016223e-01;
sigma_GOH_ini    = 4.7249662329555610e-02;
Dmax_GOH_ini     = 1;

% Evaluation mode and parameter vector
evalMode = 2;  % 1 = test, 2 = optimise

% Optimised parameters
P = [mu_matrix_ini, sigma_matrix_ini, ...
    mu_GOH_ini,    sigma_GOH_ini];


LangerAngle_perp = 0/180*pi; % Angle of Langer line with loading direction
LangerAngle_para = 90/180*pi; % Angle of Langer line with loading direction

LangerAngles = [LangerAngle_perp, LangerAngle_para];

% FEA control settings
numTimeSteps=40; %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=6; %Optimum number of iterations
max_retries=5; %Maximum number of retires
dtmin=(1/numTimeSteps)/100; %Minimum time step size
dtmax=1/numTimeSteps; %Maximum time step size

runMode='external';% 'internal' or 'external'

%Optimisation settings
maxNumberIterations=100; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-6; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';
optimisationMethod = 2; % 1= fminsearch, Nelder-Mead, 2=lsqnonlin, Levenberg-Marquart

%% LOAD EXPERIMENTAL DATA

data_perp = load(dataName_1);
data_para = load(dataName_2);

stretch_exp_perp_raw = data_perp.X1;
stress_exp_perp_raw = data_perp.Y1;

stretch_exp_para_raw = data_para.X1;
stress_exp_para_raw = data_para.Y1;


%% Different peak-based initial guesses

% Locate each experimental stress peak
[~,indexPeak_perp] = max(stress_exp_perp_raw);
[~,indexPeak_para] = max(stress_exp_para_raw);

% Stretch at each peak
lambdaPeak_perp = stretch_exp_perp_raw(indexPeak_perp);
lambdaPeak_para = stretch_exp_para_raw(indexPeak_para);

% Convert peak stretch to Lagrange strain
% Used only as an initial optimizer estimate
mu_matrix_ini = 0.5*(lambdaPeak_perp^2 - 1);
mu_GOH_ini    = 0.5*(lambdaPeak_para^2 - 1);

% Rebuild P using the new, different initial values
P = [mu_matrix_ini, sigma_matrix_ini, ...
     mu_GOH_ini,    sigma_GOH_ini];

fprintf('Initial matrix mu = %.4f\n',mu_matrix_ini);
fprintf('Initial GOH mu    = %.4f\n',mu_GOH_ini);
%%
% Set stretch levels for simulation
stretchLoad_perp = max(stretch_exp_perp_raw);
stretchLoad_para = max(stretch_exp_para_raw);

displacementMagnitude_perp=(stretchLoad_perp*sampleHeight)-sampleHeight;
displacementMagnitude_para=(stretchLoad_para*sampleHeight)-sampleHeight;

displacementMagnitudes = [displacementMagnitude_perp, displacementMagnitude_para];

%%

stretch_plot_perp = linspace(1,stretchLoad_perp,numPlotPointsGraphs);
stretch_plot_para = linspace(1,stretchLoad_para,numPlotPointsGraphs);

stress_plot_perp = interp1(stretch_exp_perp_raw,stress_exp_perp_raw,stretch_plot_perp,'pchip','extrap'); %ppval(pp_perp,stretch_plot_perp);
stress_plot_para = interp1(stretch_exp_para_raw,stress_exp_para_raw,stretch_plot_para,'pchip','extrap'); %ppval(pp_para,stretch_plot_para);

%% ---------------- STRESS vs STRETCH ----------------
hf1 = cFigure; hold on;

title('Stretch stress curves optimised','FontSize',fontSize)
xlabel('\lambda Stretch [.]','FontSize',fontSize)
ylabel('\sigma Cauchy stress [MPa]','FontSize',fontSize)

Hn(1) = plot(stretch_exp_perp_raw,stress_exp_perp_raw,'r.','MarkerSize',markerSize2);
Hn(2) = plot(stretch_plot_perp,stress_plot_perp,'r-','LineWidth',lineWidth1);

Hn(3) = plot(stretch_exp_para_raw,stress_exp_para_raw,'g.','MarkerSize',markerSize2);
Hn(4) = plot(stretch_plot_para,stress_plot_para,'g-','LineWidth',lineWidth1);

Hn(5) = plot(NaN,NaN,'ko','MarkerSize',markerSize1,'LineWidth',lineWidth2,'MarkerFaceColor','r');
Hn(6) = plot(NaN,NaN,'ko','MarkerSize',markerSize1,'LineWidth',lineWidth2,'MarkerFaceColor','g');

Hn(7) = plot(NaN,NaN,'rx','MarkerSize',markerSize1);
Hn(8) = plot(NaN,NaN,'gx','MarkerSize',markerSize1);

legend(Hn([1 2 3 4 5 6]), ...
    {'HumanBack - perp','interp - perp', ...
    'HumanBack - para','interp - para', ...
    'FEA perp','FEA para'}, ...
    'Location','northwest')
grid on
box on
axis square
axis tight
set(gca,'FontSize',fontSize,'LineWidth',1.5)
view(2)

drawnow


%% ---------------- STRAIN ENERGY vs STRETCH ----------------
hf2 = cFigure; hold on;

title('Strain energy vs stretch','FontSize',fontSize)
xlabel('\lambda Stretch','FontSize',fontSize)
ylabel('Strain energy density','FontSize',fontSize)

He(1) = plot(NaN,NaN,'r-','LineWidth',lineWidth1);
He(2) = plot(NaN,NaN,'g-','LineWidth',lineWidth1);
legend(He,{'perp','para'},'Location','northwest')

grid on
box on
axis square
axis tight
set(gca,'FontSize',fontSize,'LineWidth',1.5)

drawnow


%% ---------------- DAMAGE vs STRETCH ----------------
hf3 = cFigure; hold on;

title('Damage vs stretch','FontSize',fontSize)
xlabel('\lambda Stretch','FontSize',fontSize)
ylabel('Damage','FontSize',fontSize)

Hd(1) = plot(NaN,NaN,'r-','LineWidth',lineWidth1);
Hd(2) = plot(NaN,NaN,'g-','LineWidth',lineWidth1);
legend(Hd,{'perp','para'},'Location','northwest')

grid on
box on
axis square
axis tight
set(gca,'FontSize',fontSize,'LineWidth',1.5)

drawnow
%% CREATING MESHED BOX

%Create box 1
boxDim=[sampleWidth sampleThickness sampleHeight]; %Dimensions
boxEl=[numElementsWidth numElementsThickness numElementsHeight]; %Number of elements
[box1]=hexMeshBox(boxDim,boxEl);
E=box1.E;
V=box1.V;
Fb=box1.Fb;
faceBoundaryMarker=box1.faceBoundaryMarker;

X=V(:,1); Y=V(:,2); Z=V(:,3);
VE=[mean(X(E),2) mean(Y(E),2) mean(Z(E),2)];

elementMaterialIndices=ones(size(E,1),1);

%%

% Plotting boundary surfaces

cFigure; hold on;
title('Model surfaces','FontSize',fontSize);
gpatch(Fb,V,faceBoundaryMarker,'k',0.5);
colormap(gjet(6)); icolorbar;
axisGeom(gca,fontSize);
drawnow;

%% DEFINE BC's

%Define supported node sets
logicFace=faceBoundaryMarker==1;
Fr=Fb(logicFace,:);
bcSupportList_X=unique(Fr(:));

logicFace=faceBoundaryMarker==3;
Fr=Fb(logicFace,:);
bcSupportList_Y=unique(Fr(:));

logicFace=faceBoundaryMarker==5;
Fr=Fb(logicFace,:);
bcSupportList_Z=unique(Fr(:));

%Prescribed displacement nodes
logicPrescribe=faceBoundaryMarker==6;
Fr=Fb(logicPrescribe,:);
bcPrescribeList=unique(Fr(:));

%%
% Visualize BC's
cFigure; hold on;
title('Complete model','FontSize',fontSize);

gpatch(Fb,V,'kw','k',0.5);
plotV(V(bcSupportList_X,:),'r.','MarkerSize',markerSize1);
plotV(V(bcSupportList_Y,:),'g.','MarkerSize',markerSize1);
plotV(V(bcSupportList_Z,:),'b.','MarkerSize',markerSize1);
plotV(V(bcPrescribeList,:),'k.','MarkerSize',markerSize1);

axisGeom(gca,fontSize);
drawnow;

%% DEFINE FIBRE DIRECTIONS

R = euler2DCM([0,-LangerAngle_perp,0]);
e1 = (R*[1 0 0]')';
e1_dir_perp = e1(ones(size(E,1),1),:);

e2 = [0 1 0];
e2_dir_perp = e2(ones(size(E,1),1),:);
e3_dir_perp = cross(e1_dir_perp,e2_dir_perp,2);

R = euler2DCM([0,-LangerAngle_para,0]);
e1 = (R*[1 0 0]')';
e1_dir_para = e1(ones(size(E,1),1),:);

e2 = [0 1 0];
e2_dir_para = e2(ones(size(E,1),1),:);
e3_dir_para = cross(e1_dir_para,e2_dir_para,2);

%%
% Visualizing material directions

[VE]=patchCentre(E,V);
gamma_rad = gamma*pi/180; % Convert gamma to radians

hf=cFigure;
subplot(1,2,1);  hold on;
title('Material directions - perp','FontSize',fontSize);

gpatch(Fb,V,'kw','none',0.25);
hf(1)=quiverVec(VE,e1_dir_perp,mean(pointSpacings),'r');
hf(2)=quiverVec(VE,e2_dir_perp,mean(pointSpacings)/2,'g');
hf(3)=quiverVec(VE,e3_dir_perp,mean(pointSpacings)/2,'b');

% --- fiber directions ---
R_pos = euler2DCM([0, -LangerAngle_perp + gamma_rad, 0]);
R_neg = euler2DCM([0, -LangerAngle_perp - gamma_rad, 0]);

fiber1_dir_perp = repmat((R_pos*[1 0 0]')', size(E,1), 1);
fiber2_dir_perp = repmat((R_neg*[1 0 0]')', size(E,1), 1);

hf(4)=quiverVec(VE, fiber1_dir_perp, mean(pointSpacings), 'm'); % +gamma
hf(5)=quiverVec(VE, fiber2_dir_perp, mean(pointSpacings), 'c'); % -gamma

legend(hf,{'e1-direction','e2-direction','e3-direction','fiber +\gamma','fiber -\gamma'});
axisGeom(gca,fontSize);
camlight headlight;

subplot(1,2,2);  hold on;
title('Material directions - para','FontSize',fontSize);

gpatch(Fb,V,'kw','none',0.25);
hf(1)=quiverVec(VE,e1_dir_para,mean(pointSpacings),'r');
hf(2)=quiverVec(VE,e2_dir_para,mean(pointSpacings)/2,'g');
hf(3)=quiverVec(VE,e3_dir_para,mean(pointSpacings)/2,'b');

R_pos_para = euler2DCM([0, -LangerAngle_para + gamma_rad, 0]);
R_neg_para = euler2DCM([0, -LangerAngle_para - gamma_rad, 0]);

fiber1_dir_para = repmat((R_pos_para*[1 0 0]')', size(E,1), 1);
fiber2_dir_para = repmat((R_neg_para*[1 0 0]')', size(E,1), 1);

hf(4)=quiverVec(VE, fiber1_dir_para, mean(pointSpacings), 'm'); % +gamma
hf(5)=quiverVec(VE, fiber2_dir_para, mean(pointSpacings), 'c'); % -gamma

legend(hf,{'e1-direction (fiber)','e2-direction','e3-direction','fiber +\gamma','fiber -\gamma'});
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

%Control section
febio_spec.Control.analysis='STATIC';
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=1/numTimeSteps;
febio_spec.Control.solver.max_refs=max_refs;
febio_spec.Control.solver.qn_method.max_ups=max_ups;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;


%% FEBio material definition: elastic damage + solid mixture

materialName1 = 'Material1';

febio_spec.Material.material{1}.ATTR.id   = 1;
febio_spec.Material.material{1}.ATTR.name = materialName1;

% The outer material is now the solid mixture
febio_spec.Material.material{1}.ATTR.type = 'solid mixture';

%% Component 1: damaged Ogden ground matrix

febio_spec.Material.material{1}.solid{1}.ATTR.type = 'uncoupled elastic damage';
febio_spec.Material.material{1}.solid{1}.k = k_ini;

% Elastic Ogden material
febio_spec.Material.material{1}.solid{1}.elastic.ATTR.type = 'Ogden';
febio_spec.Material.material{1}.solid{1}.elastic.density = 1;
febio_spec.Material.material{1}.solid{1}.elastic.c1 = c1_ini;
febio_spec.Material.material{1}.solid{1}.elastic.m1 = m1_ini;
febio_spec.Material.material{1}.solid{1}.elastic.k  = k_ini;

% Matrix damage law
febio_spec.Material.material{1}.solid{1}.damage.ATTR.type = 'CDF log-normal';

febio_spec.Material.material{1}.solid{1}.damage.mu = mu_matrix_ini;

febio_spec.Material.material{1}.solid{1}.damage.sigma = sigma_matrix_ini;

febio_spec.Material.material{1}.solid{1}.damage.Dmax = Dmax_matrix_ini;

febio_spec.Material.material{1}.solid{1}.criterion.ATTR.type = 'DC max normal Lagrange strain';


%% Component 2: damaged GOH, +gamma

febio_spec.Material.material{1}.solid{2}.ATTR.type = 'elastic damage';

% Elastic GOH material
febio_spec.Material.material{1}.solid{2}.elastic.ATTR.type = 'HGO unconstrained';

febio_spec.Material.material{1}.solid{2}.elastic.c = cGOH_ini/2;

febio_spec.Material.material{1}.solid{2}.elastic.k1 = k1_ini;
febio_spec.Material.material{1}.solid{2}.elastic.k2 = k2_ini;
febio_spec.Material.material{1}.solid{2}.elastic.kappa = kappa_ini;
febio_spec.Material.material{1}.solid{2}.elastic.gamma = +gamma;
febio_spec.Material.material{1}.solid{2}.elastic.k = kGOH_ini/2;

% GOH damage law
febio_spec.Material.material{1}.solid{2}.damage.ATTR.type ='CDF log-normal';

febio_spec.Material.material{1}.solid{2}.damage.mu = mu_GOH_ini;
febio_spec.Material.material{1}.solid{2}.damage.sigma = sigma_GOH_ini;
febio_spec.Material.material{1}.solid{2}.damage.Dmax = Dmax_GOH_ini;

febio_spec.Material.material{1}.solid{2}.criterion.ATTR.type = 'DC max normal Lagrange strain';


%% Component 3: damaged GOH, -gamma

febio_spec.Material.material{1}.solid{3}.ATTR.type = 'elastic damage';

% Elastic GOH material
febio_spec.Material.material{1}.solid{3}.elastic.ATTR.type = 'HGO unconstrained';

febio_spec.Material.material{1}.solid{3}.elastic.c = cGOH_ini/2;

febio_spec.Material.material{1}.solid{3}.elastic.k1 = k1_ini;
febio_spec.Material.material{1}.solid{3}.elastic.k2 = k2_ini;
febio_spec.Material.material{1}.solid{3}.elastic.kappa = kappa_ini;
febio_spec.Material.material{1}.solid{3}.elastic.gamma = -gamma;
febio_spec.Material.material{1}.solid{3}.elastic.k = kGOH_ini/2;

% Same GOH damage parameters as component 2
febio_spec.Material.material{1}.solid{3}.damage.ATTR.type = 'CDF log-normal';

febio_spec.Material.material{1}.solid{3}.damage.mu = mu_GOH_ini;
febio_spec.Material.material{1}.solid{3}.damage.sigma = sigma_GOH_ini;
febio_spec.Material.material{1}.solid{3}.damage.Dmax = Dmax_GOH_ini;

febio_spec.Material.material{1}.solid{3}.criterion.ATTR.type = 'DC max normal Lagrange strain';

% Mesh section
% -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
partName1='Part1';
febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
febio_spec.Mesh.Elements{1}.ATTR.type='hex8'; %Element type
febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

% -> NodeSets
nodeSetName1='bcSupportList_X';
nodeSetName2='bcSupportList_Y';
nodeSetName3='bcSupportList_Z';
nodeSetName4='bcPrescribeList';

febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcSupportList_X);

febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList_Y);

febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.VAL=mrow(bcSupportList_Z);

febio_spec.Mesh.NodeSet{4}.ATTR.name=nodeSetName4;
febio_spec.Mesh.NodeSet{4}.VAL=mrow(bcPrescribeList);

%MeshDomains section
febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

%MeshData section
% -> ElementData
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
febio_spec.MeshData.ElementData{1}.ATTR.type='mat_axis';

for q=1:1:size(E,1)
    febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
    febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir_perp(q,:);
    febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir_perp(q,:);
end

%Boundary condition section
% -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.name='zero_displacement_x';
febio_spec.Boundary.bc{1}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.x_dof=1;
febio_spec.Boundary.bc{1}.y_dof=0;
febio_spec.Boundary.bc{1}.z_dof=0;

febio_spec.Boundary.bc{2}.ATTR.name='zero_displacement_y';
febio_spec.Boundary.bc{2}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName2;
febio_spec.Boundary.bc{2}.x_dof=0;
febio_spec.Boundary.bc{2}.y_dof=1;
febio_spec.Boundary.bc{2}.z_dof=0;

febio_spec.Boundary.bc{3}.ATTR.name='zero_displacement_z';
febio_spec.Boundary.bc{3}.ATTR.type='zero displacement';
febio_spec.Boundary.bc{3}.ATTR.node_set=nodeSetName3;
febio_spec.Boundary.bc{3}.x_dof=0;
febio_spec.Boundary.bc{3}.y_dof=0;
febio_spec.Boundary.bc{3}.z_dof=1;

febio_spec.Boundary.bc{4}.ATTR.name='prescibed_displacement_z';
febio_spec.Boundary.bc{4}.ATTR.type='prescribed displacement';
febio_spec.Boundary.bc{4}.ATTR.node_set=nodeSetName4;
febio_spec.Boundary.bc{4}.dof='z';
febio_spec.Boundary.bc{4}.value.ATTR.lc=1;
febio_spec.Boundary.bc{4}.value.VAL=displacementMagnitude_perp;
febio_spec.Boundary.bc{4}.relative=0;

%LoadData section
% -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
%febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='sz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

febio_spec.Output.logfile.element_data{2}.ATTR.file=febioLogFileName_stretch;
febio_spec.Output.logfile.element_data{2}.ATTR.data='Uz';
febio_spec.Output.logfile.element_data{2}.ATTR.delim=',';

febioLogFileName_energy=[febioFebFileNamePart,'_energy_out.txt'];
febioLogFileName_damage=[febioFebFileNamePart,'_damage_out.txt'];

febio_spec.Output.logfile.element_data{3}.ATTR.file=febioLogFileName_energy;
febio_spec.Output.logfile.element_data{3}.ATTR.data='sed';
febio_spec.Output.logfile.element_data{3}.ATTR.delim=',';

% febio_spec.Output.logfile.element_data{4}.ATTR.file=febioLogFileName_damage;
% febio_spec.Output.logfile.element_data{4}.ATTR.data='damage';
% febio_spec.Output.logfile.element_data{4}.ATTR.delim=',';

% febio_spec.Output.logfile.element_data{4}.ATTR.file = ...
%     febioLogFileName_damage;
%
% febio_spec.Output.logfile.element_data{4}.ATTR.data = 'D';
%
% febio_spec.Output.logfile.element_data{4}.ATTR.delim = ',';


% Plotfile section
febio_spec.Output.plotfile.compression=0;


% febio_spec.Output.plotfile.var{1}.ATTR.type='damage';
febio_spec.Output.plotfile.var{1}.ATTR.type='strain energy density';

%% Creating febio analysis structure

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

%% Create structures for optimization

%What should be known to the objective function:
objectiveStruct.stretch_exp_perp_raw=stretch_exp_perp_raw;
objectiveStruct.stress_exp_perp_raw=stress_exp_perp_raw;
objectiveStruct.stretch_exp_para_raw=stretch_exp_para_raw;
objectiveStruct.stress_exp_para_raw=stress_exp_para_raw;

objectiveStruct.febioAnalysis=febioAnalysis;
objectiveStruct.febio_spec=febio_spec;
objectiveStruct.febioFebFileName=febioFebFileName;
objectiveStruct.febioLogFileName_disp = fullfile(savePath,febioLogFileName_disp);
objectiveStruct.febioLogFileName_stress = fullfile(savePath,febioLogFileName_stress);
objectiveStruct.febioLogFileName_stretch = fullfile(savePath,febioLogFileName_stretch);
objectiveStruct.febioLogFileName_energy = fullfile(savePath,febioLogFileName_energy);
objectiveStruct.febioLogFileName_damage = fullfile(savePath,febioLogFileName_damage);

objectiveStruct.fixedMaterial = ...
    [c1_ini m1_ini k1_ini k2_ini kappa_ini];
objectiveStruct.fixedDmax = ...
    [Dmax_matrix_ini, Dmax_GOH_ini];

objectiveStruct.parNormFactors=P; %This will normalize the parameters to ones(size(P))
objectiveStruct.Pb_struct.xx_c=P; %Parameter constraining centre
objectiveStruct.Pb_struct.xxlim = [
    0.45   1.30;   % matrix mu
    0.01   1.00;   % matrix sigma
    0.10   5.00;   % GOH mu
    0.01   2.00    % GOH sigma
    ];
objectiveStruct.cGOH = cGOH_ini;
objectiveStruct.k_factor=k_factor;

objectiveStruct.LangerAngles = LangerAngles;
objectiveStruct.displacementMagnitudes = displacementMagnitudes;

objectiveStruct.hf = hf1;
objectiveStruct.h=[Hn(5),Hn(6),Hn(7),Hn(8)];

objectiveStruct.hf_energy = hf2;
objectiveStruct.he = He;

objectiveStruct.hf_damage = hf3;
objectiveStruct.hd = Hd;

objectiveStruct.method = optimisationMethod;

Pn=P./objectiveStruct.parNormFactors;

switch evalMode
    case 1
        [errorVal,simData] = objectiveFunctionIFEA(Pn, objectiveStruct);
    case 2
        %% start optimization
        switch optimisationMethod
            case 1 %fminsearch and Nelder-Mead
                OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
                OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',1e-2,...
                    'DiffMaxChange',0.5);
                [Pn_opt,OPT_out.fval,OPT_out.exitflag,OPT_out.output]= fminsearch(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,OPT_options);
            case 2 %lsqnonlin and Levenberg-Marquardt
                OPT_options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
                OPT_options = optimoptions(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
                    'MaxIter',maxNumberIterations,...
                    'TolFun',functionTolerance,...
                    'TolX',parameterTolerance,...
                    'Display',displayTypeIterations,...
                    'FinDiffRelStep',1e-2,...
                    'DiffMaxChange',0.5);
                [Pn_opt,OPT_out.resnorm,OPT_out.residual]= lsqnonlin(@(Pn) objectiveFunctionIFEA(Pn,objectiveStruct),Pn,[],[],OPT_options);
        end

        %% Unnormalize and constrain parameters

        [errorVal,simData] = objectiveFunctionIFEA(Pn_opt, objectiveStruct);

        P_opt=Pn_opt.*objectiveStruct.parNormFactors; %Scale back, undo normalization

        %Constraining parameters
        for q=1:1:numel(P_opt)
            [P_opt(q)]=boxconstrain(P_opt(q),objectiveStruct.Pb_struct.xxlim(q,1),objectiveStruct.Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));
        end

        disp_text=sprintf('%6.16e,',P_opt); disp_text=disp_text(1:end-1);
        disp(['P_opt=',disp_text]);

        paramNames = { ...
            'mu_matrix','sigma_matrix', ...
            'mu_GOH','sigma_GOH'};
        for q = 1:numel(P_opt)
            lb = objectiveStruct.Pb_struct.xxlim(q,1);
            ub = objectiveStruct.Pb_struct.xxlim(q,2);

            if abs(P_opt(q)-lb) < 1e-3*abs(ub-lb)
                fprintf('WARNING: %s is close to LOWER bound\n',paramNames{q});
            elseif abs(P_opt(q)-ub) < 1e-3*abs(ub-lb)
                fprintf('WARNING: %s is close to UPPER bound\n',paramNames{q});
            end
        end

end

%% Import FEBio results

% Displacements
N_disp_mat1=simData(1).dataStruct_disp.data; %Displacement
timeVec1=simData(1).dataStruct_disp.time; %Time
V_DEF1=N_disp_mat1+repmat(V,[1 1 size(N_disp_mat1,3)]); %Deformed coordinate set
DN_magnitude1=sqrt(sum(N_disp_mat1(:,:,end).^2,2)); %Current displacement magnitude

N_disp_mat2=simData(2).dataStruct_disp.data; %Displacement
timeVec2=simData(2).dataStruct_disp.time; %Time
V_DEF2=N_disp_mat2+repmat(V,[1 1 size(N_disp_mat2,3)]); %Deformed coordinate set
DN_magnitude2=sqrt(sum(N_disp_mat2(:,:,end).^2,2)); %Current displacement magnitude


%%
% Plotting the simulated results using |anim8| to visualize and animate
% deformations

% Create basic view and store graphics handle to initiate animation
hf=cFigure; %Open figure
subplot(1,2,1);
title('Displ. magnitude [mm]','Interpreter','Latex')
hp1=gpatch(Fb,V_DEF1(:,:,end),DN_magnitude1,'k',1,2); %Add graphics object to animate
hp1.Marker='.'; hp1.MarkerSize=markerSize2; hp1.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize);
colormap(cMap); colorbar;
caxis([0 max([max(DN_magnitude1) max(DN_magnitude2)])]); caxis manual;
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically
view(140,30);
camlight headlight;

subplot(1,2,2);
title('Displ. magnitude [mm]','Interpreter','Latex')
hp2=gpatch(Fb,V_DEF2(:,:,end),DN_magnitude2,'k',1,2); %Add graphics object to animate
hp2.Marker='.'; hp2.MarkerSize=markerSize2; hp2.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize);
colormap(cMap); colorbar;
caxis([0 max([max(DN_magnitude1) max(DN_magnitude2)])]); caxis manual;
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically
view(140,30);
camlight headlight;

% Set up animation features
animStruct.Time=timeVec1; %The time vector
for qt=1:1:size(N_disp_mat1,3) %Loop over time increments

    DN_magnitude1=sqrt(sum(N_disp_mat1(:,:,qt).^2,2)); %Current displacement magnitude
    DN_magnitude2=sqrt(sum(N_disp_mat2(:,:,qt).^2,2)); %Current displacement magnitude

    %Set entries in animation structure
    animStruct.Handles{qt}=[hp1 hp1 hp2 hp2]; %Handles of objects to animate
    animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
    animStruct.Set{qt}={V_DEF1(:,:,qt),DN_magnitude1,V_DEF2(:,:,qt),DN_magnitude2}; %Property values for to set in order to animate
end
anim8(hf,animStruct); %Initiate animation feature
drawnow;

%%

%Access data
E_stress_mat1=simData(1).dataStruct_stress.data;
E_stretch_mat1=simData(1).dataStruct_stretch.data;
[CV1]=faceToVertexMeasure(E,V,E_stress_mat1(:,:,end));

E_stress_mat2=simData(2).dataStruct_stress.data;
E_stretch_mat2=simData(2).dataStruct_stretch.data;
[CV2]=faceToVertexMeasure(E,V,E_stress_mat2(:,:,end));

%%
% Plotting the simulated results using |anim8| to visualize and animate
% deformations

% Create basic view and store graphics handle to initiate animation
hf=cFigure; %Open figure  /usr/local/MATLAB/R2020a/bin/glnxa64/jcef_helper: symbol lookup error: /lib/x86_64-linux-gnu/libpango-1.0.so.0: undefined symbol: g_ptr_array_copy
subplot(1,2,1);
title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
hp1=gpatch(Fb,V_DEF1(:,:,end),CV1,'k',1,2); %Add graphics object to animate
hp1.Marker='.'; hp1.MarkerSize=markerSize2; hp1.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize);
colormap(gca,cMap); colorbar;
caxis([min(E_stress_mat1(:)) max(E_stress_mat1(:))]);
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically
view(140,30);
camlight headlight;

subplot(1,2,2);
title('$\sigma_{zz}$ [MPa]','Interpreter','Latex')
hp2=gpatch(Fb,V_DEF2(:,:,end),CV2,'k',1,2); %Add graphics object to animate
hp2.Marker='.'; hp2.MarkerSize=markerSize2; hp2.FaceColor='interp';
gpatch(Fb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

axisGeom(gca,fontSize);
colormap(gca,cMap); colorbar;
caxis([min(E_stress_mat2(:)) max(E_stress_mat2(:))]);
axis(axisLim([V_DEF1; V_DEF2])); %Set axis limits statically
view(140,30);
camlight headlight;

% Set up animation features
animStruct.Time=timeVec1; %The time vector
for qt=1:1:size(N_disp_mat1,3) %Loop over time increments

    [CV1]=faceToVertexMeasure(E,V,E_stress_mat1(:,:,qt));
    [CV2]=faceToVertexMeasure(E,V,E_stress_mat2(:,:,qt));

    %Set entries in animation structure
    animStruct.Handles{qt}=[hp1 hp1 hp2 hp2]; %Handles of objects to animate
    animStruct.Props{qt}={'Vertices','CData','Vertices','CData'}; %Properties of objects to animate
    animStruct.Set{qt}={V_DEF1(:,:,qt),CV1,V_DEF2(:,:,qt),CV2}; %Property values for to set in order to animate
end
anim8(hf,animStruct); %Initiate animation feature
drawnow;

%%

function [errorVal,simData] = objectiveFunctionIFEA(Pn, objectiveStruct)

%% Access input structure
stretch_exp_perp_raw = objectiveStruct.stretch_exp_perp_raw;
stress_exp_perp_raw  = objectiveStruct.stress_exp_perp_raw;
stretch_exp_para_raw = objectiveStruct.stretch_exp_para_raw;
stress_exp_para_raw = objectiveStruct.stress_exp_para_raw;

febioAnalysis = objectiveStruct.febioAnalysis;
febio_spec = objectiveStruct.febio_spec;
febioFebFileName = objectiveStruct.febioFebFileName;
febioLogFileName_disp = objectiveStruct.febioLogFileName_disp;
febioLogFileName_stress = objectiveStruct.febioLogFileName_stress;
febioLogFileName_stretch = objectiveStruct.febioLogFileName_stretch;
febioLogFileName_energy = objectiveStruct.febioLogFileName_energy;
febioLogFileName_damage = objectiveStruct.febioLogFileName_damage;


parNormFactors = objectiveStruct.parNormFactors;
Pb_struct = objectiveStruct.Pb_struct;
k_factor = objectiveStruct.k_factor;

LangerAngles = objectiveStruct.LangerAngles;
displacementMagnitudes = objectiveStruct.displacementMagnitudes;

hf = objectiveStruct.hf;
h = objectiveStruct.h;
hf_energy = objectiveStruct.hf_energy;
he = objectiveStruct.he;
hf_damage = objectiveStruct.hf_damage;
hd = objectiveStruct.hd;

%% Set/get material parameters
% Unnormalize and constrain parameters
P = Pn.*parNormFactors; %Scale back, undo normalization

%Constraining parameters
for q=1:1:numel(P)
    [P(q)]=boxconstrain(P(q),Pb_struct.xxlim(q,1),Pb_struct.xxlim(q,2),objectiveStruct.Pb_struct.xx_c(q));
end

disp(['Trying   : ',sprintf(repmat('%6.8e ',[1,numel(P)]),P)]);

% Fixed intact material parameters
fixedMaterial = objectiveStruct.fixedMaterial;

c1    = fixedMaterial(1);
m1    = fixedMaterial(2);
k1    = fixedMaterial(3);
k2    = fixedMaterial(4);
kappa = fixedMaterial(5);

% Ground-matrix damage parameters
mu_matrix    = P(1);
sigma_matrix = P(2);

mu_GOH       = P(3);
sigma_GOH    = P(4);

Dmax_matrix  = objectiveStruct.fixedDmax(1);
Dmax_GOH     = objectiveStruct.fixedDmax(2);

% Derived fixed parameters
kOgden = c1*k_factor;
cGOH   = objectiveStruct.cGOH;
kGOH   = cGOH*k_factor;


%% Update damaged Ogden matrix
febio_spec.Material.material{1}.solid{1}.k = kOgden;

febio_spec.Material.material{1}.solid{1}.elastic.c1 = c1;
febio_spec.Material.material{1}.solid{1}.elastic.m1 = m1;
febio_spec.Material.material{1}.solid{1}.elastic.k  = kOgden;

febio_spec.Material.material{1}.solid{1}.damage.mu = mu_matrix;

febio_spec.Material.material{1}.solid{1}.damage.sigma = sigma_matrix;

febio_spec.Material.material{1}.solid{1}.damage.Dmax = Dmax_matrix;


%% Update damaged GOH +gamma

febio_spec.Material.material{1}.solid{2}.elastic.c = cGOH/2;
febio_spec.Material.material{1}.solid{2}.elastic.k1 = k1;
febio_spec.Material.material{1}.solid{2}.elastic.k2 = k2;
febio_spec.Material.material{1}.solid{2}.elastic.kappa = kappa;
febio_spec.Material.material{1}.solid{2}.elastic.k = kGOH/2;

febio_spec.Material.material{1}.solid{2}.damage.mu = mu_GOH;
febio_spec.Material.material{1}.solid{2}.damage.sigma = sigma_GOH;
febio_spec.Material.material{1}.solid{2}.damage.Dmax = Dmax_GOH;


%% Update damaged GOH -gamma

febio_spec.Material.material{1}.solid{3}.elastic.c = cGOH/2;
febio_spec.Material.material{1}.solid{3}.elastic.k1 = k1;
febio_spec.Material.material{1}.solid{3}.elastic.k2 = k2;
febio_spec.Material.material{1}.solid{3}.elastic.kappa = kappa;
febio_spec.Material.material{1}.solid{3}.elastic.k = kGOH/2;

febio_spec.Material.material{1}.solid{3}.damage.mu = mu_GOH;
febio_spec.Material.material{1}.solid{3}.damage.sigma = sigma_GOH;
febio_spec.Material.material{1}.solid{3}.damage.Dmax = Dmax_GOH;

nElem = size(febio_spec.Mesh.Elements{1}.elem.VAL,1);

% Set orientation data and deformation magnitude
for i = 1:1:numel(LangerAngles)

    % Create material axis vectors
    R = euler2DCM([0,-LangerAngles(i),0]); % Rotation matrix
    e1 = (R*[1 0 0]')'; % Rotated e1 vector
    e1_dir = e1(ones(nElem,1),:); % e1 copied for all elements
    e2 = [0 1 0];% (R*[0 1 0]')'; % Rotated e1 vector
    e2_dir = e2(ones(nElem,1),:); % e1 copied for all elements

    %Update the MeshData section with new material axis vectors
    for q=1:1:nElem
        febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
        febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir(q,:);
        febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir(q,:);
    end

    % Set the current displacement magnitude
    febio_spec.Boundary.bc{4}.value.VAL=displacementMagnitudes(i);

    % Export to .feb file
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

    % Run the current model
    [runFlags(i)]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

    if runFlags(i) == 1
        % Importing nodal displacements from a log file
        simData(i).dataStruct_disp=importFEBio_logfile(febioLogFileName_disp,0,1);

        % Importing element stress from a log file
        simData(i).dataStruct_stress=importFEBio_logfile(febioLogFileName_stress,0,1);

        % Importing element stretch from a log file
        simData(i).dataStruct_stretch=importFEBio_logfile(febioLogFileName_stretch,0,1);

        % Importing strain energy density
        simData(i).dataStruct_energy = importFEBio_logfile(febioLogFileName_energy,0,1);

        % simData(i).dataStruct_damage = ...
        %     importFEBio_logfile(febioLogFileName_damage,0,1);

    else
        simData(i).dataStruct_disp=NaN;
        simData(i).dataStruct_stress=NaN;
        simData(i).dataStruct_stretch=NaN;
        simData(i).dataStruct_energy=NaN;
        simData(i).dataStruct_damage = NaN;
    end
end


if all(runFlags==1)

    % simData(1) = perpendicular
    stretch_perp_sim = squeeze(mean(mean(simData(1).dataStruct_stretch.data,1),2));
    stress_perp_sim  = squeeze(mean(mean(simData(1).dataStruct_stress.data,1),2));
    energy_perp_sim  = squeeze(mean(mean(simData(1).dataStruct_energy.data,1),2));

    % simData(2) = parallel
    stretch_para_sim = squeeze(mean(mean(simData(2).dataStruct_stretch.data,1),2));
    stress_para_sim  = squeeze(mean(mean(simData(2).dataStruct_stress.data,1),2));
    energy_para_sim  = squeeze(mean(mean(simData(2).dataStruct_energy.data,1),2));

    % Damage from same CDF log-normal law
    % normalCDF = @(x) 0.5.*(1 + erf(x./sqrt(2))); % avoids needing Statistics Toolbox
    %
    % Ecrit_perp = 0.5*(stretch_perp_sim.^2 - 1);
    % Ecrit_para = 0.5*(stretch_para_sim.^2 - 1);
    %
    % damage_perp_sim = Dmax .* normalCDF((Ecrit_perp - mu)./sigma);
    % damage_para_sim = Dmax .* normalCDF((Ecrit_para - mu)./sigma);

    % % Approximate maximum Lagrange strain for this uniaxial test
    % Ecrit_perp = 0.5*(stretch_perp_sim.^2 - 1);
    % Ecrit_para = 0.5*(stretch_para_sim.^2 - 1);
    %
    % % FEBio damage is irreversible, so use the maximum previous value
    % Xhist_perp = cummax(max(Ecrit_perp,0));
    % Xhist_para = cummax(max(Ecrit_para,0));
    %
    % damage_perp_sim = zeros(size(Xhist_perp));
    % damage_para_sim = zeros(size(Xhist_para));
    %
    % indPerp = Xhist_perp > 0;
    % indPara = Xhist_para > 0;
    %
    % damage_perp_sim(indPerp) = Dmax .* 0.5 .* erfc( ...
    %     -log(Xhist_perp(indPerp)./mu) ./ (sigma*sqrt(2)));
    %
    % damage_para_sim(indPara) = Dmax .* 0.5 .* erfc( ...
    %     -log(Xhist_para(indPara)./mu) ./ (sigma*sqrt(2)));

    % Read actual damage calculated by FEBio
    % damage_perp_sim = squeeze(mean(mean( ...
    %     simData(1).dataStruct_damage.data,1),2));
    %
    % damage_para_sim = squeeze(mean(mean( ...
    %     simData(2).dataStruct_damage.data,1),2));

    % Interpolating simulation data for the experimental stretch points
    stress_perp_sim_exp = interp1(stretch_perp_sim,stress_perp_sim,stretch_exp_perp_raw,'pchip','extrap');
    stress_para_sim_exp = interp1(stretch_para_sim,stress_para_sim,stretch_exp_para_raw,'pchip','extrap');

    fprintf('Max experimental perp stretch = %.4f\n', max(stretch_exp_perp_raw));
    fprintf('Max FEA perp stretch          = %.4f\n', max(stretch_perp_sim));
    fprintf('Number of FEA perp points     = %d\n', numel(stretch_perp_sim));
    fprintf('Final FEA perp stress         = %.4f MPa\n', stress_perp_sim(end));

    fprintf('Max experimental para stretch = %.4f\n', max(stretch_exp_para_raw));
    fprintf('Max FEA para stretch          = %.4f\n', max(stretch_para_sim));
    fprintf('Number of FEA para points     = %d\n', numel(stretch_para_sim));
    fprintf('Final FEA para stress         = %.4f MPa\n', stress_para_sim(end));


    % Update stress-stretch graph
    figure(hf);
    set(h(1),'XData',stretch_perp_sim,'YData',stress_perp_sim);
    set(h(2),'XData',stretch_para_sim,'YData',stress_para_sim);

    set(h(3),'XData',stretch_exp_perp_raw,'YData',stress_perp_sim_exp);
    set(h(4),'XData',stretch_exp_para_raw,'YData',stress_para_sim_exp);

    drawnow; axis tight; drawnow;

    % Update strain-energy graph
    figure(hf_energy);
    set(he(1),'XData',stretch_perp_sim,'YData',energy_perp_sim);
    set(he(2),'XData',stretch_para_sim,'YData',energy_para_sim);
    drawnow; axis tight; drawnow;

    % Update damage graph
    % figure(objectiveStruct.hf_damage)
    % set(objectiveStruct.hd(1),'XData',stretch_perp_sim,'YData',damage_perp_sim);
    % set(objectiveStruct.hd(2),'XData',stretch_para_sim,'YData',damage_para_sim);
    % drawnow; axis tight; drawnow;

    % Compute error metrics
    diff1 = stress_exp_perp_raw(2:end) - stress_perp_sim_exp(2:end);
    diff2 = stress_exp_para_raw(2:end) - stress_para_sim_exp(2:end);

    scale1 = max(abs(stress_exp_perp_raw));
    scale2 = max(abs(stress_exp_para_raw));

    % errorVal1 = diff1 ./ scale1;
    % errorVal2 = diff2 ./ scale2;
    %
    % % Because n_perp = 24 and n_para = 22, make total directional weight equal
    % errorVal1 = errorVal1 ./ sqrt(numel(errorVal1));
    % errorVal2 = errorVal2 ./ sqrt(numel(errorVal2));
    %
    % errorVal = [errorVal1(:); errorVal2(:)];

    % Normalised residuals
    errorVal1 = diff1 ./ scale1;
    errorVal2 = diff2 ./ scale2;

    % Equal total contribution of perp and para
    errorVal1 = errorVal1 ./ sqrt(numel(errorVal1));
    errorVal2 = errorVal2 ./ sqrt(numel(errorVal2));

    RMSE_perp = sqrt(mean((diff1./scale1).^2));
    RMSE_para = sqrt(mean((diff2./scale2).^2));

    disp(['RMSE perp: ', sprintf('%6.8e',RMSE_perp)]);
    disp(['RMSE para: ', sprintf('%6.8e',RMSE_para)]);

    % Final-point error: force final FEA perp stress down
    finalPerpError = ...
        (stress_perp_sim_exp(end) - stress_exp_perp_raw(end)) ./ scale1;

    % Drop-depth error: force FEA peak-to-final drop to match experiment
    dropExp_perp = max(stress_exp_perp_raw) - stress_exp_perp_raw(end);
    dropFEA_perp = max(stress_perp_sim_exp) - stress_perp_sim_exp(end);

    dropDepthError_perp = ...
        (dropFEA_perp - dropExp_perp) ./ scale1;

    % End-slope penalty: penalise FEA if it is still rising at the end
    endSlopePenalty_perp = ...
        max(0, stress_perp_sim_exp(end) - stress_perp_sim_exp(end-1)) ./ scale1;


    % --------------------------------------------------------
    % Extra constraints to force parallel post-peak drop
    % --------------------------------------------------------

    finalParaError = ...
        (stress_para_sim_exp(end) - stress_exp_para_raw(end)) ./ scale2;

    dropExp_para = max(stress_exp_para_raw) - stress_exp_para_raw(end);
    dropFEA_para = max(stress_para_sim_exp) - stress_para_sim_exp(end);

    dropDepthError_para = ...
        (dropFEA_para - dropExp_para) ./ scale2;

    endSlopePenalty_para = ...
        max(0, stress_para_sim_exp(end) - stress_para_sim_exp(end-1)) ./ scale2;


    % Weights
    wFinal_perp = 2;
    wDrop_perp  = 2;
    wSlope_perp = 0.5;

    wFinal_para = 2;
    wDrop_para  = 2;
    wSlope_para = 0.5;

    % Directional weights for final run
    dirW_perp = 2.0;   % because RMSE_perp is currently much worse
    dirW_para = 1.0;

    errorVal = [
        dirW_perp * errorVal1(:);
        dirW_para * errorVal2(:);

        wFinal_perp * finalPerpError;
        wDrop_perp  * dropDepthError_perp;
        wSlope_perp * endSlopePenalty_perp;

        wFinal_para * finalParaError;
        wDrop_para  * dropDepthError_para;
        wSlope_para * endSlopePenalty_para
        ];
else
    errorVal = [NaN(size(stress_exp_perp_raw)); NaN(size(stress_exp_para_raw))];
end

SSE = sum(errorVal.^2);
disp(['SSE: ', sprintf('%6.8e',SSE)]);

n = numel(errorVal);
RMSE = sqrt(SSE/n);
disp(['RMSE: ', sprintf('%6.8e',RMSE)]);

if objectiveStruct.method == 1
    errorVal = SSE; % use sum of squares for Nelder-Mead
end


end

