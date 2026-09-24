
function geo = allGeometries(geometryID, pointSpacing, distRefine)
    switch geometryID
        case 1
            [F,V,indRefine,V_refine, V1, p, d, t] = geometry1(pointSpacing, distRefine);
        case 2
            [F,V,indRefine,V_refine, V1, p, d, t] = geometry2(pointSpacing, distRefine);
        case 3
            [F,V,indRefine,V_refine, V1, p, d, t] = geometry3(pointSpacing, distRefine);
        case 4
            [F,V,indRefine,V_refine, V1, p, d, t] = geometry4(pointSpacing, distRefine);
        case 5
            [F,V,indRefine,V_refine, V1, p, d, t] = geometry5(pointSpacing, distRefine);   
    end
    geo.F = F;
    geo.V = V;
    geo.indRefine = indRefine;
    geo.V_refine = V_refine;
    geo.V1 = V1;
    geo.p = p;
    geo.d = d;
    geo.t = t;


end

%% first geometry
function [F,V,indRefine,V_refine, V1, p, d, t] = geometry1(pointSpacing, distRefine)
% p=[1,1];
p=[3,3];
s=1; %cut width
c1=3.8; %cutspacing
c2=3.8;
d1=22.04;
d2=22.04;
t=3; %thickness

L1=d1-c1;
L2=d2-c2;
sc1=c1/2-s/2;
sc2=c2/2-s/2;
d=[d1,d2];

p1=[0,-s/2];
p2=[L1/2,-s/2];
p3=[L1/2,0];
p4=p3+[sc1,0];
p5=p4-[0,L2/2];
p6=p5+[s/2,0];
p7=p6-[0,sc2];
p8=p7-[L1/2,0];
p9=p8-[0,s/2];
p10=p9-[sc1,0];
p11=p10+[0,L2/2];
p12=p11-[s/2,0];


V1=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12];

V_refine = [p2;p3; p5;p6; p8;p9; p11;p12];
V_refine(:,3)=0;

interpMethod='linear';
closeLoopOpt=1; %Option for closed curve
indMust=1:length(V1);%Ensure corners are included
[V1]=evenlySpaceCurve(V1,pointSpacing,interpMethod,closeLoopOpt,indMust);%Resampling sampling for desired line spacing

%Defining a region
regionCell={V1}; %A region between V1 and V2 (V2 forms a hole inside V1)
resampleCurveOpt=0; %Option to turn on/off resampling of input boundary curves

[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);

V(:,3)=0;
[~,indRefine] = minDist(V_refine,V);
end


%% second geometry
function [F,V,indRefine,V_refine, V1, p, d, t] = geometry2(pointSpacing, distRefine)
%% input parameters
% w=60;%sample width
% pointSpacing =2;
p=[3,3];
d1 = 60/p(1);
d2 = 60/p(2);
d  = [d1 d2];
s=1; %cut width
c1=13; %cutspacing
c2=13;
% d1=29;
% d2=29;
t=3; %thickness

L1=d1-c1;
L2=d2-c2;
sc1=c1/2-s/2;
sc2=c2/2-s/2;
d=[d1,d2];

p1=[0,-s/2];
p2=[L1/2,-s/2];
p3=[L1/2,0];
p4=p3+[sc1,0];
p5=p4-[0,L2/2];
p6=p5+[s/2,0];
p7=p6-[0,sc2];
p8=p7-[L1/2,0];
p9=p8-[0,s/2];
p10=p9-[sc1,0];
p11=p10+[0,L2/2];
p12=p11-[s/2,0];


V1=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12];

V_refine = [p2;p3; p5;p6; p8;p9; p11;p12];
V_refine(:,3)=0;

interpMethod='linear';
closeLoopOpt=1; %Option for closed curve
indMust=1:length(V1);%Ensure corners are included
[V1]=evenlySpaceCurve(V1,pointSpacing,interpMethod,closeLoopOpt,indMust);%Resampling sampling for desired line spacing


%Defining a region
regionCell={V1}; %A region between V1 and V2 (V2 forms a hole inside V1)
resampleCurveOpt=0; %Option to turn on/off resampling of input boundary curves

[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);
V(:,3)=0;

[~,indRefine] = minDist(V_refine,V);
end

%% third geometry
function [F,V,indRefine,V_refine, V1, p, d, t] = geometry3(pointSpacing, distRefine)

% p=[1,1];
p=[6,3]; 
L=9.5;
N=7.5;
T=1;
s2=7.5; %I width
s3=1; %blade width
d1=10.5;
d2=21;
d=[d1,d2];
t=3; %thickness

p1=[0,0];
p2=[0,-N/2];
p3=[-s2/2,-N/2];
p4=[-s2/2,-N/2-T];
p5=[-T/2,-N/2-T];
p6=[-T/2,-N/2-T-s3-T-N/2];
p7=[-T/2-s3-s2/2,-N/2-T-s3-T-N/2];
p8=[-T/2-s3-s2/2,-T-s3-T-N/2];
p9=[-s3-T/2,-T-s3-T-N/2];
p10=[-s3-T/2,-s3-T-N/2];
p11=[-s2/2-s3,-s3-T-N/2];
p12=[-s2/2-s3,0];

V1=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12];

V_refine = [p3;p4; p5;p8; p9;p10; p11];
V_refine(:,3)=0;


%
interpMethod='linear';
closeLoopOpt=1; %Option for closed curve
indMust=1:length(V1);%Ensure corners are included
[V1]=evenlySpaceCurve(V1,pointSpacing,interpMethod,closeLoopOpt,indMust);%Resampling sampling for desired line spacing


% Defining a region
regionCell={V1}; %A region between V1 and V2 (V2 forms a hole inside V1)
resampleCurveOpt=0; %Option to turn on/off resampling of input boundary curves

[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);
V(:,3)=0;

[~,indRefine] = minDist(V_refine,V);
end

%% fourth geometry 
function [F,V,indRefine,V_refine, V1, p, d, t] = geometry4(pointSpacing, distRefine)
% p=[1,1];
p=[3,3];
h1=15.36;
h2=7.68;
n=8.60;
s1=0.66;
s2=0.46;
s3=1.42;
s4=1.82;
s6=2.38;
T=1;
d1=19.2;
d2=20;
d=[d1,d2];
t=1; %thickness

p1=[0,0];
p2=[0,-h1/2];
p3=[T/2,-h1/2];
p4=[T/2,-h1/2-s4];
p5=[-n/2+T/2,-h1/2-s4];
p6=[-n/2+T/2,-h1/2-s4+h2/2-T/2];
p7=[-n/2-T/2,-h1/2-s4+h2/2-T/2];
p8=[-n/2-T/2,-h1/2-s4];
p9=[-h1/2+T/2,-h1/2-s4];
p10=[-h1/2+T/2,-h1/2-s4-T/2];
p11=[-h1/2+T/2-s3,-h1/2-s4-T/2];
p12=[-h1/2+T/2-s3,-T/2-s4];
p13=[-h1/2-s3,-T/2-s4];
p14=[-h1/2-s3,-T/2];
p15=[-n/2-T/2,-T/2];
p16=[-n/2-T/2,-h2/2];
p17=[-n/2+T/2,-h2/2];
p18=[-n/2+T/2,-T/2];
p19=[-s3,-T/2];
p20=[-s3,0];
V1=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12;p13;p14;p15;p16;p17;p18;p19;p20];

p21=[-s2,-h2/2-s1];
p22=[-s2,-h2/2-s1-T];
p23=[-s2-h2,-h2/2-s1-T];
p24=[-s2-h2,-h2/2-s1];
V2=[p21;p22;p23;p24];


interpMethod='linear';
closeLoopOpt=1; %Option for closed curve
indMust=1:length(V1);%Ensure corners are included
[V1]=evenlySpaceCurve(V1,pointSpacing,interpMethod,closeLoopOpt,indMust);

indMust=1:length(V2);%Ensure corners are included
[V2]=evenlySpaceCurve(V2,pointSpacing,interpMethod,closeLoopOpt,indMust);


V_refine = [p2;p5;p6;p7;p8;p9;p12;p14;p15;p16;p17;p18;p19;p21;p22;p23;p24];
V_refine(:,3)=0;

% Defining a region
regionCell={V1,V2}; %A region between V1 and V2 (V2 forms a hole inside V1)
resampleCurveOpt=0; %Option to turn on/off resampling of input boundary curves
[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);
V(:,3)=0;
[~,indRefine] = minDist(V_refine,V);

V(:,1)=V(:,1)-min(V(:,1));
V(:,2)=V(:,2)-min(V(:,2));

end

%% fifth geometry
function [F,V,indRefine,V_refine, V1, p, d, t] = geometry5(pointSpacing, distRefine)
%% input parameters
% p=[1,1];
p=[5,3];
l1=1.69*4.5;
l2=3.86*4.5;
s1=0.24*4.5;
s2=0.5*4.5;
s3=0.12*4.5;
s4=0.13*4.5;
s5=0.22*4.5;
T=1;
d1=10.7750;
d2=20.53;
d=[d1,d2];
t=3; %thickness

p1=[0,0];
p2=[0,-l2/2-s1];
p3=[-s4,-l2/2-s1];
p4=[-s4,-s1];
p5=[-s4-T,-s1];
p6=[-s4-T,-s1-l2/2+T/2];
p7=[-s4-T-l1/2,-s1-l2/2+T/2];
p8=[-s4-T-l1/2,T/2];
p9=[-T-l1/2,T/2];
p10=[-T-l1/2,-l2/2+T/2];
p11=[-l1/2,-l2/2+T/2];
p12=[-l1/2,0];

V1=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11;p12];

interpMethod='linear';
closeLoopOpt=1; %Option for closed curve
indMust=1:length(V1);%Ensure corners are included
[V1]=evenlySpaceCurve(V1,pointSpacing,interpMethod,closeLoopOpt,indMust);

V_refine = [p4;p5;p6;p10;p11;p12];
V_refine(:,3)=0;

% % Defining a region
regionCell={V1}; %A region between V1 and V2 (V2 forms a hole inside V1)
resampleCurveOpt=0; %Option to turn on/off resampling of input boundary curves
[F,V]=regionTriMesh2D(regionCell,pointSpacing,resampleCurveOpt,0);

V(:,1)=V(:,1)-min(V(:,1));
V(:,2)=V(:,2)-min(V(:,2));

V(:,3)=0;

[~,indRefine] = minDist(V_refine,V);
end
