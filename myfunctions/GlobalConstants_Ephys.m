

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 0:  Define constants that will not be varied randomly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global GKS_mul VP_d VP_f GCal_mul
GKS_mul = 8;
VP_d = 12;
VP_f = 8;
GCal_mul = 1.9;

global cipa CiPA_scale_IKr CiPA_scale_IKs CiPA_scale_IK1 CiPA_scale_ICaL CiPA_scale_INaL
if cipa == 1
    CiPA_scale_IKr = 1.119;
    CiPA_scale_IKs = 1.648;
    CiPA_scale_IK1 = 1.414;
    CiPA_scale_ICaL = 1.018;
    CiPA_scale_INaL = 2.274;
else
    CiPA_scale_IKr = 1;
    CiPA_scale_IKs = 1;
    CiPA_scale_IK1 = 1;
    CiPA_scale_ICaL = 1;
    CiPA_scale_INaL = 1;
end

%extracellular ionic concentrations
global Nao Cao Ko;
Nao=140.0;
Cao=1.8;
Ko=5.4;

%physical constants
global R T F Cm;
R=8314.0;
T=310.0;
F=96485.0;
Cm=1.0;                     %uF

%cell geometry
global L rad vcell Ageo Acap vmyo vnsr vjsr vss;
L=0.01;
rad=0.0011;
vcell=1000*3.14*rad*rad*L;
Ageo=2*3.14*rad*rad+2*3.14*rad*L;
Acap=2*Ageo;
vmyo=0.68*vcell;
vnsr=0.0552*vcell;
vjsr=0.0048*vcell;
vss=0.02*vcell;
% 
%CaMK constants
global aCaMK bCaMK CaMKo;
aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;

%reverse potential constant
global PKNa;
PKNa=0.01833;

%INa constants
global Ahf Ahs;  % fraction of channels go through fast/slow inactivation
Ahf=0.99;
Ahs=1.0-Ahf;

%INaL constants  
global thL;      
thL=200.0;

%Ito constants
global delta_epi;
% % % if celltype==1
% % %     delta_epi=1.0-(0.95/(1.0+exp((V+70.0)/5.0)));
% % % else
% % %     delta_epi=1.0;
% % % end

%ICaL, ICaNa, ICaK constants
global Aff tjca Kmn k2n zca;
Aff=0.6;
tjca=75.0;
Kmn=0.002;
k2n=1000.0;
zca=2.0;

%INaCa_i constants
global kna1 kna2 kna3 kasymm wna wca wnaca kcaon kcaoff qna qca zna;
kna1=15.0;
kna2=5.0;
kna3=88.12;
kasymm=12.5;
wna=6.0e4;
wca=6.0e4;
wnaca=5.0e3;
kcaon=1.5e6;
kcaoff=5.0e3;
qna=0.5224;
qca=0.1670;
zna=1.0;

%INaCa_ss constants
global KmCaAct;
KmCaAct=150.0e-6;

%INaK constants
global k1p k1m k2p k2m k3p k3m k4p k4m Knao0 delta; % Knai0 strip off global in order to modify as PKA-P
k1p=949.5;
k1m=182.4;
k2p=687.2;
k2m=39.4;
k3p=1899.0;
k3m=79300.0;
k4p=639.0;
k4m=40.0;
Knai0=9.073;
Knao0=27.78;
delta=-0.1550;
global Kki Kko MgADP MgATP Kmgatp H eP Khp Knap Kxkur; 
Kki=0.5;
Kko=0.3582;
MgADP=0.05;
MgATP=9.8;
Kmgatp=1.698e-7;
H=1.0e-7;
eP=4.2;
Khp=1.698e-7;
Knap=224.0;
Kxkur=292.0;
global zk;
zk=1.0;

%calcium buffer constants
global cmdnmax kmcmdn trpnmax BSRmax KmBSR BSLmax KmBSL csqnmax kmcsqn; % kmtrpn strip off global in order to modify as PKA-P
cmdnmax=0.05;
if celltype==1
    cmdnmax=cmdnmax*1.3;
end
kmcmdn=0.00238;
trpnmax=0.07;
kmtrpn=0.0005;
BSRmax=0.047;
KmBSR=0.00087;
BSLmax=1.124;
KmBSL=0.0087;
csqnmax=10.0;
kmcsqn=0.8;

%jsr constants
global bt a_rel;
bt=4.75;
a_rel=0.5*bt;

global GNa GNaL Gto GKr GKs GK1 Gncx GKb GpCa PCa Pnak PNab PCab KmCaMK KmCaM;

GNa=75;

GNaL=0.0075 * CiPA_scale_INaL;
if celltype==1
    GNaL=GNaL*0.6;
end

Gto=0.02;
if celltype==1
    Gto=Gto*4.0;
elseif celltype==2
    Gto=Gto*4.0;
end

GKr=0.046 * CiPA_scale_IKr;
if celltype==1
    GKr=GKr*1.3;
elseif celltype==2
    GKr=GKr*0.8;
end

GKs=0.0034 * CiPA_scale_IKs;
GK1=0.1908 * CiPA_scale_IK1;
if celltype==1
    GKs=GKs*1.4;
end
if celltype==1
    GK1=GK1*1.2;
elseif celltype==2
    GK1=GK1*1.3;
end

Gncx=0.0008;                %GNaCa
if celltype==1
    Gncx=Gncx*1.1;
elseif celltype==2
    Gncx=Gncx*1.4;
end

GKb=0.003;
if celltype==1
    GKb=GKb*0.6;
end

GpCa=0.0005; % sarcolemmal calcium pump current IpCa
 
PCa=0.0001 * CiPA_scale_ICaL;
if celltype==1
    PCa=PCa*1.2;
elseif celltype==2
    PCa=PCa*2.5;
end

Pnak=30;
if celltype==1
    Pnak=Pnak*0.9;
elseif celltype==2
    Pnak=Pnak*0.7;
end

PNab=3.75e-10;
PCab=2.5e-8;

KmCaMK=0.15;
KmCaM=0.0015;


% Add-on factors collectively representing activity of SERCA, RyR, Leak, and nsr to jsr Ca translocation
global SERCA_total RyR_total Leak_total Trans_total
SERCA_total = 1 ;                                                                         
RyR_total = 1 ;
Leak_total = 1 ;
Trans_total = 1 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gs - conductance or Ks - maximal current or Ps - permeability
% log-normally distributed sigma = 0.2624
% n_Gs = 17
baseline_Gs = [GNa GNaL Gto GKr GKs GK1 Gncx PCa Pnak PNab PCab GKb GpCa SERCA_total RyR_total Leak_total Trans_total] ;


