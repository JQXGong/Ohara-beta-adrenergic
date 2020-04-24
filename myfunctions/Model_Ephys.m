function deriv = Model_Ephys(t,statevar,Istim,settings)

global GKS_mul VP_d VP_f GCal_mul
global Nao Cao Ko;
global R T F Cm;
global Acap vmyo vnsr vjsr vss;
global aCaMK bCaMK CaMKo;
global PKNa;
global Ahf Ahs;
global thL;
global delta_epi;
global Aff tjca Kmn k2n zca;
global kna1 kna2 kna3 kasymm wna wca wnaca kcaon kcaoff qna qca zna;
global KmCaAct;
global k1p k1m k2p k2m k3p k3m k4p k4m Knao0 delta; % Knai0 strip off global in order to modify as PKA-P
global Kki Kko MgADP MgATP Kmgatp H eP Khp Knap Kxkur;
global zk;
global cmdnmax kmcmdn trpnmax BSRmax KmBSR BSLmax KmBSL csqnmax kmcsqn; % kmtrpn strip off global in order to modify as PKA-P
global bt a_rel Jrel_inf tau_rel;
global GNa GNaL Gto GKr GKs GK1 Gncx GKb GpCa PCa Pnak PNab PCab KmCaMK KmCaM;
global celltype

global SERCA_total RyR_total Leak_total Trans_total


statevarcell = num2cell(statevar) ;

[V Nai Nass Ki Kss Cai Cass Cansr Cajsr m hf hs j hsp jp mL hL hLp a iF ...
    iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 ...
    xk1 Jrelnp Jrelp JrelBP CaMKt...% all the PKA variables following
    hPf hPs jP hBPf hBPs jBP dP fPf fPs fcaPf fcaPs fBPf fcaBPf...
    xs1P xs2P JrelP ] = deal(statevarcell{:}) ;

Id = Istim; % The stimulus is defined in the wrapper code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calculate the stimulus current (simulus current defined here)
% if t < settings.stimdelay
%     Id = 0; 
% elseif t < settings.stimdelay + settings.stimdur
%     Id = Istim;
% else
%     Id = 0; 
% end

%% 
% celltype=0; %endo = 0, epi = 1, M = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INa_P = settings.fINaP;
ICaL_P = settings.fICaLP;
IKs_P = settings.fIKsP;
INaK_P = settings.fINaKP;
RyR_P = settings.fRyRP;
SERCA_P = settings.fPLBP;
TnI_P = settings.fTnIP;
IKb_P = settings.fIKurP; % according to O'Hara 2011 Heart Rhythm paper, background pottasium modified conductance
% % not very sure how to compute the proportion of phosphorylation on IKb channel, use the IKur computed from signaling cascade for now

% % % % settings.fICaLP = min(max(((y(98) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1);
% % % % settings.fIKsP = min(max(((y(99) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);    
% % % % settings.fPLBP = min(max((y(100) - 6.591000e-001) / (9.945000e-001 - 6.591000e-001), 0), 1);    
% % % % settings.fTnIP = min(max((y(101) - 6.735188e-001) / (9.991797e-001 - 6.735188e-001), 0), 1);    
% % % % settings.fINaP = min(max((y(102) - 2.394795e-001) / (9.501431e-001 - 2.394795e-001), 0), 1);    
% % % % settings.fINaKP = min(max((y(103) - 1.263453e-001) / (9.980137e-001 - 1.263453e-001), 0), 1);    
% % % % settings.fRyRP = min(max(((y(104) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);            
% % % % settings.fIKurP = min(max((y(105) - 5.893798e-002) / (3.937470e-001 - 5.893798e-002), 0), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% Compute ionic currents
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%

%CaMK constants
% % update so that bCaMK has 90% portion sensitive to PP1, a loose link
% % between CaMKII activity and PKA phosphorylation
% % KmCaMK=0.15;

aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;
KmCaM=0.0015;
%update CaMK
CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/Cass);
CaMKa=CaMKb+CaMKt;

betaCaMKII =  bCaMK * (0.1 + (0.9  * settings.Whole_cell_PP1 / 0.1371));
 
% fraction of CaMK that has trapped calmodulin
if isfield(settings, 'PP1block')
    dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt;
else
    dCaMKt=aCaMK*CaMKb*(CaMKb+CaMKt)-betaCaMKII*CaMKt; % This equation the same as the one in myokit ydot(7)... (checked: MES)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reversal potentials
ENa=(R*T/F)*log(Nao/Nai);
EK=(R*T/F)*log(Ko/Ki);
PKNa=0.01833;
EKs=(R*T/F)*log((Ko+PKNa*Nao)/(Ki+PKNa*Nai));

%convenient shorthand calculations
vffrt=V*F*F/(R*T);
vfrt=V*F/(R*T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate INa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gating NP
mss=1.0/(1.0+exp((-(V+39.57))/9.871));
tm=1.0/(6.765*exp((V+11.64)/34.77)+8.552*exp(-(V+77.42)/5.955));
dm=(mss-m)/tm;
hss=1.0/(1+exp((V+82.90)/6.086));
thf=1.0/(1.432e-5*exp(-(V+1.196)/6.285)+6.149*exp((V+0.5096)/20.27));
ths=1.0/(0.009794*exp(-(V+17.95)/28.05)+0.3343*exp((V+5.730)/56.66));
Ahf=0.99;
Ahs=1.0-Ahf;
dhf=(hss-hf)/thf;
dhs=(hss-hs)/ths;
h=Ahf*hf+Ahs*hs;
jss=hss;
tj=2.038+1.0/(0.02136*exp(-(V+100.6)/8.281)+0.3052*exp((V+0.9941)/38.45));
dj=(jss-j)/tj;

% gating CaMK-P
hssp=1.0/(1+exp((V+89.1)/6.086));
thsp=3.0*ths;
dhsp=(hssp-hsp)/thsp;
hp=Ahf*hf+Ahs*hsp;
tjp=1.46*tj;
djp=(jss-jp)/tjp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gating PKA-P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GNa increased after PKA-P
GNaP = GNa*2.7 ;
VP_h = 5 ; % mV works for PKA-P only and Both-P, on top of the relative baseline
             
hPss =  1.0/(1+exp((V+82.90 + VP_h)/6.086));
dhPf = (hPss - hPf)/thf ;
dhPs = (hPss - hPs)/ths ;
hP = Ahf*hPf+Ahs*hPs;

jPss = hPss ; 
djP = (jPss-jP)/tj;



hBPss = 1.0/(1+exp((V+89.1+ VP_h)/6.086)) ;
thBPf = thf ; % CaMK did not modify, keep with NP
thBPs = thsp ; % stay same as CaMK

dhBPf = (hBPss - hBPf)/thBPf ;
dhBPs = (hBPss - hBPs)/thBPs ;

hBP = Ahf*hBPf+Ahs*hBPs; % fast fraction stay same as PKA, slow fraction with modified time constant
jBPss = hBPss ; % stay same as PKA
tjBP = tjp ; % stay same as CaMK
djBP = (jBPss - jBP)/tjBP ;

% % result in adding state variables
% % hBPs, jBP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Putting together the channels behavior and fraction

fINap=(1.0/(1.0+KmCaMK/CaMKa)); % CaMK-P fraction

fINa_P = INa_P ; % PKA-P fraction as assigned as input, take the value 0 or 1
fINa_BP = fINap*fINa_P ;
fINa_CaMKonly = fINap-fINa_BP ;
fINa_PKAonly = fINa_P-fINa_BP ;


INa_NP = GNa*(V-ENa)*m^3.0*h*j ; % Non-Phosphorylated 
INa_CaMK = GNa*(V-ENa)*m^3.0*hp*jp ;
INa_PKA = GNaP*(V-ENa)*m^3.0*hP*jP ;
INa_BP = GNaP*(V-ENa)*m^3.0*hBP*jBP ;

% 4 population 
INa = (1-fINa_CaMKonly-fINa_PKAonly-fINa_BP)*INa_NP + fINa_CaMKonly*INa_CaMK + fINa_PKAonly*INa_PKA + fINa_BP*INa_BP ;

% calculate INaL
mLss=1.0/(1.0+exp((-(V+42.85))/5.264));
tmL=tm;
dmL=(mLss-mL)/tmL;
hLss=1.0/(1.0+exp((V+87.61)/7.488));
thL=200.0;
dhL=(hLss-hL)/thL;
hLssp=1.0/(1.0+exp((V+93.81)/7.488));
thLp=3.0*thL;
dhLp=(hLssp-hLp)/thLp;
fINaLp=(1.0/(1.0+KmCaMK/CaMKa));

INaL=GNaL*(V-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate Ito
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ass=1.0/(1.0+exp((-(V-14.34))/14.82));
ta=1.0515/(1.0/(1.2089*(1.0+exp(-(V-18.4099)/29.3814)))+3.5/(1.0+exp((V+100.0)/29.3814)));
da=(ass-a)/ta;
iss=1.0/(1.0+exp((V+43.94)/5.711));
if celltype==1
    delta_epi=1.0-(0.95/(1.0+exp((V+70.0)/5.0)));
else
    delta_epi=1.0;
end
tiF=4.562+1/(0.3933*exp((-(V+100.0))/100.0)+0.08004*exp((V+50.0)/16.59));
tiS=23.62+1/(0.001416*exp((-(V+96.52))/59.05)+1.780e-8*exp((V+114.1)/8.079));
tiF=tiF*delta_epi;
tiS=tiS*delta_epi;
AiF=1.0/(1.0+exp((V-213.6)/151.2));
AiS=1.0-AiF;
diF=(iss-iF)/tiF;
diS=(iss-iS)/tiS;
i=AiF*iF+AiS*iS;
assp=1.0/(1.0+exp((-(V-24.34))/14.82));
dap=(assp-ap)/ta;
dti_develop=1.354+1.0e-4/(exp((V-167.4)/15.89)+exp(-(V-12.23)/0.2154));
dti_recover=1.0-0.5/(1.0+exp((V+70.0)/20.0));
tiFp=dti_develop*dti_recover*tiF;
tiSp=dti_develop*dti_recover*tiS;
diFp=(iss-iFp)/tiFp;
diSp=(iss-iSp)/tiSp;
ip=AiF*iFp+AiS*iSp;
fItop=(1.0/(1.0+KmCaMK/CaMKa));

Ito=Gto*(V-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate ICaL, ICaNa, ICaK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gating NP 
dss=1.0/(1.0+exp((-(V+3.940))/4.230));
td=0.6+1.0/(exp(-0.05*(V+6.0))+exp(0.09*(V+14.0)));
dd=(dss-d)/td;
fss=1.0/(1.0+exp((V+19.58)/3.696));
tff=7.0+1.0/(0.0045*exp(-(V+20.0)/10.0)+0.0045*exp((V+20.0)/10.0));
tfs=1000.0+1.0/(0.000035*exp(-(V+5.0)/4.0)+0.000035*exp((V+5.0)/6.0));
Aff=0.6;
Afs=1.0-Aff;
dff=(fss-ff)/tff;
dfs=(fss-fs)/tfs;
f=Aff*ff+Afs*fs;
fcass=fss; % CDI and VDI takes the same steady state, but different time constants
tfcaf=7.0+1.0/(0.04*exp(-(V-4.0)/7.0)+0.04*exp((V-4.0)/7.0));
tfcas=100.0+1.0/(0.00012*exp(-V/3.0)+0.00012*exp(V/7.0));
Afcaf=0.3+0.6/(1.0+exp((V-10.0)/10.0));
Afcas=1.0-Afcaf;
dfcaf=(fcass-fcaf)/tfcaf;
dfcas=(fcass-fcas)/tfcas;
fca=Afcaf*fcaf+Afcas*fcas;
tjca=75.0;
djca=(fcass-jca)/tjca;

% gating CaMK-P
tffp=2.5*tff; % fast component of the VDI
dffp=(fss-ffp)/tffp;
fp=Aff*ffp+Afs*fs;
tfcafp=2.5*tfcaf; % fast component of the CDI
dfcafp=(fcass-fcafp)/tfcafp;
fcap=Afcaf*fcafp+Afcas*fcas;

Kmn=0.002;
k2n=1000.0;
km2n=jca*1.0;
anca=1.0/(k2n/km2n+(1.0+Kmn/Cass)^4.0);
dnca=anca*k2n-nca*km2n;
% PhiCaL=4.0*vffrt*(Cass*exp(2.0*vfrt)-0.341*Cao)/(exp(2.0*vfrt)-1.0);
% PhiCaL=4.0*vffrt*(Cass*exp(2.0*vfrt)-Cao)/(exp(2.0*vfrt)-1.0);
PhiCaNa=1.0*vffrt*(0.75*Nass*exp(1.0*vfrt)-0.75*Nao)/(exp(1.0*vfrt)-1.0);
PhiCaK=1.0*vffrt*(0.75*Kss*exp(1.0*vfrt)-0.75*Ko)/(exp(1.0*vfrt)-1.0);
zca=2.0;
PCap=1.1*PCa;
PCaNa=0.00125*PCa;
PCaK=3.574e-4*PCa;
PCaNap=0.00125*PCap;
PCaKp=3.574e-4*PCap;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gating PKA-P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % PCaP = 2.5*PCa ;
PCaP = GCal_mul*PCa ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCaP = 1.66*PCa ;
% PCaP = 2.5*PCa ;
PCaNaP=0.00125*PCaP;
PCaKP=3.574e-4*PCaP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% as stated in O'Hara 2011 integrated model
% % VP_d = 14.0 ;% mV, d gate, activation for ICaL
% % VP_f = 8.0 ; % mV, f gate, inactivation for ICaL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PKA-P new steady state and variables
dPss=1.0/(1.0+exp((-(V+3.940+VP_d))/4.230));
ddP = (dPss-dP)/td;
fPss=1.0/(1.0+exp((V+19.58+VP_f)/3.696)); % ddP dfPf dfPs dfcaPf dfcaPs dfBPf dfcaBPf
dfPf=(fPss-fPf)/tff;
dfPs=(fPss-fPs)/tfs;
fcaPss = fPss; % steady state update, compute CDI gating again
dfcaPf=(fcaPss-fcaPf)/tfcaf;
dfcaPs=(fcaPss-fcaPs)/tfcas;
fP = Aff*fPf+Afs*fPs;
fcaP = Afcaf*fcaPf+Afcas*fcaPs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Both-P population takes on dP, for f gate, takes ss from PKA and tau from CaMK (only fast component was modified)
fBPss = fPss;
dfBPf = (fBPss-fBPf)/tffp;
fBP = Aff*fBPf+Afs*fPs; % only fast component modified, slow component same as fPs
fcaBPss = fcaPss;
dfcaBPf = (fcaBPss-fcaBPf)/tfcafp;  
fcaBP = Afcaf*fcaBPf+Afcas*fcaPs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fICaLp=(1.0/(1.0+KmCaMK/CaMKa)); % CaMK-P fraction

fICaL_P = ICaL_P; % PKA-P fraction
fICaL_BP = fICaLp*fICaL_P;
fICaL_CaMKonly = fICaLp-fICaL_BP;
fICaL_PKAonly = fICaL_P-fICaL_BP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % this PhiCaL with ceiling may only apply to PKA-P and BP ICaL, not apply
% to NP and CaMK-P channel...
% PhiCaL=4.0*vffrt*(Cass*exp(2.0*vfrt)-0.341*Cao)/(exp(2.0*vfrt)-1.0);
% will only use for PKA-P and when Cass is exceeding 0.03
Cass_ceiling = 0.025 ; %%%%%%%%%%%%%%%%%%
if Cass > Cass_ceiling
   PhiCaL=4.0*vffrt*(Cass_ceiling*exp(2.0*vfrt)-0.341*Cao)/(exp(2.0*vfrt)-1.0);
else
   PhiCaL=4.0*vffrt*(Cass*exp(2.0*vfrt)-0.341*Cao)/(exp(2.0*vfrt)-1.0);
end
% PhiCaL=4.0*vffrt*(Cass*exp(2.0*vfrt)-0.341*Cao)/(exp(2.0*vfrt)-1.0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICaL_NP=PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca);
ICaL_CaMK = PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaL_PKA = PCaP*PhiCaL*dP*(fP*(1.0-nca)+jca*fcaP*nca);
ICaL_BP = PCaP*PhiCaL*dP*(fBP*(1.0-nca)+jca*fcaBP*nca);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICaNa_NP=PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca);
ICaNa_CaMK = PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaNa_PKA = PCaNaP*PhiCaNa*dP*(fP*(1.0-nca)+jca*fcaP*nca);
ICaNa_BP = PCaNaP*PhiCaNa*dP*(fBP*(1.0-nca)+jca*fcaBP*nca);

ICaK_NP=PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca);
ICaK_CaMK = PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);
ICaK_PKA = PCaKP*PhiCaK*dP*(fP*(1.0-nca)+jca*fcaP*nca);
ICaK_BP = PCaKP*PhiCaK*dP*(fBP*(1.0-nca)+jca*fcaBP*nca);

% 4 population
ICaL = (1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaL_NP + fICaL_CaMKonly*ICaL_CaMK + fICaL_PKAonly*ICaL_PKA + fICaL_BP*ICaL_BP;
ICaNa = (1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaNa_NP + fICaL_CaMKonly*ICaNa_CaMK + fICaL_PKAonly*ICaNa_PKA + fICaL_BP*ICaNa_BP;
ICaK = (1-fICaL_CaMKonly-fICaL_PKAonly-fICaL_BP)*ICaK_NP + fICaL_CaMKonly*ICaK_CaMK + fICaL_PKAonly*ICaK_PKA + fICaL_BP*ICaK_BP;
% 


%% calculate IKr
xrss=1.0/(1.0+exp((-(V+8.337))/6.789));
txrf=12.98+1.0/(0.3652*exp((V-31.66)/3.869)+4.123e-5*exp((-(V-47.78))/20.38));
txrs=1.865+1.0/(0.06629*exp((V-34.70)/7.355)+1.128e-5*exp((-(V-29.74))/25.94));
Axrf=1.0/(1.0+exp((V+54.81)/38.21));
Axrs=1.0-Axrf;
dxrf=(xrss-xrf)/txrf;
dxrs=(xrss-xrs)/txrs;
xr=Axrf*xrf+Axrs*xrs;
rkr=1.0/(1.0+exp((V+55.0)/75.0))*1.0/(1.0+exp((V-10.0)/30.0));
IKr=GKr*sqrt(Ko/5.4)*xr*rkr*(V-EK);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate IKs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xs1ss=1.0/(1.0+exp((-(V+11.60))/8.932));
txs1=817.3+1.0/(2.326e-4*exp((V+48.28)/17.80)+0.001292*exp((-(V+210.0))/230.0));
dxs1=(xs1ss-xs1)/txs1;
xs2ss=xs1ss;
txs2=1.0/(0.01*exp((V-50.0)/20.0)+0.0193*exp((-(V+66.54))/31.0));
dxs2=(xs2ss-xs2)/txs2;
KsCa=1.0+0.6/(1.0+(3.8e-5/Cai)^1.4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gating PKA-P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing Monica version parameters 

% GKsP = 3.3*GKs ;
% VP_xs = 0 ;
% raw_taux2 = 1 ;
% raw_taux1 = -238.8611 ;
% taux1factor = 2.5865 ;  

GKsP = GKS_mul*GKs ;
VP_xs = 0 ;
raw_taux2 = 1 ;
taux1factor = 2.75 ;  
raw_taux1 = (1-taux1factor)*1.0./(2.326e-4*exp((10+48.28)./17.80)+0.001292.*exp((-(10+210.0))./230.0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xs1Pss=1.0/(1.0+exp((-(V+VP_xs+11.60))/8.932));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xs2Pss = xs1Pss; % xs1 and xs2 takes on the same steady state

txs1P = raw_taux1 + 817.3 + taux1factor*1./((2.326*10^(-4)*exp((V+48.28)/17.8)) + 0.001292*exp(-(V+210.0)/230.0)) ;
txs2P = txs2*raw_taux2 ;
dxs1P=(xs1Pss-xs1P)/txs1P;
dxs2P=(xs2Pss-xs2P)/txs2P;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% no CaMK-P in IKs
IKs_NP=GKs*KsCa*xs1*xs2*(V-EKs);
IKs_PKA = GKsP*KsCa*xs1P*xs2P*(V-EKs);
fIKs_PKA = IKs_P;
IKs = (1-fIKs_PKA)*IKs_NP + fIKs_PKA*IKs_PKA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate IK1
xk1ss=1.0/(1.0+exp(-(V+2.5538*Ko+144.59)/(1.5692*Ko+3.8115)));
txk1=122.2/(exp((-(V+127.2))/20.36)+exp((V+236.8)/69.33));
dxk1=(xk1ss-xk1)/txk1;
rk1=1.0/(1.0+exp((V+105.8-2.6*Ko)/9.493));
IK1=GK1*sqrt(Ko)*rk1*xk1*(V-EK);

%% calculate INaCa_i
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
hca=exp((qca*V*F)/(R*T));
hna=exp((qna*V*F)/(R*T));
h1=1+Nai/kna3*(1+hna);
h2=(Nai*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+Nai/kna1*(1+Nai/kna2);
h5=Nai*Nai/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+Nao/kna3*(1.0+1.0/hna);
h8=Nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+Nao/kna1*(1.0+Nao/kna2);
h11=Nao*Nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*Cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*Cai*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0/(1.0+(KmCaAct/Cai)^2.0);
zna=1.0;
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

%% calculate INaCa_ss
h1=1+Nass/kna3*(1+hna);
h2=(Nass*hna)/(kna3*h1);
h3=1.0/h1;
h4=1.0+Nass/kna1*(1+Nass/kna2);
h5=Nass*Nass/(h4*kna1*kna2);
h6=1.0/h4;
h7=1.0+Nao/kna3*(1.0+1.0/hna);
h8=Nao/(kna3*hna*h7);
h9=1.0/h7;
h10=kasymm+1.0+Nao/kna1*(1+Nao/kna2);
h11=Nao*Nao/(h10*kna1*kna2);
h12=1.0/h10;
k1=h12*Cao*kcaon;
k2=kcaoff;
k3p=h9*wca;
k3pp=h8*wnaca;
k3=k3p+k3pp;
k4p=h3*wca/hca;
k4pp=h2*wnaca;
k4=k4p+k4pp;
k5=kcaoff;
k6=h6*Cass*kcaon;
k7=h5*h2*wna;
k8=h8*h11*wna;
x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
KmCaAct=150.0e-6;
allo=1.0/(1.0+(KmCaAct/Cass)^2.0);
JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
JncxCa=E2*k2-E1*k1;
INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate INaK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not a substrate of CaMK only 2 populations, NP and PKA-P
% in this case can simply introduce 
k1p=949.5;
k1m=182.4;
k2p=687.2;
k2m=39.4;
k3p=1899.0;
k3m=79300.0;
k4p=639.0;
k4m=40.0;

Knai0_np = 9.073;
Knai0_P = Knai0_np*0.7;
Knai0=9.073;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fINaK_PKA = INaK_P;
if fINaK_PKA == 0
   Knai0=9.073;
else
   Knai0=(1-fINaK_PKA)*Knai0_np + fINaK_PKA*Knai0_P; % PKA-P effetc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Knao0=27.78;
delta=-0.1550;
Knai=Knai0*exp((delta*V*F)/(3.0*R*T));
Knao=Knao0*exp(((1.0-delta)*V*F)/(3.0*R*T));
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

P=eP/(1.0+H/Khp+Nai/Knap+Ki/Kxkur);
a1=(k1p*(Nai/Knai)^3.0)/((1.0+Nai/Knai)^3.0+(1.0+Ki/Kki)^2.0-1.0);
b1=k1m*MgADP;
a2=k2p;
b2=(k2m*(Nao/Knao)^3.0)/((1.0+Nao/Knao)^3.0+(1.0+Ko/Kko)^2.0-1.0);
a3=(k3p*(Ko/Kko)^2.0)/((1.0+Nao/Knao)^3.0+(1.0+Ko/Kko)^2.0-1.0);
b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
b4=(k4m*(Ki/Kki)^2.0)/((1.0+Nai/Knai)^3.0+(1.0+Ki/Kki)^2.0-1.0);
x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
E1=x1/(x1+x2+x3+x4);
E2=x2/(x1+x2+x3+x4);
E3=x3/(x1+x2+x3+x4);
E4=x4/(x1+x2+x3+x4);
zk=1.0;
JnakNa=3.0*(E1*a3-E2*b3);
JnakK=2.0*(E4*b1-E3*a1);
% short cut here cause INaK_P will only be given 0 or 1
INaK=Pnak*(zna*JnakNa+zk*JnakK);

%% background current and membrane pumps
%calculate IKb
xkb=1.0/(1.0+exp(-(V-14.48)/18.34));

%%%%%%%%%%%%%%%%%%%%
% added PKA effects according to O'Hara 2011
IKb_NP=GKb*xkb*(V-EK);
GKbP = GKb*2.5 ;
IKb_PKA = GKbP*xkb*(V-EK);

fIKb_P = IKb_P;
IKb = (1-fIKb_P)*IKb_NP + fIKb_P*IKb_PKA ;
%%%%%%%%%%%%%%%%%%%%

%calculate INab
INab=PNab*vffrt*(Nai*exp(vfrt)-Nao)/(exp(vfrt)-1.0);

%calculate ICab
ICab=PCab*4.0*vffrt*(Cai*exp(2.0*vfrt)-0.341*Cao)/(exp(2.0*vfrt)-1.0);

%calculate IpCa
IpCa=GpCa*Cai/(0.0005+Cai);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%update the membrane voltage
dv=-(1/Cm)*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate diffusion fluxes
JdiffNa=(Nass-Nai)/2.0;
JdiffK=(Kss-Ki)/2.0;
Jdiff=(Cass-Cai)/0.2;

%% calculate ryanodione receptor calcium induced calcium release from the jsr

% NP channel behavior
bt=4.75;
a_rel=0.5*bt; % relation between a_rel and bt always holds
Jrel_inf=a_rel*(-ICaL)/(1.0+(1.5/Cajsr)^8.0);
if celltype==2
    Jrel_inf=Jrel_inf*1.7;
end
tau_rel=bt/(1.0+0.0123/Cajsr);

if tau_rel<0.001
   tau_rel=0.001; 
end

dJrelnp=(Jrel_inf-Jrelnp)/tau_rel;

% CaMK-P channel behavior
btp=1.25*bt;
a_relp=0.5*btp; % relation between a_rel and bt always holds
Jrel_infp=a_relp*(-ICaL)/(1.0+(1.5/Cajsr)^8.0);
if celltype==2
    Jrel_infp=Jrel_infp*1.7;
end
tau_relp=btp/(1.0+0.0123/Cajsr);

if tau_relp<0.001
   tau_relp=0.001; 
end

dJrelp=(Jrel_infp-Jrelp)/tau_relp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PKA-P channel behavior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_relP = a_rel*2.5; % relation between a_rel and bt always holds

Jrel_infP=a_relP*(-ICaL)/(1.0+(1.5/Cajsr)^8.0);
if celltype==2
    Jrel_infP=Jrel_infP*1.7;
end
tau_relP=tau_rel*0.75 ;

if tau_relP<0.001
   tau_relP=0.001; 
end

dJrelP=(Jrel_infP-JrelP)/tau_relP; % a new variable JrelP introduced
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Both-P
% a_relBP = a_relp*1.75; % relation between a_rel and bt always holds
a_relBP = a_relp*2.5; % relation between a_rel and bt always holds

Jrel_infBP=a_relBP*(-ICaL)/(1.0+(1.5/Cajsr)^8.0);
if celltype==2
    Jrel_infBP=Jrel_infBP*1.7;
end
tau_relBP=tau_relp*0.75 ;

if tau_relBP<0.001
   tau_relBP=0.001; 
end

dJrelBP=(Jrel_infBP-JrelBP)/tau_relBP; % a new variable JrelBP introduced

    
% 4 population, with both-P population takes on PKA behavior    
fJrelp=(1.0/(1.0+KmCaMK/CaMKa)); % CaMK-P fraction
% fJrelp = fCaMK ;
% % fJrelp = fCaMK_RyR ;
fJrel_PKA = RyR_P;
fJrel_BP = fJrelp*fJrel_PKA;
fJrel_CaMKonly = fJrelp - fJrel_BP;
fJrel_PKAonly = fJrel_PKA - fJrel_BP;

Jrel= RyR_total* ((1.0-fJrel_CaMKonly-fJrel_PKAonly-fJrel_BP)*Jrelnp + fJrel_CaMKonly*Jrelp + fJrel_PKAonly*JrelP + fJrel_BP*JrelBP );



%% calculate serca pump, ca uptake flux
% NP SERCA
Jupnp=0.004375*Cai/(Cai+0.00092);
% CaMK-P SERCA 
Jupp=2.75*0.004375*Cai/(Cai+0.00092-0.00017);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PKA-P SERCA
% Vmax stay unchanged, KmPLB = 4.97 mM
% % Both-P population, takes on CaMK Vmax with PKA KmPLB
mul_PKA = 0.54;
JupP = 0.004375*Cai/(Cai+0.00092*mul_PKA);
Jup_BP = 2.75*0.004375*Cai/(Cai+(0.00092-0.00017)*mul_PKA);

if celltype==1
    Jupnp=Jupnp*1.3;
    Jupp=Jupp*1.3;
    JupP=JupP*1.3;
    Jup_BP=Jup_BP*1.3;
end

% 4 population
fJupp=(1.0/(1.0+KmCaMK/CaMKa)); % CaMK-P channel fraction
% fJupp = fCaMK ;
fJup_P = SERCA_P; % PKA-P channel fraction
fJup_BP = fJupp*fJup_P;
fJup_CaMKonly = fJupp - fJup_BP;
fJup_PKAonly = fJup_P - fJup_BP;


Jleak=Leak_total*(0.0039375*Cansr/15.0);
Jup= SERCA_total*((1.0-fJup_CaMKonly-fJup_PKAonly-fJup_BP)*Jupnp + fJup_CaMKonly*Jupp + fJup_PKAonly*JupP + fJup_BP*Jup_BP) - Jleak;

%calculate tranlocation flux
Jtr=Trans_total*((Cansr-Cajsr)/100.0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calcium buffer constants
% similar to INaK take the short cut here

cmdnmax=0.05;
if celltype==1
    cmdnmax=cmdnmax*1.3;
end
kmcmdn=0.00238;
trpnmax=0.07;

kmtrpn_np = 0.0005;
kmtrpn_P = kmtrpn_np*1.6 ;

kmtrpn = 0.0005;
fTnI_PKA = TnI_P;
if fTnI_PKA == 0
   kmtrpn=0.0005;
else
   kmtrpn=(1-fTnI_PKA)*kmtrpn_np + fTnI_PKA*kmtrpn_P;
end
BSRmax=0.047;
KmBSR=0.00087;
BSLmax=1.124;
KmBSL=0.0087;
csqnmax=10.0;
kmcsqn=0.8;

%update intracellular concentrations, using buffers for cai, cass, cajsr
dnai=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
dnass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;

dki=-(Ito+IKr+IKs+IK1+IKb+Id-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
dkss=-(ICaK)*Acap/(F*vss)-JdiffK;

Bcai=1.0/(1.0+cmdnmax*kmcmdn/(kmcmdn+Cai)^2.0+trpnmax*kmtrpn/(kmtrpn+Cai)^2.0);
dcai=Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);

Bcass=1.0/(1.0+BSRmax*KmBSR/(KmBSR+Cass)^2.0+BSLmax*KmBSL/(KmBSL+Cass)^2.0);
dcass=Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);

dcansr=Jup-Jtr*vjsr/vnsr;

Bcajsr=1.0/(1.0+csqnmax*kmcsqn/(kmcsqn+Cajsr)^2.0);
dcajsr=Bcajsr*(Jtr-Jrel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deriv=[dv dnai dnass dki dkss dcai dcass dcansr dcajsr dm dhf dhs ...
  dj dhsp djp dmL dhL dhLp da diF diS dap diFp diSp dd dff dfs dfcaf ...
  dfcas djca dnca dffp dfcafp dxrf dxrs dxs1 dxs2 dxk1 dJrelnp dJrelp dJrelBP dCaMKt... % all the PKA variables following
  dhPf dhPs djP dhBPf dhBPs djBP ddP dfPf dfPs dfcaPf dfcaPs dfBPf dfcaBPf...
  dxs1P dxs2P dJrelP ]';


%%  COMPARING STATE VARIABLES IN THIS CODE VS ORD CODE
% in Jingqi Code        In ORD code
%     1: dv             v=X(1);
%     2: dnai           nai=X(2);
%     3: dnass          nass=X(3);
%     4: dki            ki=X(4);
%     5: dkss           kss=X(5);
%     6: dcai           cai=X(6);
%     7: dcass          cass=X(7);
%     8: dcansr         cansr=X(8);    
%     9: dcajsr         cajsr=X(9);    
%     10: dm            m=X(10);
%     11: dhf           hf=X(11);
%     12: dhs ...       hs=X(12);    
%     13: dj            j=X(13);
%     14: dhsp          hsp=X(14);
%     15: djp           jp=X(15);
%     16: dmL           mL=X(16);
%     17: dhL           hL=X(17);
%     18: dhLp          hLp=X(18);
%     19: da            a=X(19);
%     20: diF           iF=X(20);   
%     21: diS           iS=X(21);
%     22: dap           ap=X(22);
%     23: diFp          iFp=X(23);
%     24: diSp          iSp=X(24);
%     25: dd            d=X(25);
%     26: dff           ff=X(26);
%     27: dfs           fs=X(27);
%     28: dfcaf ...     fcaf=X(28);
%     29: dfcas         fcas=X(29);    
%     30: djca          jca=X(30);
%     31: dnca          nca=X(31);
%     32: dffp          ffp=X(32);
%     33: dfcafp        fcafp=X(33);
%     34: dxrf          xrf=X(34);
%     35: dxrs          xrs=X(35);
%     36: dxs1          xs1=X(36);
%     37: dxs2          xs2=X(37);
%     38: dxk1          xk1=X(38);
%     39: dJrelnp       Jrelnp=X(39);
%     40: dJrelp        Jrelp=X(40);  
%     41: dJrelBP       --------------------------- This is a new variable introduced 
%     42: dCaMKt...     CaMKt=X(41);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% all the PKA variables following
%     43: dhPf          
%     44: dhPs          
%     45: djP           
%     46: dhBPf             
%     47: dhBPs             
%     48: djBP          
%     49: ddP           
%     50: dfPf          
%     51: dfPs          
%     52: dfcaPf            
%     53: dfcaPs            
%     54: dfBPf             
%     55: dfcaBPf...        
%     56: dxs1P             
%     57: dxs2P             
%     58: dJrelP ]';        





