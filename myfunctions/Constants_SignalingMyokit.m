%
% Constants for Heijman2011
%
function c = Constants_SignalingMyokit(iso_conc)


%%

% iso
c.iso_L = iso_conc; % concentration of isoproterenol
% IBMX
c.ibmx = 0.0;% Concentration of IBMX (in the range 0 to 100) not used in OharaBA

% phys
c.phys_T = 310.0;% Temperature
c.F = 96487.0;% Faraday's constant
c.phys_R = 8314.0;% Gas constant


%% Cell geometry
c.length = 0.01;% Cell length
c.cell_pi =  3.14159265358979312e+00;%pi
c.radius = 0.0011;% Cell radius
c.geoArea = 2.0 * c.cell_pi * c.radius * (c.radius + c.length);% Geometric membrane area
c.volume = 1000.0 * c.cell_pi * c.radius * c.radius * c.length;% Cell volume

c.v_cav = 0.02 * c.volume;% Volume of the caveolar subspace
c.v_eca = 0.04 * c.volume;% Volume of the extracaveolar subspace
c.v_cyt = c.volume * 0.678;% Volume of the Cytoplasm / Myoplasm

c.vr_cav = c.volume / c.v_cav;% Ratio of whole volume to caveolar subspace volume
c.vr_eca = c.volume / c.v_eca;% Ratio of whole volume to extracaveolar subspace volume
c.vr_cyt = c.volume / c.v_cyt;% Ratio of whole volume to cytoplasm volume


%% cAMP
c.j_cav_eca = 5e-15 * 1000000.0;% Rate of cAMP diffusion between caveolar and cytosolic compartments
c.j_cav_cyt = 7.5e-14 * 1000000.0;% Rate of cAMP diffusion between caveolar and extracaveolar compartments
c.j_eca_cyt = 9e-15 * 1000000.0;% Rate of cAMP diffusion between extracaveolar and cytosolic compartments


%% PKA (See pg 29 - 32 of Heijman supplementary document)
c.PKA_tot = 0.5;                    % Total cellular concentration of PKA holoenzyme (umol/L)
c.f_cav = 0.0388;                   % Fraction of PKA located in caveolar compartment
c.f_eca = 0.1;                      % Fraction of PKA located in extracaveolar compartment
c.f_cyt = 1.0 - c.f_cav - c.f_eca;  % Fraction of PKA located in cytosolic compartment

c.PKA_eca = c.f_eca * c.PKA_tot * c.vr_eca;% Concentration of PKA in the extracaveolar compartment
c.PKA_cav = c.f_cav * c.PKA_tot * c.vr_cav;% Concentration of PKA in the caveolar compartment
c.PKA_cyt = c.f_cyt * c.PKA_tot * c.vr_cyt;% Concentration of PKA in the cytosolic compartment

c.PKI_tot = 0.2 * c.PKA_tot;% Total cellular concentration of PKA inhibitor
c.f_pki_cav = c.f_cav;              % Fraction of PKI located in caveolar compartment
c.f_pki_eca = c.f_eca;              % Fraction of PKI located in extracaveolar compartment
c.f_pki_cyt = 1.0 - c.f_pki_cav - c.f_pki_eca;% Fraction of PKI located in cytosolic compartment

c.PKI_cav = c.f_pki_cav * c.PKI_tot * c.vr_cav;% Concentration of protein kinase inhibitor in caveolar compartment
c.PKI_eca = c.f_pki_eca * c.PKI_tot * c.vr_eca;% Concentration of protein kinase inhibitor in extracaveolar compartment
c.PKI_cyt = c.f_pki_cyt * c.PKI_tot * c.vr_cyt;% Concentration of protein kinase inhibitor in cytosolic compartment

c.f_pki = 50.0;                     % Forward rate for inhibition of C subunit by PKI
c.K_pki = 0.01 / 50.0;              % Equilibrium value for inhibition of C subunit by PKI (umol/L)
c.b_pki = c.f_pki * c.K_pki;        % Backward rate for inhibition of C subunit by PKI

% pka_cav
c.pka_cav_f1 = 100.0;% Caveolar forward rate for binding of the first cAMP to PKA
c.pka_cav_f2 = 100.0;% Caveolar forward rate for binding of the second cAMP to PKA
c.pka_cav_f3 = 100.0;% Caveolar forward rate for dissociation of C subunit

c.pka_cav_K1 = 2.4984;% Caveolar equilibrium value for the binding of the first cAMP to PKA
c.pka_cav_K2 = 11.359;% Caveolar equilibrium value for the binding of the second cAMP to PKA
c.pka_cav_K3 = 0.3755;% Caveolar equilibrium value for dissociation of C subunit

c.pka_cav_b1 = c.pka_cav_f1 * c.pka_cav_K1;% Caveolar backward rate for binding of the first cAMP to PKA
c.pka_cav_b2 = c.pka_cav_f2 * c.pka_cav_K2;% Caveolar backward rate for binding of the second cAMP to PKA
c.pka_cav_b3 = c.pka_cav_f3 * c.pka_cav_K3;% Caveolar backward rate for dissociation of C subunit

% pka_eca
c.pka_eca_f1 = c.pka_cav_f1;% Extracaveolar forward rate for binding of the first cAMP to PKA
c.pka_eca_f2 = c.pka_cav_f2;% Extracaveolar forward rate for binding of the second cAMP to PKA
c.pka_eca_f3 = c.pka_cav_f3;% Extracaveolar forward rate for dissociation of C subunit

c.pka_eca_K1 = c.pka_cav_K1;% Extracaveolar equilibrium value for the binding of the first cAMP to PKA
c.pka_eca_K2 = c.pka_cav_K2;% Extracaveolar equilibrium value for the binding of the second cAMP to PKA
c.pka_eca_K3 = c.pka_cav_K3;% Extracaveolar equilibrium value for dissociation of C subunit

c.pka_eca_b1 = c.pka_eca_f1 * c.pka_eca_K1;% Extracaveolar backward rate for binding of the first cAMP to PKA
c.pka_eca_b2 = c.pka_eca_f2 * c.pka_eca_K2;% Extracaveolar backward rate for binding of the second cAMP to PKA
c.pka_eca_b3 = c.pka_eca_f3 * c.pka_eca_K3;% Extracaveolar backward rate for dissociation of C subunit

% pka_cyt
c.pka_cyt_f1 = c.pka_cav_f1;% Cytosolic forward rate for binding of the first cAMP to PKA
c.pka_cyt_f2 = c.pka_cav_f2;% Cytosolic forward rate for binding of the second cAMP to PKA
c.pka_cyt_f3 = c.pka_cav_f3;% Cytosolic forward rate for dissociation of C subunit

c.pka_cyt_K1 = 0.1088;% Cytosolic equilibrium value for the binding of the first cAMP to PKA
c.pka_cyt_K2 = 0.4612;% Cytosolic equilibrium value for binding of the second cAMP to PKA
c.pka_cyt_K3 = 0.3755;% Cytosolic equilibrium value for dissociation of C subunit

c.pka_cyt_b1 = c.pka_cyt_f1 * c.pka_cyt_K1;% Cytosolic backward rate for binding of the first cAMP to PKA
c.pka_cyt_b2 = c.pka_cyt_f2 * c.pka_cyt_K2;% Cytosolic backward rate for binding of the second cAMP to PKA
c.pka_cyt_b3 = c.pka_cyt_f3 * c.pka_cyt_K3;% Cytosolic backward rate for dissociation of C subunit

%% PP1 (see pg 32-33 of Heijman Supplementary document)
c.f = 0.3;% Fractional increase in PP1 after Inh1 knockout
c.kp = 0.010145;% Rate of phosphorylation of inhibitor 1 by PKA
c.kdp = 0.0035731;% Rate of dephosphorylation if inhibitor 1
c.Kp = 0.001469;% Affinity of inhibitor 1 for PKA catalytic subunit
c.Kdp =  1.95259999999999991e-05;% Affinity of inhibitor 1 for PP2A

c.PP1_eca = 0.1;% PP1 concentration in the extracaveolar compartment
c.PP1_cav = 0.25;% PP1 concentration in the caveolar compartment
c.PP1_cyt = 0.2;% PP1 concentration in the cytosolic compartment

c.pp1_K = 0.001;% Affinity foor PP1 Inhibitor 1 binding

c.inhib1_tot = c.f / (1.0 - c.f) * c.pp1_K + c.f * c.PP1_cyt;% Concentration of phosphatase inhibitor 1 in the cytosolic compartment
c.PP2A = 1.0; % PP2A concentration? 


%% PDE (see pg 26 - 29 of Heijman supplementary document)
c.PDE2_tot = 0.029268;% Total cellular concentration of PDE2
c.f_pde2_cav = 0.16957;% Fraction of PDE2 located in caveolar compartment
c.f_pde2_eca =  2.12570000000000006e-04;% Fraction of PDE2 located in extracaveolar compartment
c.f_pde2_cyt = 1.0 - c.f_pde2_cav - c.f_pde2_eca;% Fraction of PDE2 located in cytosolic compartment
c.f_pde2_part = c.f_pde2_cav + c.f_pde2_eca;% Fraction of PDE2 located in the "particulate fraction" (cav + eca)

c.f_pde_part = 0.2;% Fraction of total PDE located in the "particulate fraction" (cav + eca)
c.f_pde4_part = 0.125;% Fraction of PDE4 located in the "particulate fraction" (cav + eca)

c.f_pde4_cav = 0.12481;% Fraction of PDE4 in caveolar compartment
c.f_pde4_eca = c.f_pde4_part - c.f_pde4_cav;% Fraction of PDE4 located in the extracaveolar compartment
c.f_pde4_cyt = 1.0 - c.f_pde4_part;% Fraction of PDE4 in cytosolic compartment

c.kPDE2 = 20.0;% Rate of cAMP hydrolysis by PDE2
c.kPDE3 = 2.5;% Rate of cAMP hydrolysis by PDE3
c.kPDE4 = 4.0;% Rate of cAMP hydrolysis by PDE4

c.KmPDE2 = 50.0;% Affinity of PDE2 for cAMP
c.KmPDE3 = 0.8;% Affinity of PDE3 for cAMP
c.KmPDE4 = 1.4;% Affinity of PDE4 for cAMP

c.KmIbmxPde2 = 21.58;% Affinity of IBMX for PDE2
c.h_ibmx_pde2 = 1.167;% Hill coefficient for inhibition of PDE2 by IBMX
c.h_ibmx_pde3 = 0.7629;% Hill coefficient for inhibition of PDE3 by IBMX
c.h_ibmx_pde4 = 0.9024;% Hill coefficient for inhibition of PDE4 by IBMX

c.KmIbmxPde3 = 2.642;% Affinity of IBMX for PDE3
c.KmIbmxPde4 = 11.89;% Affinity of IBMX for PDE4

c.KPDEp = 0.52218;
c.delta_k_pde34 = 3.0;% Increase in PDE3 / PDE4 activity after phosphorylation
c.ff_pde3_cyt = 0.35;% Fraction of PDE in cytosol that is of type 3
c.kfPDEp = 0.0196;% Rate of phosphorylation by PKA of PDE3 and PDE4
c.r_pde34_frac = 3.71;% Ratio of PDE3 to PDE4 in particulate fraction (78:21)
c.r_pde3_cyt = c.ff_pde3_cyt / (1.0 - c.ff_pde3_cyt);% Relative contribution of PDE3 to cytosolic PDE3
c.kbPDEp = c.KPDEp * c.kfPDEp;% Rate of dephosphorylation of PDE3 and PDE4

c.pde_PDE3_tot_alpha = c.r_pde3_cyt * (c.f_pde4_part * (1.0 + c.r_pde34_frac - c.r_pde34_frac * c.f_pde2_part - c.f_pde_part) + c.f_pde2_part * (c.f_pde_part - 1.0)) + c.r_pde34_frac * c.f_pde4_part * (c.f_pde_part - c.f_pde2_part);
c.pde_PDE3_tot_beta = c.f_pde4_part * (1.0 + c.r_pde34_frac + c.f_pde_part * (c.r_pde3_cyt - c.r_pde34_frac)) - c.f_pde_part * (1.0 + c.r_pde3_cyt);

c.PDE3_tot = c.pde_PDE3_tot_alpha / c.pde_PDE3_tot_beta * c.PDE2_tot;% Total cellular concentration of PDE3
c.PDE4_tot = ((c.f_pde_part - c.f_pde2_part) * c.PDE2_tot + c.f_pde_part * c.PDE3_tot) / ((1.0 + c.r_pde34_frac) * c.f_pde4_part - c.f_pde_part);% Total cellular concentration of PDE4

c.ibmx_h2 = c.ibmx ^ c.h_ibmx_pde2;
c.ibmx_h3 = c.ibmx ^ c.h_ibmx_pde3;
c.ibmx_h4 = c.ibmx ^ c.h_ibmx_pde4;

c.ibmx2 = (1.0 - c.ibmx_h2 / (c.KmIbmxPde2 + c.ibmx_h2)) * c.PDE2_tot;
c.ibmx3 = (1.0 - c.ibmx_h3 / (c.KmIbmxPde3 + c.ibmx_h3)) * c.PDE3_tot;
c.ibmx4 = (1.0 - c.ibmx_h4 / (c.KmIbmxPde4 + c.ibmx_h4)) * c.PDE4_tot;

c.f_pde3_cav = c.r_pde34_frac * c.f_pde4_part * c.PDE4_tot / c.PDE3_tot;
c.f_pde3_cyt = 1.0 - c.f_pde3_cav;

c.PDE2_cav = c.ibmx2 * c.f_pde2_cav * c.vr_cav;% Concentration of PDE2 in Caveolar subspace
c.PDE2_eca = c.ibmx2 * c.f_pde2_eca * c.vr_eca;% Concentration of PDE2 in Extracaveolar subspace
c.PDE2_cyt = c.ibmx2 * c.f_pde2_cyt * c.vr_cyt;% Concentration of PDE2 in Cytosolic subspace
c.PDE3_tot = c.pde_PDE3_tot_alpha / c.pde_PDE3_tot_beta * c.PDE2_tot;% Total cellular concentration of PDE3

c.PDE3_tot = c.pde_PDE3_tot_alpha / c.pde_PDE3_tot_beta * c.PDE2_tot;% Total cellular concentration of PDE3
c.PDE3_cav = c.ibmx3 * c.f_pde3_cav * c.vr_cav;% Concentration of PDE3 in Caveolar subspace
c.PDE3_cyt = c.ibmx3 * c.f_pde3_cyt * c.vr_cyt;% Concentration of PDE3 in Cytosolic subspace

c.PDE4_cav = c.ibmx4 * c.f_pde4_cav * c.vr_cav;% Concentration of PDE4 in Caveolar subspace
c.PDE4_eca = c.ibmx4 * c.f_pde4_eca * c.vr_eca;% Concentration of PDE4 in Extracaveolar subspace
c.PDE4_cyt = c.ibmx4 * c.f_pde4_cyt * c.vr_cyt;% Concentration of PDE4 in Cytosolic subspace


%% Adrenergic Receptor and G Protein Activation (pg 18 - 24 of Heijman supplementary document) 
c.beta_R_b1_tot = 0.85 * 0.025;         % Total cellular beta-1 adrenergic receptor concentration
c.beta_R_b2_tot = 0.15 * 0.025;% Total cellular beta-2 adrenergic receptor concentration
c.Gs_tot = 224.0 * c.beta_R_b1_tot;     % Total Gs protein concentration
c.Gi_tot = 0.5;                         % Total Gi protein concentration

c.f_Gs_eca = 0.5664;% Fraction of Gs proteins located in extracaveolar space
c.f_Gs_cav = 0.0011071;% Fraction of Gs proteins located in caveolar subspace
c.f_Gs_cyt = 1.0 - c.f_Gs_cav - c.f_Gs_eca;% Fraction of Gs proteins located in cytoplasm

c.f_Gi_cav = 0.85;% Fraction of Gi proteins located in caveolar subspace
c.f_Gi_eca = 1.0 - c.f_Gi_cav;% Fraction of Gi proteins located in extracaveolar space

c.f_Rb1_cav = 0.081161;% Fraction of beta-1 adrenergic receptors located in caveolar subspace
c.f_Rb1_eca = 0.48744;% Fraction of beta-1 adrenergic receptors located in extra caveolar space
c.f_Rb1_cyt = 1.0 - c.f_Rb1_cav - c.f_Rb1_eca;% Fraction of beta-1 adrenergic receptors located in cytoplasm

c.f_Rb2_cav = 0.85;% Fraction of beta-2 adrenergic receptors located in caveolar subspace
c.f_Rb2_eca = 1.0 - c.f_Rb2_cav;% Fraction of beta-2 adrenergic receptors located in extracaveolar space

c.k_b1_l = 0.567;% beta-1 receptor / ligand low affinity constant
c.k_b1_c = 2.449;% beta-1 receptor / G-protein affinity constant
c.k_b1_h = 0.062;% beta-1 receptor / ligand high affinity constant

c.k_b2_n = 1.053;% Phosph. beta-2 receptor / ligand low affinity constant
c.k_b2_h = 0.012;% beta-2 receptor / ligand high affinity constant
c.k_b2_a = 1.6655;% Phosph. beta-2 receptor / Gi-protein affinity constant
c.k_b2_c = 1.8463;% beta-2 receptor / G-protein affinity constant
c.k_b2_l = 1.053;% beta-2 receptor / ligand low affinity constant
c.k_b2_f = 0.1;% Phosph. beta-2 receptor / ligand high affinity constant

c.k_act1_Gs = 4.9054;% Activation rate for Gs by high affinity complex
c.k_act2_Gs = 0.25945;% Activation rate for Gs by low affinity complex

c.k_act1_Gi = 4.0;% Activation rate for Gi by high affinity complex
c.k_act2_Gi = 0.05;% Activation rate for Gi by low affinity complex

c.k_hydr_Gs = 0.8;% Gs GTP to GDP hydrolysis constant
c.k_hydr_Gi = c.k_hydr_Gs;% Gi GTP to GDP hydrolysis constant

c.k_reas_Gs = 1210.0;% Reassociation rate for Gs subunits
c.k_reas_Gi = c.k_reas_Gs;% Reassociation rate for Gi subunits

c.rate_bds = 0.35;

c.k_grk_dp = c.rate_bds * 0.0009833;% Rate for GRK dephosphorylation
c.k_grk_p = c.rate_bds * 0.00133;% Rate for GRK dependent receptor desensitization

c.k_pka_p = c.rate_bds * 0.0065;% Rate for (PKA dependent receptor) desensitization
c.k_pka_dp = 0.15629 * c.k_pka_p;% Rate for PKA dephosphorylation

% beta_cav
c.beta_cav_GRK = 1.0;
c.beta_cav_R_b1_tot = c.f_Rb1_cav * c.beta_R_b1_tot * c.vr_cav;% Total concentration of beta1AR in the caveolar subspace
c.beta_cav_k_GsAct_b2 = 1.0;% Ratio of beta2AR incorporated in RGs_Tot and LRGs_Tot
c.beta_cav_R_b2_tot = c.f_Rb2_cav * c.beta_R_b2_tot * c.vr_cav;% Total concentration of beta2AR in the caveolar subspace
c.beta_cav_Rb2_pka_f_a = (c.k_b2_f + c.iso_L) * (c.k_b2_n + c.iso_L) / c.k_b2_n;
c.beta_cav_Gs_f_c22 = c.k_b2_c * c.k_b2_h * c.k_b1_l * (c.k_b1_h + c.iso_L) * (c.k_b2_l + c.iso_L);
c.beta_cav_Gs_f_a = c.k_b1_l * c.k_b2_l * (c.k_b1_h + c.iso_L) * (c.k_b2_h + c.iso_L);
c.beta_cav_Gs_f_c33 = c.k_b1_c * c.k_b2_c * c.k_b1_h * c.k_b2_h * (c.k_b1_l + c.iso_L) * (c.k_b2_l + c.iso_L);
c.beta_cav_Gs_f_c11 = c.k_b1_c * c.k_b1_h * c.k_b2_l * (c.k_b2_h + c.iso_L) * (c.k_b1_l + c.iso_L);

% beta_eca
c.beta_eca_GRK = 1.0;
c.beta_eca_k_GsAct_b2 = 1.0;% Ratio of beta2AR incorporated in RGs_Tot and LRGs_Tot
c.beta_eca_R_b2_tot = c.f_Rb2_eca * c.beta_R_b2_tot * c.vr_eca;% Total concentration of beta2AR in the extracaveolar space
c.beta_eca_R_b1_tot = c.f_Rb1_eca * c.beta_R_b1_tot * c.vr_eca;% Total concentration of beta1AR in the extracaveolar space
c.beta_eca_Rb2_pka_f_a = (c.k_b2_f + c.iso_L) * (c.k_b2_n + c.iso_L) / c.k_b2_n;
c.beta_eca_Gs_f_c11 = c.k_b1_c * c.k_b1_h * c.k_b2_l * (c.k_b2_h + c.iso_L) * (c.k_b1_l + c.iso_L);
c.beta_eca_Gs_f_c33 = c.k_b1_c * c.k_b2_c * c.k_b1_h * c.k_b2_h * (c.k_b1_l + c.iso_L) * (c.k_b2_l + c.iso_L);
c.beta_eca_Gs_f_a = c.k_b1_l * c.k_b2_l * (c.k_b1_h + c.iso_L) * (c.k_b2_h + c.iso_L);
c.beta_eca_Gs_f_c22 = c.k_b2_c * c.k_b2_h * c.k_b1_l * (c.k_b1_h + c.iso_L) * (c.k_b2_l + c.iso_L);

% beta_cyt
c.beta_cyt_R_b1_tot = c.f_Rb1_cyt * c.beta_R_b1_tot * c.vr_cyt;% Total concentration of beta-1 AR in the cytoplasm
c.beta_cyt_GRK = 1.0;
c.beta_cyt_Rb1_np_f_a = (c.k_b1_h + c.iso_L) * (c.k_b1_l + c.iso_L) / c.k_b1_l;


%% AC (see pg 24-26 of Heijman supplementary document)
c.ATP = 5000.0;% Concentration of ATP
c.KmATP = 315.0;% AC affinity for ATP
c.AC_tot = 3.0 * c.beta_R_b1_tot;% Total cellular AC concentration

c.f_AC47_eca = 0.16479;% Fraction of AC47 in extracaveolar space
c.f_AC56_cav = 0.087459;% Fraction of AC56 located in caveolae
c.f_AC56_AC47 = 1.0 / (1.0 + 0.35);% Fraction of AC that is of type 5/6

c.hGsAC47 = 1.0043;% Hill coefficient for AC47 activation
c.hGsAC56 = 1.3574;% Hill coefficient for AC56 activation
c.hGsGiAC56 = 0.6623;% Hill coefficient for Gs/Gi interaction of AC56

c.KmGsAC47 = 0.031544;% AC47 affinity for Gs
c.KmGsAC56 = 0.0852;% AC56 affinity for Gs
c.KmGiAC56 = 0.0465;% AC56 affinity for inhibition by Gi
c.KmGsGiAC56 = 0.4824;% Gs-dependence of inactivation by Gi for AC56

c.basalAC47 = 0.03135;% Basal AC47 activity
c.basalAC56 = 0.037696;% Basal AC56 activity

c.afAC47 = 3.3757;% Amplification factor for AC47
c.afAC56 = 41.32;% Amplification factor for AC56

c.vGsGiAC56 = 0.8569;% Maximum reduction in Gi inhibition, by Gs

c.AC47_cyt = (1.0 - c.f_AC47_eca) * (1.0 - c.f_AC56_AC47) * c.AC_tot * c.vr_cyt;% Concentration of AC56 in the cytoplasm
c.AC56_cav = c.f_AC56_cav * c.f_AC56_AC47 * c.AC_tot * c.vr_cav;% Concentration of AC56 in the caveolar subspace

c.fATP = c.ATP / (c.KmATP + c.ATP);
c.AC47_eca = c.f_AC47_eca * (1.0 - c.f_AC56_AC47) * c.AC_tot * c.vr_eca;% Concentration of AC56 in the extracaveolar space
c.AC56_cyt = (1.0 - c.f_AC56_cav) * c.f_AC56_AC47 * c.AC_tot * c.vr_cyt;% Concentration of AC56 in the cytoplasm


%% Substrate Phosphorylation
%%  iNaK 
c.ka_inak = 0.015265;               % Rate of INaK phosphorylation by PKA
c.kp_inak = 0.092455;               % Rate of INaK dephosphorylation by phosphatases
c.Ka_inak = 0.0011001;              % Affinity of INaK for phosphorylation by PKA
c.Kp_inak = 5.7392;                 % Affinity of INaK for dephosphorylation by phosphatases

%% iNa
c.Ka_ina = 0.10988;                 % Affinity of INa for phosphorylation by PKA
c.Kp_ina = 7.8605;                  % Affinity of INa for dephosphorylation by phosphatases
c.ka_ina = 0.01368;                 % Rate of INa phosphorylation by PKA
c.kp_ina = 0.052811;                % Rate of INa dephosphorylation by phosphatases


%% iup
c.ka_plb = 0.11348;                 % Rate of PLB phosphorylation by PKA
c.kp_plb = 0.48302;                 % Rate of PLB dephosphorylation by phosphatases
c.Ka_plb =  9.88539999999999992e-04;% Affinity of PLB for phosphorylation by PKA
c.Kp_plb = 0.80737;                 % Affinity of PLB for dephosphorylation by phosphatases


%% iKur
c.ka_ikur = 0.069537;               % Rate of IKur phosphorylation by PKA
c.kp_ikur = 0.317;                  % Rate of IKur dephosphorylation by phosphatases
c.Ka_ikur = 0.27623;                % Affinity of IKur for phosphorylation by PKA
c.Kp_ikur = 0.002331;               % Affinity of IKur for dephosphorylation by phosphatases


%% iKs
c.ka_iks = 0.16305;                 % Rate of IKs channel phosphorylation by PKA
c.kp_iks = 1.0542;                  % Rate of IKs channel dephosphorylation by phosphatases
c.Ka_iks =  9.97940000000000003e-05;% Affinity of IKs channels for phosphorylation by PKA
c.Kp_iks =  1.11470000000000002e-04;% Affinity of IKs channels for dephosphorylation by phosphatases

c.M = 0.01;                 % Binding affinity between PKA and Yotiao
c.iks_sig_L = 0.0001;       % Binding affinity between PKA and Yotiao
c.iks_sig_K = 0.01;         % Binding affinity between PP1 and Yotiao
c.Yotiao = 0.025;           % Total concentration of Yotiao
c.IKs_tot = 0.025;          % Total concentration of IKs channels
c.iks_sig_PKAf_sum = 1.0 + (c.Yotiao - c.PKA_eca) / c.M;
c.iks_sig_PKAf = c.M / 2.0 * (sqrt(c.iks_sig_PKAf_sum ^ 2.0 + 4.0 * c.PKA_eca / c.M) - c.iks_sig_PKAf_sum);% Concentration of PKA not bound to AKAPs in the extracaveolar compartment
c.iks_sig_PP1f_eca_sum = 1.0 + (c.Yotiao - c.PP1_eca) / c.iks_sig_K;
c.PP1f_eca = c.iks_sig_K / 2.0 * (sqrt(c.iks_sig_PP1f_eca_sum ^ 2.0 + 4.0 * c.PP1_eca / c.iks_sig_K) - c.iks_sig_PP1f_eca_sum);% Concentration of PP1 not bound to AKAPs in the extracaveolar compartment
c.iks_sig_IKsf_sum = 1.0 + (c.Yotiao - c.IKs_tot) / c.iks_sig_L;
c.IKsf = c.iks_sig_L / 2.0 * (sqrt(c.iks_sig_IKsf_sum ^ 2.0 + 4.0 * c.IKs_tot / c.iks_sig_L) - c.iks_sig_IKsf_sum);% Concentration of IKs not bound to AKAPs in the extracaveolar compartment
c.Yotiaof = (c.Yotiao - c.IKs_tot + c.IKsf) / ((1.0 + c.PP1f_eca / c.iks_sig_K) * (1.0 + c.iks_sig_PKAf / c.M));
c.IKs_arn = c.IKsf * c.Yotiaof * c.iks_sig_PKAf / (c.iks_sig_L * c.M);% Concentration of IKs channels that have an AKAP and local PKA but no PP1 available
c.IKs_arp = c.IKs_arn * c.PP1f_eca / c.iks_sig_K;% Concentration of IKs channels that have an AKAP, local PKA and PP1 available


%% Troponin
c.ka_tni = 0.10408;                 % Rate of Troponin phosphorylation by PKA
c.kp_tni = 0.052633;                % Rate of Troponin dephosphorylation by phosphatases
c.Ka_tni =  2.71430000000000008e-05;% Affinity of Troponin for phosphorylation by PKA
c.Kp_tni = 0.26714;                 % Affinity of Troponin for dephosphorylation by phosphatases


%% RyR and iCaL  -  Substrates with A-Kinase Anchoring Protein (AKAP) 
c.RyR_tot = 0.125;          % Total concentration of RyRs
c.RyR_akap = 0.125;         % Total concentration of RyR AKAP
c.ka_ryr = 0.0025548;       % Rate of RyR phosphorylation by PKA
c.kp_ryr = 0.0038257;       % Rate of RyR dephosphorylation by phosphatases
c.Ka_ryr =  6.62979999999999944e-05;% Affinity of RyR for phosphorylation by PKA
c.Kp_ryr = 0.043003;        % Affinity of RyR for dephosphorylation by phosphatases
c.Mr = 0.01;                % Binding affinity between PKA and ICaL AKAP
c.Lr = 0.0001;              % Binding affinity between RyR and AKAP
c.Kr = 0.01;                % Binding affinity between PP1 and RyR AKAP

c.akap_sig_RyRf_sum = 1.0 + (c.RyR_akap - c.RyR_tot) / c.Lr;
c.RyRf = c.Lr / 2.0 * (sqrt(c.akap_sig_RyRf_sum ^ 2.0 + 4.0 * c.RyR_tot / c.Lr) - c.akap_sig_RyRf_sum);% Caveolar concentration of free RyR

c.ICaL_tot = 0.025;         % Total concentration of ICaL channels
c.ICaL_akap = 0.025;        % Total concentration of ICaL AKAP
c.ka_ical =  5.10090000000000044e-04;% Rate of ICaL channel phosphorylation by PKA
c.kp_ical = 0.0006903;      % Rate of ICaL channel dephosphorylation by phosphatases
c.Ka_ical =  1.27019999999999993e-06;% Affinity of ICaL channels for phosphorylation by PKA
c.Kp_ical = 0.0063064;      % Affinity of ICaL channels for dephosphorylation by phosphatases
c.Mi = 0.01;                % Binding affinity between PKA and ICaL AKAP
c.Li = 0.0001;              % Binding affinity between ICaL channel and AKAP
c.Ki = 0.01;                % Binding affinity between PP1 and ICaL AKAP
c.akap_sig_ICaLf_sum = 1.0 + (c.ICaL_akap - c.ICaL_tot) / c.Li;
c.ICaLf = c.Li / 2.0 * (sqrt(c.akap_sig_ICaLf_sum ^ 2.0 + 4.0 * c.ICaL_tot / c.Li) - c.akap_sig_ICaLf_sum);% Caveolar concentration of free ICaL

c.akap_sig_PP1f_cav_b = c.ICaL_akap + c.RyR_akap + c.Ki + c.Kr - c.PP1_cav;
c.akap_sig_PP1f_cav_c = c.ICaL_akap * c.Kr + c.RyR_akap * c.Ki + c.Ki * c.Kr - c.PP1_cav * (c.Ki + c.Kr);
c.akap_sig_PP1f_cav_d = c.PP1_cav * c.Ki * c.Kr;
c.akap_sig_PP1f_cav_rr = -c.akap_sig_PP1f_cav_d / 27.0 * c.akap_sig_PP1f_cav_b ^ 3.0 - c.akap_sig_PP1f_cav_b * c.akap_sig_PP1f_cav_b * c.akap_sig_PP1f_cav_c * c.akap_sig_PP1f_cav_c / 108.0 + c.akap_sig_PP1f_cav_b * c.akap_sig_PP1f_cav_c * c.akap_sig_PP1f_cav_d / 6.0 + c.akap_sig_PP1f_cav_c ^ 3.0 / 27.0 + c.akap_sig_PP1f_cav_d * c.akap_sig_PP1f_cav_d / 4.0;
c.akap_sig_PP1f_cav_yi = ifthenelse((c.akap_sig_PP1f_cav_rr < 0.0) , sqrt(-c.akap_sig_PP1f_cav_rr) , 0.0);
c.akap_sig_PP1f_cav_yr = ifthenelse((c.akap_sig_PP1f_cav_rr > 0.0) , sqrt(c.akap_sig_PP1f_cav_rr) , 0.0) + c.akap_sig_PP1f_cav_d / 2.0 + c.akap_sig_PP1f_cav_b * c.akap_sig_PP1f_cav_c / 6.0 - c.akap_sig_PP1f_cav_b ^ 3.0 / 27.0;
c.akap_sig_PP1f_cav_mag = (c.akap_sig_PP1f_cav_yr * c.akap_sig_PP1f_cav_yr + c.akap_sig_PP1f_cav_yi * c.akap_sig_PP1f_cav_yi) ^ (1.0 / 6.0);
c.akap_sig_PP1f_cav_arg = atan(c.akap_sig_PP1f_cav_yi / c.akap_sig_PP1f_cav_yr) / 3.0;
c.akap_sig_PP1f_cav_x = (c.akap_sig_PP1f_cav_c / 3.0 - c.akap_sig_PP1f_cav_b * c.akap_sig_PP1f_cav_b / 9.0) / (c.akap_sig_PP1f_cav_mag * c.akap_sig_PP1f_cav_mag);

c.PP1f_cav = c.akap_sig_PP1f_cav_mag * cos(c.akap_sig_PP1f_cav_arg) * (1.0 - c.akap_sig_PP1f_cav_x) - c.akap_sig_PP1f_cav_b / 3.0;% Caveolar concentration of free PP1
c.akap_sig_PKAf_d = c.PKA_cav * c.Mi * c.Mr;
c.akap_sig_PKAf_b = c.ICaL_akap + c.RyR_akap + c.Mi + c.Mr - c.PKA_cav;
c.akap_sig_PKAf_c = c.ICaL_akap * c.Mr + c.RyR_akap * c.Mi + c.Mi * c.Mr - c.PKA_cav * (c.Mi + c.Mr);
c.akap_sig_PKAf_rr = -c.akap_sig_PKAf_d / 27.0 * c.akap_sig_PKAf_b ^ 3.0 - c.akap_sig_PKAf_b * c.akap_sig_PKAf_b * c.akap_sig_PKAf_c * c.akap_sig_PKAf_c / 108.0 + c.akap_sig_PKAf_b * c.akap_sig_PKAf_c * c.akap_sig_PKAf_d / 6.0 + c.akap_sig_PKAf_c ^ 3.0 / 27.0 + c.akap_sig_PKAf_d * c.akap_sig_PKAf_d / 4.0;
c.akap_sig_PKAf_yr = ifthenelse((c.akap_sig_PKAf_rr > 0.0) , sqrt(c.akap_sig_PKAf_rr) , 0.0) + c.akap_sig_PKAf_d / 2.0 + c.akap_sig_PKAf_b * c.akap_sig_PKAf_c / 6.0 - c.akap_sig_PKAf_b ^ 3.0 / 27.0;
c.akap_sig_PKAf_yi = ifthenelse((c.akap_sig_PKAf_rr < 0.0) , sqrt(-c.akap_sig_PKAf_rr) , 0.0);
c.akap_sig_PKAf_mag = (c.akap_sig_PKAf_yr * c.akap_sig_PKAf_yr + c.akap_sig_PKAf_yi * c.akap_sig_PKAf_yi) ^ (1.0 / 6.0);
c.akap_sig_PKAf_arg = atan(c.akap_sig_PKAf_yi / c.akap_sig_PKAf_yr) / 3.0;
c.akap_sig_PKAf_x = (c.akap_sig_PKAf_c / 3.0 - c.akap_sig_PKAf_b * c.akap_sig_PKAf_b / 9.0) / (c.akap_sig_PKAf_mag * c.akap_sig_PKAf_mag);

c.akap_sig_PKAf = c.akap_sig_PKAf_mag * cos(c.akap_sig_PKAf_arg) * (1.0 - c.akap_sig_PKAf_x) - c.akap_sig_PKAf_b / 3.0;% Caveolar concentration of free PKA

c.RyR_akapf = (c.RyR_akap - c.RyR_tot + c.RyRf) / ((c.PP1f_cav / c.Kr + 1.0) * (c.akap_sig_PKAf / c.Mr + 1.0));% Caveolar concentration of free RyR AKAP
c.RyR_arn = c.RyRf * c.RyR_akapf * c.akap_sig_PKAf / (c.Lr * c.Mr);% Concentration of RyR that have an AKAP, local PKA but no PP1 available
c.RyR_arp = c.RyR_arn * c.PP1f_cav / c.Kr;% Concentration of RyR that have an AKAP, local PKA and PP1

c.ICaL_akapf = (c.ICaL_akap - c.ICaL_tot + c.ICaLf) / ((c.PP1f_cav / c.Ki + 1.0) * (c.akap_sig_PKAf / c.Mi + 1.0));% Caveolar concentration of free ICaL AKAP
c.ICaL_arn = c.ICaLf * c.ICaL_akapf * c.akap_sig_PKAf / (c.Li * c.Mi);% Concentration of ICaL channels that have an AKAP, local PKA but no PP1 available
c.ICaL_arp = c.ICaL_arn * c.PP1f_cav / c.Ki;% Concentration of ICaL channels that have an AKAP, local PKA and PP1


%% These are for electrophysiology.. check these to make sure that the one for O'Hara is consisten
% irel
c.beta_0 = 0.6667 * 4.75;
c.irel_fhat_ratio = 0.0329 + c.RyR_arn / c.RyR_tot;

% ical
c.ical_f_hat_ratio = 0.0269 + c.ICaL_arn / c.ICaL_tot;

%iks
c.iks_f_hat_ratio = 0.0306 + c.IKs_arn / c.IKs_tot;
end
