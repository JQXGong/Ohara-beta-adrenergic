function ydot = Model_SignalingMyokit(y, c)  % Remove t and pace because they are not used.... 



% Create derivatives vector
ydot = zeros(size(y,1), size(y,2));


%% Start of the Signaling state variables
beta_cav_Gs_aGTP = y(1);   %  1: Gs_aGTP_CAV
beta_eca_Gs_aGTP = y(2);   %  2: Gs_aGTP_ECAV
beta_cyt_Gs_aGTP = y(3);   %  3: Gs_a_GTP_CYT
beta_cav_Gs_bg = y(4);     %  4: Gs_bg_CAV
beta_eca_Gs_bg = y(5);     %  5: Gs_bg_ECAV
beta_cyt_Gs_bg = y(6);     %  6: Gs_bg_CYT
beta_cav_Gs_aGDP = y(7);   %  7: Gs_aGDP_CAV
beta_eca_Gs_aGDP = y(8);   %  8: Gs_aGDP_ECAV
beta_cyt_Gs_aGDP = y(9);   %  9: Gs_aGDP_CYT

cAMP_cav = y(10);   % 10: cAMP_CAVVV
cAMP_eca = y(11);   % 11: cAMP_ECAV
cAMP_cyt = y(12);  % 12: cAMP_CYT

beta_cav_Rb1_pka_tot = y(13);  % 13: R_pkap_tot_CAV
beta_eca_Rb1_pka_tot = y(14);  % 14: R_pkap_tot_ECAV
beta_cyt_Rb1_pka_tot = y(15);  % 15: R_pkap_tot_CYT
beta_cav_Rb1_grk_tot = y(16);  % 16: R_grkp_tot_CAV
beta_eca_Rb1_grk_tot = y(17);  % 17: R_grkp_tot_ECAV
beta_cyt_Rb1_grk_tot = y(18);  % 18: R_grkp_tot_CYT

pka_cav_ARC = y(19);   % 19: RLC_CAV
pka_cav_A2RC = y(20);  % 20: L2RC_CAV
pka_cav_A2R = y(21);   % 21: L2R_CAV
pka_cav_C = y(22);     % 22: C_CAV
pka_cav_PKIC = y(23);  % 23: PKI_CAV
pka_eca_ARC = y(24);   % 24: RLC_ECAV
pka_eca_A2RC = y(25);  % 25: L2RC_ECAV
pka_eca_A2R = y(26);   % 26: L2R_ECAV
pka_eca_C = y(27);     % 27: C_ECAV
pka_eca_PKIC = y(28);  % 28: PKI_ECAV
pka_cyt_ARC = y(29);   % 29: RLC_CYT
pka_cyt_A2RC = y(30);  % 30: L2RC_CYT
pka_cyt_A2R = y(31);   % 31: L2R_CYT
pka_cyt_C = y(32);     % 32: C_CYT
pka_cyt_PKIC = y(33);  % 33: PKI_CYT

PDE3_P_cav = y(34);    %34   34: PDE3_P_CAV
PDE3_P_cyt = y(35);    %35   35: PDE3_P_CYT
PDE4_P_cav = y(36);    %36   36: PDE4_P_CAV
PDE4_P_eca = y(37);    %37   37: PDE4_P_ECAV
PDE4_P_cyt = y(38);    %38   38: PDE4_P_CYT

inhib1_p = y(39);  %39       39: Inhib1_P_CYT

ICaLp = y(40);     % 40: fLCC_P
IKsp = y(41);      % 41: fIKS_P
iup_f_plb = y(42); % 42: fPLB_P
f_tni = y(43);     % 43: fTnI_P
ina_f_ina = y(44); % 44: fINa_P
f_inak = y(45);    % 45: fINaK_P
RyRp = y(46);      % 46: fRyR_P
f_ikur = y(47);    % 47: fIKur_P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% update here to make the original PKA fractions bound between [0,1]
%%%% use 0.0001 and 0.9999
%%%% for AKAP related 3 channels implement the constraint in EffectiveFraction function
%%% ICaLp IKsp RyRp

%%% iup_f_plb
iup_f_plb = ifthenelse((iup_f_plb < 0.0) ,0.0001 , iup_f_plb);
iup_f_plb = ifthenelse((iup_f_plb > 1.0) ,0.9999 , iup_f_plb);
%%% f_tni
f_tni = ifthenelse((f_tni < 0.0) ,0.0001 , f_tni);
f_tni = ifthenelse((f_tni > 1.0) ,0.9999 , f_tni);
%%% ina_f_ina
ina_f_ina = ifthenelse((ina_f_ina < 0.0) ,0.0001 , ina_f_ina);
ina_f_ina = ifthenelse((ina_f_ina > 1.0) ,0.9999 , ina_f_ina);
%%% f_inak
f_inak = ifthenelse((f_inak < 0.0) ,0.0001 , f_inak);
f_inak = ifthenelse((f_inak > 1.0) ,0.9999 , f_inak);
%%% f_ikur
f_ikur = ifthenelse((f_ikur < 0.0) ,0.0001 , f_ikur);
f_ikur = ifthenelse((f_ikur > 1.0) ,0.9999 , f_ikur);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta_cav_Rb2_pka_tot = y(48);  % 48: Rb2_pkap_tot_CAV
beta_cav_Rb2_grk_tot = y(49);  % 49: Rb2_grkp_tot_CAV
beta_cav_Gi_aGTP = y(50);      % 50: Gi_aGTP_CAV
beta_cav_Gi_bg = y(51);        % 51: Gi_bg_CAV
beta_cav_Gi_aGDP = y(52);      % 52: Gi_aGDP_CAV
beta_eca_Rb2_pka_tot = y(53);  % 53: Rb2_pkap_tot_ECAV
beta_eca_Rb2_grk_tot = y(54);  % 54: Rb2_grkp_tot_ECAV
beta_eca_Gi_aGTP = y(55);      % 55: Gi_aGTP_ECAV
beta_eca_Gi_bg = y(56);        % 56: Gi_bg_ECAV
beta_eca_Gi_aGDP = y(57);      % 57: Gi_aGDP_ECAV

%% 
%% CAVEOLAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total concentration of non-phosphorylated B1AR in the caveolar subspace
beta_cav_Rb1_np_tot = c.beta_cav_R_b1_tot - beta_cav_Rb1_pka_tot - beta_cav_Rb1_grk_tot;
% Total concentration of non-phosphorylated B2AR in the caveolar subspace
beta_cav_Rb2_np_tot = c.beta_cav_R_b2_tot - beta_cav_Rb2_pka_tot - beta_cav_Rb2_grk_tot;
% Concentration of Gi holoenzyme in the caveolar subspace
beta_cav_Gi_abg = c.f_Gi_cav * c.Gi_tot * c.vr_cav - beta_cav_Gi_aGTP - beta_cav_Gi_aGDP;
% Concentration of Gs holoenzyme in the caveolar subspace
beta_cav_Gs_abg = c.f_Gs_cav * c.Gs_tot * c.vr_cav - beta_cav_Gs_aGTP - beta_cav_Gs_aGDP;

beta_cav_Gs_f_d = beta_cav_Gs_abg * c.beta_cav_Gs_f_c33 / c.beta_cav_Gs_f_a;
beta_cav_Gs_f_b = (c.beta_cav_Gs_f_c11 + c.beta_cav_Gs_f_c22) / c.beta_cav_Gs_f_a + beta_cav_Rb1_np_tot + beta_cav_Rb2_np_tot - beta_cav_Gs_abg;
beta_cav_Gs_f_c = (c.beta_cav_Gs_f_c22 * (beta_cav_Rb1_np_tot - beta_cav_Gs_abg) + c.beta_cav_Gs_f_c11 * (beta_cav_Rb2_np_tot - beta_cav_Gs_abg) + c.beta_cav_Gs_f_c33) / c.beta_cav_Gs_f_a;
beta_cav_Gs_f_rr = -beta_cav_Gs_f_d / 27.0 * beta_cav_Gs_f_b ^ 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_c / 108.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c * beta_cav_Gs_f_d / 6.0 + beta_cav_Gs_f_c ^ 3.0 / 27.0 + beta_cav_Gs_f_d * beta_cav_Gs_f_d / 4.0;
beta_cav_Gs_f_yr = ifthenelse((beta_cav_Gs_f_rr > 0.0) , sqrt(beta_cav_Gs_f_rr) , 0.0) + beta_cav_Gs_f_d / 2.0 + beta_cav_Gs_f_b * beta_cav_Gs_f_c / 6.0 - beta_cav_Gs_f_b ^ 3.0 / 27.0;
beta_cav_Gs_f_yi = ifthenelse((beta_cav_Gs_f_rr < 0.0) , sqrt(-beta_cav_Gs_f_rr) , 0.0);
beta_cav_Gs_f_mag = (beta_cav_Gs_f_yr * beta_cav_Gs_f_yr + beta_cav_Gs_f_yi * beta_cav_Gs_f_yi) ^ (1.0 / 6.0);
beta_cav_Gs_f_arg = atan(beta_cav_Gs_f_yi / beta_cav_Gs_f_yr) / 3.0;
beta_cav_Gs_f_x = (beta_cav_Gs_f_c / 3.0 - beta_cav_Gs_f_b * beta_cav_Gs_f_b / 9.0) / (beta_cav_Gs_f_mag * beta_cav_Gs_f_mag);
beta_cav_Gs_f_r = beta_cav_Gs_f_mag * cos(beta_cav_Gs_f_arg) * (1.0 - beta_cav_Gs_f_x) - beta_cav_Gs_f_b / 3.0;
beta_cav_Gs_f_i = beta_cav_Gs_f_mag * sin(beta_cav_Gs_f_arg) * (1.0 + beta_cav_Gs_f_x);

% Concentration of free Gs in the caveolar subspace
beta_cav_Gs_f = sqrt(beta_cav_Gs_f_r * beta_cav_Gs_f_r + beta_cav_Gs_f_i * beta_cav_Gs_f_i);
% Concentration of free non-phosphorylated beta1AR in caveolar subspace
beta_cav_Rb1_f = beta_cav_Rb1_np_tot / (1.0 + c.iso_L / c.k_b1_l + beta_cav_Gs_f * (c.k_b1_h + c.iso_L) / (c.k_b1_c * c.k_b1_h));
% Concentration of non-phosphorylated Ligand / Receptor1 complexes in caveolar subspace
beta_cav_LRb1 = c.iso_L * beta_cav_Rb1_f / c.k_b1_l;
% Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in caveolar subspace
beta_cav_LRb1Gs = c.iso_L * beta_cav_Rb1_f * beta_cav_Gs_f / (c.k_b1_c * c.k_b1_h);
% Concentration of free non-phosphorylated beta2AR in caveolar subspace
beta_cav_Rb2_f = beta_cav_Rb2_np_tot / (1.0 + c.iso_L / c.k_b2_l + beta_cav_Gs_f * (c.k_b2_h + c.iso_L) / (c.k_b2_c * c.k_b2_h));
% Concentration of non-phosphorylated Ligand / Receptor2 complexes in caveolar subspace
beta_cav_LRb2 = c.iso_L * beta_cav_Rb2_f / c.k_b2_l;

% Concentration of non-phosphorylated Receptor2 / G-protein complexes in caveolar subspace
beta_cav_Rb2Gs = beta_cav_Rb2_f * beta_cav_Gs_f / c.k_b2_c;
% Concentration of non-phosphorylated Receptor1 / G-protein complexes in caveolar subspace
beta_cav_Rb1Gs = beta_cav_Rb1_f * beta_cav_Gs_f / c.k_b1_c;
% Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in caveolar subspace
beta_cav_LRb2Gs = c.iso_L * beta_cav_Rb2_f * beta_cav_Gs_f / (c.k_b2_c * c.k_b2_h);

% Concentration of total PKA-phosphorylated beta1 receptors
ydot(13) = 0.001 * (c.k_pka_p * pka_cav_C * beta_cav_Rb1_np_tot - c.k_pka_dp * beta_cav_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated beta1 receptors
ydot(16) = 0.001 * (c.k_grk_p * c.beta_cav_GRK * (beta_cav_LRb1 + beta_cav_LRb1Gs) - c.k_grk_dp * beta_cav_Rb1_grk_tot);
% Concentration of total PKA-phosphorylated beta2 receptors
ydot(48) = 0.001 * (c.k_pka_p * pka_cav_C * beta_cav_Rb2_np_tot - c.k_pka_dp * beta_cav_Rb2_pka_tot);
% Concentration of total GRK-phosphorylated beta2 receptors
ydot(49) = 0.001 * (c.k_grk_p * c.beta_cav_GRK * (beta_cav_LRb2 + beta_cav_LRb2Gs) - c.k_grk_dp * beta_cav_Rb2_grk_tot);

%% EXTRACAVEOLAR %%%%%%%%%%%%%%%%%%
%
% beta_eca
%
% Concentration of Gs holoenzyme in the extracaveolar space
beta_eca_Gs_abg = c.f_Gs_eca * c.Gs_tot * c.vr_eca - beta_eca_Gs_aGTP - beta_eca_Gs_aGDP;
% Total concentration of non-phosphorylated B2AR in the extracaveolar space
beta_eca_Rb2_np_tot = c.beta_eca_R_b2_tot - beta_eca_Rb2_pka_tot - beta_eca_Rb2_grk_tot;
% Total concentration of non-phosphorylated B1AR in the extracaveolar space
beta_eca_Rb1_np_tot = c.beta_eca_R_b1_tot - beta_eca_Rb1_pka_tot - beta_eca_Rb1_grk_tot;
beta_eca_Gs_f_d = beta_eca_Gs_abg * c.beta_eca_Gs_f_c33 / c.beta_eca_Gs_f_a;
beta_eca_Gs_f_c = (c.beta_eca_Gs_f_c22 * (beta_eca_Rb1_np_tot - beta_eca_Gs_abg) + c.beta_eca_Gs_f_c11 * (beta_eca_Rb2_np_tot - beta_eca_Gs_abg) + c.beta_eca_Gs_f_c33) / c.beta_eca_Gs_f_a;
beta_eca_Gs_f_b = (c.beta_eca_Gs_f_c11 + c.beta_eca_Gs_f_c22) / c.beta_eca_Gs_f_a + beta_eca_Rb1_np_tot + beta_eca_Rb2_np_tot - beta_eca_Gs_abg;
beta_eca_Gs_f_rr = -beta_eca_Gs_f_d / 27.0 * beta_eca_Gs_f_b ^ 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_c / 108.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c * beta_eca_Gs_f_d / 6.0 + beta_eca_Gs_f_c ^ 3.0 / 27.0 + beta_eca_Gs_f_d * beta_eca_Gs_f_d / 4.0;
beta_eca_Gs_f_yi = ifthenelse((beta_eca_Gs_f_rr < 0.0) , sqrt(-beta_eca_Gs_f_rr) , 0.0);
beta_eca_Gs_f_yr = ifthenelse((beta_eca_Gs_f_rr > 0.0) , sqrt(beta_eca_Gs_f_rr) , 0.0) + beta_eca_Gs_f_d / 2.0 + beta_eca_Gs_f_b * beta_eca_Gs_f_c / 6.0 - beta_eca_Gs_f_b ^ 3.0 / 27.0;
beta_eca_Gs_f_mag = (beta_eca_Gs_f_yr * beta_eca_Gs_f_yr + beta_eca_Gs_f_yi * beta_eca_Gs_f_yi) ^ (1.0 / 6.0);
beta_eca_Gs_f_arg = atan(beta_eca_Gs_f_yi / beta_eca_Gs_f_yr) / 3.0;
beta_eca_Gs_f_x = (beta_eca_Gs_f_c / 3.0 - beta_eca_Gs_f_b * beta_eca_Gs_f_b / 9.0) / (beta_eca_Gs_f_mag * beta_eca_Gs_f_mag);
beta_eca_Gs_f_i = beta_eca_Gs_f_mag * sin(beta_eca_Gs_f_arg) * (1.0 + beta_eca_Gs_f_x);
beta_eca_Gs_f_r = beta_eca_Gs_f_mag * cos(beta_eca_Gs_f_arg) * (1.0 - beta_eca_Gs_f_x) - beta_eca_Gs_f_b / 3.0;
% Concentration of free Gs in the caveolar subspace
beta_eca_Gs_f = sqrt(beta_eca_Gs_f_r * beta_eca_Gs_f_r + beta_eca_Gs_f_i * beta_eca_Gs_f_i);
% Concentration of free non-phosphorylated beta1AR in the extracaveolar space
beta_eca_Rb1_f = beta_eca_Rb1_np_tot / (1.0 + c.iso_L / c.k_b1_l + beta_eca_Gs_f * (c.k_b1_h + c.iso_L) / (c.k_b1_c * c.k_b1_h));
% Concentration of free non-phosphorylated beta2AR in the extracaveolar space
beta_eca_Rb2_f = beta_eca_Rb2_np_tot / (1.0 + c.iso_L / c.k_b2_l + beta_eca_Gs_f * (c.k_b2_h + c.iso_L) / (c.k_b2_c * c.k_b2_h));
% Concentration of non-phosphorylated Ligand / Receptor1 complexes in the extracaveolar space
beta_eca_LRb1 = c.iso_L * beta_eca_Rb1_f / c.k_b1_l;
% Concentration of non-phosphorylated Ligand / Receptor2 complexes in the extracaveolar space
beta_eca_LRb2 = c.iso_L * beta_eca_Rb2_f / c.k_b2_l;
% Concentration of non-phosphorylated Ligand / Receptor2 / G-protein complexes in the extracaveolar space
beta_eca_LRb2Gs = c.iso_L * beta_eca_Rb2_f * beta_eca_Gs_f / (c.k_b2_c * c.k_b2_h);
% Concentration of non-phosphorylated Ligand / Receptor1 / G-protein complexes in the extracaveolar space
beta_eca_LRb1Gs = c.iso_L * beta_eca_Rb1_f * beta_eca_Gs_f / (c.k_b1_c * c.k_b1_h);
% Concentration of non-phosphorylated Receptor2 / G-protein complexes in the extracaveolar space
beta_eca_Rb2Gs = beta_eca_Rb2_f * beta_eca_Gs_f / c.k_b2_c;
% Concentration of non-phosphorylated Receptor1 / G-protein complexes in the extracaveolar space
beta_eca_Rb1Gs = beta_eca_Rb1_f * beta_eca_Gs_f / c.k_b1_c;
beta_eca_RGs_tot = beta_eca_Rb1Gs + c.beta_eca_k_GsAct_b2 * beta_eca_Rb2Gs;
beta_eca_LRGs_tot = beta_eca_LRb1Gs + c.beta_eca_k_GsAct_b2 * beta_eca_LRb2Gs;
% Concentration of Gi holoenzyme in the extracaveolar space
beta_eca_Gi_abg = c.f_Gi_eca * c.Gi_tot * c.vr_eca - beta_eca_Gi_aGTP - beta_eca_Gi_aGDP;


beta_eca_Rb2_pka_f_c = -beta_eca_Rb2_pka_tot * c.k_b2_a * c.k_b2_f;
beta_eca_Rb2_pka_f_b = beta_eca_Gi_abg * (c.iso_L + c.k_b2_f) - beta_eca_Rb2_pka_tot * (c.k_b2_f + c.iso_L) + c.k_b2_a * c.k_b2_f * (1.0 + c.iso_L / c.k_b2_n);
% Concentration of free PKA-phosphorylated beta2AR in the extracaveolar space
beta_eca_Rb2_pka_f = (-beta_eca_Rb2_pka_f_b + sqrt(beta_eca_Rb2_pka_f_b * beta_eca_Rb2_pka_f_b - 4.0 * c.beta_eca_Rb2_pka_f_a * beta_eca_Rb2_pka_f_c)) / (2.0 * c.beta_eca_Rb2_pka_f_a);

% Concentration of total PKA-phosphorylated beta1 receptors in the extracaveolar space
ydot(14) = 0.001 * (c.k_pka_p * pka_eca_C * beta_eca_Rb1_np_tot - c.k_pka_dp * beta_eca_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated beta1 receptors in the extracaveolar space
ydot(17) = 0.001 * (c.k_grk_p * c.beta_eca_GRK * (beta_eca_LRb1 + beta_eca_LRb1Gs) - c.k_grk_dp * beta_eca_Rb1_grk_tot);
% Concentration of total PKA-phosphorylated beta2 receptors in the extracaveolar space
ydot(53) = 0.001 * (c.k_pka_p * pka_eca_C * beta_eca_Rb2_np_tot - c.k_pka_dp * beta_eca_Rb2_pka_tot);
% Concentration of total GRK-phosphorylated beta2 receptors in the extracaveolar space
ydot(54) = 0.001 * (c.k_grk_p * c.beta_eca_GRK * (beta_eca_LRb2 + beta_eca_LRb2Gs) - c.k_grk_dp * beta_eca_Rb2_grk_tot);


%% CYTOPLASM %%%%%%%%%%%
% Concentration of Gs holoenzyme in the cytoplasm
beta_cyt_Gs_abg = c.f_Gs_cyt * c.Gs_tot * c.vr_cyt - beta_cyt_Gs_aGTP - beta_cyt_Gs_aGDP;
% Total concentration of non-phosphorylated beta-1 AR in the cytoplasm
beta_cyt_Rb1_np_tot = c.beta_cyt_R_b1_tot - beta_cyt_Rb1_pka_tot - beta_cyt_Rb1_grk_tot;

beta_cyt_Rb1_np_f_b = beta_cyt_Gs_abg * (c.k_b1_h + c.iso_L) - beta_cyt_Rb1_np_tot * (c.k_b1_h + c.iso_L) + c.k_b1_c * c.k_b1_h * (1.0 + c.iso_L / c.k_b1_l);
beta_cyt_Rb1_np_f_c = -beta_cyt_Rb1_np_tot * c.k_b1_h * c.k_b1_c;
% Concentration of free non-phosphorylated beta-1 AR in the cytoplasm
Rb1_np_f = (-beta_cyt_Rb1_np_f_b + sqrt(beta_cyt_Rb1_np_f_b * beta_cyt_Rb1_np_f_b - 4.0 * c.beta_cyt_Rb1_np_f_a * beta_cyt_Rb1_np_f_c)) / (2.0 * c.beta_cyt_Rb1_np_f_a);
% Concentration of free (non-complexed) Gi in the cytoplasm
beta_cyt_Gs_f = beta_cyt_Gs_abg / (1.0 + Rb1_np_f / c.k_b1_c * (1.0 + c.iso_L / c.k_b1_h));
% Concentration of non-phosphorylated ligand / beta-1 AR complexes in the cytoplasm
LRb1_np = c.iso_L * Rb1_np_f / c.k_b1_l;
% Concentration of non-phosphorylated ligand / beta-1 AR / Gs complexes in the cytoplasm
LRb1Gs_np = c.iso_L * Rb1_np_f * beta_cyt_Gs_f / (c.k_b1_c * c.k_b1_h);
% Concentration of non-phosphorylated beta-1 AR / Gs complexes in the cytoplasm
Rb1Gs_np = beta_cyt_Gs_f * Rb1_np_f / c.k_b1_c;

% Concentration of total PKA-phosphorylated receptors in the cytoplasm
ydot(15) = 0.001 * (c.k_pka_p * pka_cyt_C * beta_cyt_Rb1_np_tot - c.k_pka_dp * beta_cyt_Rb1_pka_tot);
% Concentration of total GRK-phosphorylated receptors in the cytoplasm
ydot(18) = 0.001 * (c.k_grk_p * c.beta_cyt_GRK * (LRb1_np + LRb1Gs_np) - c.k_grk_dp * beta_cyt_Rb1_grk_tot);


%% function Mod_GprotAct() in the other version %%%

beta_cav_RGs_tot = beta_cav_Rb1Gs + c.beta_cav_k_GsAct_b2 * beta_cav_Rb2Gs;
beta_cav_LRGs_tot = beta_cav_LRb1Gs + c.beta_cav_k_GsAct_b2 * beta_cav_LRb2Gs;
beta_cav_Rb2_pka_f_c = -beta_cav_Rb2_pka_tot * c.k_b2_a * c.k_b2_f;
beta_cav_Rb2_pka_f_b = beta_cav_Gi_abg * (c.iso_L + c.k_b2_f) - beta_cav_Rb2_pka_tot * (c.k_b2_f + c.iso_L) + c.k_b2_a * c.k_b2_f * (1.0 + c.iso_L / c.k_b2_n);
% Concentration of free PKA-phosphorylated beta2AR in the caveolar subspace
beta_cav_Rb2_pka_f = (-beta_cav_Rb2_pka_f_b + sqrt(beta_cav_Rb2_pka_f_b * beta_cav_Rb2_pka_f_b - 4.0 * c.beta_cav_Rb2_pka_f_a * beta_cav_Rb2_pka_f_c)) / (2.0 * c.beta_cav_Rb2_pka_f_a);
% Concentration of free (non-complexed) Gi in the caveolar subspace
beta_cav_Gi_f = beta_cav_Gi_abg / (1.0 + beta_cav_Rb2_pka_f / c.k_b2_a * (1.0 + c.iso_L / c.k_b2_f));
% Concentration of PKA phosphorylated b2AR/Gi complexes in caveolar subspace
beta_cav_Rb2Gi = beta_cav_Rb2_pka_f * beta_cav_Gi_f / c.k_b2_a;
% Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in caveolar subspace
beta_cav_LRb2Gi = beta_cav_Rb2Gi * c.iso_L / c.k_b2_f;
% Concentration of free (non-complexed) Gi in the extracaveolar space
beta_eca_Gi_f = beta_eca_Gi_abg / (1.0 + beta_eca_Rb2_pka_f / c.k_b2_a * (1.0 + c.iso_L / c.k_b2_f));
% Concentration of PKA phosphorylated b2AR/Gi complexes in extracaveolar space
beta_eca_Rb2Gi = beta_eca_Rb2_pka_f * beta_eca_Gi_f / c.k_b2_a;
% Concentration of PKA phosphorylated Ligand/b2AR/Gi complexes in extracaveolar space
beta_eca_LRb2Gi = c.iso_L / c.k_b2_f * beta_eca_Rb2Gi;

% Concentration of active Gs alpha subunit in caveolar subspace (tri-phosphate)
ydot(1) = 0.001 * (c.k_act2_Gs * beta_cav_RGs_tot + c.k_act1_Gs * beta_cav_LRGs_tot - c.k_hydr_Gs * beta_cav_Gs_aGTP);
% Concentration of active Gi alpha subunit in caveolar subspace
ydot(50) = 0.001 * (c.k_act2_Gi * beta_cav_Rb2Gi + c.k_act1_Gi * beta_cav_LRb2Gi - c.k_hydr_Gi * beta_cav_Gi_aGTP);
% Concentration of active Gs alpha subunit in the extracaveolar space (tri-phosphate)
ydot(2) = 0.001 * (c.k_act2_Gs * beta_eca_RGs_tot + c.k_act1_Gs * beta_eca_LRGs_tot - c.k_hydr_Gs * beta_eca_Gs_aGTP);
% Concentration of active Gi alpha subunit in the extracaveolar space
ydot(55) = 0.001 * (c.k_act2_Gi * beta_eca_Rb2Gi + c.k_act1_Gi * beta_eca_LRb2Gi - c.k_hydr_Gi * beta_eca_Gi_aGTP);
% Concentration of active Gs alpha subunit in cytoplasm
ydot(3) = 0.001 * (c.k_act2_Gs * Rb1Gs_np + c.k_act1_Gs * LRb1Gs_np - c.k_hydr_Gs * beta_cyt_Gs_aGTP);

% Concentration of active Gs beta-gamma subunit in caveolar subspace
ydot(4) = 0.001 * (c.k_act2_Gs * beta_cav_RGs_tot + c.k_act1_Gs * beta_cav_LRGs_tot - c.k_reas_Gs * beta_cav_Gs_bg * beta_cav_Gs_aGDP);
% Concentration of active Gi beta-gamma subunit in caveolar subspace
ydot(51) = 0.001 * (c.k_act2_Gi * beta_cav_Rb2Gi + c.k_act1_Gi * beta_cav_LRb2Gi - c.k_reas_Gi * beta_cav_Gi_bg * beta_cav_Gi_aGDP);
% Concentration of active Gs beta-gamma subunit in the extracaveolar space
ydot(5) = 0.001 * (c.k_act2_Gs * beta_eca_RGs_tot + c.k_act1_Gs * beta_eca_LRGs_tot - c.k_reas_Gs * beta_eca_Gs_bg * beta_eca_Gs_aGDP);
% Concentration of active Gi beta-gamma subunit in the extracaveolar space
ydot(56) = 0.001 * (c.k_act2_Gi * beta_eca_Rb2Gi + c.k_act1_Gi * beta_eca_LRb2Gi - c.k_reas_Gi * beta_eca_Gi_bg * beta_eca_Gi_aGDP);
% Concentration of active Gs beta-gamma subunit in cytoplasm
ydot(6) = 0.001 * (c.k_act2_Gs * Rb1Gs_np + c.k_act1_Gs * LRb1Gs_np - c.k_reas_Gs * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP);

% Concentration of inactive Gs alpha subunit in caveolar subspace (di-phosphate)
ydot(7) = 0.001 * (c.k_hydr_Gs * beta_cav_Gs_aGTP - c.k_reas_Gs * beta_cav_Gs_bg * beta_cav_Gs_aGDP);
% Concentration of inactive Gi alpha subunit in caveolar subspace
ydot(52) = 0.001 * (c.k_hydr_Gi * beta_cav_Gi_aGTP - c.k_reas_Gi * beta_cav_Gi_bg * beta_cav_Gi_aGDP);
% Concentration of inactive Gs alpha subunit in the extracaveolar space (di-phosphate)
ydot(8) = 0.001 * (c.k_hydr_Gs * beta_eca_Gs_aGTP - c.k_reas_Gs * beta_eca_Gs_bg * beta_eca_Gs_aGDP);
% Concentration of inactive Gi alpha subunit in the extracaveolar space
ydot(57) = 0.001 * (c.k_hydr_Gi * beta_eca_Gi_aGTP - c.k_reas_Gi * beta_eca_Gi_bg * beta_eca_Gi_aGDP);
% Concentration of inactive Gs alpha subunit in cytoplasm
ydot(9) = 0.001 * (c.k_hydr_Gs * beta_cyt_Gs_aGTP - c.k_reas_Gs * beta_cyt_Gs_bg * beta_cyt_Gs_aGDP);

%% function Mod_AC() 

% only calculating constants, but does not compute derivative of state
% variables

%% function Mod_PKA() in the other version %%%%%%%%%%
% Concentration of free PKA RC subunits in the caveolar compartment
pka_cav_RCf = c.PKA_cav - pka_cav_ARC - pka_cav_A2RC - pka_cav_A2R;

% Caveolar concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(19) = 0.001 * (c.pka_cav_f1 * pka_cav_RCf * cAMP_cav - c.pka_cav_b1 * pka_cav_ARC - c.pka_cav_f2 * pka_cav_ARC * cAMP_cav + c.pka_cav_b2 * pka_cav_A2RC);
% Caveolar concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(20) = 0.001 * (c.pka_cav_f2 * pka_cav_ARC * cAMP_cav - (c.pka_cav_b2 + c.pka_cav_f3) * pka_cav_A2RC + c.pka_cav_b3 * pka_cav_A2R * pka_cav_C);
% Caveolar concentration of PKA R subunit with 2 cAMP molecules bound
ydot(21) = 0.001 * (c.pka_cav_f3 * pka_cav_A2RC - c.pka_cav_b3 * pka_cav_A2R * pka_cav_C);
% Caveolar concentration of free PKA catalytic subunit
ydot(22) = 0.001 * (c.pka_cav_f3 * pka_cav_A2RC - c.pka_cav_b3 * pka_cav_A2R * pka_cav_C + c.b_pki * pka_cav_PKIC - c.f_pki * (c.PKI_cav - pka_cav_PKIC) * pka_cav_C);
% Caveolar concentration of free PKI inactivated PKA C subunit
ydot(23) = 0.001 * (c.f_pki * (c.PKI_cav - pka_cav_PKIC) * pka_cav_C - c.b_pki * pka_cav_PKIC);

% Concentration of free PKA RC subunits in the Extracaveolar compartment
pka_eca_RCf = c.PKA_eca - pka_eca_ARC - pka_eca_A2RC - pka_eca_A2R;
% Extracaveolar rate of change in free cAMP through binding by PKA
pka_eca_dcAMP = -c.pka_eca_f1 * pka_eca_RCf * cAMP_eca + c.pka_eca_b1 * pka_eca_ARC - c.pka_eca_f2 * pka_eca_ARC * cAMP_eca + c.pka_eca_b2 * pka_eca_A2RC;

% Extracaveolar concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(24) = 0.001 * (c.pka_eca_f1 * pka_eca_RCf * cAMP_eca - c.pka_eca_b1 * pka_eca_ARC - c.pka_eca_f2 * pka_eca_ARC * cAMP_eca + c.pka_eca_b2 * pka_eca_A2RC);
% Extracaveolar concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(25) = 0.001 * (c.pka_eca_f2 * pka_eca_ARC * cAMP_eca - (c.pka_eca_b2 + c.pka_eca_f3) * pka_eca_A2RC + c.pka_eca_b3 * pka_eca_A2R * pka_eca_C);
% Extracaveolar concentration of PKA R subunit with 2 cAMP molecules bound
ydot(26) = 0.001 * (c.pka_eca_f3 * pka_eca_A2RC - c.pka_eca_b3 * pka_eca_A2R * pka_eca_C);
% Extracaveolar concentration of free PKA catalytic subunit
ydot(27) = 0.001 * (c.pka_eca_f3 * pka_eca_A2RC - c.pka_eca_b3 * pka_eca_A2R * pka_eca_C + c.b_pki * pka_eca_PKIC - c.f_pki * (c.PKI_eca - pka_eca_PKIC) * pka_eca_C);
% Extracaveolar concentration of free PKI inactivated PKA C subunit
ydot(28) = 0.001 * (c.f_pki * (c.PKI_eca - pka_eca_PKIC) * pka_eca_C - c.b_pki * pka_eca_PKIC);

% Concentration of free PKA RC subunits in the Cytosolic compartment
pka_cyt_RCf = c.PKA_cyt - pka_cyt_ARC - pka_cyt_A2RC - pka_cyt_A2R;

% Cytosolic concentration of PKA RC dimer with 1 cAMP molecule bound
ydot(29) = 0.001 * (c.pka_cyt_f1 * pka_cyt_RCf * cAMP_cyt - c.pka_cyt_b1 * pka_cyt_ARC - c.pka_cyt_f2 * pka_cyt_ARC * cAMP_cyt + c.pka_cyt_b2 * pka_cyt_A2RC);
% Cytosolic concentration of PKA RC dimer with 2 cAMP molecules bound
ydot(30) = 0.001 * (c.pka_cyt_f2 * pka_cyt_ARC * cAMP_cyt - (c.pka_cyt_b2 + c.pka_cyt_f3) * pka_cyt_A2RC + c.pka_cyt_b3 * pka_cyt_A2R * pka_cyt_C);
% Cytosolic concentration of PKA R subunit with 2 cAMP molecules bound
ydot(31) = 0.001 * (c.pka_cyt_f3 * pka_cyt_A2RC - c.pka_cyt_b3 * pka_cyt_A2R * pka_cyt_C);
% Cytosolic concentration of free PKA catalytic subunit
ydot(32) = 0.001 * (c.pka_cyt_f3 * pka_cyt_A2RC - c.pka_cyt_b3 * pka_cyt_A2R * pka_cyt_C + c.b_pki * pka_cyt_PKIC - c.f_pki * (c.PKI_cyt - pka_cyt_PKIC) * pka_cyt_C);
% Cytosolic concentration of free PKI inactivated PKA C subunit
ydot(33) = 0.001 * (c.f_pki * (c.PKI_cyt - pka_cyt_PKIC) * pka_cyt_C - c.b_pki * pka_cyt_PKIC);

%% function Mod_cAMP() in the other version


% Caveolar rate of change in free cAMP through binding by PKA
pka_cav_dcAMP = -c.pka_cav_f1 * pka_cav_RCf * cAMP_cav + c.pka_cav_b1 * pka_cav_ARC - c.pka_cav_f2 * pka_cav_ARC * cAMP_cav + c.pka_cav_b2 * pka_cav_A2RC;

% Cytosolic rate of change in free cAMP through binding by PKA
pka_cyt_dcAMP = -c.pka_cyt_f1 * pka_cyt_RCf * cAMP_cyt + c.pka_cyt_b1 * pka_cyt_ARC - c.pka_cyt_f2 * pka_cyt_ARC * cAMP_cyt + c.pka_cyt_b2 * pka_cyt_A2RC;

%PDE
% Rate of cAMP degradation by PDE2 in cytosolic subspace
dcAMP_PDE2_cyt = c.PDE2_cyt * c.kPDE2 / (1.0 + c.KmPDE2 / cAMP_cyt);
% Rate of cAMP degradation by PDE2 in extracaveolar subspace
dcAMP_PDE2_eca = c.PDE2_eca * c.kPDE2 / (1.0 + c.KmPDE2 / cAMP_eca);
% Rate of cAMP degradation by PDE2 in caveolar subspace
dcAMP_PDE2_cav = c.PDE2_cav * c.kPDE2 / (1.0 + c.KmPDE2 / cAMP_cav);
% Rate of cAMP degradation by PDE3 in caveolar subspace
dcAMP_PDE3_cav = (c.PDE3_cav + (c.delta_k_pde34 - 1.0) * PDE3_P_cav) * c.kPDE3 / (1.0 + c.KmPDE3 / cAMP_cav);
% Rate of cAMP degradation by PDE4 in cytosolic subspace
dcAMP_PDE4_cyt = (c.PDE4_cyt + (c.delta_k_pde34 - 1.0) * PDE4_P_cyt) * c.kPDE4 / (1.0 + c.KmPDE4 / cAMP_cyt);
% Rate of cAMP degradation by PDE4 in extracaveolar subspace
dcAMP_PDE4_eca = (c.PDE4_eca + (c.delta_k_pde34 - 1.0) * PDE4_P_eca) * c.kPDE4 / (1.0 + c.KmPDE4 / cAMP_eca);
% Rate of cAMP degradation by PDE4 in caveolar subspace
dcAMP_PDE4_cav = (c.PDE4_cav + (c.delta_k_pde34 - 1.0) * PDE4_P_cav) * c.kPDE4 / (1.0 + c.KmPDE4 / cAMP_cav);
% Rate of cAMP degradation by PDE3 in cytosolic subspace
dcAMP_PDE3_cyt = (c.PDE3_cyt + (c.delta_k_pde34 - 1.0) * PDE3_P_cyt) * c.kPDE3 / (1.0 + c.KmPDE3 / cAMP_cyt);

camp_cAMP_cyt_pde = dcAMP_PDE2_cyt + dcAMP_PDE3_cyt + dcAMP_PDE4_cyt;
camp_cAMP_cyt_j1 = c.j_cav_cyt * (cAMP_cav - cAMP_cyt) / c.v_cyt;
camp_cAMP_cyt_j2 = c.j_eca_cyt * (cAMP_eca - cAMP_cyt) / c.v_cyt;

camp_cAMP_eca_pde = dcAMP_PDE2_eca + dcAMP_PDE4_eca;
camp_cAMP_eca_j2 = c.j_eca_cyt * (cAMP_eca - cAMP_cyt) / c.v_eca;
camp_cAMP_eca_j1 = c.j_cav_eca * (cAMP_cav - cAMP_eca) / c.v_eca;

camp_cAMP_cav_pde = dcAMP_PDE2_cav + dcAMP_PDE3_cav + dcAMP_PDE4_cav;
camp_cAMP_cav_j2 = c.j_cav_cyt * (cAMP_cav - cAMP_cyt) / c.v_cav;
camp_cAMP_cav_j1 = c.j_cav_eca * (cAMP_cav - cAMP_eca) / c.v_cav;
% ac
%
ac_kAC47_cyt_gsa = beta_cyt_Gs_aGTP ^ c.hGsAC47;
kAC47_cyt = c.afAC47 * (c.basalAC47 + ac_kAC47_cyt_gsa / (c.KmGsAC47 + ac_kAC47_cyt_gsa));
ac_kAC56_cav_gsa = beta_cav_Gs_aGTP ^ c.hGsAC56;
gsi = beta_cav_Gs_aGTP ^ c.hGsGiAC56;
kAC56_cav = c.afAC56 * (c.basalAC56 + ac_kAC56_cav_gsa / (c.KmGsAC56 + ac_kAC56_cav_gsa)) * (1.0 - (1.0 - c.vGsGiAC56 * gsi / (c.KmGsGiAC56 + gsi)) * beta_cav_Gi_bg / (c.KmGiAC56 + beta_cav_Gi_bg));
ac_kAC47_eca_gsa = beta_eca_Gs_aGTP ^ c.hGsAC47;
kAC47_eca = c.afAC47 * (c.basalAC47 + ac_kAC47_eca_gsa / (c.KmGsAC47 + ac_kAC47_eca_gsa));
ac_kAC56_cyt_gsa = beta_cyt_Gs_aGTP ^ c.hGsAC56;
kAC56_cyt = c.afAC56 * (c.basalAC56 + ac_kAC56_cyt_gsa / (c.KmGsAC56 + ac_kAC56_cyt_gsa));

% Rate of cAMP production by AC type 4/7 in cytoplasm
dcAMP_AC47_cyt = kAC47_cyt * c.AC47_cyt * c.fATP;
% Rate of cAMP production by AC type 5/6 in cytoplasm
dcAMP_AC56_cyt = kAC56_cyt * c.AC56_cyt * c.fATP;
% Rate of cAMP production by AC type 5/6 in caveolar subspace
dcAMP_AC56_cav = kAC56_cav * c.AC56_cav * c.fATP;
% Rate of cAMP production by AC type 4/7 in extracaveolar subspace
dcAMP_AC47_eca = kAC47_eca * c.AC47_eca * c.fATP;


% Caveolar concentration of cAMP
ydot(10) = 0.001 * (pka_cav_dcAMP + dcAMP_AC56_cav - camp_cAMP_cav_pde - camp_cAMP_cav_j1 - camp_cAMP_cav_j2);
% Extracaveolar concentration of cAMP
ydot(11) = 0.001 * (pka_eca_dcAMP + dcAMP_AC47_eca - camp_cAMP_eca_pde + camp_cAMP_eca_j1 - camp_cAMP_eca_j2);
% Cytosolic concentration of cAMP
ydot(12) = 0.001 * (pka_cyt_dcAMP + dcAMP_AC47_cyt + dcAMP_AC56_cyt - camp_cAMP_cyt_pde + camp_cAMP_cyt_j1 + camp_cAMP_cyt_j2);


%% Mod_PDE_Phosphorylation() function

% Concentration of phosphorylated PDE3 in the caveolar subspace
ydot(34) = 0.001 * (c.kfPDEp * pka_cav_C * (c.PDE3_cav - PDE3_P_cav) - c.kbPDEp * PDE3_P_cav);
% Concentration of phosphorylated PDE3 in the cytosolic subspace
ydot(35) = 0.001 * (c.kfPDEp * pka_cyt_C * (c.PDE3_cyt - PDE3_P_cyt) - c.kbPDEp * PDE3_P_cyt);
% Concentration of phosphorylated PDE4 in the caveolar subspace
ydot(36) = 0.001 * (c.kfPDEp * pka_cav_C * (c.PDE4_cav - PDE4_P_cav) - c.kbPDEp * PDE4_P_cav);
% Concentration of phosphorylated PDE4 in the extracaveolar subspace
ydot(37) = 0.001 * (c.kfPDEp * pka_eca_C * (c.PDE4_eca - PDE4_P_eca) - c.kbPDEp * PDE4_P_eca);
% Concentration of phosphorylated PDE4 in the cytosolic subspace
ydot(38) = 0.001 * (c.kfPDEp * pka_cyt_C * (c.PDE4_cyt - PDE4_P_cyt) - c.kbPDEp * PDE4_P_cyt);

%% Mod_PP1_Inhibition()

pp1_PP1f_cyt_sum = c.pp1_K - c.PP1_cyt + inhib1_p;
% Concentration of uninhibited PP1 in the cytosolic compartment
PP1f_cyt = 0.5 * (sqrt(pp1_PP1f_cyt_sum ^ 2.0 + 4.0 * c.pp1_K * c.PP1_cyt) - pp1_PP1f_cyt_sum);
di = c.inhib1_tot - inhib1_p;
% Concentration of phosphorylated PP1 inhibitor 1 (cytoplasmic)
ydot(39) = 0.001 * (c.kp * pka_cyt_C * di / (c.Kp + di) - c.kdp * c.PP2A * inhib1_p / (c.Kdp + inhib1_p));

%% Mod_Channel_Phosphorylation()      

% Substrates without AKAP
% Fraction of phosphorylated PLB
ydot(42) = 0.001 * (c.ka_plb * pka_cyt_C * (1.0 - iup_f_plb) / (c.Ka_plb + 1.0 - iup_f_plb) - c.kp_plb * PP1f_cyt * iup_f_plb / (c.Kp_plb + iup_f_plb));
% Fraction of phosphorylated Troponin
ydot(43) = 0.001 * (c.ka_tni * pka_cyt_C * (1.0 - f_tni) / (c.Ka_tni + 1.0 - f_tni) - c.kp_tni * c.PP2A * f_tni / (c.Kp_tni + f_tni));
% Fraction of phosphorylated INa channels
ydot(44) = 0.001 * (c.ka_ina * pka_cav_C * (1.0 - ina_f_ina) / (c.Ka_ina + 1.0 - ina_f_ina) - c.kp_ina * c.PP1_cav * ina_f_ina / (c.Kp_ina + ina_f_ina));
% Fraction of phosphorylated INaK
ydot(45) = 0.001 * (c.ka_inak * pka_cav_C * (1.0 - f_inak) / (c.Ka_inak + 1.0 - f_inak) - c.kp_inak * c.PP1_cav * f_inak / (c.Kp_inak + f_inak));
% Fraction of phosphorylated IKur channels
ydot(47) = 0.001 * (c.ka_ikur * pka_eca_C * (1.0 - f_ikur) / (c.Ka_ikur + 1.0 - f_ikur) - c.kp_ikur * c.PP1_eca * f_ikur / (c.Kp_ikur + f_ikur));


%Substrates with AKAP
iks_sig_IKsp_dif = c.IKs_arp - IKsp;
% Concentration of phosphorylated IKs channels
ydot(41) = 0.001 * (c.ka_iks * pka_eca_C * iks_sig_IKsp_dif / (c.Ka_iks + iks_sig_IKsp_dif) - c.kp_iks * c.PP1_eca * IKsp / (c.Kp_iks + IKsp));
akap_sig_RyRp_dif = c.RyR_arp - RyRp;
ydot(46) = 0.001 * (c.ka_ryr * pka_cav_C * akap_sig_RyRp_dif / (c.Ka_ryr + akap_sig_RyRp_dif) - c.kp_ryr * c.PP1_cav * RyRp / (c.Kp_ryr + RyRp));
akap_sig_ICaLp_dif = c.ICaL_arp - ICaLp;
% Concentration of phosphorylated L-type Calcium channels
ydot(40) = 0.001 * (c.ka_ical * pka_cav_C * akap_sig_ICaLp_dif / (c.Ka_ical + akap_sig_ICaLp_dif) - c.kp_ical * c.PP1_cav * ICaLp / (c.Kp_ical + ICaLp));

% if ICaLp > c.ICaL_arp; ydot(40) = -0.001; end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reporting is ydot ever become complex, which one is the first to have a
% non-zero imaginary part
% % size(ydot) 145 x 1
% % % if isreal(ydot) ~= 1
% % %    image_part = imag(ydot) ;
% % %    [index,image_value] = find(image_part~=0)
% % %    
% % %    if isempty(index)~=1
% % %       disp(['first complex derivative is y(',num2str(index(1)),')'])
% % %       disp(['No. complex derivatives ',num2str(length(index))])
% % %       save FirstComplexIdentified.mat ydot index image_value
% % %       
% % %       ydot = [] ;
% % %    end
% % % end












%% 








% %
% % ina_pka
% %
% ina_pka_j_alpha = ifthenelse((V >= -40.0) , 0.0 , (-127140.0 * exp(0.2444 * (V + c.ina_pka_dVIn)) - 6.948e-05 * exp(-0.04391 * (V + c.ina_pka_dVIn))) * (V + c.ina_pka_dVIn + 37.78) / (1.0 + exp(0.311 * (V + c.ina_pka_dVIn + 79.23))));
% ina_pka_j_beta = ifthenelse((V >= -40.0) , 0.3 * exp(-2.535e-07 * (V + c.ina_pka_dVIn)) / (1.0 + exp(-0.1 * (V + c.ina_pka_dVIn + 32.0))) , 0.1212 * exp(-0.01052 * (V + c.ina_pka_dVIn)) / (1.0 + exp(-0.1378 * (V + c.ina_pka_dVIn + 40.14))));
% ydot(76) = ina_pka_j_alpha * (1.0 - ina_pka_j) - ina_pka_j_beta * ina_pka_j;
% ina_pka_h_beta = ifthenelse((V >= -40.0) , 1.0 / 0.13 / (1.0 + exp((V + c.ina_pka_dVIn + 27.4034) / -11.1)) , 3.56 * exp(0.079 * (V + c.ina_pka_dVIn + 7.0)) + 310000.0 * exp(0.35 * (V + c.ina_pka_dVIn + 7.0)));
% ina_pka_h_alpha = ifthenelse((V >= -40.0) , 0.0 , 0.135 * exp((87.0 + V + c.ina_pka_dVIn) / -6.8));
% ydot(74) = ina_pka_h_alpha * (1.0 - ina_pka_h) - ina_pka_h_beta * ina_pka_h;
% ina_pka_m_beta = 0.08 * exp((V + c.dVAc - 13.7299) / -11.0);
% ina_pka_m_alpha = 0.32 * (V + c.dVAc + 58.4729) / (1.0 - exp(-0.1 * (V + c.dVAc + 58.4729)));
% ydot(75) = ina_pka_m_alpha * (1.0 - ina_pka_m) - ina_pka_m_beta * ina_pka_m;


% %
% % iso
% %
% %
% % engine
% %
% time = t;
% engine_pace = pace;
% %
% % ina_np
% %
% ina_np_h_alpha = ifthenelse((V >= -40.0) , 0.0 , 0.135 * exp((87.0 + V) / -6.8));
% ina_np_h_beta = ifthenelse((V >= -40.0) , 1.0 / (0.13 * (1.0 + exp((V + 27.4034) / -11.1))) , 3.56 * exp(0.079 * (V + 7.0)) + 310000.0 * exp(0.35 * (V + 7.0)));
% ydot(14) = ina_np_h_alpha * (1.0 - ina_np_h) - ina_np_h_beta * ina_np_h;
% ina_np_j_alpha = ifthenelse((V >= -40.0) , 0.0 , (-127140.0 * exp(0.2444 * V) - 6.948e-05 * exp(-0.04391 * V)) * (V + 37.78) / (1.0 + exp(0.311 * (V + 79.23))));
% ina_np_j_beta = ifthenelse((V >= -40.0) , 0.3 * exp(-2.535e-07 * V) / (1.0 + exp(-0.1 * (V + 32.0))) , 0.1212 * exp(-0.01052 * V) / (1.0 + exp(-0.1378 * (V + 40.14))));
% ydot(16) = ina_np_j_alpha * (1.0 - ina_np_j) - ina_np_j_beta * ina_np_j;
% ina_np_m_alpha = 0.32 * (V + 58.4729) / (1.0 - exp(-0.1 * (V + 58.4729)));
% ina_np_m_beta = 0.08 * exp((13.7299 - V) / 11.0);
% ydot(15) = ina_np_m_alpha * (1.0 - ina_np_m) - ina_np_m_beta * ina_np_m;

% %
% % ina_camk
% %
% ina_camk_m_beta = 0.08 * exp((13.7299 - V) / 11.0);
% ina_camk_m_alpha = 0.32 * (V + 58.4729) / (1.0 - exp(-0.1 * (V + 58.4729)));
% ydot(77) = ina_camk_m_alpha * (1.0 - ina_camk_m) - ina_camk_m_beta * ina_camk_m;
% ina_camk_j_alpha = ifthenelse((V + c.ina_camk_dVIn >= -40.0) , 0.0 , (-127140.0 * exp(0.2444 * (V + c.ina_camk_dVIn)) + -6.948e-05 * exp(-0.04391 * (V + c.ina_camk_dVIn))) * (V + c.ina_camk_dVIn + 37.78) / (1.0 + exp(0.311 * (V + c.ina_camk_dVIn + 79.23))));
% ina_camk_j_beta = ifthenelse((V + c.ina_camk_dVIn >= -40.0) , 0.3 * exp(-2.535e-07 * (V + c.ina_camk_dVIn)) / (1.0 + exp(-0.1 * (V + c.ina_camk_dVIn + 32.0))) , 0.1212 * exp(-0.01052 * (V + c.ina_camk_dVIn)) / (1.0 + exp(-0.1378 * (V + c.ina_camk_dVIn + 40.14))));
% ydot(79) = ina_camk_j_alpha * (1.0 - ina_camk_j) - ina_camk_j_beta * ina_camk_j;
% ina_camk_h_alpha = ifthenelse((V + c.ina_camk_dVIn >= -40.0) , 0.0 , 0.135 * exp((87.0 + V + c.ina_camk_dVIn) / -6.8));
% ina_camk_h_beta = ifthenelse((V + c.ina_camk_dVIn >= -40.0) , 1.0 / (0.13 * (1.0 + exp((V + c.ina_camk_dVIn + 27.4034) / -11.1))) , 3.56 * exp(0.079 * (V + c.ina_camk_dVIn + 7.0)) + 310000.0 * exp(0.35 * (V + c.ina_camk_dVIn + 7.0)));
% ydot(78) = ina_camk_h_alpha * (1.0 - ina_camk_h) - ina_camk_h_beta * ina_camk_h;
% %
% % extra
% %
% %
% % beta
% %
% %
% % phys
% %
% %
% % inab
% %
% inab_INab_phi = V * c.FRT;
% ePhi = exp(inab_INab_phi);
% % Background Sodium current
% INab = c.P * c.F * inab_INab_phi * (Na * ePhi - c.Nao) / (ePhi - 1.0);
% %
% % stimulus
% %
% % The stimulus current
% i_stim = engine_pace * c.amplitude;
% %
% % nernst
% %
% % Reversal potential for Sodium
% ENa = c.RTF * log(c.Nao / Na);
% % Reversal potential for Chloride
% ECl = -c.RTF * log(c.Clo / Cl);
% % Reversal potential for Potassium
% EK = c.RTF * log(c.Ko / potassium_K);
% % Reversal potential for IKs
% EKs = c.RTF * log((c.Ko + c.PNaK * c.Nao) / (potassium_K + c.PNaK * Na));
% %
% % iks_np
% %
% iks_np_e = 0.0031124 + (0.02833 - 0.0031124) / (1.0 + exp(c.FRT * (V + 0.05166) / 1.5522));
% iks_np_g = 0.38839 / (1.0 + exp(c.FRT * (V + 0.15019) / -0.60693));
% iks_np_o =  4.41980000000000020e-04 * exp(-1.2022 * V * c.FRT);
% iks_np_p =  4.01729999999999991e-04 * exp( 2.08729999999999989e-04 * V * c.FRT);
% iks_np_b = 0.0056992 / (1.0 + exp(c.FRT * (V - 0.04152) / 1.3489));
% iks_np_a = 0.007399 / (1.0 + exp(c.FRT * (V - 0.031196) / -0.80019));
% iks_np_d = 0.090654 * exp(-0.11157 * V * c.FRT);
% ydot(38) = 3.0 * iks_np_a * iks_np_C6 + 2.0 * iks_np_b * iks_np_C8 + 2.0 * iks_np_g * iks_np_C3 + 2.0 * iks_np_d * iks_np_C10 - iks_np_C7 * (2.0 * iks_np_a + iks_np_b + iks_np_g + iks_np_d);
% ydot(36) = iks_np_a * iks_np_C4 + iks_np_d * iks_np_C9 - iks_np_C5 * (4.0 * iks_np_b + 4.0 * iks_np_g);
% ydot(45) = iks_np_a * iks_np_C13 + 2.0 * iks_np_g * iks_np_C12 + 4.0 * iks_np_d * iks_np_C15 - iks_np_C14 * (iks_np_b + iks_np_g + 3.0 * iks_np_d);
% ydot(44) = iks_np_b * iks_np_C14 + iks_np_g * iks_np_C11 - iks_np_C13 * (iks_np_a + 3.0 * iks_np_d);
% ydot(43) = iks_np_a * iks_np_C11 + 3.0 * iks_np_g * iks_np_C9 + 3.0 * iks_np_d * iks_np_C14 - iks_np_C12 * (2.0 * iks_np_b + 2.0 * iks_np_g + 2.0 * iks_np_d);
% ydot(48) = iks_np_p * iks_np_O1 - iks_np_o * iks_np_O2;
% ydot(47) = -(iks_np_e + iks_np_p) * iks_np_O1 + iks_np_o * iks_np_O2 + c.iks_np_t * iks_np_C15;
% ydot(34) = 3.0 * iks_np_a * iks_np_C2 + 3.0 * iks_np_b * iks_np_C4 + iks_np_d * iks_np_C7 - iks_np_C3 * (2.0 * iks_np_a + 2.0 * iks_np_b + 2.0 * iks_np_g);
% ydot(35) = 2.0 * iks_np_a * iks_np_C3 + 4.0 * iks_np_b * iks_np_C5 + iks_np_d * iks_np_C8 - iks_np_C4 * (iks_np_a + 3.0 * iks_np_b + 3.0 * iks_np_g);
% ydot(40) = iks_np_a * iks_np_C8 + 4.0 * iks_np_g * iks_np_C5 + 2.0 * iks_np_d * iks_np_C12 - iks_np_C9 * (3.0 * iks_np_b + 3.0 * iks_np_g + iks_np_d);
% ydot(39) = 2.0 * iks_np_a * iks_np_C7 + 3.0 * iks_np_b * iks_np_C9 + 3.0 * iks_np_g * iks_np_C4 + 2.0 * iks_np_d * iks_np_C11 - iks_np_C8 * (iks_np_a + 2.0 * iks_np_b + 2.0 * iks_np_g + iks_np_d);
% ydot(46) = iks_np_g * iks_np_C14 - iks_np_C15 * (4.0 * iks_np_d + c.iks_np_t) + iks_np_e * iks_np_O1;
% ydot(41) = iks_np_b * iks_np_C11 + iks_np_g * iks_np_C7 - iks_np_C10 * (2.0 * iks_np_a + 2.0 * iks_np_d);
% ydot(32) = iks_np_b * iks_np_C2 - iks_np_C1 * (4.0 * iks_np_a);
% ydot(42) = 2.0 * iks_np_a * iks_np_C10 + 2.0 * iks_np_b * iks_np_C12 + 2.0 * iks_np_g * iks_np_C8 + 3.0 * iks_np_d * iks_np_C13 - iks_np_C11 * (iks_np_a + iks_np_b + iks_np_g + 2.0 * iks_np_d);
% ydot(33) = 4.0 * iks_np_a * iks_np_C1 + 2.0 * iks_np_b * iks_np_C3 + iks_np_d * iks_np_C6 - iks_np_C2 * (3.0 * iks_np_a + iks_np_b + iks_np_g);
% ydot(37) = iks_np_b * iks_np_C7 + iks_np_g * iks_np_C2 - iks_np_C6 * (3.0 * iks_np_a + iks_np_d);
% %
% % cell
% %
% %

% inak_fhat_val = (f_inak - 0.1263453) / (0.9980137 - 0.1263453);
% % Effective fraction of phosphorylated INaK pumps
% inak_fhat = ifthenelse((inak_fhat_val < 0.0) , 0.0 , inak_fhat_val);
% inak_phi = c.ibar * c.pk / (1.0 + exp(-(V + 92.0) * c.FRT));
% % Phosphorylated component of INaK
% INaK_p = inak_phi * (Na / (Na + c.km_p)) ^ 3.0;
% % Non-phosphorylated component of INaK
% INaK_np = inak_phi * (Na / (Na + c.km_np)) ^ 3.0;
% % Sodium-Potassium pump current
% INaK = (1.0 - inak_fhat) * INaK_np + inak_fhat * INaK_p;
% %





% % iks_pka
% %
% iks_pka_e =  3.85250000000000013e-04 + (0.012406 -  3.85250000000000013e-04) / (1.0 + exp(c.FRT * (V + 0.064118) / 0.77992));
% iks_pka_b = 0.0033201 / (1.0 + exp(c.FRT * (V - 0.094217) / 0.95364));
% iks_pka_o = 0.0002373 * exp(-1.9742 * V * c.FRT);
% iks_pka_g = 0.56356 / (1.0 + exp(c.FRT * (V + 0.17986) / -0.58381));
% iks_pka_a = 0.0099415 / (1.0 + exp(c.FRT * (V - 0.044809) / -0.58172));
% iks_pka_d = 0.0657 * exp(-0.11899 * V * c.FRT);
% iks_pka_p =  2.26519999999999993e-04 * exp( 2.46889999999999989e-04 * V * c.FRT);
% ydot(70) = iks_pka_a * iks_pka_C13 + 2.0 * iks_pka_g * iks_pka_C12 + 4.0 * iks_pka_d * iks_pka_C15 - iks_pka_C14 * (iks_pka_b + iks_pka_g + 3.0 * iks_pka_d);
% ydot(60) = 2.0 * iks_pka_a * iks_pka_C3 + 4.0 * iks_pka_b * iks_pka_C5 + iks_pka_d * iks_pka_C8 - iks_pka_C4 * (iks_pka_a + 3.0 * iks_pka_b + 3.0 * iks_pka_g);
% ydot(72) = -(iks_pka_e + iks_pka_p) * iks_pka_O1 + iks_pka_o * iks_pka_O2 + c.iks_pka_t * iks_pka_C15;
% ydot(59) = 3.0 * iks_pka_a * iks_pka_C2 + 3.0 * iks_pka_b * iks_pka_C4 + iks_pka_d * iks_pka_C7 - iks_pka_C3 * (2.0 * iks_pka_a + 2.0 * iks_pka_b + 2.0 * iks_pka_g);
% ydot(69) = iks_pka_b * iks_pka_C14 + iks_pka_g * iks_pka_C11 - iks_pka_C13 * (iks_pka_a + 3.0 * iks_pka_d);
% ydot(62) = iks_pka_b * iks_pka_C7 + iks_pka_g * iks_pka_C2 - iks_pka_C6 * (3.0 * iks_pka_a + iks_pka_d);
% ydot(71) = iks_pka_g * iks_pka_C14 - iks_pka_C15 * (4.0 * iks_pka_d + c.iks_pka_t) + iks_pka_e * iks_pka_O1;
% ydot(65) = iks_pka_a * iks_pka_C8 + 4.0 * iks_pka_g * iks_pka_C5 + 2.0 * iks_pka_d * iks_pka_C12 - iks_pka_C9 * (3.0 * iks_pka_b + 3.0 * iks_pka_g + iks_pka_d);
% ydot(66) = iks_pka_b * iks_pka_C11 + iks_pka_g * iks_pka_C7 - iks_pka_C10 * (2.0 * iks_pka_a + 2.0 * iks_pka_d);
% ydot(64) = 2.0 * iks_pka_a * iks_pka_C7 + 3.0 * iks_pka_b * iks_pka_C9 + 3.0 * iks_pka_g * iks_pka_C4 + 2.0 * iks_pka_d * iks_pka_C11 - iks_pka_C8 * (iks_pka_a + 2.0 * iks_pka_b + 2.0 * iks_pka_g + iks_pka_d);
% ydot(58) = 4.0 * iks_pka_a * iks_pka_C1 + 2.0 * iks_pka_b * iks_pka_C3 + iks_pka_d * iks_pka_C6 - iks_pka_C2 * (3.0 * iks_pka_a + iks_pka_b + iks_pka_g);
% ydot(67) = 2.0 * iks_pka_a * iks_pka_C10 + 2.0 * iks_pka_b * iks_pka_C12 + 2.0 * iks_pka_g * iks_pka_C8 + 3.0 * iks_pka_d * iks_pka_C13 - iks_pka_C11 * (iks_pka_a + iks_pka_b + iks_pka_g + 2.0 * iks_pka_d);
% ydot(57) = iks_pka_b * iks_pka_C2 - iks_pka_C1 * (4.0 * iks_pka_a);
% ydot(68) = iks_pka_a * iks_pka_C11 + 3.0 * iks_pka_g * iks_pka_C9 + 3.0 * iks_pka_d * iks_pka_C14 - iks_pka_C12 * (2.0 * iks_pka_b + 2.0 * iks_pka_g + 2.0 * iks_pka_d);
% ydot(61) = iks_pka_a * iks_pka_C4 + iks_pka_d * iks_pka_C9 - iks_pka_C5 * (4.0 * iks_pka_b + 4.0 * iks_pka_g);
% ydot(73) = iks_pka_p * iks_pka_O1 - iks_pka_o * iks_pka_O2;
% ydot(63) = 3.0 * iks_pka_a * iks_pka_C6 + 2.0 * iks_pka_b * iks_pka_C8 + 2.0 * iks_pka_g * iks_pka_C3 + 2.0 * iks_pka_d * iks_pka_C10 - iks_pka_C7 * (2.0 * iks_pka_a + iks_pka_b + iks_pka_g + iks_pka_d);
% %
% % pka
% %
% %
% % ik1
% %
% ik1_IK1_np_vv = V - EK;
% ik1_IK1_np_alpha = 1.02 / (1.0 + exp(0.2385 * (ik1_IK1_np_vv - 59.215)));
% ik1_IK1_np_beta = (0.49124 * exp(0.08032 * (ik1_IK1_np_vv + 5.476)) + exp(0.06175 * (ik1_IK1_np_vv - 594.31))) / (1.0 + exp(-0.5143 * (ik1_IK1_np_vv + 4.753)));
% % Non-phosphorylated component of IK1
% IK1_np = c.ik1_Gbar * (ik1_IK1_np_alpha / (ik1_IK1_np_alpha + ik1_IK1_np_beta)) * ik1_IK1_np_vv;
% % CaMKII phosphorylated component of IK1
% IK1_camk = IK1_np * 1.2;
% % Inward rectifier Potassium current
% IK1 = (1.0 - f_ik1) * IK1_np + f_ik1 * IK1_camk;
% %
% % ito
% %
% alph_if = 1.0 / (1.0 + exp((V + 58.0) / 5.0));
% beta_i = 1.0 / (1.0 + exp((V + 19.0) / -9.0)) / 0.5 / 9.7953;
% ito_R = exp(V / 550.0);
% alph_is = 1.0 / (1.0 + exp((V + 60.0) / 5.0)) / 250.0;
% ito_x = c.ito_Gbar * a_np ^ 3.0 * ito_R * (V - EK);
% ito_is_camk_alpha = 2.46 * alph_is;
% ydot(82) = ito_is_camk_alpha * (1.0 - is_camk) - beta_i * is_camk;
% ito_a_np_alpha = 1.0 / (1.0 + exp((V - 18.4099) / -29.3814));
% a_inf = 1.0 / (1.0 + exp((V + 9.437) / -7.133));
% ito_a_np_beta = 1.0 / (1.0 + exp((V + 100.0) / 29.3814));
% a_tau = 1.0 / (ito_a_np_alpha / 1.2089 + 3.5 * ito_a_np_beta);
% ydot(19) = (a_inf - a_np) / a_tau;
% % CaMKII-phosphorylated component of ITo
% ITo_camk = ito_x * (0.7356 * if_camk + 0.2644 * is_camk);
% ito_if_np_alpha = 0.02144 * alph_if;
% ydot(20) = ito_if_np_alpha * (1.0 - if_np) - beta_i * if_np;
% ito_is_np_alpha = 0.56034 * alph_is;
% ydot(21) = ito_is_np_alpha * (1.0 - is_np) - beta_i * is_np;
% % Non-phosphorylated component of ITo
% ITo_np = ito_x * (0.7356 * if_np + 0.2644 * is_np);
% ito_if_camk_alpha = 0.04796 * alph_if;
% ydot(81) = ito_if_camk_alpha * (1.0 - if_camk) - beta_i * if_camk;
% % Transient outward Potassium current
% ITo = (1.0 - f_ito) * ITo_np + f_ito * ITo_camk;
% %
% % ikr
% %
% % Inactivation gate of IKr
% inx = 1.0 / (1.0 + exp((V + 10.0) / 15.4));
% % Rapid delayed rectifier Potassium current
% IKr = c.GKr * ikr_ac * inx * (V - EK);
% inf = 1.0 / (1.0 + exp((V + 10.085) / -4.25));
% ikr_ac_tau = 1.0 / (0.0006 * (V - 1.7384) / (1.0 - exp(-0.136 * (V - 1.7384))) - 0.0003 * (V + 38.3608) / (1.0 - exp(0.1522 * (V + 38.3608))));
% % Activation gate of IKr
% ydot(22) = (inf - ikr_ac) / ikr_ac_tau;


% INa_both = c.gNaBar * ina_pka_m ^ 3.0 * ina_camk_h * ina_camk_j * (V - ENa);
% INa_camk = c.gNaBar * ina_camk_m ^ 3.0 * ina_camk_h * ina_camk_j * (V - ENa);
% INa_pka = c.gNaBar * ina_pka_m ^ 3.0 * ina_pka_h * ina_pka_j * (V - ENa) * 1.25;
% INa_np = c.gNaBar * ina_np_m ^ 3.0 * ina_np_h * ina_np_j * (V - ENa);
% ina_f_pka_val = (ina_f_ina - 0.2394795) / (0.9501431 - 0.2394795);
% % Effective fraction of phosphorylated INa channels
% ina_f_pka = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);
% ina_f_both = ina_f_pka * camk_f_ina;
% ina_f_camk_only = camk_f_ina - ina_f_both;
% ina_f_pka_only = ina_f_pka - ina_f_both;
% ina_f_np = 1.0 - ina_f_pka_only - ina_f_camk_only - ina_f_both;
% INa = ina_f_np * INa_np + ina_f_pka_only * INa_pka + ina_f_camk_only * INa_camk + ina_f_both * INa_both;


% %
% % ctkcl
% %
% ctkcl_CTKCl_z1 = EK - ECl;
% % K+Cl- cotransporter
% CTKCl = c.KClBar * ctkcl_CTKCl_z1 / (ctkcl_CTKCl_z1 + c.ctkcl_CTKCl_z2);


% % inal
% %
% conductance = inal_m ^ 3.0 * inal_h * (V - ENa);
% % Non-phosphorylated component of INaL
% INaL_np = 0.0065 * conductance;
% h_inf = 1.0 / (1.0 + exp((V + 91.0) / 6.1));
% ydot(18) = (h_inf - inal_h) / c.tau_h;
% % CaMKII-phosphorylated component of INaL
% INaL_camk = 0.016 * conductance;
% inal_m_alpha = 0.32 * (V + 47.13) / (1.0 - exp(-0.1 * (V + 47.13)));
% inal_m_beta = 0.08 * exp(V / -11.0);
% ydot(17) = inal_m_alpha * (1.0 - inal_m) - inal_m_beta * inal_m;
% % Late Sodium current
% INaL = (1.0 - camk_f_ina) * INaL_np + camk_f_ina * INaL_camk;



% %
% % ikur
% %
% ikur_fhat_val = (f_ikur -  5.89379800000000009e-02) / (0.393747 -  5.89379800000000009e-02);
% % Effective fraction of phosphorylated IKur channels
% ikur_fhat = ifthenelse((ikur_fhat_val < 0.0) , 0.0 , ikur_fhat_val);
% IKur_p = c.gbar_np * (V - EK) / (1.0 + exp((36.0 - V) / 17.0)) * 3.62;
% IKur_np = c.gbar_np * (V - EK) / (1.0 + exp((15.0 - V) / 17.0));

% % Ultrarapid plateau Potassium current
% IKur = (1.0 - ikur_fhat) * IKur_np + ikur_fhat * IKur_p;
% %
% % ctnacl
% %
% ctnacl_CTNaCl_z1 = (ENa - ECl) ^ 4.0;
% % Na+Cl- cotransporter
% CTNaCl = c.NaClBar * ctnacl_CTNaCl_z1 / (ctnacl_CTNaCl_z1 + c.ctnacl_CTNaCl_z2);
% %
% % iclb
% %
% % Background chloride current
% IClb = c.iclb_Gbar * (V - ECl);
% %
% % iks_sig
% %
% % Fraction of phosphorylated IKs
% fp_iks = (IKsp + c.IKs_arn) / c.IKs_tot;

% % Fraction of phosphorylated ICaL channels
% fp_ICaL = (ICaLp + c.ICaL_arn) / c.ICaL_tot;
% % Fraction of phosphorylated RyR channels
% fp_RyR = (RyRp + c.RyR_arn) / c.RyR_tot;




% %
% % inaca
% %
% Na_ss3 = Na_sr ^ 3.0;
% Na_i3 = Na ^ 3.0;
% exp1 = exp(c.eta * V * c.FRT);
% exp2 = exp((c.eta - 1.0) * V * c.FRT);




% %
% % ical_np
% %
% ical_np_delta_tau = 5.0 * f_ical;
% % Permeability of non-phosphorylated channels
% ical_np_PCa = 0.0001552 * (1.0 + 0.4 * f_ical);
% ical_np_in_inf = 1.0 / (1.0 + exp((17.5 + V) / 3.0));
% in_a = 1.0 / (70.0 * (1.0 - 0.5 * f_ical) * (1.0 + exp((V + 49.1) / 10.349)));
% % Steady state value of activation
% ical_np_ac_inf = 1.0 / ((1.0 + exp((13.56 - V) / 9.45)) * (1.0 + exp((25.0 + V) / -5.0)));
% in_b = 1.0 / (1.0 + exp((V + 0.213) / -10.807));
% % Time constant of activation
% ac_tau = 0.59 + 0.8 * exp(0.052 * (V + 13.0)) / (1.0 + exp(0.132 * (V + 13.0)));
% % Time constant for inactivation at low Ca2+ levels
% ical_np_in_lo_tau = 1.0 / (in_a + in_b / 26.553);
% % Steady state value for inactivation at high Ca2+ levels
% ical_np_in_hi_inf = (0.001 + ical_np_in_inf) / 1.001;
% % Steady state value for inactivation at low Ca2+ levels
% ical_np_in_lo_inf = (0.2474 + ical_np_in_inf) / 1.2474;
% % Transition rate from open to closed
% ical_np_beta = (1.0 - ical_np_ac_inf) / ac_tau;
% % Transition rate from closed to open
% ical_np_alpha = ical_np_ac_inf / ac_tau;
% % Transition rate from inactive to active, at low Ca2+ levels
% ical_np_x = ical_np_in_lo_inf / ical_np_in_lo_tau;
% % Transition rate from active to inactive, at low Ca2+ levels
% ical_np_y = (1.0 - ical_np_in_lo_inf) / ical_np_in_lo_tau;
% %
% % diff
% %
% Idiff_Na = (Na_sr - Na) / c.diff_tau;
% Idiff_Cl = (Cl_sr - Cl) / c.diff_tau;
% %
% % sodium
% %
% %
% % iclca
% %
% iclca_i2_alpha = 0.025 / (1.0 + exp((V + 58.0) / 5.0));
% iclca_i2_beta = 0.2 / (1.0 + exp((V + 19.0) / -9.0));
% % IClCa (ITo2) inactivation gate
% ydot(23) = (iclca_i2_alpha / (iclca_i2_alpha + iclca_i2_beta) - i2) / c.iclca_tau;
% vexp = exp(V * c.FRT);
% IClCa_bar = c.PCl * V * c.FFRT * (Cl - c.Clo * vexp) / (1.0 - vexp);
% %
% % chloride
% %
% % Intracullular Chloride concentration
% ydot(11) = c.chloride_Cl_r1 * IClb + CTNaCl + CTKCl + c.chloride_Cl_r2 * Idiff_Cl;





% %
% % iks
% %
% iks_f_hat_val = (fp_iks - c.iks_f_hat_ratio) / (0.785 - c.iks_f_hat_ratio);
% % Effective fraction of phosphorylated IKs channels
% iks_f_hat = ifthenelse((iks_f_hat_val < 0.0) , 0.0 , iks_f_hat_val);
% %
% % camk
% %
% %
% % potassium
% %
% %
% % membrane
% %
% %
% % irel
% %
% beta_p = c.beta_0 * (1.0 + 0.0 * f_ryr);
% beta_np = c.beta_0 * (1.0 + 2.0 * f_ryr);
% irel_fhat_val = (fp_RyR - c.irel_fhat_ratio) / (0.9586 - c.irel_fhat_ratio);
% % Effective fraction of phosphorylated ryr channels
% irel_fhat = ifthenelse((irel_fhat_val < 0.0) , 0.0 , irel_fhat_val);
% alpha_p = 0.1125 * beta_p;
% alpha_np = 0.1125 * beta_np;
% Irel_pure = (1.0 - irel_fhat) * Irel_np + irel_fhat * Irel_p;
% %
% % ical
% %
% ical_f_hat_val = (fp_ICaL - c.ical_f_hat_ratio) / (0.9273 - c.ical_f_hat_ratio);
% % Effective fraction of phosphorylated ICaL channels
% ical_f_hat = ifthenelse((ical_f_hat_val < 0.0) , 0.0 , ical_f_hat_val);
% %
% % iup
% %
% iup_f_pka_val = (iup_f_plb - 0.6591) / (0.9945 - 0.6591);
% iup_f_pka = ifthenelse((iup_f_pka_val < 0.0) ,0.0 , iup_f_pka_val);



%% %%%%%%%%%%%%%%%%%%%%%%%%
% iup_f_both = iup_f_pka * camk_f_plb;
% iup_f_pka_only = iup_f_pka - iup_f_both;
% iup_f_camk_only = camk_f_plb - iup_f_both;
% iup_f_np = 1.0 - iup_f_pka_only - iup_f_camk_only - iup_f_both;
% % SERCA / PLB, Neta Ca2+ affinity
% Km_up = iup_f_np * c.Km_np + iup_f_pka_only * c.Km_pka + iup_f_camk_only * c.Km_camk + iup_f_both * c.Km_both;
% %
% % icab
% %
% %
% % ipca
% %
% %
% % calcium
% %
% calcium_fhat_val = (f_tni - 0.6735188) / (0.9991797 - 0.6735188);
% % Effective fraction of phosphorylated Troponin
% calcium_fhat = ifthenelse((calcium_fhat_val < 0.0) , 0.0 , calcium_fhat_val);

% % Combined affinity of phosphorylated and non-phosphorylated Troponin for Ca2+
% kt = (1.0 - calcium_fhat) * c.ktn + calcium_fhat * c.ktp;
% calcium_Ca_jsr_c = uCa_jsr * c.csqn_km;
% calcium_Ca_jsr_b = c.csqn_bar + c.csqn_km - uCa_jsr;
% % Concentration of Ca in the JSR subspace
% Ca_jsr = (sqrt(calcium_Ca_jsr_b * calcium_Ca_jsr_b + 4.0 * calcium_Ca_jsr_c) - calcium_Ca_jsr_b) / 2.0;
% calcium_Ca_CaL_d = -c.km_pro * uCa_CaL;
% calcium_Ca_CaL_c = c.ss_pro - uCa_CaL * c.km_sum;
% calcium_Ca_CaL_b = c.ss_sum - uCa_CaL;
% % Calcium concentration in the CaL subspace
% Ca_CaL = -calcium_Ca_CaL_b / 3.0 + 2.0 / 3.0 * sqrt(calcium_Ca_CaL_b * calcium_Ca_CaL_b - 3.0 * calcium_Ca_CaL_c) * cos(acos((9.0 * calcium_Ca_CaL_b * calcium_Ca_CaL_c - 2.0 * calcium_Ca_CaL_b * calcium_Ca_CaL_b * calcium_Ca_CaL_b - 27.0 * calcium_Ca_CaL_d) / (2.0 * (calcium_Ca_CaL_b * calcium_Ca_CaL_b - 3.0 * calcium_Ca_CaL_c) ^ 1.5)) / 3.0);
% calcium_Ca_sr_d = -c.km_pro * uCa_sr;
% calcium_Ca_sr_c = c.ss_pro - uCa_sr * c.km_sum;
% calcium_Ca_sr_b = c.ss_sum - uCa_sr;
% % Calcium concentration in the SR subspace
% Ca_sr = -calcium_Ca_sr_b / 3.0 + 2.0 / 3.0 * sqrt(calcium_Ca_sr_b * calcium_Ca_sr_b - 3.0 * calcium_Ca_sr_c) * cos(acos((9.0 * calcium_Ca_sr_b * calcium_Ca_sr_c - 2.0 * calcium_Ca_sr_b * calcium_Ca_sr_b * calcium_Ca_sr_b - 27.0 * calcium_Ca_sr_d) / (2.0 * (calcium_Ca_sr_b * calcium_Ca_sr_b - 3.0 * calcium_Ca_sr_c) ^ 1.5)) / 3.0);
% ksum = kt + c.kc;
% kpro = kt * c.kc;
% calcium_Ca_b = ksum - uCa + c.cbar + c.tbar;
% calcium_Ca_d = -kpro * uCa;
% calcium_Ca_c = kpro - uCa * ksum + c.tbar * c.kc + c.cbar * kt;
% % Intracellular calcium
% Ca = -calcium_Ca_b / 3.0 + 2.0 / 3.0 * sqrt(calcium_Ca_b * calcium_Ca_b - 3.0 * calcium_Ca_c) * cos(acos((9.0 * calcium_Ca_b * calcium_Ca_c - 2.0 * calcium_Ca_b * calcium_Ca_b * calcium_Ca_b - 27.0 * calcium_Ca_d) / (2.0 * (calcium_Ca_b * calcium_Ca_b - 3.0 * calcium_Ca_c) ^ 1.5)) / 3.0);
% %
% % ical_camk  ICaL Markov States
% %
% ical_camk_ss_cal_4 = 1.0 + (0.002 / Ca_sr) ^ 4.0;
% ical_camk_ss_cal_10 = 1.0 + (0.01 / Ca_sr) ^ 10.0;
% % Steady state value of activation
% ical_camk_ac_inf = 1.0 / ((1.0 + exp((4.798 + V) / -7.5699)) * (1.0 + exp((25.0 + V) / -5.0)));
% ical_camk_delta_tau = 0.1 * f_ical;
% ical_camk_in_inf = 1.0 / (1.0 + exp((29.979 + V) / 3.1775));
% % Permeability of non-phosphorylated channels
% ical_camk_PCa = 0.0002579 * (1.0 + 0.1 * f_ical);
% % Time constant for inactivation at low Ca2+ levels
% ical_camk_in_lo_tau = 1.0 / (in_a + in_b / 38.494);
% % Transition rate from high [Ca2+] model to low [Ca2+] model, active channel
% ical_camk_delta = 6.0 / ical_camk_ss_cal_4;
% ical_camk_IBar_vv = exp(2.0 * V * c.FRT);
% % Maximum current through non-phosphorylated channels
% ical_camk_IBar = ical_camk_PCa * 4.0 * V * c.FFRT * (Ca_sr * ical_camk_IBar_vv - 0.341 * c.Cao) / (ical_camk_IBar_vv - 1.0);
% % Transition rate from open to closed
% ical_camk_beta = (1.0 - ical_camk_ac_inf) / ac_tau;
% % Transition rate from closed to open
% ical_camk_alpha = ical_camk_ac_inf / ac_tau;
% % Steady state value for inactivation at high Ca2+ levels
% ical_camk_in_hi_inf = (0.0001 + ical_camk_in_inf) / 1.0001;
% % Ca2+ dependent component of inactivation constant for ICaL non-phosphorylated  channels and modulation by CaMKII
% ical_camk_inca = 32.5 - (18.0 - ical_camk_delta_tau) / ical_camk_ss_cal_4 - 10.0 / ical_camk_ss_cal_10;
% % Steady state value for inactivation at low Ca2+ levels
% ical_camk_in_lo_inf = (0.1 + ical_camk_in_inf) / 1.1;
% % Transition rate from active to inactive, at low Ca2+ levels
% ical_camk_y = (1.0 - ical_camk_in_lo_inf) / ical_camk_in_lo_tau;
% % Transition rate from inactive to active, at low Ca2+ levels
% ical_camk_x = ical_camk_in_lo_inf / ical_camk_in_lo_tau;
% % ICaL current density through non-phosphorylated channels
% ical_camk_ICaL = ical_camk_IBar * (ical_camk_O + ical_camk_Os);
% % Time constant for inactivation at high Ca2+ levels
% ical_camk_in_hi_tau = 1.0 / (in_a + in_b / ical_camk_inca);
% % Transition rate from active to inactive, at high Ca2+ levels
% ical_camk_ys = (1.0 - ical_camk_in_hi_inf) / ical_camk_in_hi_tau;
% ydot(49) = -(ical_camk_alpha + ical_camk_delta + ical_camk_y) * ical_camk_C + ical_camk_beta * ical_camk_O + c.ical_camk_theta * ical_camk_Cs + ical_camk_x * ical_camk_CI;
% ydot(50) = -(ical_camk_beta + ical_camk_delta + ical_camk_y) * ical_camk_O + ical_camk_alpha * ical_camk_C + c.ical_camk_theta * ical_camk_Os + ical_camk_x * ical_camk_OI;
% % Transition rate from inactive to active, at high Ca2+ levels
% ical_camk_xs = ical_camk_in_hi_inf / ical_camk_in_hi_tau;
% ydot(51) = -(ical_camk_alpha + c.ical_camk_theta + ical_camk_ys) * ical_camk_Cs + ical_camk_delta * ical_camk_C + ical_camk_beta * ical_camk_Os + ical_camk_xs * ical_camk_CIs;
% ydot(52) = -(ical_camk_beta + c.ical_camk_theta + ical_camk_ys) * ical_camk_Os + ical_camk_delta * ical_camk_O + ical_camk_alpha * ical_camk_Cs + ical_camk_xs * ical_camk_OIs;
% ical_camk_delta1_y_cor = ifthenelse((abs(ical_camk_y) < 1e-12) , 1e-12 , ical_camk_y);
% ical_camk_delta1_xs_cor = ifthenelse((abs(ical_camk_xs) < 1e-12) , 1e-12 , ical_camk_xs);
% % Transition rate from high [Ca2+] model to low [Ca2+] model, inactive channel
% ical_camk_delta1 = c.ical_camk_theta1 * (ical_camk_x * ical_camk_ys * ical_camk_delta) / (ical_camk_delta1_y_cor * ical_camk_delta1_xs_cor * c.ical_camk_theta);
% ydot(53) = -(ical_camk_alpha + ical_camk_delta1 + ical_camk_x) * ical_camk_CI + ical_camk_y * ical_camk_C + c.ical_camk_theta1 * ical_camk_CIs + ical_camk_beta * ical_camk_OI;
% ydot(54) = -(ical_camk_beta + ical_camk_delta1 + ical_camk_x) * ical_camk_OI + ical_camk_y * ical_camk_O + c.ical_camk_theta1 * ical_camk_OIs + ical_camk_alpha * ical_camk_CI;
% ydot(56) = -(ical_camk_beta + c.ical_camk_theta1 + ical_camk_xs) * ical_camk_OIs + ical_camk_ys * ical_camk_Os + ical_camk_delta1 * ical_camk_OI + ical_camk_alpha * ical_camk_CIs;
% ydot(55) = -(ical_camk_alpha + c.ical_camk_theta1 + ical_camk_xs) * ical_camk_CIs + ical_camk_ys * ical_camk_Cs + ical_camk_delta1 * ical_camk_CI + ical_camk_beta * ical_camk_OIs;
% 
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % *remaining*
% %
% % Sarcolemmal Calcium pump current
% IpCa = c.IpCa_bar * Ca / (c.Km_pCa + Ca);
% Idiff_Ca = (Ca_sr - Ca) / c.diff_tau;
% ical_np_ss_cal_10 = 1.0 + (0.012 / Ca_sr) ^ 10.0;
% Ileak_ryr_np = c.k_ryr_leak_np * exp(Ca_jsr / c.Km_ryr_leak_np) * (Ca_jsr - Ca_sr);
% ical_np_IBar_vv = exp(2.0 * V * c.FRT);
% % Maximum current through non-phosphorylated channels
% ical_np_IBar = ical_np_PCa * 4.0 * V * c.FFRT * (Ca_sr * ical_np_IBar_vv - 0.341 * c.Cao) / (ical_np_IBar_vv - 1.0);
% ical_np_ss_cal_4 = 1.0 + (0.0011 / Ca_sr) ^ 4.0;
% Itr = (Ca_nsr - Ca_jsr) / c.tau_tr;
% irel_x = 1.0 + 0.0123 / Ca_jsr;
% Ileak_ryr_p = c.k_ryr_leak_p * exp(Ca_jsr / c.Km_ryr_leak_p) * (Ca_jsr - Ca_sr);
% KClCa = 1.0 - 1.0 / (1.0 + (Irel_pure / c.kCaCl) ^ 2.0);
% % Fraction of CaMKII that has ...?
% bound = c.CaMK0 * (1.0 - trap) / (1.0 + c.Km / Ca_sr);
% vfrt = 2.0 * V * c.FRT;
% efrt = exp(vfrt);
% ICab = 1.995e-07 * 2.0 * c.F * vfrt * (Ca * efrt - 0.341 * c.Cao) / (efrt - 1.0);
% Idiff_sr = (Ca_sr - Ca_CaL) / c.tau_sr;
% % Maximum conductivity of IKs
% G = 0.19561 * (1.0 + 0.6 / (1.0 + (3.8e-05 / Ca) ^ 1.4));
% % Phosphorylated component of IKs
% IKs_pka = G * (iks_pka_O1 + iks_pka_O2) * (V - EKs);
% Ileak_ryr = (1.0 - irel_fhat) * Ileak_ryr_np + irel_fhat * Ileak_ryr_p;
% % Transition rate from high [Ca2+] model to low [Ca2+] model, active channel
% ical_np_delta = 14.9186 / ical_np_ss_cal_4;
% % Calcium dependent chloride current (aka ITo2)
% IClCa = IClCa_bar * KClCa * i2;
% % ICaL current density through non-phosphorylated channels
% ical_np_ICaL = ical_np_IBar * (ical_np_O + ical_np_Os);
% inaca_INaCaSR_denom2 = 1.0 + c.kSat * exp2;
% inaca_INaCaSR_denom3 = c.Km_Cao * Na_ss3 + c.KmNao3 * Ca_sr + c.KmNai3 * c.Cao * (1.0 + Ca_sr / c.Km_Cai);
% inaca_INaCaSR_num = 0.2 * c.vMax * (Na_ss3 * c.Cao * exp1 - c.Na_o3 * Ca_sr * exp2);
% inaca_INaCaSR_denom4 = c.Km_Cai * c.Na_o3 * (1.0 + Na_ss3 / c.KmNai3) + Na_ss3 * c.Cao + c.Na_o3 * Ca_sr;
% inaca_INaCaSR_denom1 = 1.0 + (c.Km_Ca / Ca_sr) ^ 2.0;
% % Sodium-Calcium exchange current into the SR subspace (SS,SR)
% INaCaSR = inaca_INaCaSR_num / (inaca_INaCaSR_denom1 * inaca_INaCaSR_denom2 * (inaca_INaCaSR_denom3 + inaca_INaCaSR_denom4));
% irel_tau_np = beta_np / irel_x;
% irel_tau_p = 0.5357 * beta_p / irel_x;
% % Non-phosphorylated component of IKs
% IKs_np = G * (iks_np_O1 + iks_np_O2) * (V - EKs);
% % Ca2+ dependent component of inactivation constant for ICaL non-phosphorylated  channels and modulation by CaMKII
% ical_np_inca = 13.825 - (6.3836 - ical_np_delta_tau) / ical_np_ss_cal_4 - 3.3696 / ical_np_ss_cal_10;
% % Fraction of active CaMKII
% active = bound + trap;
% inaca_INaCa_denom2 = 1.0 + c.kSat * exp2;
% inaca_INaCa_num = 0.8 * c.vMax * (Na_i3 * c.Cao * exp1 - c.Na_o3 * Ca * exp2);
% inaca_INaCa_denom3 = c.Km_Cao * Na_i3 + c.KmNao3 * Ca + c.KmNai3 * c.Cao * (1.0 + Ca / c.Km_Cai);
% inaca_INaCa_denom1 = 1.0 + (c.Km_Ca / Ca) ^ 2.0;
% inaca_INaCa_denom4 = c.Km_Cai * c.Na_o3 * (1.0 + Na_i3 / c.KmNai3) + Na_i3 * c.Cao + c.Na_o3 * Ca;
% % Sodium-Calcium exchange current
% INaCa = inaca_INaCa_num / (inaca_INaCa_denom1 * inaca_INaCa_denom2 * (inaca_INaCa_denom3 + inaca_INaCa_denom4));
% ydot(24) = -(ical_np_alpha + ical_np_delta + ical_np_y) * ical_np_C + ical_np_beta * ical_np_O + c.ical_np_theta * ical_np_Cs + ical_np_x * ical_np_CI;
% f_SERCA2a = 1.0 / (1.0 + (0.03 / active) ^ 2.0);
% camk_c = active / (active + c.camk_K);
% % Intracullular Chloride concentration in the SR subspace (SS,SR)
% ydot(12) = c.chloride_Cl_sr_r1 * IClCa - Idiff_Cl;
% PP1_tot = c.PP1_cav / c.vr_cav + c.PP1_eca / c.vr_eca + PP1f_cyt / c.vr_cyt;
% ydot(7) = c.camk_trap_alpha * bound * active - c.camk_trap_beta * trap * (0.1 + 0.9 * PP1_tot / 0.1371);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Irel = Ileak_ryr + Irel_pure;
% INa_cyt = INa + INab + INaL + INaK * 3.0 + INaCa * 3.0;
% ydot(25) = -(ical_np_beta + ical_np_delta + ical_np_y) * ical_np_O + ical_np_alpha * ical_np_C + c.ical_np_theta * ical_np_Os + ical_np_x * ical_np_OI;

% % Total Chloride current through the membrane
% ICl_tot = IClb + IClCa;
% % Total ICaL current
% ical_ICaL = (1.0 - ical_f_hat) * ical_np_ICaL + ical_f_hat * ical_camk_ICaL;
% % Slowly activating delayed rectifier Potassium current
% IKs = iks_f_hat * IKs_pka + (1.0 - iks_f_hat) * IKs_np;
% % Time constant for inactivation at high Ca2+ levels
% ical_np_in_hi_tau = 1.0 / (in_a + in_b / ical_np_inca);
% % Sodium concentration in the SR subspace
% ydot(9) = -(c.sodium_Na_sr_r1 * INaCaSR + Idiff_Na);
% Imax = (1.0 - f_SERCA2a) * c.iupmax + f_SERCA2a * c.iupmaxCAMK;
% % Unbuffered Calcium concentration in the CaL subspace
% ydot(4) = -c.calcium_uCa_CaL_r1 * ical_ICaL + c.calcium_uCa_CaL_r2 * Idiff_sr;
% ydot(88) = (camk_c - f_ical) / c.tau_cal;
% 
% % Total Na+ current through the membrane
% INa_tot = INa_cyt + INaCaSR * 3.0;
% camk_f_ryr_d = 1.0 / (1.0 + (c.camk_K / active) ^ 2.0);
% ydot(84) = (camk_f_ryr_d - f_ryr) / c.tau_ryr;
% irel_y = ical_ICaL / (1.0 + (1.0 / Ca_jsr) ^ 8.0);
% % Unbuffered Calcium concentration in the SR subspace
% ydot(3) = -(c.calcium_uCa_sr_r1 * INaCaSR + c.calcium_uCa_sr_r2 * Irel + Idiff_Ca + Idiff_sr);
% % Total Ca2+ current through the membrane
% ICa_tot = ical_ICaL + ICab + IpCa - INaCa * 2.0 - INaCaSR * 2.0;
% % Intracellular Sodium concentration
% ydot(8) = c.sodium_Na_r1 * INa_cyt + c.sodium_Na_r2 * Idiff_Na + CTNaCl;
% % Unbuffered concentration of Ca2+ in the JSR subspace
% ydot(5) = Itr - Irel;
% ydot(85) = (camk_c - f_ito) / c.tau_ito;
% % Transition rate from active to inactive, at high Ca2+ levels
% ical_np_ys = (1.0 - ical_np_in_hi_inf) / ical_np_in_hi_tau;
% % Total Potassium current through the membrane
% IK_tot = IK1 + IKr + IKs + IKur + ITo - 2.0 * INaK;
% ydot(87) = (camk_c - f_ik1) / c.tau_ik1;
% ydot(86) = (camk_c - camk_f_ina) / c.tau_ina;
% ydot(83) = (camk_c - camk_f_plb) / c.tau_plb;
% 
% % Transition rate from inactive to active, at high Ca2+ levels
% ical_np_xs = ical_np_in_hi_inf / ical_np_in_hi_tau;
% % Intracellular Potassium
% ydot(10) = c.potassium_K_r1 * (IK_tot + i_stim) + CTKCl;
% irel_inf_p = 1.9925 * alpha_p * irel_y;
% ydot(26) = -(ical_np_alpha + c.ical_np_theta + ical_np_ys) * ical_np_Cs + ical_np_delta * ical_np_C + ical_np_beta * ical_np_Os + ical_np_xs * ical_np_CIs;
% i_ion = INa_tot + IK_tot + ICa_tot + ICl_tot;
% ydot(27) = -(ical_np_beta + c.ical_np_theta + ical_np_ys) * ical_np_Os + ical_np_delta * ical_np_O + ical_np_alpha * ical_np_Cs + ical_np_xs * ical_np_OIs;
% irel_inf_np = alpha_np * irel_y;
% ydot(80) = -(irel_inf_p + Irel_p) / irel_tau_p;
% % The membrane potential
% ydot(1) = -(i_ion + i_stim);
% ydot(13) = -(irel_inf_np + Irel_np) / irel_tau_np;
% 
% ical_np_delta1_y_cor = ifthenelse((abs(ical_np_y) < 1e-12) , 1e-12 , ical_np_y);
% ical_np_delta1_xs_cor = ifthenelse((abs(ical_np_xs) < 1e-12) , 1e-12 , ical_np_xs);
% % Transition rate from high [Ca2+] model to low [Ca2+] model, inactive channel
% ical_np_delta1 = c.ical_np_theta1 * (ical_np_x * ical_np_ys * ical_np_delta) / (ical_np_delta1_y_cor * ical_np_delta1_xs_cor * c.ical_np_theta);
% uptake = Imax * Ca / (Ca + Km_up);
% leak = Imax * Ca_nsr / c.nsrmax;
% % Ca2+ uptake
% Iup = uptake - leak;
% % Concentration of Ca2+ in the NSR subspace
% ydot(6) = Iup - c.calcium_Ca_nsr_r1 * Itr;
% ydot(30) = -(ical_np_alpha + c.ical_np_theta1 + ical_np_xs) * ical_np_CIs + ical_np_ys * ical_np_Cs + ical_np_delta1 * ical_np_CI + ical_np_beta * ical_np_OIs;
% ydot(28) = -(ical_np_alpha + ical_np_delta1 + ical_np_x) * ical_np_CI + ical_np_y * ical_np_C + c.ical_np_theta1 * ical_np_CIs + ical_np_beta * ical_np_OI;
% ydot(31) = -(ical_np_beta + c.ical_np_theta1 + ical_np_xs) * ical_np_OIs + ical_np_ys * ical_np_Os + ical_np_delta1 * ical_np_OI + ical_np_alpha * ical_np_CIs;
% 
% % Unbuffered intracellular Calcium
% ydot(2) = -c.calcium_uCa_r1 * (ICab + IpCa - 2.0 * INaCa) - c.calcium_uCa_r2 * Iup + c.r3 * Idiff_Ca;
% ydot(29) = -(ical_np_beta + ical_np_delta1 + ical_np_x) * ical_np_OI + ical_np_y * ical_np_O + c.ical_np_theta1 * ical_np_OIs + ical_np_alpha * ical_np_CI;

end
