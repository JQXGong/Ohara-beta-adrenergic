function [settings] = EffectiveFraction(y, settings)
% Calculating effective fraction of phosphorylated substrates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ICaL
ICaLp = y(40+58);
fp_ICaL = (ICaLp + settings.const_signaling.ICaL_arn) / settings.const_signaling.ICaL_tot; % Fraction of phosphorylated ICaL channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_ICaL = ifthenelse((fp_ICaL < 0.0) ,0.0001 , fp_ICaL);
fp_ICaL = ifthenelse((fp_ICaL > 1.0) ,0.9999 , fp_ICaL);

ical_f_hat_val = (fp_ICaL - settings.const_signaling.ical_f_hat_ratio) / (0.9273 - settings.const_signaling.ical_f_hat_ratio);
ical_f_hat_val = ifthenelse((ical_f_hat_val < 0.0) , 0.0 , ical_f_hat_val); % Effective fraction of phosphorylated ICaL channels
ical_f_hat = ifthenelse((ical_f_hat_val > 1.0) , 1.0 , ical_f_hat_val); % Effective fraction of phosphorylated ICaL channels
%% IKs
IKsp = y(41+58);
fp_iks = (IKsp + settings.const_signaling.IKs_arn) / settings.const_signaling.IKs_tot;% Fraction of phosphorylated IKs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_iks = ifthenelse((fp_iks < 0.0) ,0.0001 , fp_iks);
fp_iks = ifthenelse((fp_iks > 1.0) ,0.9999 , fp_iks);

iks_f_hat_val = (fp_iks - settings.const_signaling.iks_f_hat_ratio) / (0.785 - settings.const_signaling.iks_f_hat_ratio);
iks_f_hat_val = ifthenelse((iks_f_hat_val < 0.0) , 0.0 , iks_f_hat_val);% Effective fraction of phosphorylated IKs channels
iks_f_hat = ifthenelse((iks_f_hat_val > 1.0) , 1.0 , iks_f_hat_val);% Effective fraction of phosphorylated IKs channels
%% Iup (PLB)
iup_f_plb = y(42+58);
iup_f_pka_val = (iup_f_plb - 0.6591) / (0.9945 - 0.6591);
iup_f_pka_val = ifthenelse((iup_f_pka_val < 0.0) ,0.0 , iup_f_pka_val);
iup_f_pka = ifthenelse((iup_f_pka_val > 1.0) ,1.0 , iup_f_pka_val);
%% Tni
f_tni = y(43+58);
calcium_fhat_val = (f_tni - 0.6735188) / (0.9991797 - 0.6735188);
calcium_fhat_val = ifthenelse((calcium_fhat_val < 0.0) , 0.0 , calcium_fhat_val);% Effective fraction of phosphorylated Troponin
calcium_fhat = ifthenelse((calcium_fhat_val > 1.0) , 1.0 , calcium_fhat_val);% Effective fraction of phosphorylated Troponin
%% INa
ina_f_ina = y(44+58);
ina_f_pka_val = (ina_f_ina - 0.2394795) / (0.9501431 - 0.2394795);
% % ina_f_pka = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
ina_f_pka_val = ifthenelse((ina_f_pka_val < 0.0) , 0.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
ina_f_pka = ifthenelse((ina_f_pka_val > 1.0) , 1.0 , ina_f_pka_val);% Effective fraction of phosphorylated INa channels
%% INaK
f_inak = y(45+58);
inak_fhat_val = (f_inak - 0.1263453) / (0.9980137 - 0.1263453);
inak_fhat_val = ifthenelse((inak_fhat_val < 0.0) , 0.0 , inak_fhat_val);% Effective fraction of phosphorylated INaK pumps
inak_fhat = ifthenelse((inak_fhat_val > 1.0) , 1.0 , inak_fhat_val);% Effective fraction of phosphorylated INaK pumps
%% RyR
RyRp = y(46+58);
fp_RyR = (RyRp + settings.const_signaling.RyR_arn) / settings.const_signaling.RyR_tot;% Fraction of phosphorylated RyR channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp_RyR = ifthenelse((fp_RyR < 0.0) ,0.0001 , fp_RyR);
fp_RyR = ifthenelse((fp_RyR > 1.0) ,0.9999 , fp_RyR);

irel_fhat_val = (fp_RyR - settings.const_signaling.irel_fhat_ratio) / (0.9586 - settings.const_signaling.irel_fhat_ratio);
irel_fhat_val = ifthenelse((irel_fhat_val < 0.0) , 0.0 , irel_fhat_val);% Effective fraction of phosphorylated ryr channels
irel_fhat = ifthenelse((irel_fhat_val > 1.0) , 1.0 , irel_fhat_val);% Effective fraction of phosphorylated ryr channels
%% IKur
f_ikur = y(47+58);
ikur_fhat_val = (f_ikur -  5.89379800000000009e-02) / (0.393747 -  5.89379800000000009e-02);
ikur_fhat_val = ifthenelse((ikur_fhat_val < 0.0) , 0.0 , ikur_fhat_val);% Effective fraction of phosphorylated IKur channels
ikur_fhat = ifthenelse((ikur_fhat_val > 1.0) , 1.0 , ikur_fhat_val);% Effective fraction of phosphorylated IKur channels

% Applying the effective fraction of phosphorylated substrates into
% settings struct...  Note: the commented equations below are from the
% code downloaded from Rudy website... The ones above are from Myokit

settings.fICaLP = ical_f_hat;%min(max(((y(98) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1);
settings.fIKsP = iks_f_hat;%min(max(((y(99) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);
settings.fPLBP = iup_f_pka;%min(max((y(100) - 6.591000e-001) / (9.945000e-001 - 6.591000e-001), 0), 1);
settings.fTnIP = calcium_fhat;%min(max((y(101) - 6.735188e-001) / (9.991797e-001 - 6.735188e-001), 0), 1);
settings.fINaP = ina_f_pka;%min(max((y(102) - 2.394795e-001) / (9.501431e-001 - 2.394795e-001), 0), 1);
settings.fINaKP = inak_fhat;%min(max((y(103) - 1.263453e-001) / (9.980137e-001 - 1.263453e-001), 0), 1);
settings.fRyRP = irel_fhat;%min(max(((y(104) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);
settings.fIKurP = ikur_fhat;%min(max((y(105) - 5.893798e-002) / (3.937470e-001 - 5.893798e-002), 0), 1);



%% Here is the code from Rudy website download to calculate effective fraction: 
% settings.fICaLP = ical_f_hat;%min(max(((y(98) + sigdata.ICaL_AKAP_PKA) / sigdata.ICaL_tot - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)) / (0.9273 - (0.0269 + sigdata.ICaL_AKAP_PKA / sigdata.ICaL_tot)), 0), 1);
% settings.fIKsP = iks_f_hat;%min(max(((y(99) + sigdata.IKs_AKAP_PKA) / sigdata.IKs_tot - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)) / (0.7850 - (0.0306 + sigdata.IKs_AKAP_PKA / sigdata.IKs_tot)), 0), 1);
% settings.fPLBP = iup_f_pka;%min(max((y(100) - 6.591000e-001) / (9.945000e-001 - 6.591000e-001), 0), 1);
% settings.fTnIP = calcium_fhat;%min(max((y(101) - 6.735188e-001) / (9.991797e-001 - 6.735188e-001), 0), 1);
% settings.fINaP = ina_f_pka;%min(max((y(102) - 2.394795e-001) / (9.501431e-001 - 2.394795e-001), 0), 1);
% settings.fINaKP = inak_fhat;%min(max((y(103) - 1.263453e-001) / (9.980137e-001 - 1.263453e-001), 0), 1);
% settings.fRyRP = irel_fhat;%min(max(((y(104) + sigdata.RyR_AKAP_PKA) / sigdata.RyR_tot - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)) / (0.9586 - (0.0329 + sigdata.RyR_AKAP_PKA / sigdata.RyR_tot)), 0), 1);
% settings.fIKurP = ikur_fhat;%min(max((y(105) - 5.893798e-002) / (3.937470e-001 - 5.893798e-002), 0), 1);


end