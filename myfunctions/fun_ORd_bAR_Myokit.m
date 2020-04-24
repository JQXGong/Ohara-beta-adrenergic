function dy = fun_ORd_bAR_Myokit(t, y, Istim, settings)

dy_bAR = zeros(57,1);
if settings.runSignalingPathway == 1
    dy_bAR = Model_SignalingMyokit(y(59:end), settings.const_signaling);
    [settings] = EffectiveFraction(y, settings);  % Calculate the effective fraction of phosphorylated substrates, and save it in settings. These will be read by the Ephys function
else
    [settings] = EffectiveFraction(y, settings);  % Calculate the effective fraction of phosphorylated substrates, and save it in settings. These will be read by the Ephys function
end

dy_ephys = zeros(58,1); 
if settings.runElectrophysiol == 1
    settings.Whole_cell_PP1 = 0.1371;  
    if settings.runSignalingPathway == 1
        pp1_PP1f_cyt_sum = settings.const_signaling.pp1_K - settings.const_signaling.PP1_cyt + y(97);
        % Concentration of uninhibited PP1 in the cytosolic compartment
        PP1f_cyt = 0.5 * (sqrt(pp1_PP1f_cyt_sum ^ 2.0 + 4.0 * settings.const_signaling.pp1_K * settings.const_signaling.PP1_cyt) - pp1_PP1f_cyt_sum);
        PP1_tot = settings.const_signaling.PP1_cav / settings.const_signaling.vr_cav + settings.const_signaling.PP1_eca / settings.const_signaling.vr_eca + PP1f_cyt / settings.const_signaling.vr_cyt;
        settings.Whole_cell_PP1 = PP1_tot;  %This is to calculate the fraction of CaMK that has trapped calmodulin or dCaMKt or ydot(42)

        dy_ephys = Model_Ephys(t, y(1:58), Istim, settings);
        
    end
end

dy = [dy_ephys; dy_bAR]; 


end