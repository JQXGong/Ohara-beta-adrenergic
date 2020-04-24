
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Quantitative analysis of variability in an integrated model of 
%%% human ventricular electrophysiology and ?-adrenergic signaling
%%% by Jingqi Q.X. Gong, Monica E. Susilo, Anna Sher, Cynthia J. Musante, Eric A. Sobie

%%% Published:April 21, 2020 DOI:https://doi.org/10.1016/j.yjmcc.2020.04.009

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions used: 
% 1. setInitialValues
% 2. GlobalConstants_Ephys
% 3. Constants_SignalingMyokit
% 4. EffectiveFraction
% 5. fun_ORd_bAR_Myokit
% 6.    Model_SignalingMyokit
% 7.    Model_Ephys
% 8. findAPD
% 9. ifthenelse

clear
clc
% % close all

format long


addpath('myfunctions')

iso_conc = 0.1 ;  %ISO concentration range from 0 to 1 uM
 
for PCL = 1000 % pacing cycle length in ms
    
global cipa
cipa = 1; % Dutta et al. ORd updated conductances
global celltype
celltype=1; %endo = 0, epi = 1, M = 2

if iso_conc == 0
    filename_isoval = ['ISO',num2str(iso_conc)];
elseif iso_conc == 1
    filename_isoval = ['ISO',num2str(iso_conc)];
elseif iso_conc < 0.1
    filename_isoval = ['ISO0p0',num2str(iso_conc*100)];
else
    filename_isoval = ['ISO0p',num2str(iso_conc*10)];
end

%% Set initial values

y0 = setInitialValues()'; 
y = y0 ;

myfilenamesave = ['OharaBA_PCL',num2str(PCL),filename_isoval,'.mat'];
    
%% Creating the settings variable and other constants
GlobalConstants_Ephys % This generates the constants for ORd model that are global variables

% Extracting signaling constants
settings = struct('runSignalingPathway',1,'runElectrophysiol',1);
settings.const_signaling = Constants_SignalingMyokit(iso_conc); 


settings.pcl = PCL;%Pacing cycle length in ms
settings.nbeats = round(1000/(PCL/1000)); %number of beats to run simulations
settings.storeLast = 5; %number of beats to store
settings.Istim = -80; %amplitude of stimulus current
settings.stimdur = 0.5; %duration of stimulus in ms
settings.stimdelay = 10; %simulus delay within each beat in ms


% Save the ISO concentration in settings
settings.ISO = iso_conc; % % 
settings.const_signaling.iso_L = settings.ISO;


%% Start the simulation

countstore = 0; 
APD90 = zeros(1,settings.nbeats) ;
CaT_max = zeros(1,settings.nbeats) ;

for countbeat = 1:1:settings.nbeats
    
    ticstartbeat = tic; %tic toc is to record the elapsed time every loop
    
    t_startbeat = settings.pcl*(countbeat - 1);
    t_endbeat = t_startbeat + settings.pcl;
    clear t_curbeat Y_curbeat; %Clear the last beat

    for count_stimint = 1:1:3  %Three simulus current interval: t = 0 to delay time, delay time to end of stimulus, end of stimulus to end of PCL
        if count_stimint == 1
            cur_interval = [t_startbeat, t_startbeat + settings.stimdelay];
            Id = 0; 
        elseif count_stimint == 2
            cur_interval = [t_startbeat +  settings.stimdelay, t_startbeat + settings.stimdelay + settings.stimdur];
            Id = settings.Istim; 
        elseif count_stimint == 3
            cur_interval = [t_startbeat + settings.stimdelay + settings.stimdur, t_endbeat];
            Id = 0; 
        end
        [t_int, Y_int] = ode15s(@(t, y) fun_ORd_bAR_Myokit(t, y, Id, settings),cur_interval,y) ;
        
        y = Y_int(end,:); %rewrite the initial values for the ode15s
        
        %Save the current beat
        if count_stimint == 1
            Y_curbeat = Y_int; 
            t_curbeat = t_int; 
        else
            Y_curbeat = [Y_curbeat; Y_int]; 
            t_curbeat = [t_curbeat; t_int];
        end
    end
    
    y = Y_curbeat(end,:); %re-write the initial values for the next beat for ode15s
    try
        APD90(countbeat) = findAPD(t_curbeat,Y_curbeat(:,1),0.90);
        CaT_max(countbeat) = max(Y_curbeat(:,6));
    catch
        APD90(countbeat) = 0;%findAPD(t_curbeat,Y_curbeat(:,1),0.90);
        CaT_max(countbeat) = 0; %max(Y_curbeat(:,6));
    end
    
    % Save the last settings.storeLast beats
    if countbeat == 1
        
        Y_store = Y_curbeat;
        t_store = t_curbeat;
        t_reset_store = t_curbeat - t_startbeat;  %This is the time that is reset every PCL
    elseif countbeat > settings.storeLast
        startbeat_idx = find(t_reset_store == 0); 
        Y_store = Y_store(startbeat_idx(2):end,:); 
        Y_store = [Y_store; Y_curbeat];
        t_store = t_store(startbeat_idx(2):end); 
        t_store = [t_store; t_curbeat];
        t_reset_store = t_reset_store(startbeat_idx(2):end);
        t_reset_store = [t_reset_store; t_curbeat - t_startbeat];  %This is the time that is reset every PCL
    else
        Y_store = [Y_store; Y_curbeat];
        t_store = [t_store; t_curbeat];
        t_reset_store = [t_reset_store; t_curbeat - t_startbeat];  %This is the time that is reset every PCL
        
    end
    

    
    fprintf('Beat# %6.0f: APD90 = %6.2f, CaTmax = %6.8f \n', countbeat, APD90(countbeat), CaT_max(countbeat))
    telapsedbeat(countbeat) = toc(ticstartbeat); %tic toc is to record the elapsed time every loop

end

avg_calctime = mean(telapsedbeat)


% 
save(myfilenamesave);

clearvars -except iso_conc

end
