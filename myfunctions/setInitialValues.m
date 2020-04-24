
% % set up the initial values for oharaBA_Ephys and oharaBA_Signaling

function initalvalue = setInitialValues()


% % oharaBA_Electrophysiology, first 1-58 state variables
    %initial conditions for state variables
    V=-87;
    Nai=7;
    Nass=Nai;
    Ki=145;
    Kss=Ki;
    Cai=1.0e-4;
    Cass=Cai;
    Cansr=1.2;
    Cajsr=Cansr;
    m=0;
    hf=1;
    hs=1;
    j=1;
    hsp=1;
    jp=1;
    mL=0;
    hL=1;
    hLp=1;
    a=0;
    iF=1;
    iS=1;
    ap=0;
    iFp=1;
    iSp=1;
    d=0;
    ff=1;
    fs=1;
    fcaf=1;
    fcas=1;
    jca=1;
    nca=0;
    ffp=1;
    fcafp=1;
    xrf=0;
    xrs=0;
    xs1=0;
    xs2=0;
    xk1=1;
    Jrelnp=0;
    Jrelp=0;
    CaMKt=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INaF PKA gates and 2 more for Both-P
    hPf = 1 ;
    hPs = 1 ;
    jP = 1 ;
    hBPf = 1 ;
    hBPs = 1 ;
    jBP = 1 ;
    % ICaL PKA gates and 2 more for Both-P
    dP = 0;
    fPf = 1;
    fPs = 1;
    fcaPf = 1;
    fcaPs = 1;
    fBPf = 1;
    fcaBPf = 1;
    % IKs PKA gates
    xs1P = 0;
    xs2P = 0;
    % INaK takes a short cut here
    % % % % % % % %
    % RyR PKA variable
    JrelP = 0;
    JrelBP = 0 ;
    
    %statevar_i is the vector for initial conditions for state variables
    initial_Ephys = [V Nai Nass Ki Kss Cai Cass Cansr Cajsr m hf hs j hsp jp ...
        mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp ...
        xrf xrs xs1 xs2 xk1 Jrelnp Jrelp JrelBP CaMKt...% all the PKA variables following
        hPf hPs jP hBPf hBPs jBP dP fPf fPs fcaPf fcaPs fBPf fcaBPf...
        xs1P xs2P JrelP ];
    
    % % then 59:end (in total 57 statevars) for signaling cascade
    % % copy from Heijman code
    
   initial_Signaling = [0.00685041638458665,0.0184627603007976,0.000731420577213056,0.00745773247314215,0.0191017987408719,0.00115141243826747,0.000607316088556676,0.000639038440072507,0.000419991861054322,0.347102959606005,9.62359241535767,0.474081735738211,0.0149041813757831,0.203016833596288,0.00944463350378086,2.49592854373432e-10,1.18055788874765e-09,7.07824478944671e-11,0.0904820284659604,0.00276490711096605,0.225475702283053,0.0326565916584703,0.192819110624505,0.205444874210056,0.174057375932567,0.817161796756964,0.567249910261073,0.249911886495890,0.0646928309115710,0.0664997605558791,0.489063888619456,0.362113356111496,0.126950532507959,0.0236821659448037,0.0128402905095188,0.00637363047239019,4.29171113639322e-05,0.00917039986149184,0.0282662056977524,0.000673713947839317,0.000765988420110534,0.592167467082831,0.673518785672382,0.239479458960528,0.126345311579566,0.00410693810508171,0.0589379755147718,0.0275455839709412,8.91799633266019e-10,0.00159632196638178,0.00209911481235842,0.000502792845976641,0.0110248953370551,1.13428924662652e-10,0.000364315164237569,0.000705656306851923,0.000341341142614041];
   initalvalue = [initial_Ephys, initial_Signaling]; 
    
