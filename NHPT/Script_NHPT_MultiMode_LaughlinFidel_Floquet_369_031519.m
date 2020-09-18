%% Non-Hermitian Perturbation Theory for CRyP in a Multimode Cavity
% see FewPhotonTheory.lyx

% % % CHECK EXACTLY HOW DARK COUNTS COME INTO g2 (be able to predict
% % % minimum g2 for a given dark count/signal rate in each SPCM)

% % % 09/20/18
% % % % Corrected the Laughlin fidelity defintions for higher n. Requires
% % % % accounting for the different contributions from |1,0,1| and
% % % % |0,2,0|, as a function of the # image copies around the cone tip,
% % % % i.e. modenum_step. The Laughlin state at the tip is predominantly
% % % % composed of |1,0,1|, which naively minimizes interactions by
% % % % already keeping the photons far away frome ach other, and only
% % % % includes a small amount of |0,2,0| to minimize interactions by 
% % % % using interference to eliminate the amplitude where they might have
% % % % overlapped


%% MAKE A RANDOM ATOMIC DISTRIBUTION (FIXED, for easy parameter dependence;
% just repeat this twice and check if it matters, if you're smart)

% atom cloud parameters; pretend homogeneous and symmetric (cube)
Nat =800; % 300 gives decent trends, but not enough for consistent fidelities
                % i.e. sometimes draw a funny distribution
atom_side = 230; % microns, cube side length - effectively a diameter
% atom_side_z = 14; % microns; longitudinal [i.e. from "slicing"] - effectively a diameter
Z_RMS = 20; % RMS width of gaussian distribution longitudinally
lambda_long = 1; 

XY = (rand(Nat, 2)-0.5)*atom_side; % homogeneous square
% Z = (rand(Nat, 1)-0.5)*atom_side_z; % z axis, generally thinner
Z = normrnd(0, Z_RMS, [Nat,1]); % GAUSSIAN; for flat, (rand(Nat, 1)-0.5)*atom_side_z; 
XYZ = [XY, Z]; % helper for later
X = XY(:,1); Y = XY(:,2);




%% Laughlin fidelity vs modulation amplitude (sideband/carrier ratio)
% NOTE: draw a random atomic distribution using the top cell


mult=1;

laughfidel = [];
rateout = [];    
pop02 = []; 
pop11 = []; 
popOther = []; 
phase11vs02 = [];  
rateSecond = [];  

% GammaRVals = [0.01:0.01:0.05];
% for GammaR = GammaRVals % MHz, Rydberg linewidth
%     


% % % % % % % % POSSIBLE RAMP PARAMETERS % % % % % % % % 
% rampparm = [0.4];
% XLABEL = 'Ratio g_{0,2}/g_1    |    \Omega\propto g_1';

% rampparm = -10:-10:-200;
% XLABEL = 'C_6/10^6 [MHz*microns^6]';

% rampparm = 0:-1:-6;
% XLABEL = 'P-State Detuning [MHz]';

% rampparm = 2.5:0.5:5;
% XLABEL = 'Omega [MHz]';

rampparm = [1];
XLABEL = 'g_{02}/g_1 Ratio     |    \Omega\propto g_1';

% rampparm = 13;
% XLABEL = 'g_{total} [MHz]';

% rampparm = 0.01:0.01:0.05;
% XLABEL = '\Gamma_R [MHz]';

% rampparm = 0:0.2:2;
% XLABEL = 'Detuning of L=6 relative to L=0,3 [MHz]';

% rampparm = 1.3;
% XLABEL = 'Sideband Effective Cav Linewidth [MHz]';

% rampparm = 0;
% XLABEL = 'Probe detuning from single-photon resonance [MHz]';

% rampparm = 3;
% XLABEL = 'Angular Momentum of Lowest Mode';


% rampparm = 15; % [1,2,3,4,5,6,7:2:25];
% XLABEL = 'Cloud Longitudinal RMS Thickness [microns]';4


% rampparm = 1;
% XLABEL = 'Coupling Strength Multiplier at Fixed Rotation Angle';


% rampparm = [1]; % cloud thickness multiplier, ALSO modifies g appropriately
% XLABEL = 'Thinning of ELAT, with g modified to compensate assuming constant density';


for kx= 1:length(rampparm) % MHz, Rydberg linewidth  
    MULT=1;
    GammaR = 0.05;
    ratio = rampparm(kx);
    C6 = -150*10^6; 
    DE = 0; % NOTE: nonzero P-state detuning
    % detunings
    D2 = 0; % two photon detuning
    DP = 0; % probe detuning for two-photon experiments
    thickness = 15; % overwrites z_RMS, microns
    bwaist = 90; % blue waist, microns
    Blue_gaussian = false;

    kappa = [1.3, 1.3, 1.3]; % MHz   % OLD CAVITY WAS 1.5 MHz, NEW CAVITY IS EXPECTED 1.3 MHz
    DC = [0,0,0]; % detuning of cavity modes, MANUAL 
    
    
    g_modifier = [ratio, 1, ratio]; % MULTIPLIER FOR G ON EACH MODE, i.e. Floquet style
    g_modifier = g_modifier/sqrt(sum(abs(g_modifier).^2));
%     g_modifier(1) = 0; g_modifier(3)=0; % TEMPORARY HACK
    disp('REMOVE TEMPORARY HACK');
    
    % key parameters:
    g_per_atom = 15*MULT; %5.3/sqrt(Nat); % 2
    Omega_per_atom = 2.4*MULT*g_modifier(2); % maximum
    modenum_step = 3; % i.e. first HG mode is TEM0modenum_first, second mode is TEM0(modenum_first+modenum_step)
    modenum_first = 3;
    
    % % % % PARAMETERS THAT YOU MIGHT WANT TO SCAN:
    Blue_waist_x = bwaist; Blue_waist_y = bwaist;
    % van der Waal's interaction strength
    intres = false; % true -> use C3; false -> use C6; CAN DO THIS BETTER! CAN JUST USE TRUE 2-PARTICLE POT
    C3 = 3200; % MHZ*microns^3
    
    

    Ncav =3;

    ModeType = 'LG'; % 'HG'
    FLAG_FB = false; % if false, forward only; if true, include forward and backward modes (within Ncav) 
    mode_waist = 20; %  microns; same along both axes, for now
    E_step = 0; % energy increase of each cavity mode from the previous one


    % possible Blue-beam shapes
    Blue_center_x = 0;

    

    % state lifetimes
    GammaP = 6.0001; % MHz, excited state lifetime
    



    % gaussian atomic distribution
    % length scale in microns

    % input mode:
    C=zeros(1, Ncav);
    probeon = [2];
    for ix = 1:length(probeon)
        C(probeon(ix))=1/sqrt(length(probeon));
    end
%     C=[0.2i, 1, 0.2i];
    C = C/sqrt(sum(abs(C.^2))); % ensure proper normalization
    if length(C) > Ncav
        disp('WARNING! Trying to probe on non-existent modes');
        C = C(1:Ncav); % don't allow probing on non-existent modes =P
        find = probeon<=Ncav;
        probeon = probeon(find);
    end

    % DC = DC-DC(probeon); % make the probed mode resonant
    % C(2)=0.01;% input mode: amplitude overlap with each cavity mode

    % check parameters for reasonableness
    if FLAG_FB && mod(Ncav, 2) == 1
        disp('Warning! for Forward-Backward, Ncav should be even; setting to Ncav+1');
        Ncav = Ncav+1;
        C(end+1)=0;
        DC(end+1)=DC(Ncav/2);
    end

    if FLAG_FB % if forward-backward, ensure that F and B are degenerate
       DC(Ncav/2+1:Ncav)=DC(1:Ncav/2); 
    end

%     params = struct('Ncav', Ncav, 'Nat', Nat, 'atom_side', atom_side, 'atom_side_z', atom_side_z, 'ModeType', ModeType, ...
%         'modenum_first', modenum_first, 'modenum_step', modenum_step, 'mode_waist', mode_waist, ...
%         'E_step', E_step, 'g_per_atom', g_per_atom, 'Omega_per_atom', Omega_per_atom, 'DC', DC, ...
%         'DE', DE, 'D2', D2, 'kappa', kappa, 'GammaP', GammaP, 'GammaR', GammaR, 'C6', C6, 'C', C);

    % calculate g based on cavity modes
    g = zeros(Nat, Ncav);


    for ix = 1:Ncav
        if FLAG_FB && ix>Ncav/2 % make the second half of the modes the "backward" modes
           modenum = modenum_first+(ix-1-Ncav/2)*modenum_step;
        else
           modenum = modenum_first+(ix-1)*modenum_step;
        end
    %     modenum = modenum_first+(ix-1)*modenum_step;

        switch ModeType
            case 'HG'
                field = HermiteGauss(XY, [mode_waist,mode_waist], [0, modenum]); % at each atom, for this cavity mode
            case 'LG'
                field = LaguerreGauss(XY, [mode_waist,mode_waist], modenum); % at each atom, for this cavity mode
            otherwise
                disp('Unknown mode type!');
        end
        

        % account for longitudinal phase winding, if forward-backward
        if FLAG_FB && ix<=Ncav/2  
            field = field.*exp(2*pi*1i*Z/lambda_long);
        elseif FLAG_FB && ix>Ncav/2 
            field = field.*exp(-2*pi*1i*Z/lambda_long);
        end

        g(:, ix) = g_per_atom*g_modifier(ix)*field;

        fieldsave{ix}=field/sqrt(sum(abs(field).^2)); % for checking how
                                                   % orthogonal the collective modes are
    end

    % if nothing else, blue is homogeneous:
    Omega = Omega_per_atom*ones(Nat, 1); % blue Rabi frequencies, for each atom
    % make Blue Gaussian shape
    if Blue_gaussian
        Omega = Omega_per_atom*sqrt(exp(-2*((X+Blue_center_x).^2/Blue_waist_x^2+Y.^2/Blue_waist_y^2)));
    end

    % find distances between each atom, calculate U
    dist = sqrt((repmat(X, [1, Nat])-repmat(X', [Nat, 1])).^2+(repmat(Y, [1, Nat])-repmat(Y', [Nat, 1])).^2+...
        +(repmat((thickness/Z_RMS)*Z, [1, Nat])-repmat((thickness/Z_RMS)*Z', [Nat, 1])).^2);
    
    
    if ~intres
        U = C6./dist.^6; % Rydberg Rydberg interactions
                     % for now, diagonal elements are Inf; I think that is
                     % acceptable... since no rydberg atom should be
                     % interacting with itself anyway
    else
        U= C3./dist.^3; % resonant dipole-dipole interactions 
    end


    % %% Test parameters at linear order by solving 1st order transmission (fast)

    probeMin = -50*mult;
    probeMax = 50*mult;
    dp = 0.2*mult;

    dprobes_trans =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)     

    trans = zeros(length(dprobes_trans),1);

    for dix = 1:length(dprobes_trans)
        dprobe = dprobes_trans(dix);

        dc =  DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime

        [Ac] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', false);

        for tix = 1:length(probeon)
            trans(dix, tix) = abs(Ac(probeon(tix)))^2/2.367;
        end

    end


    figure(5123); close(5123); f=figure(5123);
    set(f, 'color', 'w');
    colors = 'rbkm';
    legvals = {};
    for tix = 1:size(trans,2) % for each probed mode...
        plot(dprobes_trans , trans(:, tix), ['.', colors(tix)]);
        hold on;

        % fit transmission, to provide "effective" parameters at linear level
        pguess = [DC(1), DE, D2, GammaR, g_per_atom*sqrt(Nat)*10/atom_side/sqrt(3), Omega_per_atom, 1,0];
        pout = fitEIT1D(dprobes_trans, trans(:, tix), pguess);
        plotEIT1D(dprobes_trans, pout, ['-', colors(tix)], 'LineWidth', 2); % show result
        legvals{(tix-1)*2+1} = ['Mode ', num2str(tix)];
        legvals{(tix-1)*2+2} = ['Mode ', num2str(tix)];

        disp(['==================== Mode ', num2str(tix), ' ====================']);
        disp(['g=', num2str(pout(5))]);    
        disp(['\Omega=', num2str(pout(6))]);
        disp(['\Gamma_R=', num2str(pout(4))]);
        disp('');
    end
    xlabel('\delta_c [MHz]');
    ylabel('Transmission [frac bare cav]');
    hold off;
    title(['Effective g=', num2str(pout(5))]);
    set(f, 'Position', [50, 600, 600, 400])
    legend(legvals)


    % text(5, max(trans(tix))/2, ['g=', num2str(pout(5))], 'Interpreter', 'latex', 'FontSize', 14);    
    % text(5, max(trans(tix))/1.7, ['$\Omega=$', num2str(pout(6))], 'Interpreter', 'latex', 'FontSize', 14);
    % text(5, max(trans(tix))/1.5, ['$\Gamma_R=$', num2str(pout(4))], 'Interpreter', 'latex', 'FontSize', 14);

    % %% Run Two-Photon Non-Hermitian Perturbation Theory Calculation (slower)

    probeMin = 0+DP; %-0.4;
    probeMax = 0+DP; % 0.4;
    dp = 0.1;

    % at the moment, pairSwap results seem to depend a lot on Nat...

    probeStrength = 0.01; % perturbation parameter... shouldn't really matter, keep it small
    dprobes =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)



    twoToOne = zeros(length(dprobes),1);
    pairSwap = zeros(length(dprobes),1);
    pairColl = zeros(length(dprobes),1);
    g2_allmode = zeros(length(dprobes),1);
    Acs = []; Bccs = [];

    % parfor dix = 1:length(dprobes)
    for dix = 1:length(dprobes)
        dprobe = dprobes(dix);

        dc = DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime


        % *** NO probeStrength below.. need to add it? does it matter? **
        [Ac, Bcc,M2,V2] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', true, 'ShowTiming', false);

        twoToOne(dix) = abs(Bcc(1))^2*2/abs(Ac(1))^4;
        if Ncav>1
            pairSwap(dix) = abs(Bcc(3))^2/abs(Bcc(1))^2;
            pairColl(dix) = abs(Bcc(3))^2;
        end
        Acs(:, dix) = Ac;
        Bccs(:, dix) = Bcc;
    %     g2_allmode(dix) = sum(abs(Bcc).^2)*2/sum(abs(Ac).^2)^2;
    end

    g2_allmode
    Bccs

    % %% Visualize output
    % 
    % figure(6000); close(6000); f=figure(6000);
    % set(f, 'Position', [50, 50, 600, 400]);
    % plot(dprobes , twoToOne, 'k');
    % xlabel('\delta_c [MHz]');
    % ylabel('g_2(0)');
    % ylim([0, 1.0])

    % Visualize various auto/cross correlations, ON PROBE MODES
    corr_all = ones(Ncav, Ncav);
    ctrate = zeros(Ncav, Ncav);
    BccMod = [];
    for ix = 1:Ncav
        for jx = ix:Ncav
            totix=(2*Ncav-ix)*(ix-1)/2+jx;

            if ix==jx
               corr_all(ix, jx) = 2*abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix)*sqrt(2);
            else
               corr_all(ix, jx) = abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix);
            end
            ctrate(ix, jx) = abs(Ac(ix)^2*Ac(jx)^2);
        end
    end
    ctrate = ctrate/sum(sum(ctrate));

    f=figure(171); imagesc(corr_all)
    set(f, 'Position', [500, 500, 500, 475])
    caxis([0, 1]);
    densitycbar
    axis equal;
    xlabel('Mode j');
    ylabel('Mode i');
    title('Correlations between modes i and j');
    set(f, 'color', 'w')
    set(gca, 'XTick', 1:Ncav);
    set(gca, 'YTick', 1:Ncav);

    for ix = 1:Ncav
        for jx = ix:Ncav
            if corr_all(ix, jx)<100
                t=text(0.7+(jx-1),ix-0.2, ['g_2=',num2str(corr_all(ix, jx), '%0.2f')]);
            else
                t=text(0.7+(jx-1),ix-0.2, ['g_2>100']);
            end
            t2=text(0.55+(jx-1),ix+0.2, ['Rate~', num2str(ctrate(ix, jx), '%0.3f')]);
            set(t, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
            set(t2, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end


    % % % LAUGHLIN FIDELITY
    % % % for NCav arbitrary, 'LG', EVERY MODE INCLUDED, STATE ON-CENTER
    if Ncav>2 % && modenum_step == 1
        laugh = zeros(size(Bccs));
        laugh(3)=1;
        laugh(Ncav+1)=-1/sqrt(factorial(2*modenum_step+modenum_first)*factorial(modenum_first)/2/(factorial(modenum_step+modenum_first))^2); % See Script_InteractionFactors, LaughlinConeTip.lyx
        laugh = laugh/sqrt(sum(abs(laugh).^2));
        amp2p = Bccs/sqrt(sum(abs(Bccs).^2));
        fidel=abs(amp2p'*laugh)^2;
%         disp(['Laughlin fidelity = ', num2str(round(100*fidel)), '%']);
        rate2p = sum(abs(BccMod).^2)/2.367^2;
        rateSecond(end+1) = 100*rate2p/(sum(abs(Ac).^2)/2.367);
    end
   
    disp(['Fidelity = ', num2str(100*fidel), '%']);
    disp(['Two photon rate relative to empty cavity = ', num2str(100*rate2p), '%']);
    laughfidel(end+1)=100*fidel;
    rateout(end+1) =100*rate2p;

    
    pop02(end+1) = abs(Bccs(3).^2);
    pop11(end+1) = abs(Bccs(4).^2);
    popOther(end+1) = sum(abs(Bccs([1,2,5,6]).^2));
    phase11vs02(end+1) = angle(Bccs(4)/Bccs(3));
end

% figure; plot(GammaRVals, laughfidel, 'ob');
% hold on; plot(GammaRVals, rateout, 'or');
% legend({'Laughlin Fidelity', 'Two photon rate'});
% ylabel('%');
% xlabel('\Gamma_R')
% %%

figure(951); hold off;
plot(rampparm, laughfidel, '-ob');
hold on; plot(rampparm, rateSecond, '-or');
plot(rampparm, 10*rateSecond, '-or', 'MarkerFaceColor', 'r');
legend({'Laughlin Fidelity', 'SECOND photon rate vs. empty cav', '10x SECOND photon rate'});
ylabel('%');
xlabel(XLABEL)
title('Fidelity and data rate vs. Floquet sideband to carrier ratio');
ylim([0, 100]);
set(gcf, 'color', 'w')


figure(111951); hold off;
plot(rampparm, pop02, '-ok', 'MarkerFaceColor', 'm'); hold on;
plot(rampparm, pop11, '-ok', 'MarkerFaceColor', 'g');
plot(rampparm, popOther, '-ok');
legend({'2\gamma Rate in 02', '2\gamma Rate in 11', '2\gamma Rate other'})
set(gcf, 'color', 'w')
title('Comparison of two photon populations in relevant pairs of modes');
xlabel(XLABEL)


figure(111952); hold off;
phase11vs02(phase11vs02<0)=phase11vs02(phase11vs02<0)+2*pi;
plot(rampparm, phase11vs02, '-ok', 'MarkerFaceColor', 'k');
hold on; plot([min(rampparm), max(rampparm)], [pi, pi], '-r');
ylim([0, 2*pi]);
legend({'Phase between the 11 and 02 states'})
xlabel(XLABEL)
ylabel('Radians');


set(gcf, 'color', 'w')















%% Laughlin fidelity for ratio=0.7, vs. P-state detuning (LOOK PARTICULARLY AT PHASE!!)
% NOTE: draw a random atomic distribution using the top cell

waist = 90; % blue waist, microns

mult=1/sqrt(3);

laughfidel = [];
rateout = [];    
pop02 = [];
pop11 = [];
popOther = [];
phase11vs02 = [];

% GammaRVals = [0.01:0.01:0.05];
% for GammaR = GammaRVals % MHz, Rydberg linewidth
%     

ratios = [0.7];
DEs = [-14:1:0];


GammaR = 0.1;
for kx= 1:length(DEs) % MHz, Rydberg linewidth  
    ratio = ratios(1);
    g_modifier = [ratio, 1, ratio]; % MULTIPLIER FOR G ON EACH MODE, i.e. Floquet style
    g_modifier = g_modifier/sqrt(sum(g_modifier.^2))*sqrt(3);
    
    
    % key parameters:
    g_per_atom = 1.5*mult*40/sqrt(Nat); %5.3/sqrt(Nat); % 2
    Omega_per_atom = mult*6.0*g_modifier(2); % maximum
    modenum_step = 1; % i.e. first HG mode is TEM0modenum_first, second mode is TEM0(modenum_first+modenum_step)

    
    % % % % PARAMETERS THAT YOU MIGHT WANT TO SCAN:
    Blue_waist_x = waist; Blue_waist_y = waist;
    % van der Waal's interaction strength
    intres = false; % true -> use C3; false -> use C6; CAN DO THIS BETTER! CAN JUST USE TRUE 2-PARTICLE POT
    C6 = -6*10^6;  % 56 *10^6 for 100S; 8*56*10^6 for 121S; 5.6*10^6 for 81S
    C3 = 3200; % MHZ*microns^3
    
    
    DE = mult*DEs(kx); % NOTE: nonzero P-state detuning

    Ncav =3;

    ModeType = 'LG'; % 'HG'
    FLAG_FB = false; % if false, forward only; if true, include forward and backward modes (within Ncav)
    modenum_first = 0; 
    mode_waist = 20; %  microns; same along both axes, for now
    E_step = 0; % energy increase of each cavity mode from the previous one


    % possible Blue-beam shapes
    Blue_gaussian = false;
    Blue_center_x = 0;

    % detunings
    DC =((1:Ncav) -1)*E_step;
    D2 = 0;

    % state lifetimes

    GammaP = 6.0001; % MHz, excited state lifetime
    



    % gaussian atomic distribution
    % length scale in microns

    % input mode:
    C=zeros(1, Ncav);
    probeon = [2];
    for ix = 1:length(probeon)
        C(probeon(ix))=1/sqrt(length(probeon));
    end
    if length(C) > Ncav
        disp('WARNING! Trying to probe on non-existent modes');
        C = C(1:Ncav); % don't allow probing on non-existent modes =P
        find = probeon<=Ncav;
        probeon = probeon(find);
    end

    % DC = DC-DC(probeon); % make the probed mode resonant
    % C(2)=0.01;% input mode: amplitude overlap with each cavity mode

    % check parameters for reasonableness
    if FLAG_FB && mod(Ncav, 2) == 1
        disp('Warning! for Forward-Backward, Ncav should be even; setting to Ncav+1');
        Ncav = Ncav+1;
        C(end+1)=0;
        DC(end+1)=DC(Ncav/2);
    end

    if FLAG_FB % if forward-backward, ensure that F and B are degenerate
       DC(Ncav/2+1:Ncav)=DC(1:Ncav/2); 
    end

    params = struct('Ncav', Ncav, 'Nat', Nat, 'atom_side', atom_side, 'atom_side_z', atom_side_z, 'ModeType', ModeType, ...
        'modenum_first', modenum_first, 'modenum_step', modenum_step, 'mode_waist', mode_waist, ...
        'E_step', E_step, 'g_per_atom', g_per_atom, 'Omega_per_atom', Omega_per_atom, 'DC', DC, ...
        'DE', DE, 'D2', D2, 'kappa', kappa, 'GammaP', GammaP, 'GammaR', GammaR, 'C6', C6, 'C', C);

    % calculate g based on cavity modes
    g = zeros(Nat, Ncav);


    for ix = 1:Ncav
        if FLAG_FB && ix>Ncav/2 % make the second half of the modes the "backward" modes
           modenum = modenum_first+(ix-1-Ncav/2)*modenum_step;
        else
           modenum = modenum_first+(ix-1)*modenum_step;
        end
    %     modenum = modenum_first+(ix-1)*modenum_step;

        switch ModeType
            case 'HG'
                field = HermiteGauss(XY, [mode_waist,mode_waist], [0, modenum]); % at each atom, for this cavity mode
            case 'LG'
                field = LaguerreGauss(XY, [mode_waist,mode_waist], modenum); % at each atom, for this cavity mode
            otherwise
                disp('Unknown mode type!');
        end
        

        % account for longitudinal phase winding, if forward-backward
        if FLAG_FB && ix<=Ncav/2  
            field = field.*exp(2*pi*1i*Z/lambda_long);
        elseif FLAG_FB && ix>Ncav/2 
            field = field.*exp(-2*pi*1i*Z/lambda_long);
        end

        g(:, ix) = g_per_atom*g_modifier(ix)*field;

        fieldsave{ix}=field/sqrt(sum(abs(field).^2)); % for checking how
                                                   % orthogonal the collective modes are
    end

    % if nothing else, blue is homogeneous:
    Omega = Omega_per_atom*ones(Nat, 1); % blue Rabi frequencies, for each atom
    % make Blue Gaussian shape
    if Blue_gaussian
        Omega = Omega_per_atom*sqrt(exp(-2*((X+Blue_center_x).^2/Blue_waist_x^2+Y.^2/Blue_waist_y^2)));
    end

    % find distances between each atom, calculate U
    dist = sqrt((repmat(X, [1, Nat])-repmat(X', [Nat, 1])).^2+(repmat(Y, [1, Nat])-repmat(Y', [Nat, 1])).^2+...
        +(repmat(Z, [1, Nat])-repmat(Z', [Nat, 1])).^2);
    
    
    if ~intres
        U = C6./dist.^6; % Rydberg Rydberg interactions
                     % for now, diagonal elements are Inf; I think that is
                     % acceptable... since no rydberg atom should be
                     % interacting with itself anyway
    else
        U= C3./dist.^3; % resonant dipole-dipole interactions 
    end


    % %% Test parameters at linear order by solving 1st order transmission (fast)

    probeMin = -40*mult;
    probeMax = 40*mult;
    dp = 0.2*mult;

    dprobes_trans =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)     

    trans = zeros(length(dprobes_trans),1);

    for dix = 1:length(dprobes_trans)
        dprobe = dprobes_trans(dix);

        dc =  DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime

        [Ac] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', false);

        for tix = 1:length(probeon)
            trans(dix, tix) = abs(Ac(probeon(tix)))^2/2.367;
        end

    end


    figure(5123); close(5123); f=figure(5123);
    set(f, 'color', 'w');
    colors = 'rbkm';
    legvals = {};
    for tix = 1:size(trans,2) % for each probed mode...
        plot(dprobes_trans , trans(:, tix), ['.', colors(tix)]);
        hold on;

        % fit transmission, to provide "effective" parameters at linear level
        pguess = [DC(1), DE, D2, GammaR, g_per_atom*sqrt(Nat)/9, Omega_per_atom, 1,0];
        pout = fitEIT1D(dprobes_trans, trans(:, tix), pguess);
        plotEIT1D(dprobes_trans, pout, ['-', colors(tix)], 'LineWidth', 2); % show result
        legvals{(tix-1)*2+1} = ['Mode ', num2str(tix)];
        legvals{(tix-1)*2+2} = ['Mode ', num2str(tix)];

        disp(['==================== Mode ', num2str(tix), ' ====================']);
        disp(['g=', num2str(pout(5))]);    
        disp(['\Omega=', num2str(pout(6))]);
        disp(['\Gamma_R=', num2str(pout(4))]);
        disp('');
    end
    xlabel('\delta_c [MHz]');
    ylabel('Transmission [frac bare cav]');
    hold off;
    title(['Effective g=', num2str(pout(5))]);
    set(f, 'Position', [50, 600, 600, 400])
    legend(legvals)


    % text(5, max(trans(tix))/2, ['g=', num2str(pout(5))], 'Interpreter', 'latex', 'FontSize', 14);    
    % text(5, max(trans(tix))/1.7, ['$\Omega=$', num2str(pout(6))], 'Interpreter', 'latex', 'FontSize', 14);
    % text(5, max(trans(tix))/1.5, ['$\Gamma_R=$', num2str(pout(4))], 'Interpreter', 'latex', 'FontSize', 14);

    % %% Run Two-Photon Non-Hermitian Perturbation Theory Calculation (slower)

    probeMin = 0; %-0.4;
    probeMax = 0; % 0.4;
    dp = 0.1;

    % at the moment, pairSwap results seem to depend a lot on Nat...

    probeStrength = 0.01; % perturbation parameter... shouldn't really matter, keep it small
    dprobes =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)



    twoToOne = zeros(length(dprobes),1);
    pairSwap = zeros(length(dprobes),1);
    pairColl = zeros(length(dprobes),1);
    g2_allmode = zeros(length(dprobes),1);
    Acs = []; Bccs = [];

    % parfor dix = 1:length(dprobes)
    for dix = 1:length(dprobes)
        dprobe = dprobes(dix);

        dc = DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime


        % *** NO probeStrength below.. need to add it? does it matter? **
        [Ac, Bcc,M2,V2] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', true, 'ShowTiming', false);

        twoToOne(dix) = abs(Bcc(1))^2*2/abs(Ac(1))^4;
        if Ncav>1
            pairSwap(dix) = abs(Bcc(3))^2/abs(Bcc(1))^2;
            pairColl(dix) = abs(Bcc(3))^2;
        end
        Acs(:, dix) = Ac;
        Bccs(:, dix) = Bcc;
    %     g2_allmode(dix) = sum(abs(Bcc).^2)*2/sum(abs(Ac).^2)^2;
    end

    g2_allmode
    Bccs

    % %% Visualize output
    % 
    % figure(6000); close(6000); f=figure(6000);
    % set(f, 'Position', [50, 50, 600, 400]);
    % plot(dprobes , twoToOne, 'k');
    % xlabel('\delta_c [MHz]');
    % ylabel('g_2(0)');
    % ylim([0, 1.0])

    % Visualize various auto/cross correlations, ON PROBE MODES
    corr_all = ones(Ncav, Ncav);
    ctrate = zeros(Ncav, Ncav);
    BccMod = [];
    for ix = 1:Ncav
        for jx = ix:Ncav
            totix=(2*Ncav-ix)*(ix-1)/2+jx;

            if ix==jx
               corr_all(ix, jx) = 2*abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix)*sqrt(2);
            else
               corr_all(ix, jx) = abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix);
            end
            ctrate(ix, jx) = abs(Ac(ix)^2*Ac(jx)^2);
        end
    end
    ctrate = ctrate/sum(sum(ctrate));

    f=figure(171); imagesc(corr_all)
    set(f, 'Position', [500, 500, 500, 475])
    caxis([0, 1]);
    densitycbar
    axis equal
    xlabel('Mode j');
    ylabel('Mode i');
    title('Correlations between modes i and j');
    set(f, 'color', 'w')
    set(gca, 'XTick', 1:Ncav);
    set(gca, 'YTick', 1:Ncav);

    for ix = 1:Ncav
        for jx = ix:Ncav
            if corr_all(ix, jx)<100
                t=text(0.7+(jx-1),ix-0.2, ['g_2=',num2str(corr_all(ix, jx), '%0.2f')]);
            else
                t=text(0.7+(jx-1),ix-0.2, ['g_2>100']);
            end
            t2=text(0.55+(jx-1),ix+0.2, ['Rate~', num2str(ctrate(ix, jx), '%0.3f')]);
            set(t, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
            set(t2, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end


    % % % LAUGHLIN FIDELITY
    % % % for NCav arbitrary, 'LG', EVERY MODE INCLUDED, STATE ON-CENTER
    if Ncav>2 % && modenum_step == 1
        laugh = zeros(size(Bccs));
        laugh(3)=1;
        laugh(Ncav+1)=-1/sqrt(factorial(2*modenum_step)/2/(factorial(modenum_step))^2); % See Script_InteractionFactors, LaughlinConeTip.lyx
        laugh = laugh/sqrt(sum(abs(laugh).^2));
        amp2p = Bccs/sqrt(sum(abs(Bccs).^2));
        fidel=abs(amp2p'*laugh)^2;
%         disp(['Laughlin fidelity = ', num2str(round(100*fidel)), '%']);
        rate2p = sum(abs(BccMod).^2)/2.367^2;
    end
   
    disp(['Fidelity = ', num2str(100*fidel), '%']);
    disp(['Two photon rate relative to empty cavity = ', num2str(100*rate2p), '%']);
    laughfidel(end+1)=100*fidel;
    rateout(end+1) =100*rate2p;

    
    pop02(end+1) = abs(Bccs(3).^2);
    pop11(end+1) = abs(Bccs(4).^2);
    popOther(end+1) = sum(abs(Bccs([1,2,5,6]).^2));
    phase11vs02(end+1) = angle(Bccs(4)/Bccs(3));
end

% figure; plot(GammaRVals, laughfidel, 'ob');
% hold on; plot(GammaRVals, rateout, 'or');
% legend({'Laughlin Fidelity', 'Two photon rate'});
% ylabel('%');
% xlabel('\Gamma_R')
%%

figure(951); hold off;
plot(DEs, laughfidel, '-ob');
hold on; plot(DEs, rateout, '-or');
legend({'Laughlin Fidelity', '2\gamma rate vs. empty cav'});
ylabel('%');
xlabel('P-state detuning [MHz]')
title('Fidelity and data rate vs. Floquet sideband to carrier ratio');
ylim([0, 100]);
set(gcf, 'color', 'w')


figure(111951); hold off;
plot(DEs, pop02, '-ok', 'MarkerFaceColor', 'm'); hold on;
plot(DEs, pop11, '-ok', 'MarkerFaceColor', 'g');
plot(DEs, popOther, '-ok');
legend({'2\gamma Rate in 02', '2\gamma Rate in 11', '2\gamma Rate other'})
set(gcf, 'color', 'w')
title('Comparison of two photon populations in relevant pairs of modes');
xlabel('P-state detuning [MHz]')


figure(111952); hold off;
phase11vs02(phase11vs02<0)=phase11vs02(phase11vs02<0)+2*pi;
plot(DEs, phase11vs02, '-ok', 'MarkerFaceColor', 'k');
hold on; plot([min(DEs), max(DEs)], [pi, pi], '--r')
ylim([0, 2*pi]);
legend({'Phase between the 11 and 02 states', '\pi'})
xlabel('P-state detuning [MHz]')
ylabel('Radians');


set(gcf, 'color', 'w')






%% Laughlin fidelity vs Floquet asymmetry (-1 to +1 band ratio)
% NOTE: draw a random atomic distribution using the top cell

waist = 90; % blue waist, microns

mult=1/sqrt(3);

laughfidel = [];
rateout = [];    
pop02 = [];
pop11 = [];
popOther = [];
phase11vs02 = [];

% GammaRVals = [0.01:0.01:0.05];
% for GammaR = GammaRVals % MHz, Rydberg linewidth
%     

ratios = [1]; %0.65
asymmetries = [1]; % 0.3:0.1:3;

GammaR = 0.1;
for kx= 1:length(asymmetries) % MHz, Rydberg linewidth  
    asymmetry = asymmetries(kx);
    ratio = 1;
    g_modifier = [ratio/asymmetry, 1, ratio*asymmetry]; % MULTIPLIER FOR G ON EACH MODE, i.e. Floquet style
    g_modifier = g_modifier/sqrt(sum(g_modifier.^2))*sqrt(3);
    
    
    % key parameters:
    g_per_atom = 1.5*mult*40/sqrt(Nat); %5.3/sqrt(Nat); % 2
    Omega_per_atom = mult*6.0*g_modifier(2); % maximum
    modenum_step = 1; % i.e. first HG mode is TEM0modenum_first, second mode is TEM0(modenum_first+modenum_step)

    
    % % % % PARAMETERS THAT YOU MIGHT WANT TO SCAN:
    Blue_waist_x = waist; Blue_waist_y = waist;
    % van der Waal's interaction strength
    intres = false; % true -> use C3; false -> use C6; CAN DO THIS BETTER! CAN JUST USE TRUE 2-PARTICLE POT
    C6 = -6*10^6;  % 56 *10^6 for 100S; 8*56*10^6 for 121S; 5.6*10^6 for 81S
    C3 = 3200; % MHZ*microns^3
    
    
    DE = mult*0; % NOTE: nonzero P-state detuning

    Ncav =3;

    ModeType = 'LG'; % 'HG'
    FLAG_FB = false; % if false, forward only; if true, include forward and backward modes (within Ncav)
    modenum_first = 0; 
    mode_waist = 20; %  microns; same along both axes, for now
    E_step = 0; % energy increase of each cavity mode from the previous one


    % possible Blue-beam shapes
    Blue_gaussian = false;
    Blue_center_x = 0;

    % detunings
    DC =((1:Ncav) -1)*E_step;
    D2 = 0;

    % state lifetimes
    kappa = [1.3]; % MHz   % OLD CAVITY WAS 1.5 MHz, NEW CAVITY IS EXPECTED 1.3 MHz
    GammaP = 6.0001; % MHz, excited state lifetime
    



    % gaussian atomic distribution
    % length scale in microns

    % input mode:
    C=zeros(1, Ncav);
    probeon = [2];
    for ix = 1:length(probeon)
        C(probeon(ix))=1/sqrt(length(probeon));
    end
    if length(C) > Ncav
        disp('WARNING! Trying to probe on non-existent modes');
        C = C(1:Ncav); % don't allow probing on non-existent modes =P
        find = probeon<=Ncav;
        probeon = probeon(find);
    end

    % DC = DC-DC(probeon); % make the probed mode resonant
    % C(2)=0.01;% input mode: amplitude overlap with each cavity mode

    % check parameters for reasonableness
    if FLAG_FB && mod(Ncav, 2) == 1
        disp('Warning! for Forward-Backward, Ncav should be even; setting to Ncav+1');
        Ncav = Ncav+1;
        C(end+1)=0;
        DC(end+1)=DC(Ncav/2);
    end

    if FLAG_FB % if forward-backward, ensure that F and B are degenerate
       DC(Ncav/2+1:Ncav)=DC(1:Ncav/2); 
    end

    params = struct('Ncav', Ncav, 'Nat', Nat, 'atom_side', atom_side, 'atom_side_z', atom_side_z, 'ModeType', ModeType, ...
        'modenum_first', modenum_first, 'modenum_step', modenum_step, 'mode_waist', mode_waist, ...
        'E_step', E_step, 'g_per_atom', g_per_atom, 'Omega_per_atom', Omega_per_atom, 'DC', DC, ...
        'DE', DE, 'D2', D2, 'kappa', kappa, 'GammaP', GammaP, 'GammaR', GammaR, 'C6', C6, 'C', C);

    % calculate g based on cavity modes
    g = zeros(Nat, Ncav);


    for ix = 1:Ncav
        if FLAG_FB && ix>Ncav/2 % make the second half of the modes the "backward" modes
           modenum = modenum_first+(ix-1-Ncav/2)*modenum_step;
        else
           modenum = modenum_first+(ix-1)*modenum_step;
        end
    %     modenum = modenum_first+(ix-1)*modenum_step;

        switch ModeType
            case 'HG'
                field = HermiteGauss(XY, [mode_waist,mode_waist], [0, modenum]); % at each atom, for this cavity mode
            case 'LG'
                field = LaguerreGauss(XY, [mode_waist,mode_waist], modenum); % at each atom, for this cavity mode
            otherwise
                disp('Unknown mode type!');
        end
        

        % account for longitudinal phase winding, if forward-backward
        if FLAG_FB && ix<=Ncav/2  
            field = field.*exp(2*pi*1i*Z/lambda_long);
        elseif FLAG_FB && ix>Ncav/2 
            field = field.*exp(-2*pi*1i*Z/lambda_long);
        end

        g(:, ix) = g_per_atom*g_modifier(ix)*field;

        fieldsave{ix}=field/sqrt(sum(abs(field).^2)); % for checking how
                                                   % orthogonal the collective modes are
    end

    % if nothing else, blue is homogeneous:
    Omega = Omega_per_atom*ones(Nat, 1); % blue Rabi frequencies, for each atom
    % make Blue Gaussian shape
    if Blue_gaussian
        Omega = Omega_per_atom*sqrt(exp(-2*((X+Blue_center_x).^2/Blue_waist_x^2+Y.^2/Blue_waist_y^2)));
    end

    % find distances between each atom, calculate U
    dist = sqrt((repmat(X, [1, Nat])-repmat(X', [Nat, 1])).^2+(repmat(Y, [1, Nat])-repmat(Y', [Nat, 1])).^2+...
        +(repmat(Z, [1, Nat])-repmat(Z', [Nat, 1])).^2);
    
    
    if ~intres
        U = C6./dist.^6; % Rydberg Rydberg interactions
                     % for now, diagonal elements are Inf; I think that is
                     % acceptable... since no rydberg atom should be
                     % interacting with itself anyway
    else
        U= C3./dist.^3; % resonant dipole-dipole interactions 
    end


    % %% Test parameters at linear order by solving 1st order transmission (fast)

    probeMin = -40*mult;
    probeMax = 40*mult;
    dp = 0.2*mult;

    dprobes_trans =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)     

    trans = zeros(length(dprobes_trans),1);

    for dix = 1:length(dprobes_trans)
        dprobe = dprobes_trans(dix);

        dc =  DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime

        [Ac] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', false);

        for tix = 1:length(probeon)
            trans(dix, tix) = abs(Ac(probeon(tix)))^2/2.367;
        end

    end


    figure(5123); close(5123); f=figure(5123);
    set(f, 'color', 'w');
    colors = 'rbkm';
    legvals = {};
    for tix = 1:size(trans,2) % for each probed mode...
        plot(dprobes_trans , trans(:, tix), ['.', colors(tix)]);
        hold on;

        % fit transmission, to provide "effective" parameters at linear level
        pguess = [DC(1), DE, D2, GammaR, g_per_atom*sqrt(Nat)/9, Omega_per_atom, 1,0];
        pout = fitEIT1D(dprobes_trans, trans(:, tix), pguess);
        plotEIT1D(dprobes_trans, pout, ['-', colors(tix)], 'LineWidth', 2); % show result
        legvals{(tix-1)*2+1} = ['Mode ', num2str(tix)];
        legvals{(tix-1)*2+2} = ['Mode ', num2str(tix)];

        disp(['==================== Mode ', num2str(tix), ' ====================']);
        disp(['g=', num2str(pout(5))]);    
        disp(['\Omega=', num2str(pout(6))]);
        disp(['\Gamma_R=', num2str(pout(4))]);
        disp('');
    end
    xlabel('\delta_c [MHz]');
    ylabel('Transmission [frac bare cav]');
    hold off;
    title(['Effective g=', num2str(pout(5))]);
    set(f, 'Position', [50, 600, 600, 400])
    legend(legvals)


    % text(5, max(trans(tix))/2, ['g=', num2str(pout(5))], 'Interpreter', 'latex', 'FontSize', 14);    
    % text(5, max(trans(tix))/1.7, ['$\Omega=$', num2str(pout(6))], 'Interpreter', 'latex', 'FontSize', 14);
    % text(5, max(trans(tix))/1.5, ['$\Gamma_R=$', num2str(pout(4))], 'Interpreter', 'latex', 'FontSize', 14);

    % %% Run Two-Photon Non-Hermitian Perturbation Theory Calculation (slower)

    probeMin = 0; %-0.4;
    probeMax = 0; % 0.4;
    dp = 0.1;

    % at the moment, pairSwap results seem to depend a lot on Nat...

    probeStrength = 0.01; % perturbation parameter... shouldn't really matter, keep it small
    dprobes =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)



    twoToOne = zeros(length(dprobes),1);
    pairSwap = zeros(length(dprobes),1);
    pairColl = zeros(length(dprobes),1);
    g2_allmode = zeros(length(dprobes),1);
    Acs = []; Bccs = [];

    % parfor dix = 1:length(dprobes)
    for dix = 1:length(dprobes)
        dprobe = dprobes(dix);

        dc = DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime


        % *** NO probeStrength below.. need to add it? does it matter? **
        [Ac, Bcc,M2,V2] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', true, 'ShowTiming', false);

        twoToOne(dix) = abs(Bcc(1))^2*2/abs(Ac(1))^4;
        if Ncav>1
            pairSwap(dix) = abs(Bcc(3))^2/abs(Bcc(1))^2;
            pairColl(dix) = abs(Bcc(3))^2;
        end
        Acs(:, dix) = Ac;
        Bccs(:, dix) = Bcc;
    %     g2_allmode(dix) = sum(abs(Bcc).^2)*2/sum(abs(Ac).^2)^2;
    end

    g2_allmode
    Bccs

    % %% Visualize output
    % 
    % figure(6000); close(6000); f=figure(6000);
    % set(f, 'Position', [50, 50, 600, 400]);
    % plot(dprobes , twoToOne, 'k');
    % xlabel('\delta_c [MHz]');
    % ylabel('g_2(0)');
    % ylim([0, 1.0])

    % Visualize various auto/cross correlations, ON PROBE MODES
    corr_all = ones(Ncav, Ncav);
    ctrate = zeros(Ncav, Ncav);
    BccMod = [];
    for ix = 1:Ncav
        for jx = ix:Ncav
            totix=(2*Ncav-ix)*(ix-1)/2+jx;

            if ix==jx
               corr_all(ix, jx) = 2*abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix)*sqrt(2);
            else
               corr_all(ix, jx) = abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix);
            end
            ctrate(ix, jx) = abs(Ac(ix)^2*Ac(jx)^2);
        end
    end
    ctrate = ctrate/sum(sum(ctrate));

    f=figure(171); imagesc(corr_all)
    set(f, 'Position', [500, 500, 500, 475])
    caxis([0, 1]);
    densitycbar
    axis equal
    xlabel('Mode j');
    ylabel('Mode i');
    title('Correlations between modes i and j');
    set(f, 'color', 'w')
    set(gca, 'XTick', 1:Ncav);
    set(gca, 'YTick', 1:Ncav);

    for ix = 1:Ncav
        for jx = ix:Ncav
            if corr_all(ix, jx)<100
                t=text(0.7+(jx-1),ix-0.2, ['g_2=',num2str(corr_all(ix, jx), '%0.2f')]);
            else
                t=text(0.7+(jx-1),ix-0.2, ['g_2>100']);
            end
            t2=text(0.55+(jx-1),ix+0.2, ['Rate~', num2str(ctrate(ix, jx), '%0.3f')]);
            set(t, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
            set(t2, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end


    % % % LAUGHLIN FIDELITY
    % % % for NCav arbitrary, 'LG', EVERY MODE INCLUDED, STATE ON-CENTER
    if Ncav>2 % && modenum_step == 1
        laugh = zeros(size(Bccs));
        laugh(3)=1;
        laugh(Ncav+1)=-1/sqrt(factorial(2*modenum_step)/2/(factorial(modenum_step))^2); % See Script_InteractionFactors, LaughlinConeTip.lyx
        laugh = laugh/sqrt(sum(abs(laugh).^2));
        amp2p = Bccs/sqrt(sum(abs(Bccs).^2));
        fidel=abs(amp2p'*laugh)^2;
%         disp(['Laughlin fidelity = ', num2str(round(100*fidel)), '%']);
        rate2p = sum(abs(BccMod).^2)/2.367^2;
    end
   
    disp(['Fidelity = ', num2str(100*fidel), '%']);
    disp(['Two photon rate relative to empty cavity = ', num2str(100*rate2p), '%']);
    laughfidel(end+1)=100*fidel;
    rateout(end+1) =100*rate2p;

    
    pop02(end+1) = abs(Bccs(3).^2);
    pop11(end+1) = abs(Bccs(4).^2);
    popOther(end+1) = sum(abs(Bccs([1,2,5,6]).^2));
    phase11vs02(end+1) = angle(Bccs(4)/Bccs(3));
end

% figure; plot(GammaRVals, laughfidel, 'ob');
% hold on; plot(GammaRVals, rateout, 'or');
% legend({'Laughlin Fidelity', 'Two photon rate'});
% ylabel('%');
% xlabel('\Gamma_R')

figure(951); hold off;
plot(asymmetries, laughfidel, '-ob');
hold on; plot(asymmetries, rateout, '-or');
legend({'Laughlin Fidelity', '2\gamma rate vs. empty cav'});
ylabel('%');
xlabel('Asymmetry between 0 and 2 states (-1 and +1 bands)')
title('Fidelity and data rate vs. Floquet band asymmetry');
ylim([0, 100]);
set(gcf, 'color', 'w')


figure(111951); hold off;
plot(asymmetries, pop02, '-ok', 'MarkerFaceColor', 'm'); hold on;
plot(asymmetries, pop11, '-ok', 'MarkerFaceColor', 'g');
plot(asymmetries, popOther, '-ok');
legend({'2\gamma Rate in 02', '2\gamma Rate in 11', '2\gamma Rate other'})
set(gcf, 'color', 'w')
title('Comparison of two photon populations in relevant pairs of modes');
xlabel('Asymmetry between 0 and 2 states (-1 and +1 bands)')


figure(111952); hold off;
phase11vs02(phase11vs02<0)=phase11vs02(phase11vs02<0)+2*pi;
plot(asymmetries, phase11vs02, '-ok', 'MarkerFaceColor', 'k');
hold on; plot([min(asymmetries), max(asymmetries)], [pi, pi], '-r');
ylim([0, 2*pi]);
legend({'Phase between the 11 and 02 states'})
xlabel('Asymmetry between 0 and 2 states (-1 and +1 bands)')
ylabel('Radians');


set(gcf, 'color', 'w')






















































%% Laughlin fidelity vs g for fixed rotation angle, fixed interactions 09/19/18
% NOTE: draw a random atomic distribution using the top cell

laughfidel = [];
rateout = [];    
gout = [];

% GammaRVals = [0.01:0.01:0.05];
% for GammaR = GammaRVals % MHz, Rydberg linewidth
%     
  
mults = 0.5:1:5;

GammaR = 0.000001;
for kx= 1:length(mults) % MHz, Rydberg linewidth  
    mult = mults(kx);
    % key parameters:
    g_per_atom = 4*mult*40/sqrt(Nat); %5.3/sqrt(Nat); % 2
    Omega_per_atom = mult*1; % maximum
    modenum_step = 1; % i.e. first HG mode is TEM0modenum_first, second mode is TEM0(modenum_first+modenum_step)

    
    % % % % PARAMETERS THAT YOU MIGHT WANT TO SCAN:
    Blue_waist_x = waist; Blue_waist_y = waist;
    % van der Waal's interaction strength
    intres = false; % true -> use C3; false -> use C6; CAN DO THIS BETTER! CAN JUST USE TRUE 2-PARTICLE POT
    C6 = 10^-4*10^6;  % 56 *10^6 for 100S; 8*56*10^6 for 121S; 5.6*10^6 for 81S
    C3 = 3200; % MHZ*microns^3
    
    
    DE = mult*0; % NOTE: nonzero P-state detuning

    Ncav =3;

    ModeType = 'LG'; % 'HG'
    FLAG_FB = false; % if false, forward only; if true, include forward and backward modes (within Ncav)
    modenum_first = 0; 
    mode_waist = 20; %  microns; same along both axes, for now
    E_step = 0; % energy increase of each cavity mode from the previous one


    % possible Blue-beam shapes
    Blue_gaussian = false;
    Blue_center_x = 0;

    % detunings
    DC =((1:Ncav) -1)*E_step;
    D2 = 0;

    % state lifetimes
    kappa = [1.3]; % MHz   % OLD CAVITY WAS 1.5 MHz, NEW CAVITY IS EXPECTED 1.3 MHz
    GammaP = 6.0001; % MHz, excited state lifetime
    



    % gaussian atomic distribution
    % length scale in microns

    % input mode:
    C=zeros(1, Ncav);
    probeon = [2];
    for ix = 1:length(probeon)
        C(probeon(ix))=1/sqrt(length(probeon));
    end
    if length(C) > Ncav
        disp('WARNING! Trying to probe on non-existent modes');
        C = C(1:Ncav); % don't allow probing on non-existent modes =P
        find = probeon<=Ncav;
        probeon = probeon(find);
    end

    % DC = DC-DC(probeon); % make the probed mode resonant
    % C(2)=0.01;% input mode: amplitude overlap with each cavity mode

    % check parameters for reasonableness
    if FLAG_FB && mod(Ncav, 2) == 1
        disp('Warning! for Forward-Backward, Ncav should be even; setting to Ncav+1');
        Ncav = Ncav+1;
        C(end+1)=0;
        DC(end+1)=DC(Ncav/2);
    end

    if FLAG_FB % if forward-backward, ensure that F and B are degenerate
       DC(Ncav/2+1:Ncav)=DC(1:Ncav/2); 
    end

    params = struct('Ncav', Ncav, 'Nat', Nat, 'atom_side', atom_side, 'atom_side_z', atom_side_z, 'ModeType', ModeType, ...
        'modenum_first', modenum_first, 'modenum_step', modenum_step, 'mode_waist', mode_waist, ...
        'E_step', E_step, 'g_per_atom', g_per_atom, 'Omega_per_atom', Omega_per_atom, 'DC', DC, ...
        'DE', DE, 'D2', D2, 'kappa', kappa, 'GammaP', GammaP, 'GammaR', GammaR, 'C6', C6, 'C', C);

    % calculate g based on cavity modes
    g = zeros(Nat, Ncav);


    for ix = 1:Ncav
        if FLAG_FB && ix>Ncav/2 % make the second half of the modes the "backward" modes
           modenum = modenum_first+(ix-1-Ncav/2)*modenum_step;
        else
           modenum = modenum_first+(ix-1)*modenum_step;
        end
    %     modenum = modenum_first+(ix-1)*modenum_step;

        switch ModeType
            case 'HG'
                field = HermiteGauss(XY, [mode_waist,mode_waist], [0, modenum]); % at each atom, for this cavity mode
            case 'LG'
                field = LaguerreGauss(XY, [mode_waist,mode_waist], modenum); % at each atom, for this cavity mode
            otherwise
                disp('Unknown mode type!');
        end

        % account for longitudinal phase winding, if forward-backward
        if FLAG_FB && ix<=Ncav/2  
            field = field.*exp(2*pi*1i*Z/lambda_long);
        elseif FLAG_FB && ix>Ncav/2 
            field = field.*exp(-2*pi*1i*Z/lambda_long);
        end

        g(:, ix) = g_per_atom*field;

        fieldsave{ix}=field/sqrt(sum(abs(field).^2)); % for checking how
                                                   % orthogonal the collective modes are
    end

    % if nothing else, blue is homogeneous:
    Omega = Omega_per_atom*ones(Nat, 1); % blue Rabi frequencies, for each atom
    % make Blue Gaussian shape
    if Blue_gaussian
        Omega = Omega_per_atom*sqrt(exp(-2*((X+Blue_center_x).^2/Blue_waist_x^2+Y.^2/Blue_waist_y^2)));
    end

    % find distances between each atom, calculate U
    dist = sqrt((repmat(X, [1, Nat])-repmat(X', [Nat, 1])).^2+(repmat(Y, [1, Nat])-repmat(Y', [Nat, 1])).^2+...
        +(repmat(Z, [1, Nat])-repmat(Z', [Nat, 1])).^2);
    
    
    if ~intres
        U = C6./dist.^6; % Rydberg Rydberg interactions
                     % for now, diagonal elements are Inf; I think that is
                     % acceptable... since no rydberg atom should be
                     % interacting with itself anyway
    else
        U= C3./dist.^3; % resonant dipole-dipole interactions 
    end


    % %% Test parameters at linear order by solving 1st order transmission (fast)

    probeMin = -20*mult;
    probeMax = 20*mult;
    dp = 0.1*mult;

    dprobes_trans =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)     

    trans = zeros(length(dprobes_trans),1);

    for dix = 1:length(dprobes_trans)
        dprobe = dprobes_trans(dix);

        dc =  DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime

        [Ac] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', false);

        for tix = 1:length(probeon)
            trans(dix, tix) = abs(Ac(probeon(tix)))^2;
        end

    end


    figure(5123); close(5123); f=figure(5123);
    set(f, 'color', 'w');
    colors = 'rbkm';
    legvals = {};
    for tix = 1:size(trans,2) % for each probed mode...
        plot(dprobes_trans , trans(:, tix), ['.', colors(tix)]);
        hold on;

        % fit transmission, to provide "effective" parameters at linear level
        pguess = [DC(1), DE, D2, GammaR, g_per_atom*sqrt(Nat)/9, Omega_per_atom, 1,0];
        pout = fitEIT1D(dprobes_trans, trans(:, tix), pguess);
        plotEIT1D(dprobes_trans, pout, ['-', colors(tix)], 'LineWidth', 2); % show result
        legvals{(tix-1)*2+1} = ['Mode ', num2str(tix)];
        legvals{(tix-1)*2+2} = ['Mode ', num2str(tix)];

        disp(['==================== Mode ', num2str(tix), ' ====================']);
        disp(['g=', num2str(pout(5))]);    
        disp(['\Omega=', num2str(pout(6))]);
        disp(['\Gamma_R=', num2str(pout(4))]);
        disp('');
    end
    xlabel('\delta_c [MHz]');
    ylabel('Transmission [arb]');
    hold off;
    title(['Effective g=', num2str(pout(5))]);
    set(f, 'Position', [50, 600, 600, 400])
    legend(legvals)


    % text(5, max(trans(tix))/2, ['g=', num2str(pout(5))], 'Interpreter', 'latex', 'FontSize', 14);    
    % text(5, max(trans(tix))/1.7, ['$\Omega=$', num2str(pout(6))], 'Interpreter', 'latex', 'FontSize', 14);
    % text(5, max(trans(tix))/1.5, ['$\Gamma_R=$', num2str(pout(4))], 'Interpreter', 'latex', 'FontSize', 14);

    % %% Run Two-Photon Non-Hermitian Perturbation Theory Calculation (slower)

    probeMin = 0; %-0.4;
    probeMax = 0; % 0.4;
    dp = 0.1;

    % at the moment, pairSwap results seem to depend a lot on Nat...

    probeStrength = 0.01; % perturbation parameter... shouldn't really matter, keep it small
    dprobes =  (probeMin:dp:probeMax); % probe frequency, gets subtracted from all of the other detunings 
                % (which are in the rotating frame of the probe)



    twoToOne = zeros(length(dprobes),1);
    pairSwap = zeros(length(dprobes),1);
    pairColl = zeros(length(dprobes),1);
    g2_allmode = zeros(length(dprobes),1);
    Acs = []; Bccs = [];

    % parfor dix = 1:length(dprobes)
    for dix = 1:length(dprobes)
        dprobe = dprobes(dix);

        dc = DC+1i* kappa/2-dprobe; % list of cavity mode detunings
        de =  DE+1i* GammaP/2-dprobe; % single-photon detuning
        d2 =  D2+1i* GammaR/2-dprobe; % two-photon detuning, Rydberg lifetime


        % *** NO probeStrength below.. need to add it? does it matter? **
        [Ac, Bcc,M2,V2] = NHPT_MM(C, dc, de, d2, g, Omega, U, 'TwoPhoton', true, 'ShowTiming', false);

        twoToOne(dix) = abs(Bcc(1))^2*2/abs(Ac(1))^4;
        if Ncav>1
            pairSwap(dix) = abs(Bcc(3))^2/abs(Bcc(1))^2;
            pairColl(dix) = abs(Bcc(3))^2;
        end
        Acs(:, dix) = Ac;
        Bccs(:, dix) = Bcc;
    %     g2_allmode(dix) = sum(abs(Bcc).^2)*2/sum(abs(Ac).^2)^2;
    end

    g2_allmode
    Bccs

    % %% Visualize output
    % 
    % figure(6000); close(6000); f=figure(6000);
    % set(f, 'Position', [50, 50, 600, 400]);
    % plot(dprobes , twoToOne, 'k');
    % xlabel('\delta_c [MHz]');
    % ylabel('g_2(0)');
    % ylim([0, 1.0])

    % Visualize various auto/cross correlations, ON PROBE MODES
    corr_all = ones(Ncav, Ncav);
    ctrate = zeros(Ncav, Ncav);
    BccMod = [];
    for ix = 1:Ncav
        for jx = ix:Ncav
            totix=(2*Ncav-ix)*(ix-1)/2+jx;

            if ix==jx
               corr_all(ix, jx) = 2*abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix)*sqrt(2);
            else
               corr_all(ix, jx) = abs(Bcc(totix))^2/abs(Ac(ix)^2*Ac(jx)^2);
               BccMod(totix) = Bcc(totix);
            end
            ctrate(ix, jx) = abs(Ac(ix)^2*Ac(jx)^2);
        end
    end
    ctrate = ctrate/sum(sum(ctrate));

    f=figure(171); imagesc(corr_all)
    set(f, 'Position', [500, 500, 500, 475])
    caxis([0, 1]);
    densitycbar
    axis equal
    xlabel('Mode j');
    ylabel('Mode i');
    title('Correlations between modes i and j');
    set(f, 'color', 'w')
    set(gca, 'XTick', 1:Ncav);
    set(gca, 'YTick', 1:Ncav);

    for ix = 1:Ncav
        for jx = ix:Ncav
            if corr_all(ix, jx)<100
                t=text(0.7+(jx-1),ix-0.2, ['g_2=',num2str(corr_all(ix, jx), '%0.2f')]);
            else
                t=text(0.7+(jx-1),ix-0.2, ['g_2>100']);
            end
            t2=text(0.55+(jx-1),ix+0.2, ['Rate~', num2str(ctrate(ix, jx), '%0.3f')]);
            set(t, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
            set(t2, 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        end
    end


    % % % LAUGHLIN FIDELITY
    % % % for NCav arbitrary, 'LG', EVERY MODE INCLUDED, STATE ON-CENTER
    if Ncav>2 % && modenum_step == 1
        laugh = zeros(size(Bccs));
        laugh(3)=1;
        laugh(Ncav+1)=-1/sqrt(factorial(2*modenum_step)/2/(factorial(modenum_step))^2); % See Script_InteractionFactors, LaughlinConeTip.lyx
        laugh = laugh/sqrt(sum(abs(laugh).^2));
        amp2p = Bccs/sqrt(sum(abs(Bccs).^2));
        fidel=abs(amp2p'*laugh)^2;
%         disp(['Laughlin fidelity = ', num2str(round(100*fidel)), '%']);
        rate2p = sum(abs(BccMod).^2)/2.367^2;
    end
   
    disp(['Fidelity = ', num2str(100*fidel), '%']);
    disp(['Two photon rate relative to empty cavity = ', num2str(100*rate2p), '%']);
    laughfidel(end+1)=100*fidel;
    rateout(end+1) =100*rate2p;
    gout(end+1) = pout(5);
end

% figure; plot(GammaRVals, laughfidel, 'ob');
% hold on; plot(GammaRVals, rateout, 'or');
% legend({'Laughlin Fidelity', 'Two photon rate'});
% ylabel('%');
% xlabel('\Gamma_R')

figure; plot(gout, laughfidel, 'ob');
hold on; plot(gout, rateout, 'or');
legend({'Laughlin Fidelity', '2\gamma rate vs. empty cav'});
ylabel('%');
xlabel('g [MHz]')
ylim([0,100]);










































