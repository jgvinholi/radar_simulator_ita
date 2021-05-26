%% Clearing workspace
clear
close all
clc
%% Defining Parameters
% Given Parameters
Pt = 1e5; % Potência de Transmissão
fp = 3e9; % Frequência de Transmissão
Tp = 2e-6; % Duração do Pulso
R = 6e4; % Distância do Alvo ao Radar
L = 6; % Perdas do Radar em dB
Gt = 1; Gr = Gt; % Linear !!!!!!!(Ganho de Transmissão e Recepção)
Sigma = 1; % (Radar Cross Section)
% Defining Missing Parameters
c = physconst('lightspeed');
fs = 2/Tp; % Sample Rate
PRI = Tp/0.004; % Pulse Repetition Interval for 0,4% Duty Cicle
Lamb = c/fp; % Wavelength
Nofig = 0; % Noise Figure in dB
% keep memory requirements low
Npulse = 1e3;
Npulsebuffsize = 10000;
% Defining Noise Power Based on Desired SNR
dSNR = [0 5 10 15 20];
No = Pt.*Gt.*Gr.*Lamb.^2.*Sigma./(((4.*pi).^3.*R.^4.*db2pow(L)).*db2pow(dSNR));
Notemp = No./(physconst('Boltzmann').*fs.*db2pow(Nofig));
Notemp = Notemp*sqrt(2);
%% Transmition
% Waveform
hwav = phased.RectangularWaveform('SampleRate',fs,'PulseWidth',Tp,...
    'PRF',1/PRI,'OutputFormat','Pulses','NumPulses',1);
% Antenna
hant = phased.IsotropicAntennaElement('FrequencyRange',...
    [1e9 10e9]);
% Antenna Plataform
htxplat = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
% Transmitter
htx = phased.Transmitter('PeakPower',Pt,'Gain',pow2db(Gt),...
    'LossFactor',L,'InUseOutputPort',true,...
    'CoherentOnTransmit',true);
% Radiator
hrad = phased.Radiator('Sensor',hant,...
    'PropagationSpeed',c,...
    'OperatingFrequency',fp);
%% Target
% Target Model
htgt{1} = phased.RadarTarget('Model','Nonfluctuating',...
    'MeanRCS',Sigma,'PropagationSpeed',c,...
    'OperatingFrequency',fp);
% Target Plataform
htgtplat{1} = phased.Platform('InitialPosition',[R; 0; 0],...
    'Velocity',[0;0;0]);
% Target Model (noise)
htgt{2} = phased.RadarTarget('Model','Nonfluctuating',...
    'MeanRCS',0,'PropagationSpeed',c,...
    'OperatingFrequency',fp);
% Target Plataform (noise)
htgtplat{2} = phased.Platform('InitialPosition',[R; 0; 0],...
    'Velocity',[0;0;0]);
%% Propagation
% Channel 1 (target)
hspace{1} = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fp,'TwoWayPropagation',true,...
    'SampleRate',fs);
% Channel 2 (noise)
hspace{2} = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fp,'TwoWayPropagation',true,...
    'SampleRate',fs);
%% Reception
% Waveform Collection
hcol = phased.Collector('Sensor',hant,...
    'PropagationSpeed',c,...
    'Wavefront','Plane','OperatingFrequency',fp);
% Receiver
hrec = phased.ReceiverPreamp('Gain',pow2db(Gr),'NoiseMethod','Noise temperature',...
        'ReferenceTemperature',Notemp(1),'SampleRate',fs,...
        'NoiseFigure',Nofig,'EnableInputPort',true);
% Matched Filter
hmf = phased.MatchedFilter(...
   'Coefficients',getMatchedFilter(hwav),...
   'GainOutputPort',false);
% Get group delay of matched filter
Gd = length(hmf.Coefficients)-1;

%% Implementing the Radar Model
% Time steps and grid
Tstp = PRI;
Tgrid = unigrid(0,1/fs,Tstp,'[)');
rangegates = c*Tgrid/2;
% for fix Tx
txpos = htxplat.InitialPosition;
txvel = htxplat.Velocity;
% for fix target
tgtpos = htgtplat{1}.InitialPosition;
tgtvel = htgtplat{1}.Velocity;
% Range and Angle to Target
[tgtrng,tgtang] = rangeangle(htgtplat{1}.InitialPosition,...
    htxplat.InitialPosition);
rangeidx = val2ind(tgtrng,rangegates(2)-rangegates(1),rangegates(1))-1;
% for fix PRF waveform
sigwav = step(hwav);
% for coherent transmition
[sigtx,txstatus] = step(htx,sigwav);
% for unchanging radiation angle
sigrad = step(hrad,sigtx,tgtang);
% Allocate array for received echoes
rxsig = zeros(fs*Tstp,Npulsebuffsize); 
mfsig = zeros(fs*Tstp,Npulsebuffsize); 
% Allocate arrays for hypoteses
h1 = zeros(1,Npulse);
h0 = zeros(1,Npulse);
h1a = zeros(length(No),Npulse);
h0a = zeros(length(No),Npulse);
% buferize
Nbuff = floor(Npulse/Npulsebuffsize);
Nrem = Npulse - Nbuff*Npulsebuffsize;
for k=1:length(No)
    release(hrec)
    hrec.ReferenceTemperature = Notemp(k);
    for r = 1:2
        % for unchanging range to target
        sigbounce = step(hspace{r},sigrad,txpos,tgtpos,txvel,tgtvel);
        % for nonfluctuating target
        sigtgt = step(htgt{r},sigbounce);
        % for unchanging angle to target
        sigcol = step(hcol,sigtgt,tgtang);
        for m = 1:Nbuff
            for n = 1:Npulsebuffsize
                % Update the Tx position
                % [txpos,txvel] = step(htxplat,Tstp);
                % Update the target position
                % [tgtpos,tgtvel] = step(htgtplat,Tstp);
                % Get the range and angle to the target
                % [tgtrng,tgtang] = rangeangle(tgtpos,txpos);
                % Generate the pulse
                % sigwav = step(hwav);
                % Transmit the pulse. Output transmitter status
                % [sigtx,txstatus] = step(htx,sigwav);
                % Radiate the pulse toward the target
                % sigrad = step(hrad,sigtx,tgtang);
                % Propagate the pulse to the target and back in free space
                % sigbounce = step(hspace,sigrad,txpos,tgtpos,txvel,tgtvel);
                % Reflect the pulse off the target
                % sigtgt = step(htgt,sigbounce);
                % Collect the echo from the incident angle at the antenna
                % sigcol = step(hcol,sigtgt,tgtang);
                % Receive the echo at the antenna when not transmitting
                rxsig(:,n) = (step(hrec,sigcol,~txstatus));
            end
            % Apply Matched Filter
            mfsig = real(exp(1i*4*pi/Lamb*R)*step(hmf,rxsig));
            % Shift the matched filter output
            mfsig=[mfsig(Gd+1:end,:); mfsig(1:Gd,:)];
            if r == 1
                h1(1,(1:Npulsebuffsize) + (m-1)*Npulsebuffsize) = mfsig(rangeidx,:);
            else
                h0(1,(1:Npulsebuffsize) + (m-1)*Npulsebuffsize) = mfsig(rangeidx,:);
            end
        end
        if (Nrem > 0)
            for n = 1:Nrem
                % Update the Tx position
                % [txpos,txvel] = step(htxplat,Tstp);
                % Update the target position
                % [tgtpos,tgtvel] = step(htgtplat,Tstp);
                % Get the range and angle to the target
                % [tgtrng,tgtang] = rangeangle(tgtpos,txpos);
                % Generate the pulse
                % sigwav = step(hwav);
                % Transmit the pulse. Output transmitter status
                % [sigtx,txstatus] = step(htx,sigwav);
                % Radiate the pulse toward the target
                % sigrad = step(hrad,sigtx,tgtang);
                % Propagate the pulse to the target and back in free space
                % sigbounce = step(hspace,sigrad,txpos,tgtpos,txvel,tgtvel);
                % Reflect the pulse off the target
                % sigtgt = step(htgt,sigbounce);
                % Collect the echo from the incident angle at the antenna
                % sigcol = step(hcol,sigtgt,tgtang);
                % Receive the echo at the antenna when not transmitting
                rxsig(:,n) = step(hrec,sigcol,~txstatus);
            end
            % Apply Matched Filter
            mfsig = real(exp(1i*4*pi/Lamb*R)*step(hmf,rxsig));
            % Shift the matched filter output
            mfsig=[mfsig(Gd+1:end,:); mfsig(1:Gd,:)];
            if r == 1
                h1(1,(1:Nrem) + Nbuff*Npulsebuffsize) = mfsig(rangeidx,1:Nrem);
            else
                h0(1,(1:Nrem) + Nbuff*Npulsebuffsize) = mfsig(rangeidx,1:Nrem);
            end
        end
        if r == 1
            % Noncoherently integrate the received echoes
            %rxsig(:,k) = pulsint(rxsig(:,:,k),'noncoherent');
            % Plot rangegates power
            figure(k)
            plot(rangegates,(mfsig(:,1))); hold on;
            xlabel('Meters'); ylabel('Amplitude');
            title('Matched Signal Amplitude per Distance')
            plot([R,R],[0 max(abs(mfsig(:,1)))],'r--');
        end
    end
    range_estimates
    %% Create Histogram of Outputs
    h1a(k,:) = (h1(1,:));
    h0a(k,:) = (h0(1,:));
    thresh_low = min([h1a(k,:), h0a(k,:)]);
    thresh_hi  = max([h1a(k,:), h0a(k,:)])^2;
    nbins = 100;
    binedges = linspace(thresh_low,thresh_hi,nbins);
    figure(length(No)+k)
    histogram(h0a(k,:),binedges)
    hold on
    histogram(h1a(k,:),binedges)
    hold off
    title('Target-Absent Vs Target-Present Histograms')
    xlabel('Signal Amplitude'); ylabel('Number of Signals');
    legend('Target Absent','Target Present');
    %% Compare Simulated and Theoretical Pd and Pfa
    nbins = 1000;
    thresh_steps = linspace(thresh_low,thresh_hi,nbins);
    sim_pd = zeros(1,nbins);
    sim_pfa = zeros(1,nbins);
    for r = 1:nbins
        thresh = thresh_steps(r);
        sim_pd(r) = sum(h1a(k,:) >= thresh);
        sim_pfa(r) = sum(h0a(k,:) >= thresh);
    end
    sim_pd = sim_pd/(Npulse);
    sim_pfa = sim_pfa/(Npulse);
    pfa_diff = diff(sim_pfa);
    idx = (pfa_diff == 0);
    sim_pfa(idx) = [];
    sim_pd(idx) = [];
    minpfa = 1e-6;
    N = sum(sim_pfa >= minpfa);
    sim_pfa = fliplr(sim_pfa(1:N)).';
    sim_pd = fliplr(sim_pd(1:N)).';
    % plot the theoretical Pfa curve
    figure(2*length(No)+k)
    rocsnr(dSNR(k),'SignalType',...
        'NonfluctuatingNoncoherent',...
        'MinPfa',minpfa,'NumPoints',N,'NumPulses',1);
    hold on
    semilogx(sim_pfa,sim_pd,'r.')
    title('Simulated and Theoretical ROC Curves')
    xlabel('Pfa')
    ylabel('Pd'), ylim([-inf 1])
    grid on
    legend('Theoretical','Simulated','Location','SE');
end
