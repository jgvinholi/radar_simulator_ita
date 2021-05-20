%% Clearing workspace
clear
close all
clc
%% Defining Given Parameters
Pt = 1e5; % Potência de Transmissão
fp = 3e9; % Frequência de Transmissão
Tp = 2e-6; % Duração do Pulso
R = 6e4; % Distância do Alvo ao Radar
L = 6; % Perdas do Radar em dB
Gt = 1; Gr = Gt; % (Ganho de Transmissão e Recepção)
Sigma = 1; % (Radar Cross Section)
%% Defining Missing Parameters
Srate = 1e6; % Sample Rate
PRI = Tp/0.004; % Pulse Repetition Interval for 0,4% Duty Cicle
Lamb = physconst('lightspeed')/fp; % Wavelength
numpulses = 1e1; % Number of pulses to simulate
%% Defining Noise Power Based on Desired SNR
dSNR = [0 5 10 15 20];
No = pow2db(Pt*Gt*Gr*Lamb^2*Sigma./(((4*pi)^3*R^4*db2pow(L))*db2pow(dSNR)));
%% Waveform
hwav = phased.RectangularWaveform('SampleRate',Srate,'PulseWidth',Tp,...
    'PRF',1/PRI,'OutputFormat','Pulses','NumPulses',1);
%% Antenna
hant = phased.IsotropicAntennaElement('FrequencyRange',...
    [1e9 10e9]);
%% Target Model
htgt = phased.RadarTarget('Model','Nonfluctuating',...
    'MeanRCS',Sigma,'PropagationSpeed',physconst('LightSpeed'),...
    'OperatingFrequency',fp);
%% Antenna and Target Plataforms
htxplat = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
htgtplat = phased.Platform('InitialPosition',[R; 0; 0],...
    'Velocity',[0;0;0]);
%% Determining Range and Angle to Target
[tgtrng,tgtang] = rangeangle(htgtplat.InitialPosition,...
    htxplat.InitialPosition);
%% Transmitter
htx = phased.Transmitter('PeakPower',Pt,'Gain',Gt,...
    'LossFactor',L,'InUseOutputPort',true,...
    'CoherentOnTransmit',true);
%% Waveform Radiation and Collection
hrad = phased.Radiator('Sensor',hant,...
    'PropagationSpeed',physconst('LightSpeed'),...
    'OperatingFrequency',fp);
hcol = phased.Collector('Sensor',hant,...
    'PropagationSpeed',physconst('LightSpeed'),...
    'Wavefront','Plane','OperatingFrequency',fp);
%% Propagation
hspace = phased.FreeSpace(...
    'PropagationSpeed',physconst('LightSpeed'),...
    'OperatingFrequency',fp,'TwoWayPropagation',false,...
    'SampleRate',Srate);
%% Receiver
for k=1:length(No)
    hrec = phased.ReceiverPreamp('Gain',Gr,'NoiseMethod','Noise power',...
        'NoisePower',No(k),'SampleRate',Srate,...
        'EnableInputPort',true,'SeedSource','Property','Seed',42);
    % Implementing the Radar Model
    T = 1/hwav.PRF; % Time step between pulses
    txpos = htxplat.InitialPosition; % Get antenna position
    if k==1
        rxsig = zeros(hwav.SampleRate*T,numpulses); % Allocate array for received echoes
        rxsig(:,:,length(No)) = 0;
    end
    for n = 1:numpulses
        % Update the target position
        [tgtpos,tgtvel] = step(htgtplat,T);
        % Get the range and angle to the target
        [tgtrng,tgtang] = rangeangle(tgtpos,txpos);
        % Generate the pulse
        sig = step(hwav);
        % Transmit the pulse. Output transmitter status
        [sig,txstatus] = step(htx,sig);
        % Radiate the pulse toward the target
        sig = step(hrad,sig,tgtang);
        % Propagate the pulse to the target in free space
        sig = step(hspace,sig,txpos,tgtpos,[0;0;0],tgtvel);
        % Reflect the pulse off the target
        sig = step(htgt,sig);
        % Propagate the echo to the antenna in free space
        sig = step(hspace,sig,tgtpos,txpos,tgtvel,[0;0;0]);
        % Collect the echo from the incident angle at the antenna
        sig = step(hcol,sig,tgtang);
        % Receive the echo at the antenna when not transmitting
        rxsig(:,n,k) = step(hrec,sig,~txstatus);
    end
    release(hrec)
    % Noncoherently integrate the received echoes
    rxsig(:,:,k) = pulsint(rxsig(:,:,k),'noncoherent');
    t = unigrid(0,1/hrec.SampleRate,T,'[)');
    rangegates = (physconst('LightSpeed')*t)/2;
    figure(k)
    plot(rangegates,rxsig(:,:,k)); hold on;
    xlabel('Meters'); ylabel('Power');
    ylim = get(gca,'YLim');
    plot([tgtrng,tgtrng],[0 ylim(2)],'r');
end
