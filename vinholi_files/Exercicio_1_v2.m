%% Clearing workspace
clear
close all
clc
%% Defining Given Parameters
Pt = 1e5; % Pot�ncia de Transmiss�o
fp = 3e9; % Frequ�ncia de Transmiss�o
Tp = 2e-6; % Dura��o do Pulso
R = 6e4; % Dist�ncia do Alvo ao Radar
L = 6; % Perdas do Radar em dB
Gt = 1; Gr = Gt; % (Ganho de Transmiss�o e Recep��o)
Sigma = 1; % (Radar Cross Section)
%% Defining Missing Parameters
Srate = 1e6; % Sample Rate
PRI = Tp/0.004; % Pulse Repetition Interval for 0,4% Duty Cicle
Lamb = physconst('lightspeed')/fp; % Wavelength
numpulses = 1e3; % Number of pulses to simulate
%% Defining Noise Power Based on Desired SNR
dSNR = [40];
No = Pt*Gt*Gr*Lamb^2*Sigma./(((4*pi)^3*R^4*db2pow(L))*db2pow(dSNR));
%% Waveform
hwav = phased.RectangularWaveform('SampleRate',Srate,'PulseWidth',Tp,...
    'PRF',1/PRI,'OutputFormat','Pulses','NumPulses',1);
%% Antenna
% hant = phased.IsotropicAntennaElement('FrequencyRange',...
%     [1e9 10e9]);
hant = phased.CrossedDipoleAntennaElement('FrequencyRange',[1/10,2500]*fp);
%% Target Model
htgt = phased.RadarTarget('Model','Nonfluctuating',...
    'MeanRCS',Sigma,'PropagationSpeed',physconst('LightSpeed'),...
    'OperatingFrequency',fp);
%% Antenna and Target Plataforms
htxplat = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
htgtplat = phased.Platform('InitialPosition',[0; R; 0],...
    'Velocity',[0;1e4;0]);
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
        'EnableInputPort',true);
    %% Implementing the Radar Model
    T = 1/hwav.PRF; % Time step between pulses
    txpos = htxplat.InitialPosition; % Get antenna position
    if k==1
        rxsig = zeros(hwav.SampleRate*T,numpulses); % Allocate array for received echoes
        rxsig(:,:,length(No)) = 0;
    end
    range = zeros(1,numpulses);
    angle = zeros(2,numpulses);
    for n = 1:numpulses
        % Update the target position
        [tgtpos,tgtvel] = step(htgtplat,T);
        % Get the range and angle to the target
        [tgtrng,tgtang] = rangeangle(tgtpos,txpos);
        range(n) = tgtrng;
        angle(:, n) = tgtang;
        
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
    %rxsig(:,k) = pulsint(rxsig(:,:,k),'noncoherent');
    t = unigrid(0,1/hrec.SampleRate,T,'[)');
    rangegates = (physconst('LightSpeed')*t)/2;
    figure(k)
    xlabel('Meters'); ylabel('Power');
    ylim = get(gca,'YLim');
    hold on;
    
    for i=1:100:numpulses
        plot(rangegates,abs(rxsig(:,i,k)));
    end
    
    fig_phase = figure;
    figure(fig_phase)
    hold on;
    for i=1:100:numpulses
        plot(rangegates,phase(rxsig(:,i,k)));
    end
    
%     for i=1:floor(length(range)/10):length(range)
%       range(i)
%       plot ( [range(i), range(i)],[0, ylim(2)],'r');
%     end
    
%     %% Determining P_d
%     if k ==1
%         x = zeros(length(No),numpulses);
%     end
%     x(k,:) = abs(rxsig(402,:,k));
%     Th = linspace(0.5*sqrt(No(1,k)),4*sqrt(No(1,k)),10);
%     if k==1
%         Pd = zeros(length(No),length(Th));
%     end
%     for r = 1:length(Th)
%         Pd(k,r) = sum(x(k,:)>Th(1,r))/numpulses;
%     end
%     %% Determining P_fa
%     x(k,:) = abs(rxsig(2,:,k));
%     Th = linspace(0.5*sqrt(No(1,k)),4*sqrt(No(1,k)),10);
%     if k==1
%         Pfa = zeros(length(No),length(Th));
%     end
%     for r = 1:length(Th)
%         Pfa(k,r) = sum(x(k,:)>Th(1,r))/numpulses;
%     end
end
% %% Ploting ROC
% figure(k+1)
% Labels = cell(1,length(No));
% for i = 1:length(No)
%     semilogx(Pfa(i,:),Pd(i,:)), hold on
%     Labels{1,i} = strcat(num2str(dSNR(i)),' SNR');
% end
% title('ROC'), xlabel('P_{fa}'), ylabel('P_d')
% xlim([1e-4 0.5])
% legend(Labels)