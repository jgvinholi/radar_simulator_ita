%% Clearing workspace
clear
close all
clc
%% Defining Parameters
% Parameters for testing
SNR = 20:1:25; 
NTrial = 1e3;
NPulse = 10;
Pfa = 1e-3;
% Given Parameters
Pt = 1e5; % Potência de Transmissão
fc = 3e9; % Frequência de Transmissão
Tau = 2e-6; % Duração do Pulso
L = 6; % Perdas do Radar em dB
Gt = 0; % Ganho de Transmissão em dB
Gr = 0; % Ganho de Recepção em dB
Rt = 6e4; % Base Range for SNR calculation
Sigma = 1; % RCS
% Complementry Parameters
c = 3e8; % Lightspeed aproximation
Lambda = c/fc; % Wavelenght
fs = 2/Tau; % Sample Rate
fpr = 0.004/Tau; % Pulse Repetition Interval for 0,4% Duty Cicle
NoFig = 0; %Noise Figure
% Required Noise Temperature given desired SNRs
No = (Pt*db2pow(Gr+Gt)*Lambda^2*Sigma)./((4*pi)^3*Rt^4*db2pow(L+SNR));
NoTemp = (No)/(physconst('Boltzmann')*db2pow(NoFig)*1/Tau);
% keep memory requirements low
BuffSize = NPulse*1000;
%% Transmition
% Waveform
Waveform = phased.RectangularWaveform('SampleRate',fs,'PulseWidth',Tau,...
    'PRF',fpr,'OutputFormat','Pulses','NumPulses',1);
% Antenna
Antenna = phased.IsotropicAntennaElement('FrequencyRange',...
    [1e9 10e9]);
% Antenna Plataform
AntennaPlatform = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
% Transmitter
Transmitter = phased.Transmitter('PeakPower',Pt,'Gain',Gt,...
    'LossFactor',L,'InUseOutputPort',true,...
    'CoherentOnTransmit',true);
% Radiator
Radiator = phased.Radiator('Sensor',Antenna,...
    'PropagationSpeed',c,...
    'OperatingFrequency',fc);
%% Target
% Target Model
Target = phased.RadarTarget('Model','Nonfluctuating',...
    'MeanRCS',Sigma,'PropagationSpeed',c,...
    'OperatingFrequency',fc);
% Target Plataform
TargetPlataform = phased.Platform('InitialPosition',[Rt; 0; 0],...
    'Velocity',[0;0;0]);
%% Propagation
Channel = phased.FreeSpace('PropagationSpeed',c,...
    'OperatingFrequency',fc,'TwoWayPropagation',true,...
    'SampleRate',fs);
%% Reception
% Waveform Collection
Collector = phased.Collector('Sensor',Antenna,...
    'PropagationSpeed',c,...
    'Wavefront','Plane','OperatingFrequency',fc);
% Receiver
Receiver = phased.ReceiverPreamp('Gain',Gr,'NoiseMethod','Noise temperature',...
        'ReferenceTemperature',NoTemp(1),'SampleRate',fs,...
        'NoiseFigure',NoFig,'EnableInputPort',true);
%% Signal Processing
% Time steps and grid
FastTimeGrid = unigrid(0,1/fs,1/fpr,'[)');
RangeGates = c*FastTimeGrid/2;
% Matched Filter
MatchedFilter = phased.MatchedFilter(...
   'Coefficients',getMatchedFilter(Waveform),...
   'GainOutputPort',true);
% Get group delay of matched filter
GroupDelay = length(MatchedFilter.Coefficients)-1;
%% Implementing the Radar Model
% for fix Tx
AntennaPosition = AntennaPlatform.InitialPosition;
AntennaVelocity = AntennaPlatform.Velocity;
% for fix target
TargetPosition = TargetPlataform.InitialPosition;
TargetVelocity = TargetPlataform.Velocity;
% Range and Angle to Target
[TargetRange,TargetAngle] = rangeangle(TargetPosition,AntennaPosition);
RangeIndex = val2ind(TargetRange,RangeGates(2)-RangeGates(1),RangeGates(1));
% for fix PRF waveform
SigWav = step(Waveform);
% for coherent transmition
[SigTx,TxStatus] = step(Transmitter,SigWav);
% for unchanging radiation angle
SigRad = step(Radiator,SigTx,TargetAngle);
% for unchanging range to target
SigProp = step(Channel,SigRad,AntennaPosition,TargetPosition,...
    AntennaVelocity,TargetVelocity);
% for nonfluctuating target
SigTgt = step(Target,SigProp);
% for unchanging angle to target
SigCol = step(Collector,SigTgt,TargetAngle);
% Allocate arrays
SigRx = zeros(size(RangeGates,2),BuffSize);
SigNonInt = zeros(size(RangeGates,2),BuffSize/NPulse);
SigCoInt = zeros(size(RangeGates,2),BuffSize/NPulse);
SigNcInt = zeros(size(RangeGates,2),BuffSize/NPulse);
SigBi = zeros(size(RangeGates,2),BuffSize/NPulse);
SigRec = zeros(size(RangeGates,2),BuffSize/NPulse);
K = zeros(NPulse,1);
H1Co = zeros(1,NTrial);
H0Co = zeros(1,NTrial);
H1Nc = zeros(1,NTrial);
H0Nc = zeros(1,NTrial);
H1CoInt = zeros(1,NTrial);
H0CoInt = zeros(1,NTrial);
H1NcInt = zeros(1,NTrial);
H0NcInt = zeros(1,NTrial);
H1Bi = zeros(1,NTrial);
H0Bi = zeros(1,NTrial);
H1K = zeros(1,NTrial);
H0K = zeros(1,NTrial);
SimPdCo = zeros(1,length(SNR));
SimPfaCo = zeros(1,length(SNR));
SimPdNc = zeros(1,length(SNR));
SimPfaNc = zeros(1,length(SNR));
SimPdCoInt = zeros(1,length(SNR));
SimPfaCoInt = zeros(1,length(SNR));
SimPdNcInt = zeros(1,length(SNR));
SimPfaNcInt = zeros(1,length(SNR));
SimPdBi = zeros(1,length(SNR));
SimPfaBi = zeros(1,length(SNR));
SimPdK = zeros(1,length(SNR));
SimPfaK = zeros(1,length(SNR));
% buferize
NBuff = floor((NPulse*NTrial)/BuffSize);
NRem = (NPulse*NTrial) - NBuff*BuffSize;
PowCoInt = zeros(NBuff,100);
PowNcInt = zeros(NBuff,100);
for k=1:length(SNR)
    release(Receiver)
    Receiver.ReferenceTemperature = NoTemp(k);
    for m = 1:NBuff
        for n = 1:BuffSize
            % Receive the echo at the antenna when not transmitting
            SigRx(:,n) = step(Receiver,SigCol,~(TxStatus>0));  
        end
        % Apply Matched Filter
        [SigMf, MfGain] = step(MatchedFilter,SigRx);
        % Shift the matched filter output
        SigMf = [SigMf(GroupDelay+1:end,:); SigMf(1:GroupDelay,:)];
        % Defining the Threshold for Noncoherent Detection
        ThresholdNc = (2*No(k)*db2pow(npwgnthresh(Pfa,1,'noncoherent')))*db2pow(MfGain);
        % Condition pulses for integration
        for r = 1:(BuffSize/NPulse)
            SigToInt = SigMf(:,(1:NPulse) + (r-1)*NPulse);
            % Pass on a single pulse
            SigNonInt(:,r) = SigToInt(:,1);
            % Coherentily Integrate NPulses
            SigCoInt(:,r) = pulsint(SigToInt,'coherent');
            % Noncoherentily Integrate NPulses
            Angs = angle(SigToInt) + 2*pi*randn(size(RangeGates,2),NPulse);
            SigToIntNc = abs(SigToInt).*exp(1j*Angs);
            SigNcInt(:,r) = pulsint(SigToIntNc,'noncoherent');
            % Binaraly Integrate NPulses
            SigToIntBi = abs(SigToInt);
            SigToIntBi(SigToIntBi>sqrt(ThresholdNc/4)) = 1;
            SigToIntBi(SigToIntBi<sqrt(ThresholdNc/4)) = 0;
            SigBi(:,r) = sum(SigToIntBi,2);
            % Recursively Integrate NPulses
            SigToIntK = abs(SigToInt);
            SigToIntK(SigToIntK>sqrt(ThresholdNc/4)) = 1;
            SigToIntK(SigToIntK<sqrt(ThresholdNc/4)) = 0;
            for w = 1:NPulse
               K(w,1) = 0.85^(NPulse-w);
            end
            SigRec(:,r) = SigToIntK*K;
        end
        % Extract desired information
        H1Co((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(RangeIndex,:);
        H0Co((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(2,:);
        H1Nc((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(RangeIndex,:);
        H0Nc((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(2,:);
        H1CoInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigCoInt(RangeIndex,:);
        H0CoInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigCoInt(2,:);
        H1NcInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNcInt(RangeIndex,:);
        H0NcInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNcInt(2,:);
        H1Bi((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigBi(RangeIndex,:);
        H0Bi((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigBi(2,:);
        H1K((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigRec(RangeIndex,:);
        H0K((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigRec(2,:);
        if SNR(k) == 0
            for z = 1:10
                SigCoIntVar = zeros(size(RangeGates,2),floor((BuffSize/z)));
                SigNcIntVar = zeros(size(RangeGates,2),floor((BuffSize/z)));
                for r = 1:floor((BuffSize/z))
                SigToIntVar = SigMf(:,(1:z) + (r-1)*z);
                % Coherentily Integrate NPulses
                SigCoIntVar(:,r) = pulsint(SigToIntVar,'coherent');
                % Noncoherentily Integrate NPulses
                Angs = angle(SigToIntVar) + 2*pi*randn(size(RangeGates,2),z);
                SigToIntVarNc = abs(SigToIntVar).*exp(1j*Angs);
                SigNcIntVar(:,r) = pulsint(SigToIntVarNc,'noncoherent');
                end
                EcoCoInt = SigCoIntVar(RangeIndex,:);
                EcoNcInt = SigNcIntVar(RangeIndex,:);
                NoCoInt = [SigCoIntVar(1:RangeIndex-1,:); SigCoIntVar(RangeIndex+1:end,:)];
                NoNcInt = [SigNcIntVar(1:RangeIndex-1,:); SigNcIntVar(RangeIndex+1:end,:)];
                PowCoInt(m,z) = mean(real(EcoCoInt.^2))/mean(var(NoCoInt));
                PowNcInt(m,z) = mean(EcoNcInt.^2)/(mean(mean(NoNcInt.^2,1)));
                AAA(z) = mean(EcoNcInt.^2);
                BBB(z) = mean(real(EcoCoInt)).^2;
                CCC(z) = mean(mean(NoNcInt.^2,1));
                DDD(z) = mean(var(NoCoInt));
            end
            SNRCoInt = pow2db(mean(PowCoInt,1));
            SNRNcInt = pow2db(mean(PowNcInt,1));
        end
    end
    if (NRem > 0)
        parfor n = 1:NRem
            % Receive the echo at the antenna when not transmitting
            SigRx(:,n) = step(Receiver,SigCol,~TxStatus);
        end
        % Apply Matched Filter
        [SigMf, MfGain] = step(MatchedFilter,SigRx);
        % Shift the matched filter output
        SigMf = [SigMf(GroupDelay+1:end,:); SigMf(1:GroupDelay,:)];
        % Defining the Threshold for Noncoherent Detection
        ThresholdNc = (2*No(k)*db2pow(npwgnthresh(Pfa,1,'noncoherent')))*db2pow(MfGain);
        % Condition pulses for integration
        for r = 1:(BuffSize/NPulse)
            SigToInt = SigMf(:,(1:NPulse) + (r-1)*NPulse);
            % Pass on a single pulse
            SigNonInt(:,r) = SigToInt(:,1);
            % Coherentily Integrate NPulses
            SigCoInt(:,r) = pulsint(SigToInt,'coherent');
            % Noncoherentily Integrate NPulses
            Angs = angle(SigToInt) + 2*pi*randn(size(RangeGates,2),NPulse);
            SigToIntNc = abs(SigToInt).*exp(1j*Angs);
            SigNcInt(:,r) = pulsint(SigToIntNc,'noncoherent');
            % Binaraly Integrate NPulses
            SigToIntBi = abs(SigToInt);
            SigToIntBi(SigToIntBi>sqrt(ThresholdNc/4)) = 1;
            SigToIntBi(SigToIntBi<sqrt(ThresholdNc/4)) = 0;
            SigBi(:,r) = sum(SigToIntBi,2);
            % Recursively Integrate NPulses
            SigToIntK = abs(SigToInt);
            SigToIntK(SigToIntK>sqrt(ThresholdNc/4)) = 1;
            SigToIntK(SigToIntK<sqrt(ThresholdNc/4)) = 0;
            for w = 1:NPulse
               K(w,1) = 0.85^(NPulse-w);
            end
            SigRec(:,r) = SigToIntK*K;
        end
        % Extract desired information
        H1Co((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(RangeIndex,:);
        H0Co((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(2,:);
        H1Nc((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(RangeIndex,:);
        H0Nc((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNonInt(2,:);
        H1CoInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigCoInt(RangeIndex,:);
        H0CoInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigCoInt(2,:);
        H1NcInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNcInt(RangeIndex,:);
        H0NcInt((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigNcInt(2,:);
        H1Bi((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigBi(RangeIndex,:);
        H0Bi((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigBi(2,:);
        H1K((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigRec(RangeIndex,:);
        H0K((1:BuffSize/NPulse) + (m-1)*BuffSize/NPulse) = SigRec(2,:);
    end
    %% Compare Simulated and Theoretical Pd and Pfa
    % Coherent Detector
    H1Co = real(H1Co);%*exp(-1j*4*pi*Rt/Lambda));
    H0Co = real(H0Co);%*exp(-1j*4*pi*Rt/Lambda));
    H1CoInt = real(H1CoInt);%*exp(-1j*4*pi*Rt/Lambda));
    H0CoInt = real(H0CoInt);%*exp(-1j*4*pi*Rt/Lambda));
    % Noncoherent Detector
    H1Nc = abs(H1Nc);
    H0Nc = abs(H0Nc);
    H1NcInt = abs(H1NcInt);
    H0NcInt = abs(H0NcInt);
    % Defining the Remaining Threshold
    ThresholdNcInt = (2*No(k)*db2pow(npwgnthresh(Pfa,NPulse,'noncoherent')))...
        *db2pow(MfGain);
    ThresholdCo = (2*No(k)*db2pow(npwgnthresh(Pfa,1,'coherent')))*db2pow(MfGain);
    ThresholdCoInt = (2*No(k)*db2pow(npwgnthresh(Pfa,NPulse,'coherent')))...
        *db2pow(MfGain);
    % Estimate Pd
    CFAR_ON = 0;
    K0 = 2;
    gwndw = 6;
    rwndw = 30;
    if CFAR_ON
        SimPdNc(k) = sum( cfar_filter(H1Nc, K0, gwndw, rwndw) )/(NTrial);
        SimPfaNc(k) = sum(cfar_filter(H0Nc, K0, gwndw, rwndw))/(NTrial);
        SimPdCo(k) = sum( cfar_filter(H1Co, K0, gwndw, rwndw) )/(NTrial);
        SimPfaCo(k) = sum(cfar_filter(H0Co, K0, gwndw, rwndw))/(NTrial);
        SimPdCoInt(k) = sum( cfar_filter(H1CoInt, K0, gwndw, rwndw) )/(NTrial);
        SimPfaCoInt(k) = sum(cfar_filter(H0CoInt, K0, gwndw, rwndw))/(NTrial);
        SimPdNcInt(k) = sum(cfar_filter(H1NcInt, K0, gwndw, rwndw) )/(NTrial);
        SimPfaNcInt(k) = sum( cfar_filter(H0NcInt, K0, gwndw, rwndw) )/(NTrial);
%         SimPdBi(k) = sum(H1Bi >= 0.955*NPulse^(0.8))/(NTrial);
%         SimPfaBi(k) = sum(H0Bi >= 0.955*NPulse^(0.8))/(NTrial);
%         SimPdK(k) = sum(H1K >= 0.955*(sum(K))^(0.8))/(NTrial);
%         SimPfaK(k) = sum(H0K >= 0.955*(sum(K))^(0.8))/(NTrial);
    else
        SimPdNc(k) = sum(H1Nc >= sqrt(ThresholdNc))/(NTrial);
        SimPfaNc(k) = sum(H0Nc >= sqrt(ThresholdNc))/(NTrial);
        SimPdCo(k) = sum(H1Co >= sqrt(ThresholdCo))/(NTrial);
        SimPfaCo(k) = sum(H0Co >= sqrt(ThresholdCo))/(NTrial);
        SimPdCoInt(k) = sum(H1CoInt >= sqrt(ThresholdCoInt))/(NTrial);
        SimPfaCoInt(k) = sum(H0CoInt >= sqrt(ThresholdCoInt))/(NTrial);
        SimPdNcInt(k) = sum(H1NcInt >= sqrt(ThresholdNcInt))/(NTrial);
        SimPfaNcInt(k) = sum(H0NcInt >= sqrt(ThresholdNcInt))/(NTrial);
        SimPdBi(k) = sum(H1Bi >= 0.955*NPulse^(0.8))/(NTrial);
        SimPfaBi(k) = sum(H0Bi >= 0.955*NPulse^(0.8))/(NTrial);
        SimPdK(k) = sum(H1K >= 0.955*(sum(K))^(0.8))/(NTrial);
        SimPfaK(k) = sum(H0K >= 0.955*(sum(K))^(0.8))/(NTrial);
    end
    
end
% Plot the ROC curve
figure(1), hold on, grid on
plot(SNR,SimPdCo)
plot(SNR,SimPdCoInt)
plot(SNR,SimPdNc)
plot(SNR,SimPdNcInt)
plot(SNR,SimPdBi)
plot(SNR,SimPdK)
title(['ROC Curves for ' num2str(Pfa) ' P_{fa}'])
xlabel('SNR (dB)')
ylabel('P_{d}')
legend('Single Pulse Coherent',[num2str(NPulse) ' Pulse Coherent'],...
    'Single Pulse Noncoherent',[num2str(NPulse) ' Pulse Noncoherent'],...
    [num2str(NPulse) ' Pulse Binary'],[num2str(NPulse) ' Pulse Recursive'],...
    'Location','southeast');
% Plot SNR gain x number of pulses integrated
% figure(2), hold on, grid on
% plot(1:100,SNRCoInt)
% plot(1:100,SNRNcInt)
% title('SNR gain per number of integratade pulses')
% xlabel('Number of pulses')
% ylabel('SNR gain')
% legend('Coherent Integration','Noncoherent Integration','Location','southeast')

