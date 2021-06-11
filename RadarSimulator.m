classdef RadarSimulator < handle

    properties
        % Default Parameters
        antenna_rpm = 12; % rpm
        beam_aperture = 3; % degrees
        c = 3e8;
        dcycle = 0.004;
        dSNR = 20;
        dvdt = 0; % Aceleração.
        fp = 3e9; % Frequência de Transmissão
        Gt = 1;
        Gr = 1; % Linear !!!!!!!(Ganho de Transmissão e Recepção)
        L = 6; % Perdas do Radar em dB
        Nofig = 0; % Noise Figure in dB
        Pt = 1e5; % Potência de Transmissão
        R = 6e4; % Distância do Alvo ao Radar
        Sigma = 1; % (Radar Cross Section)
        SimDuration = 1; % Simulation duration in seconds.
        Tp = 2e-6; % Duração do Pulso
        V_0 = 20; % Velocidade do alvo 20 m/s.
        buffer_length = 1000; % Number of pulses to receive before prediction. 
    end
    
    properties (SetAccess = private)
        fs
        gd
        hant
        hcol
        hmf
        hrad
        hrec
        hspace
        htgt
        htgtplat
        htx
        htxplat
        hwav
        lamb
        No
        Notemp
        npulse
        number_of_predictions
        PRF
        PRI
        rangegates
        rxsig
        sigrad
        sigtx
        sigwav
        txpos
        txstatus
        txvel
    end
    methods
        function obj = RadarSimulator % radar_params is a struct with parameters
            update_parameters(obj)
        end
        function obj = update_parameters(obj)
            
            obj.fs = 2/obj.Tp; % Sample Rate
            obj.PRI = obj.Tp/obj.dcycle; % Pulse Repetition Interval for 0,5 % Duty Cicle
            obj.PRF = 1/obj.PRI;
            obj.lamb = obj.c/obj.fp; % Wavelength
            obj.npulse = ceil(obj.SimDuration/obj.PRI)
%             obj.Npulse_per_rev = floor( obj.beam_aperture*60/(360*obj.antenna_rpm*obj.PRI) ); % Number of pulses per revolution. Rotating Radar.
            obj.number_of_predictions = ceil( obj.buffer_length/obj.npulse );
            % keep memory requirements low
            % Npulse = 1000;
            % Defining Noise Power Based on Desired SNR
            obj.No = obj.Pt.*obj.Gt.*obj.Gr.*obj.lamb.^2.*obj.Sigma./(((4.*pi).^3.*obj.R.^4.*db2pow(obj.L)).*db2pow(obj.dSNR));
            obj.Notemp = obj.No./(physconst('Boltzmann').*1/obj.Tp.*db2pow(obj.Nofig));

            %% Transmition
            % Waveform
            obj.hwav = phased.RectangularWaveform('SampleRate',obj.fs,'PulseWidth',obj.Tp,...
                'PRF',1/obj.PRI,'OutputFormat','Pulses','NumPulses',1);
            % Antenna
            obj.hant = phased.IsotropicAntennaElement('FrequencyRange',...
                [1e9 10e9]);
            % Antenna Plataform
            obj.htxplat = phased.Platform('InitialPosition',[0;0;0],'Velocity',[0;0;0]);
            % Transmitter
            obj.htx = phased.Transmitter('PeakPower',obj.Pt,'Gain',pow2db(obj.Gt),...
                'LossFactor',obj.L,'InUseOutputPort',true,...
                'CoherentOnTransmit',true);
            % Radiator
            obj.hrad = phased.Radiator('Sensor',obj.hant,...
                'PropagationSpeed',obj.c,...
                'OperatingFrequency',obj.fp);
            %% Target
            % Target Model
            obj.htgt = phased.RadarTarget('Model','Nonfluctuating',...
                'MeanRCS',obj.Sigma,'PropagationSpeed',obj.c,...
                'OperatingFrequency',obj.fp);
            % Target Plataform
            obj.htgtplat = phased.Platform('InitialPosition',[obj.R; 0; 0],...
                'Velocity',[obj.V_0;0;0]);
            %% Propagation
            % Channel 1 (target)
            obj.hspace = phased.FreeSpace('PropagationSpeed',obj.c,...
                'OperatingFrequency',obj.fp,'TwoWayPropagation',true,...
                'SampleRate',obj.fs);
            %% Reception
            % Waveform Collection
            obj.hcol = phased.Collector('Sensor',obj.hant,...
                'PropagationSpeed',obj.c,...
                'Wavefront','Plane','OperatingFrequency',obj.fp);
            % Receiver
            obj.hrec = phased.ReceiverPreamp('Gain',pow2db(obj.Gr),'NoiseMethod','Noise temperature',...
                    'ReferenceTemperature',obj.Notemp,'SampleRate',obj.fs,...
                    'NoiseFigure',obj.Nofig,'EnableInputPort',true);
            % Matched Filter
            obj.hmf = phased.MatchedFilter(...
               'Coefficients',getMatchedFilter(obj.hwav),...
               'GainOutputPort',false);
            % Get group delay of matched filter
            obj.gd = length(obj.hmf.Coefficients)-1;

            %% Implementing the Radar Model
            % Time steps and grid
            Tstp = obj.PRI;
            Tgrid = unigrid(0,1/obj.fs,Tstp, '[)');
            obj.rangegates = obj.c*Tgrid/2;
            % for fix Tx
            obj.txpos = obj.htxplat.InitialPosition;
            obj.txvel = obj.htxplat.Velocity;
%             obj.rangeidx = val2ind(tgtrng,rangegates(2)-rangegates(1),rangegates(1))-1;
            obj.rxsig = zeros( round(obj.fs*Tstp), obj.npulse); 
            obj.hrec.ReferenceTemperature = obj.Notemp;
            obj.sigwav = step(obj.hwav);
            % for coherent transmition
            [obj.sigtx, obj.txstatus] = step(obj.htx, obj.sigwav);
        end
        function runSim(obj)
            for n = 1:obj.npulse
                % Update velocity and position.
                [tgtpos, tgtvel] = step(obj.htgtplat, obj.PRI);
                [tgtrng,tgtang] = rangeangle(tgtpos, obj.htxplat.InitialPosition);
                % Update radiated signal with the angle
                obj.sigrad = step(obj.hrad, obj.sigtx, tgtang);
                % Reflect wave on target.
                sigbounce = step(obj.hspace, obj.sigrad, obj.txpos, tgtpos, obj.txvel, tgtvel);
                % for nonfluctuating target
                sigtgt = step(obj.htgt, sigbounce);
                % angle to target
                sigcol = step(obj.hcol, sigtgt, tgtang);

                obj.rxsig(:,n) = step(obj.hrec, sigcol, ~obj.txstatus);
            end
            % Apply Matched Filter
            %             mfsig = real(exp(1i*4*pi/Lamb*R)*step(hmf,rxsig)); % coherent
            mfsig = step(obj.hmf,obj.rxsig);

            % Shift the matched filter output
            mfsig=[mfsig(obj.gd+1:end, :); mfsig(1:obj.gd, :)];
            %  h1(1,(1:Nrem) + Nbuff*Npulsebuffsize) = mfsig(rangeidx,1:Nrem);
            %     h0(1,(1:Nrem) + Nbuff*Npulsebuffsize) = mfsig(rangeidx,1:Nrem);


            % Noncoherently integrate the received echoes
            %rxsig(:,k) = pulsint(rxsig(:,:,k),'noncoherent');
            % Plot rangegates power
            figure
            plot(obj.rangegates, abs( mfsig(:,1)) ); hold on;
            %             plot(rangegates,( pulsint(mfsig(:,1:10),'noncoherent') ) ); hold on;
            xlabel('Meters'); 
            ylabel('Amplitude');
            title('Matched Signal Amplitude per Distance')
            plot([obj.R,obj.R],[0 max(abs(mfsig(:,1)))],'r--');
            [~, range_estimates_index] = findpeaks( abs(mfsig(:,1)), 'MinPeakHeight', max(abs( mfsig(:,1) ) )*0.8 );
            range_estimates = obj.rangegates(range_estimates_index);

            %% Doppler Estimation
            max_speed = dop2speed(obj.PRF/2,obj.lamb)/2;
            speed_res = 2*max_speed/obj.npulse;
            for i=1:numel(range_estimates)
                [Pxx, Fx] = periodogram((mfsig(range_estimates_index(i),:)).',[],1024,obj.PRF, 'power','centered'); % Calculate power spectrum of rxsig along pulses.
                speed_vec = dop2speed(Fx,obj.lamb)/2; % Translate to speed each frequency detected.
                Pxx_norm = Pxx/max(Pxx); % Normalize Sxx to have its peak at 0 dB.
                [~,detected_index] = findpeaks(pow2db(Pxx_norm),'MinPeakHeight',-5); % Select peaks whose power is at least 5dB below the max power.
                speed_target = speed_vec(detected_index);
                speed_target
            end
        end
    end

end

% %% Create Histogram of Outputs
% h1a(k,:) = abs(h1(1,:));
% h0a(k,:) = abs(h0(1,:));
% thresh_low = min([h1a(k,:), h0a(k,:)]);
% thresh_hi  = max([h1a(k,:), h0a(k,:)]);
% nbins = 100;
% binedges = linspace(thresh_low,thresh_hi,nbins);
% figure(length(No)+k)
% histogram(h0a(k,:),binedges)
% hold on
% histogram(h1a(k,:),binedges)
% hold off
% title('Target-Absent Vs Target-Present Histograms')
% xlabel('Signal Amplitude'); ylabel('Number of Signals');
% legend('Target Absent','Target Present');
% %% Compare Simulated and Theoretical Pd and Pfa
% nbins = 1000;
% thresh_steps = linspace(thresh_low,thresh_hi,nbins);
% sim_pd = zeros(1,nbins);
% sim_pfa = zeros(1,nbins);
% for r = 1:nbins
%     thresh = thresh_steps(r);
%     sim_pd(r) = sum(h1a(k,:) >= thresh);
%     sim_pfa(r) = sum(h0a(k,:) >= thresh);
% end
% sim_pd = sim_pd/(Npulse);
% sim_pfa = sim_pfa/(Npulse);
% pfa_diff = diff(sim_pfa);
% idx = (pfa_diff == 0);
% sim_pfa(idx) = [];
% sim_pd(idx) = [];
% minpfa = 1e-6;
% N = sum(sim_pfa >= minpfa);
% sim_pfa = fliplr(sim_pfa(1:N)).';
% sim_pd = fliplr(sim_pd(1:N)).';
% % plot the theoretical Pfa curve
% figure(2*length(No)+k)
% rocsnr(dSNR(k),'SignalType',...
%     'NonfluctuatingNoncoherent',...
%     'MinPfa',minpfa,'NumPoints',N,'NumPulses',1);
% hold on
% semilogx(sim_pfa,sim_pd,'r.')
% title('Simulated and Theoretical ROC Curves')
% xlabel('Pfa')
% ylabel('Pd'), ylim([-inf 1])
% grid on
% legend('Theoretical','Simulated','Location','SE');
