classdef RadarSimulator < handle

    properties
        % Default Parameters
        antenna_rpm = 12; % rpm
        beam_aperture = 3; % degrees
        c = 3e8;
        dcycle = 0.004;
        dSNR = 20;
        dvdt = 0; % Aceleração.
        fp = 1e8; % Frequência de Transmissão
        Gt = 1;
        Gr = 1; % Linear !!!!!!!(Ganho de Transmissão e Recepção)
        L = 6; % Perdas do Radar em dB
        Nofig = 0; % Noise Figure in dB
        Pt = 1e5; % Potência de Transmissão
        R = 600; % Distância do Alvo ao Radar
        Sigma = 1; % (Radar Cross Section)
        SimDuration = 2; % Simulation duration in seconds.
        Tp = 1e-6; % Duração do Pulso
        V_0 = 100; % Velocidade inicial do alvo 100 m/s.
        acceleration = 30; % A = 30m/s^2
        npulses_per_prediction = 100; % Number of pulses to receive before prediction. 
        rangegates
        range_res
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
        max_range
        max_speed
        No
        Notemp
        npulse
        n_predictions
        PRF
        PRI
        rxsig
        sigrad
        sigtx
        sigwav
        speed_res
        tgtpos
        tgtvel
        tgtang
        txpos
        txstatus
        txvel
    end
    methods
        function obj = RadarSimulator % radar_params is a struct with parameters
            update_parameters(obj)
        end
        function obj = update_parameters(obj)
            
            obj.fs = 5*1/obj.Tp; % Sample Rate 10x a frequência do pulsador.
            obj.PRI = obj.Tp/obj.dcycle; % Pulse Repetition Interval for 0,5 % Duty Cicle
            obj.PRF = 1/obj.PRI;
            obj.lamb = obj.c/obj.fp; % Wavelength
            
            obj.max_speed = dop2speed(obj.PRF/2,obj.lamb)/2;
            obj.speed_res = 2*obj.max_speed/obj.npulses_per_prediction;
            obj.range_res = obj.c*(1/obj.fs)/2;
            obj.max_range = obj.c*obj.PRI/2;
            % Defining Noise Power Based on Desired SNR
            obj.No = obj.Pt.*obj.Gt.*obj.Gr.*obj.lamb.^2.*obj.Sigma./(((4.*pi).^3.*obj.R.^4.*db2pow(obj.L)).*db2pow(obj.dSNR));
            obj.Notemp = obj.No./(physconst('Boltzmann').*1/obj.Tp.*db2pow(obj.Nofig));

            %% Transmition
            % Waveform
            obj.hwav = phased.RectangularWaveform('SampleRate',obj.fs,'PulseWidth',obj.Tp,...
                'PRF',1/obj.PRI,'OutputFormat','Pulses','NumPulses',1);
            % Antenna
            obj.hant = phased.IsotropicAntennaElement('FrequencyRange',...
                [obj.fp/10 obj.fp*10]);
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
                'InitialVelocity',[obj.V_0;0;0], 'Acceleration', [obj.acceleration; 0; 0], 'MotionModel', 'Acceleration');
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
            Tgrid = unigrid(0,1/obj.fs,obj.PRI, '[)');
            obj.rangegates = obj.c*Tgrid/2;
            % for fix Tx
            obj.txpos = obj.htxplat.InitialPosition;
            obj.txvel = obj.htxplat.Velocity;
            obj.tgtpos = obj.htgtplat.InitialPosition;
            obj.tgtvel = obj.htgtplat.InitialVelocity;
            [~, obj.tgtang] = rangeangle(obj.tgtpos, obj.txpos);
%             obj.rangeidx = val2ind(tgtrng,rangegates(2)-rangegates(1),rangegates(1))-1;
            obj.rxsig = zeros( round(obj.fs*obj.PRI), obj.npulses_per_prediction, obj.n_predictions); 
            obj.hrec.ReferenceTemperature = obj.Notemp;
            obj.sigwav = step(obj.hwav);
            % for coherent transmition
            [obj.sigtx, obj.txstatus] = step(obj.htx, obj.sigwav);
        end
        
        
        function result_struct = runSim(obj, SimDuration)
            pred_range = cell(obj.n_predictions, 1);
            pred_speed = pred_range;
            real_range = pred_range;
            real_speed = pred_range;
            error_range = pred_range;
            error_speed = pred_speed;
            obj.npulse = ceil(SimDuration/obj.PRI);
%             obj.Npulse_per_rev = floor( obj.beam_aperture*60/(360*obj.antenna_rpm*obj.PRI) ); % Number of pulses per revolution. Rotating Radar.
            obj.n_predictions = ceil( obj.npulse/obj.npulses_per_prediction); % Number of predictions to make, based on the number of pulses to use for each prediction and the number of total pulses to send.
            for m = 1:obj.n_predictions
                for n = 1:obj.npulses_per_prediction
                    % Update velocity and position.
                    % Update radiated signal with the angle
                    obj.sigrad = step(obj.hrad, obj.sigtx, obj.tgtang);
                    % Reflect wave on target.
                    sigbounce = step(obj.hspace, obj.sigrad, obj.txpos, obj.tgtpos, obj.txvel, obj.tgtvel);
                    % for nonfluctuating target
                    sigtgt = step(obj.htgt, sigbounce);
                    % angle to target
                    sigcol = step(obj.hcol, sigtgt, obj.tgtang);
                    obj.rxsig(:,n, m) = step(obj.hrec, sigcol, ~obj.txstatus);
                    [obj.tgtpos, obj.tgtvel] = step(obj.htgtplat, obj.PRI);
                    [~,obj.tgtang] = rangeangle(obj.tgtpos, obj.htxplat.InitialPosition);
                end
                
                real_range{m} = norm(obj.tgtpos);
                real_speed{m} = norm(obj.tgtvel);
                % Apply Matched Filter
                %             mfsig = real(exp(1i*4*pi/Lamb*R)*step(hmf,rxsig)); % coherent
                mfsig = step(obj.hmf,obj.rxsig(:,:,m) );

                % Shift the matched filter output
                mfsig=[mfsig(obj.gd+1:end, :); mfsig(1:obj.gd, :)];
                %  h1(1,(1:Nrem) + Nbuff*Npulsebuffsize) = mfsig(rangeidx,1:Nrem);
                %     h0(1,(1:Nrem) + Nbuff*Npulsebuffsize) = mfsig(rangeidx,1:Nrem);


                % Noncoherently integrate the received echoes
                %rxsig(:,k) = pulsint(rxsig(:,:,k),'noncoherent');
                % Plot rangegates power
%                 figure
%                 plot(obj.rangegates, abs( mfsig(:,1)) ); hold on;
                %             plot(rangegates,( pulsint(mfsig(:,1:10),'noncoherent') ) ); hold on;
%                 xlabel('Meters'); 
%                 ylabel('Amplitude');
%                 title('Matched Signal Amplitude per Distance')
%                 plot([obj.R,obj.R],[0 max(abs(mfsig(:,1)))],'r--');

% Find peaks based on threshold:
%                 [~, range_estimates_index] = findpeaks( abs(mfsig(:,1)), 'MinPeakHeight', max(abs( mfsig(:,1) ) )*0.8 );
% Find the max value (assuming only one target is present).
                [~, range_estimates_index] = max( pulsint(mfsig, 'noncoherent' ) );
                pred_range{m} = obj.rangegates(range_estimates_index);
                error_range{m} = pred_range{m} - real_range{m};
                %% Doppler Estimation
                numel(pred_range{m});
                speed_target = {} ;
                for i=1:numel(pred_range{m})
                    [Pxx, Fx] = periodogram((mfsig(range_estimates_index(i),:)).',[],512,obj.PRF, 'power','centered'); % Calculate power spectrum of rxsig along pulses.
                    speed_vec = dop2speed(Fx,obj.lamb)/2; % Translate to speed each frequency detected.
                    Pxx_norm = Pxx/max(Pxx); % Normalize Sxx to have its peak at 0 dB.
%                     [~,detected_index] = findpeaks(pow2db(Pxx_norm),'MinPeakHeight',-5); % Select peaks whose power is at least 5dB below the max power.
                    [~,detected_index] = max(pow2db(Pxx_norm) );
                    prediction_sp = speed_vec(detected_index); % Speed whose sign follows the doppler frequency sign. That is, if the target is getting away from the platform, the sign is negative. 
                    speed_target{end+1} = -prediction_sp;
                end
                pred_speed{m} = cell2mat(speed_target); % The minus sign corrects the sign of the prediction.
                error_speed{m} = pred_speed{m} - real_speed{m};
            end
            mse_range = sqrt( mean( vertcat(error_range{:}).^2 ) );
            mse_speed = sqrt( mean( vertcat(error_speed{:}).^2 ) );
%             result_struct = struct( 'pred_range', pred_range, 'pred_speed', pred_speed, ...
%                 'real_range', real_range, 'real_speed', real_speed, 'error_range', error_range, ...
%                 'error_speed', error_speed, 'mse_range', mse_range, 'mse_speed', mse_speed);
            result_struct.pred_range = pred_range;
            result_struct.pred_speed = pred_speed;
            result_struct.real_range = real_range;
            result_struct.real_speed = real_speed;
            result_struct.error_range = error_range;
            result_struct.error_speed = error_speed;
            result_struct.mse_range = mse_range;
            result_struct.mse_speed = mse_speed;
        end
        
    end
    methods(Static)
       function result = flatten_result(result)
            result.pred_range = vertcat(result.pred_range{:});
            result.pred_speed = vertcat(result.pred_speed{:});
            result.real_range = vertcat(result.real_range{:});
            result.real_speed = vertcat(result.real_speed{:});
            result.error_range = vertcat(result.error_range{:});
            result.error_speed = vertcat(result.error_speed{:});
        end

        function plot_range_speed_pred(result_struct, radobj)
            linewidth = 1.2;
            result = RadarSimulator.flatten_result(result_struct);
            time_vec = linspace(0, radobj.SimDuration, length(result.pred_range));
            SNR = radobj.dSNR;
            %% Plot range estimate and compare with real value.
                figure
                title([ 'Range Estimation with SNR = ', num2str(SNR), ' dB'] )
                xlabel('Time [s]')
                ylabel('Range [m]')
                hold on
                plot(time_vec, result.pred_range, 'LineWidth',linewidth);
                plot(time_vec, result.real_range, '--', 'LineWidth', linewidth);
                legend({'Prediction','Real'},'Location','southeast')
                hold off

                %% Plot speed estimate and comapre with real value.
                figure
                title([ 'Speed Estimation with SNR = ', num2str(SNR), ' dB'] )
                xlabel('Time [s]')
                ylabel('Speed [m/s]')
                hold on
                plot(time_vec, result.pred_speed, 'LineWidth', linewidth);
                plot(time_vec, result.real_speed, '--', 'LineWidth', linewidth);
                legend({'Prediction','Real'},'Location','southeast')
                hold off
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
