classdef RadarSimulator < handle
% Parameters for testing
% SNR = -15:1:15; 
% NTrial = 100 ;
% NPulse = 10;
% Pfa = 1e-3;
% % Given Parameters
% Pt = 1e5; % Pot?ncia de Transmiss?o
% fc = 3e9; % Frequ?ncia de Transmiss?o
% Tau = 2e-6; % Dura??o do Pulso
% L = 6; % Perdas do Radar em dB
% Gt = 0; % Ganho de Transmiss?o em dB
% Gr = 0; % Ganho de Recep??o em dB
% Rt = 6e4; % Base Range for SNR calculation
% Sigma = 1; % RCS
% % Complementry Parameters
% c = 3e8; % Lightspeed aproximation
% Lambda = c/fc; % Wavelenght
% fs = 2/Tau; % Sample Rate
% fpr = 0.004/Tau; % Pulse Repetition Interval for 0,4% Duty Cicle
% NoFig = 0; %Noise Figure
% % Required Noise Temperature given desired SNRs
% No = (Pt*db2pow(Gr+Gt)*Lambda^2*Sigma)./((4*pi)^3*Rt^4*db2pow(L+SNR));
% NoTemp = (No)/(physconst('Boltzmann')*db2pow(NoFig)*1/Tau);
% % keep memory requirements low
% BuffSize = NPulse*1000;
    properties
        % Default Parameters
        antenna_rpm = 12; % rpm
        beam_aperture = 3; % degrees
        c = 3e8;
        dcycle = 0.001;
        dSNR = 0;
        dSCR = 10;
        dvdt = 0; % Acelera??o.
        fp = 3e9; % Frequ?ncia de Transmiss?o
        Gt = 1;
        Gr = 1; % Linear !!!!!!!(Ganho de Transmiss?o e Recep??o)
        L = 6; % Perdas do Radar em dB
        Nofig = 0; % Noise Figure in dB
        Pt = 1; % Pot?ncia de Transmiss?o
        Pr = 1;
        R = 60000; % Dist?ncia do Alvo ao Radar
        Sigma = 1; % (Radar Cross Section)
        SimDuration = 2; % Simulation duration in seconds.
        Tp = 2e-6; % Dura??o do Pulso
        V_0 = 0; % Velocidade inicial do alvo 100 m/s.
        acceleration = 0; % A = 30m/s^2
        npulses_per_prediction = 10; % Number of pulses to receive before prediction. 
        Integrator = 'noncoherent';
        threshold = 0.95;
        binary_integ_percent_detected = 0.499; % If target is detected in more than x% of the pulses, it will be assumed to be a detection.
        rangegates
        range_res
        PFA = 1e-3;
        No
        Co
        CFAR_ON = 1;
        CLUTTER_ON = 1;
        PRED_VEL = 0;
        K0 = 2;
        rwndw = 40;
        gwndw = 10;
        cluttype = 'rayleigh';
        thresh
        signal_cfar
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
            
            obj.fs = 2*1/obj.Tp; % Sample Rate 10x a frequ?ncia do pulsador.
            obj.PRI = obj.Tp/obj.dcycle; % Pulse Repetition Interval for 0,5 % Duty Cicle
            obj.PRF = 1/obj.PRI;
            obj.lamb = obj.c/obj.fp; % Wavelength
            
            obj.max_speed = dop2speed(obj.PRF/2,obj.lamb)/2;
            obj.speed_res = 2*obj.max_speed/obj.npulses_per_prediction;
            obj.range_res = obj.c*obj.Tp/2;
            obj.max_range = obj.c*obj.PRI/2;
            % Defining Noise Power Based on Desired SNR
%             obj.Pr = obj.Pt*obj.Gt.*obj.Gr.*obj.lamb.^2.*obj.Sigma./( (4.*pi).^3.*obj.R.^4.*db2pow(obj.L) );
            obj.No = obj.Pr/db2pow(obj.dSNR);
            obj.Co = obj.Pr/db2pow(obj.dSCR);
%             obj.Notemp = obj.No./(physconst('Boltzmann').*1/obj.Tp.*db2pow(obj.Nofig));
%             obj.No = obj.Pt/10^(obj.dSNR/10);
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
%             obj.hrec = phased.ReceiverPreamp('Gain',pow2db(obj.Gr),'NoiseMethod','Noise temperature',...
%                     'ReferenceTemperature',obj.Notemp,'SampleRate',obj.fs,...
%                     'NoiseFigure',obj.Nofig,'EnableInputPort',true);
            obj.hrec = phased.ReceiverPreamp('Gain',pow2db(obj.Gr),'NoiseMethod','Noise power',...
            'SampleRate',obj.fs ,'EnableInputPort',true, 'NoisePower', obj.No);
            % Matched Filter
            obj.hmf = phased.MatchedFilter(...
               'Coefficients',getMatchedFilter(obj.hwav),...
               'GainOutputPort',true);
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
%             obj.hrec.ReferenceTemperature = obj.Notemp;
            obj.sigwav = step(obj.hwav);
            [obj.sigtx, obj.txstatus] = step(obj.htx, obj.sigwav);
            
%             npower = noisepow(noise_bw,receiver.NoiseFigure,...
%             receiver.ReferenceTemperature);
            
        end
        
        
        function result_struct = runSim(obj, SimDuration)
            pred_range = cell(obj.n_predictions, 1);
            obj.thresh = cell(obj.n_predictions, 1);
            obj.signal_cfar = cell(obj.n_predictions, 1);
            pred_speed = pred_range;
            real_range = pred_range;
            real_speed = pred_range;
            error_range = pred_range;
            error_speed = pred_speed;
            obj.npulse = ceil(SimDuration/obj.PRI);
%             obj.Npulse_per_rev = floor( obj.beam_aperture*60/(360*obj.antenna_rpm*obj.PRI) ); % Number of pulses per revolution. Rotating Radar.
            obj.n_predictions = ceil( obj.npulse/obj.npulses_per_prediction); % Number of predictions to make, based on the number of pulses to use for each prediction and the number of total pulses to send.
            obj.rxsig = zeros( round(obj.fs*obj.PRI), obj.npulses_per_prediction, obj.n_predictions); 
            for m = 1:obj.n_predictions
                obj.thresh{m} = zeros(size( obj.rxsig(:,1,1) ));
                obj.signal_cfar{m} =  zeros(size( obj.rxsig(:,1,1) ));
%                 if mod(m, 10) == 1
                    clut_real = repmat( obj.clutter_gen( size( obj.rxsig(:,1,1) ) , obj.Co/2, 'rayleigh'), 1, obj.npulses_per_prediction);
                    clut_imag = repmat( obj.clutter_gen( size( obj.rxsig(:,1,1) ) , obj.Co/2, 'rayleigh'), 1, obj.npulses_per_prediction);
%                 end
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
                    sigcol( sigcol ~= 0) = sqrt(obj.Pr)*sigcol( sigcol ~= 0)./abs(sigcol( sigcol ~= 0));
%                     obj.rxsig(:,n, m) = step(obj.hrec, sigcol, ~obj.txstatus);
                    obj.rxsig(:,n,m) = (sigcol + sqrt(obj.No/2)*( randn( size(sigcol) ) + 1i*randn( size(sigcol) )));
%                     var(obj.rxsig(24:end,n, m))
                    [obj.tgtpos, obj.tgtvel] = step(obj.htgtplat, obj.PRI);
                    [~,obj.tgtang] = rangeangle(obj.tgtpos, obj.htxplat.InitialPosition);
                end
                if obj.CLUTTER_ON
                    obj.rxsig(:,:,m) = obj.rxsig(:,:,m) + clut_real + 1i*clut_imag;
                end
                real_range{m} = norm(obj.tgtpos);
                real_speed{m} = norm(obj.tgtvel);
                % Apply Matched Filter
                [mfsig, mfgain] = step(obj.hmf,obj.rxsig(:, :, m) );
%                 mfgain
                % Shift the matched filter output
                mfsig=[mfsig(obj.gd+1:end, :); mfsig(1:obj.gd, :)];
                if ~obj.CFAR_ON
                    if strcmp(obj.Integrator, 'coherent') ||  strcmp(obj.Integrator, 'noncoherent') 
                        obj.threshold = sqrt( 2*obj.No* db2pow(npwgnthresh(obj.PFA, obj.npulses_per_prediction, obj.Integrator) )*db2pow(mfgain) );
    %                 elseif strcmp(obj.Integrator, 'binary  ')
    %                    obj.threshold = sqrt( obj.No * db2pow(npwgnthresh(obj.PFA, 1, 'noncoherent') )*db2pow(mfgain) ); 
                    elseif strcmp(obj.Integrator, 'none')
                        obj.threshold = sqrt( 2*obj.No * db2pow(npwgnthresh(obj.PFA, 1, 'noncoherent') )*db2pow(mfgain) ); 
                    end
                    if strcmp( obj.Integrator, 'noncoherent')
                        int_sig = pulsint( mfsig, 'noncoherent' );
                        range_estimates_index = find( int_sig > obj.threshold);
                    elseif strcmp( obj.Integrator, 'coherent')
                        int_sig =  real( pulsint(mfsig, 'coherent' )); 
                        range_estimates_index = find( int_sig > obj.threshold);
                    elseif strcmp( obj.Integrator, 'none')
                        mfsig = abs( mfsig(:, 1) );
                        range_estimates_index = find( mfsig > obj.threshold);
                    end
                else
                    if strcmp( obj.Integrator, 'noncoherent')
                        int_sig = pulsint( mfsig, 'noncoherent' );
                    elseif strcmp( obj.Integrator, 'coherent')
                        int_sig =  real( pulsint(mfsig, 'coherent' )); 
                    elseif strcmp( obj.Integrator, 'none')
                        int_sig = abs( mfsig(:, 1) );
                    end
                    [sigbin, signal_cfar, thresh] = obj.cfar_filter(int_sig, obj.K0, obj.gwndw, obj.rwndw, obj.cluttype);
                    range_estimates_index = find(sigbin);
                end
                obj.signal_cfar{m} = signal_cfar;
                obj.thresh{m} = thresh;
%                 [~, range_estimates_index] = 
%                 if length(range_estimates_index) > 1
%                     range_estimates_index = range_estimates_index(1);
%                 end
                pred_range{m} = obj.rangegates(range_estimates_index);
                error_range{m} = ( abs(pred_range{m} - real_range{m}) );
                %% Doppler Estimation
                if obj.PRED_VEL
                    speed_target = {} ;
                    for i=1:numel(pred_range{m})
                        [Pxx, Fx] = periodogram((obj.rxsig(range_estimates_index(i),:)).',[],512,obj.PRF, 'power','centered'); % Calculate power spectrum of rxsig along pulses.
                        speed_vec = dop2speed(Fx,obj.lamb)/2; % Translate to speed each frequency detected.
                        Pxx_norm = Pxx/max(Pxx); % Normalize Sxx to have its peak at 0 dB.
    %                     [~,detected_index] = findpeaks(pow2db(Pxx_norm),'MinPeakHeight',-5); % Select peaks whose power is at least 5dB below the max power.
                        [~,detected_index] = max(pow2db(Pxx_norm) );
                        prediction_sp = speed_vec(detected_index); % Speed whose sign follows the doppler frequency sign. That is, if the target is getting away from the platform, the sign is negative. 
                        speed_target{end+1} = -prediction_sp;
                    end
                    pred_speed{m} = cell2mat(speed_target); % The minus sign corrects the sign of the prediction.
                    error_speed{m} = abs(pred_speed{m} - real_speed{m});
                end
            end
%             mse_range = sqrt( mean( vertcat(error_range{:}).^2 ) );
%             mse_speed = sqrt( mean( vertcat(error_speed{:}).^2 ) );
            result_struct.pred_range = pred_range;
            result_struct.pred_speed = pred_speed;
            result_struct.real_range = real_range;
            result_struct.real_speed = real_speed;
            result_struct.error_range = error_range;
            result_struct.error_speed = error_speed;
%             result_struct = obj.flatten_result(result_struct);
            result_struct.detected_target = zeros(length(result_struct.error_range), 1);
            result_struct.false_alarms = zeros(length(result_struct.error_range), 1);
            for i=1:length(result_struct.error_range)
                result_struct.detected_target(i) = any(cell2mat( result_struct.error_range(i) ) <= obj.range_res/2);
                result_struct.false_alarms(i) = any(cell2mat( result_struct.error_range(i) ) >= obj.range_res*4);
            end
            result_struct.pd = sum(result_struct.detected_target)/length(result_struct.detected_target);
            result_struct.pfa = sum(result_struct.false_alarms)/length(result_struct.false_alarms);
            %             result_struct.pd = any( result_struct.error_range <= obj.range_res )/length( result_struct.error_range );
%             if isempty(result_struct.pd)
%                 result_struct.pd = 0;
%             end
        %             result_struct.mse_range = mse_range;
%             result_struct.mse_speed = mse_speed;
        end
        
        function range_estimates_index = binary_integration(obj, mfsig)
            mfsig = abs(mfsig) > obj.threshold; 
            mfsig = sum(mfsig, 2)./size(mfsig,2);
            range_estimates_index = find(mfsig > obj.binary_integ_percent_detected);
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
        
        function clut_sig = clutter_gen(vec_size, power, type)
            if strcmp(type, 'rayleigh')
                scale = sqrt(power/2); % https://www.sciencedirect.com/topics/engineering/rayleigh-distribution
                clut_sig = raylrnd(scale, vec_size(1), vec_size(2) );
            elseif strcmp(type, 'weibull')
                k = 1.5; % Shape
                lamb = sqrt(power/gamma(2/k + 1)); % Scale
                clut_sig = wblrnd(lamb, k, vec_size(1), vec_size(2));
            end
        end
        function [sigbin, signal, thresh] = cfar_filter(signal,K0, gwndw, rwndw, type)
            % K0 -> CFAR threshold constant.
            % gwndw -> size of guard window. 
            % rwndw -> size of reference window.
            sigbin = zeros(length(signal), 1);
            thresh = sigbin;
            if mod(rwndw, 2) ~= 0
               rwndw = rwndw + 1;
            end
            if mod(gwndw, 2) ~= 0
               gwndw = gwndw + 1;
            end
            signal = 999*(signal - min(signal))/( max(signal) - min(signal ) ) + 1;
            signal = 10*log10(signal);
            for i=1:length(signal)
                intv0 = max([1, i - (gwndw/2 + rwndw/2)]);
                intv1 = min( [ length(signal), i + (gwndw/2 + rwndw/2)]);

                intv_grd_0 = max([1, i - gwndw/2]);
                intv_grd_1 = min( [length(signal), i + gwndw/2] );
%                 nel_avg =  (intv1-intv0) - (intv_grd_1-intv_grd_0) ;
                samples_ref = signal( [intv0:(intv_grd_0-1) (intv_grd_1+1):intv1  ] );
                
%                 averager = ( sum( signal(intv0:intv1) ) - sum( signal( intv_grd_0:intv_grd_1 ) ) )/nel_avg;
                averager = mean(samples_ref);
                if strcmp(type, 'weibull')
                    padronizer = std(samples_ref);
                    thresh(i) = averager + K0*padronizer;
                elseif strcmp(type, 'rayleigh')
                    thresh(i) = averager + K0;
                end
                sigbin(i) = signal(i) > thresh(i);
            end
        end
    end
end


