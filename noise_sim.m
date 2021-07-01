dbstop if error
obj = RadarSimulator;
simtime = 1; % seconds
obj.npulses_per_prediction = 10;
SNR = -10:5:15; % db
Pd_ncoh = zeros(length(SNR), 1);
Pd_coh = zeros(length(SNR), 1);
Pd_bin = zeros(length(SNR), 1);
Pd_none = zeros(length(SNR), 1);
pfa_ncoh = zeros(length(SNR), 1);
pfa_coh = zeros(length(SNR), 1);
pfa_bin = zeros(length(SNR), 1);
pfa_none = zeros(length(SNR), 1);



for i=1:length(SNR)
    obj = RadarSimulator;
    obj.K0 = 3.5;
    obj.dSNR = SNR(i);
    
%     obj.K0 = K0*sqrt(obj.npulses_per_prediction);
    obj.Integrator = 'noncoherent';
    obj = update_parameters(obj);
    result_struct_ncoh = runSim(obj, simtime);
    Pd_ncoh(i) = result_struct_ncoh.pd;
    pfa_ncoh(i) = result_struct_ncoh.pfa;
    
    obj.Integrator = 'coherent';
    obj = update_parameters(obj);
    result_struct_coh = runSim(obj, simtime);
    Pd_coh(i) = result_struct_coh.pd;
    pfa_coh(i) = result_struct_coh.pfa;
    
    
%     obj.K0 = K0;
    obj.Integrator = 'none';
    obj = update_parameters(obj);
    result_struct_none = runSim(obj, simtime);
    
    Pd_none(i) = result_struct_none.pd;
    pfa_none(i) = result_struct_none.pfa;
    
%     obj.Integrator = 'binary';
%     obj = update_parameters(obj);
%     result_struct_bin = runSim(obj, simtime);
%     Pd_bin(i) = result_struct_bin.pd;
end

figure
hold on;
grid on

plot(SNR, Pd_bin)
plot(SNR, Pd_coh)
plot(SNR, Pd_ncoh)
plot(SNR, Pd_none)
title('Pd')
legend('Binary', 'Coherent', 'Noncoherent', 'None')


figure
hold on;
grid on

plot(SNR, pfa_bin)
plot(SNR, pfa_coh)
plot(SNR, pfa_ncoh)
plot(SNR, pfa_none)
title('PFA')
legend('Binary', 'Coherent', 'Noncoherent', 'None')

