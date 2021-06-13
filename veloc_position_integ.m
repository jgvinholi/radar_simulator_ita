%% Simulating radar with noncoherent integration of 100 pulses.
%% Initialize Radar Object.
radobj = RadarSimulator;
sim_duration = 5; % time in seconds to simulate.

%% SNR = -5 dB.
radobj.dSNR = -5; % dB
radobj = update_parameters(radobj);
result_struct_0 = runSim(radobj, sim_duration);
RadarSimulator.plot_range_speed_pred(result_struct_0, radobj)


%% SNR = 0 dB.
radobj.dSNR = 0; % dB
radobj = update_parameters(radobj);
result_struct_1 = runSim(radobj, sim_duration);
RadarSimulator.plot_range_speed_pred(result_struct_1, radobj)

%% SNR = 5 dB.
radobj.dSNR = 5; % dB
radobj = update_parameters(radobj);
result_struct_2 = runSim(radobj, sim_duration);
RadarSimulator.plot_range_speed_pred(result_struct_2, radobj)


%% SNR = 10 dB.
radobj.dSNR = 10; % dB
radobj = update_parameters(radobj);
result_struct_3 = runSim(radobj, sim_duration);
RadarSimulator.plot_range_speed_pred(result_struct_3, radobj)

%% SNR = 15 dB.
radobj.dSNR = 15; % dB
radobj = update_parameters(radobj);
result_struct_4 = runSim(radobj, sim_duration);
RadarSimulator.plot_range_speed_pred(result_struct_4, radobj)


%% SNR = 20 dB.
radobj.dSNR = 20; % dB
radobj = update_parameters(radobj);
result_struct_5 = runSim(radobj, sim_duration);
RadarSimulator.plot_range_speed_pred(result_struct_5, radobj)


