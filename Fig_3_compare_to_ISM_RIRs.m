%% complementary code for publication 
% Scḧafer, M., Prawda, K., Rabenstein, R., Schlecht, S.J., "Distribution
% of Modal Damping in Absorptive Shoebox Rooms", in Proc. WASPAA 2023, New Paltz, NY, USA, Oct 22--25, 2023

% Compare RIRs obatined with proposed method and ISM 
% and plot Figure 3

% created by K. P.

%% housekeeping 
clear; clc; close all;

%% read audiofiles
[h.ism, fs] = audioread('ism_rir.wav');
h.density = audioread('simple_rir.wav');

time = 0:1/fs:(length(h.ism)-1)/fs; % time axis

%% plot Figure 3
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

f = figure(1);clf; hold on; grid on;
plot(time, h.ism,'LineWidth',1.5);
plot(time, h.density,'r:','LineWidth',1);
xlim([0,0.04])
ylim(0.04*[-1 1])
xlabel('Time (seconds)', 'Interpreter','latex', 'FontSize',12)
ylabel('Amplitude (linear)', 'Interpreter','latex', 'FontSize',12)
legend('ISM','Proposed','location', 'southeast', 'Interpreter','latex', 'FontSize',12);
box on
set(gca, 'Ytick', [-0.04, 0, 0.04], 'FontSize',12)
f.Position(end) = 210;