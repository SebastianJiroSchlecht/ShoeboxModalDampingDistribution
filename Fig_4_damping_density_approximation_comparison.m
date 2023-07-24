%% complementary code for publication 
% Scḧafer, M., Prawda, K., Rabenstein, R., Schlecht, S.J., "Distribution
% of Modal Damping in Absorptive Shoebox Rooms", in Proc. WASPAA 2023, New Paltz, NY, USA, Oct 22--25, 2023

% Compare different methods for shoebox room damping density approximation 
% and plot Figure 4
% uses shoebox2modes() and histwv() 

% created by Sebastian J. Schlecht, Saturday, 15 April 2023
% updated by K. P.

%% housekeeping 
clear; clc; close all;

%% set simulation parameters 
fs = 8000;            % Sample frequency (samples/s)

% room dimensions
Lx = 8.9;   % in [m]
Ly = 6.3;   % in [m]
Lz = 3.6;   % in [m]
L = [Lx, Ly, Lz];

c = 343;    % speed of sound in [m/s]

rng('default');
lower_bound = 0.5*ones(4,3); %m source and receiver positions need to be at least this far away from a surface

recs = (L-lower_bound).*rand(4,3) + lower_bound; % Receiver position [x y z] (m)
srcs = (L-lower_bound).*rand(4,3) + lower_bound; % Source position [x y z] (m)

limitsTime = 1.3;               % lenght of impulse response (seconds)
time = (1:limitsTime*fs).'/fs;  % time axis seconds

bandpassEdges = [2000 4000];     % Hz

numberOfCases = 4; % compare 4 approximation methods 

% reflection coefficients 
r = 0.97 * ones(1,6);
r(1) = 0.85;
r(3) = 0.9;
r(6) = 0.75;

%% run simulations
for it = 1:numberOfCases

    rec = recs(it,:);
    src = srcs(it,:);
    switch it
        case {1, 2}
            useMeanAmplitude = false;
            useMonteCarlo = false;
        case {3}
            useMeanAmplitude = true;
            useMonteCarlo = false;
        case {4}
            useMeanAmplitude = true;
            useMonteCarlo = true;
    end

    [smu{it}, A{it}] = ...
        shoebox2modes(L,c, src, rec, r, bandpassEdges, 'simple', 'dirac');
    
    amplitude = A{it};
    decay = real(smu{it});

    if useMonteCarlo
        fraction = 0.001;
        ind = randperm(numel(amplitude),floor(numel(amplitude) * fraction));
        amplitude = amplitude(ind);
        decay = decay(ind);
    end

    powerDecay = decay*2;
    powerAmplitude = abs(amplitude).^2 / numel(decay);
 
    if useMeanAmplitude
        powerAmplitude = mean(powerAmplitude).*ones(size(amplitude));
    end
    
    [density(:,it),damping(:,it)] = histwv(powerDecay,powerAmplitude, -40, 0, 100);

end

%% colors
c1 =0*[116, 156, 117]./255;
c2 = [46, 134, 171]./255;
c3 = [255, 203, 119]./255;
c4 = [239, 45, 86]./255;

cMap = [c1; c2; c3; c4];
%% plot comparison
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
f = figure(1); hold on; grid on;
mylines = plot(damping,density,'LineWidth',2);
set(mylines,{'Color'},num2cell(cMap, 2))
colororder(cMap)
set(gca, 'FontSize', 12)
xlabel('Damping $\sigma$ (linear)', 'FontSize',12, 'Interpreter','latex')
ylabel('Density $H(\sigma)$ (linear)', 'FontSize',12, 'Interpreter','latex')
legend({'Source Receiver Position 1','Source Receiver Position 2','Mean Amplitude','Monte Carlo Selection'}, 'FontSize',12, 'Interpreter','latex')

f.Position(4) = 270;
box on
mylines(1).LineWidth = 4;