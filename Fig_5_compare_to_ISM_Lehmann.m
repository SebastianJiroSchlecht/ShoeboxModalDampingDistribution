%% complementary code for publication 
% Scḧafer, M., Prawda, K., Rabenstein, R., Schlecht, S.J., "Distribution
% of Modal Damping in Absorptive Shoebox Rooms", in Proc. WASPAA 2023, New Paltz, NY, USA, Oct 22--25, 2023

% Compare proposed method to ISM and Lehmann methods
% and plot Figure 5
% uses:
% * shoebox2modes() 
% * histwv() 
% * modalSynthesis() 
% * zeroPhaseBandpass()
% * compute_echograms_mics()    (in './shoebox-roomsim')
% * render_mic_rirs()           (in './shoebox-roomsim')
% * rir_generator()             (in './RIR-Generator/')
% * fast_ISM_RoomResp()         (in './fastISM')
% * edc()

% Sebastian J. Schlecht, Thursday, 13 April 2023
% shoebox to damping density comparison to ISM
%% housekeeping 
clear; clc; close all;
addpath('./shoebox-roomsim')
addpath('./fastISM')

%% set simulation parameters
fs = 8000;          % Sample frequency (samples/s)

% room dimensions
Lx = 8.9;   % in [m]
Ly = 6.3;   % in [m]
Lz = 3.6;   % in [m]
L = [Lx, Ly, Lz];

c = 343;    % speed of sound in [m/s]

rec = [3.3 1.5 2.2];            % Receiver position [x y z] (m)
src = [0.8 4.5 1.1];            % Source position [x y z] (m)
limitsTime = 1.3;               % lenght of impulse response (seconds)
time = (1:limitsTime*fs).'/fs;  % seconds

bandpassEdges = [2000 4000];     % Hz

numberOfCases = 1;%3;

%% simulation
for it = 1:numberOfCases
    % reflection coefficients
    r = 0.97 * ones(1,6);
    r(3*it) = 0.9;
    r(it) = 0.85;
    r(6*it) = 0.75;

    [smu{it}, A{it}] = ...
        shoebox2modes(L,c, src, rec, r, bandpassEdges, 'simple', 'dirac');
    
    amplitude = A{it};
    decay = real(smu{it});

    powerDecay = decay*2;
    powerAmplitude = abs(amplitude).^2 / numel(decay);

    useMeanAmplitude = true;
    if useMeanAmplitude
        powerAmplitude = mean(powerAmplitude).*ones(size(amplitude));
    end
    
    [density(:,it),damping(:,it)] = histwv(powerDecay,powerAmplitude, -60, 0, 400);

    z = 0*density(:,it);
    h.density(:,it) = modalSynthesis(damping(:,it),z,density(:,it),z,limitsTime*fs,fs); 
    h.density(:,it) = sqrt(h.density(:,it));

    % image sources
    [abs_echograms, rec_echograms, echograms] = compute_echograms_mics(L, src, rec, 1-r.^2, limitsTime);
    h_temp = render_mic_rirs(abs_echograms, [], fs);
    h.ism(:,it) = h_temp(1:limitsTime*fs,:);
    
    % Lehmann
    h_temp = fast_ISM_RoomResp(fs,r,'t60',limitsTime,src,rec,L); % ,'Diffuse_dB',13
    h.lehmann(:,it) = h_temp;%(1:limitsTime*fs,:);
       
end

h.ism = zeroPhaseBandpass(h.ism,bandpassEdges,fs);
h.lehmann = zeroPhaseBandpass(h.lehmann,bandpassEdges,fs); 

%% colors for plotting
c2 =[239, 45, 86]./255;
c1 = [46, 134, 171]./255;
c3 = [255, 203, 119]./255;
%% plot Figure 5
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
cuton = 0.2*fs;
f = figure(1); hold on; grid on;

lines(3) = plot(time(1:length(h.lehmann)),edc(h.lehmann,cuton),'-','color', c3,'LineWidth',2);
lines(1) = plot(time,edc(h.ism,cuton),'-','color', c1, 'LineWidth',2);
lines(2) = plot(time,edc(h.density,cuton),'--','color', c2,'LineWidth',2);

set(gca, 'FontSize',12)
ylim([-80,40]);
xlabel('Time (seconds)', 'Interpreter','latex', 'FontSize',12)
ylabel('Energy Decay Curve (dB)', 'Interpreter','latex', 'FontSize',12);

box on
legend(lines, {'ISM','Proposed','Lehmann'}, 'Interpreter','latex', 'FontSize',12)
f.Position(end) = 270;


