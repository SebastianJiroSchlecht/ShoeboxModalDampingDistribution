function [smu, amplitude, lambda] = shoebox2modes(L,c, src, rec, r, bandpassEdges, evcase, sourcetype)
%% Modal Synthesis for the lossless case 
%  based on the derivations of Rudi (from 28.03.2032) for all reflection factors 
%  are equal to 1. 
%  
%  Room is excited by a point source 
%
% Input parameters L = [Lx Ly Lz] room dimensions in meters
% c speed of sound in meters/second
% src - source pos
% rec - receiver pos
% r - reflection factors
% bandpassEdges - lower and upper limit of frequency

%% Parameters
[r10, r11, r20, r21, r30, r31] = feval(@(x)x{:}, num2cell(r));

rhox = log(r10*r11); 
rhoy = log(r20*r21); 
rhoz = log(r30*r31); 
rho = [rhox, rhoy, rhoz];

T_ = L./c; %travel time
%% Calculate eigenvalues (either with grid or monte carlo) 

% roughly set the necessary mode indices to cover all frequencies up to the
% limit
K = ceil(bandpassEdges(2) .* (2*L)/c);

[I1,I2,I3] = ndgrid(0:K(1)-1,0:K(2)-1,0:K(3)-1);
index = [I1(:), I2(:), I3(:)];
index = index(2:end,:); % remove [0 0 0]

M = size(index,1);

smu = zeros(1,M);

if(mean(r) == 1)
    evcase = 'lossless';
end

switch evcase 
    case 'lossless'
        mu = 1:M;
        k = index(mu, :);
        s = 1j.*k*c*pi./L;
        smu(mu) = sqrt(sum(s.^2,2));
        
    case 'simple'
        mu = 1:M;
        k = index(mu, :);
        s = 1./(2*T_).*rho + 1j*pi./T_.*k;
        smu(mu) = -sqrt(sum(s.^2,2)); %sqrt(sx^2 + sy^2 + sz^2);
        

    case 'numerical'
        [smu, ] = lossyEVs_numerical_v2(L, r, c, index);
        
    case 'numerical_approx'
        [smu] = lossyEVs_numerical_approxKutt(index, T_, r);
    case 'feedback'
        % to be added eigenvalue calculation via feedback loops
end

%% order for frequencies
[freqs, reari] = sort(imag(smu)/(2*pi));
% cut at frequency limit
reari = reari(freqs < bandpassEdges(2));

smu = smu(reari);
index = index(reari,:,:);
M = size(index,1);

%% Find mode types 
% This simplifies the search for mode types later to use the correct
% weigthing factors 

% indices of axial modes 
mx = find(index(:,1) ~= 0 & index(:,2) == 0 & index(:,3) == 0);
my = find(index(:,1) == 0 & index(:,2) ~= 0 & index(:,3) == 0);
mz = find(index(:,1) == 0 & index(:,2) == 0 & index(:,3) ~= 0);

% index tangential modes
mt1 = find(index(:,1) ~= 0 & index(:,2) ~= 0 & index(:,3) == 0);
mt2 = find(index(:,1) == 0 & index(:,2) ~= 0 & index(:,3) ~= 0);
mt3 = find(index(:,1) ~= 0 & index(:,2) == 0 & index(:,3) ~= 0);
% oblique modes 
mo = find(index(:,1) ~= 0 & index(:,2) ~= 0 & index(:,3) ~= 0);

%% Eigenfunctions for the losless case 

% we use here already the higher level counting by mu = [kx, ky, kz] 
% e.g., kx = index(1,mu); 

phi = @(x,y,z,mu) prod(cos(index(mu,:)*pi.*[x y z]./L),2);

%% Modal scaling factors 
% The factor K(0,mu) from eq (26) is not included because we know that this
% factor will cancel out in the synthesis and is not further determined 

V = prod(L); %Volume in m3, Lx*Ly*Lz; 

% Asign lambda according to eq. (27) 
lambda = zeros(1,M); 
% kx = ky = kz = 0; 
lambda(1) = 1; 
% all axial modes 
lambda(mx) = 2; 
lambda(my) = 2; 
lambda(mz) = 2; 
% all tangential modes 
lambda(mt1) = 4; 
lambda(mt2) = 4; 
lambda(mt3) = 4; 
% all oblique modess
lambda(mo) = 8; 

% calculate lambda according to eq. (26) 
nmu = V./lambda; 

%% Excitation 
% we asume a dirac at position src = [xe, ye, ze] at t = 0

% amplitude of the dirac excitation 
fs = 1;

% excitation weights 
switch sourcetype
    case 'dirac' 
        mu = 1:M; 
        phi_xs = phi(src(1),src(2),src(3),mu).';
    case 'distributed'
        x0 = 0.1; % width of the raised cosines 
        xe = src(1); 
        ye = src(2); 
        ze = src(3);
        
        Lx = L(1); 
        Ly = L(2); 
        Lz = L(3);
        % perform each calculation individually 
        func = @(var,k,ext,wid,len) cos(k*pi*var/len).*1/wid.*(1 + ...
            cos(2*pi/wid*(var - ext)));
        
        phi_xs = zeros(1,M); 
        for n = 1:M 
            xcomp = integral(@(var)func(var, index(n,1), xe, x0, Lx), xe - x0/2, xe + x0/2);
            ycomp = integral(@(var)func(var, index(n,2), ye, x0, Ly), ye - x0/2, ye + x0/2);
            zcomp = integral(@(var)func(var, index(n,3), ze, x0, Lz), ze - x0/2, ze + x0/2);
            phi_xs(n) = xcomp*ycomp*zcomp;
        end
end


%% Observation position 
% observation of the sound pressure at rec = [xr, yr, zr]
% basis functions for synthesis 
mu = 1:M; 
phi_xo = phi(rec(1),rec(2),rec(3),mu);

%% Modal Amplitude
amplitude = 2*fs/V*lambda.*phi_xs.*phi_xo.';
% amplitude = amplitude ./ imag(smu); % magic scaling factor
% amplitude = amplitude * 8; % more random factors
% amplitude = amplitude / (4 * pi);


% %% just take one type of modes
% smu = smu(mo);

%% Filter modes
frequency = imag(smu)/(2*pi);
id = bandpassEdges(1) < frequency & frequency < bandpassEdges(2);

%% Format
smu = smu(id);
amplitude = amplitude(id);

end