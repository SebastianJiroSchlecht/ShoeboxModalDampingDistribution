function h = modalSynthesis(decay, frequency, amplitude, phase, length, fs)

h = zeros(length,1);
numModes = numel(decay);
time = (0:length-1).'/fs; % seconds

for it = 1:numModes
    mode = amplitude(it) .* exp(time .* decay(it) ) .* cos(2*pi*time*frequency(it) + phase(it));
    
    % sum to impulse response
    h = h + mode;
end

% h = h / sqrt(numModes);