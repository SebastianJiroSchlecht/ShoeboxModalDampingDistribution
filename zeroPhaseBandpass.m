function y = zeroPhaseBandpass(x,cutoff,fs)

bpFilt = designfilt('bandpassiir','FilterOrder',16, ...
         'HalfPowerFrequency1',cutoff(1),'HalfPowerFrequency2',cutoff(2), ...
         'SampleRate',fs);

% zero initial conditions
len = fs;
numCh = size(x,2);
xx = [zeros(len,numCh); x; zeros(len,numCh)];

yy = filtfilt(bpFilt,xx);

y = yy(len+1:(end-len),:);