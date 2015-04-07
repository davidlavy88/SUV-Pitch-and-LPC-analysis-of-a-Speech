%% Complete Project - David Lavy : Ptch Frequency Contour
%% Part 2: Pitch Detector
% Method used: Autocorrelation-Based Pitch Estiamtor
clear all, clc, close all
% [y,Fs,N]=wavread('we were away a year ago_lrr.wav');
[y,Fs,N]=wavread('H.1.wav');
cl = 60; % Clipper parameter for center clipping level (in percentage)
smooth = 5; %5 point median smoother

% Specifications of the window
wintype = 'hamming';
winlen = 400; % Length segments of 30msec
winshift = 100; % Length shift of 10msec

%% Step 1. Sampling Rate Conversion
Fres = 10000; % 10kHz
m = lcm(Fs,Fres); %Least common multiple of Fsin, Fsout
fs = m/Fs;      % inputs for the function 
fres = m/Fres;  % resample

Y = resample(y,fs,fres); % Resample the signal to 10kHz
N = length(Y);

%% Step 2. Highpass filter the signal to eliminate DC offset and hum
% Creating a high pass filter with low freq cutoff of 200Hz
Fcut = 200;
% Fstop = Fcut-25;
% Fpass = Fcut+25;
Fstop = 80;
Fpass = 150;
d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
  'PassbandFrequency',Fpass,'SampleRate',Fres,'DesignMethod','equiripple');
D = mean(grpdelay(d)) % filter delay in samples

Yf1 = filter(d,[Y; zeros(D,1)]); % Append D zeros to the input data
Yf1 = Yf1(D+1:end);              % Shift data to compensate for delay

%% Step 3. Lowpass filter the signal at a freq cutoff of 900Hz
Fcut = 900;
Fstop = 920;
Fpass = 880;
dl = designfilt('lowpassfir','StopbandFrequency',Fstop, ...
  'PassbandFrequency',Fpass,'SampleRate',Fres,'DesignMethod','equiripple');
D = mean(grpdelay(dl)) % filter delay in samples

Yf = filter(dl,[Yf1; zeros(D,1)]); % Append D zeros to the input data
Yf = Yf(D+1:end);                  % Shift data to compensate for delay

pitchp_low=75;
pitchp_high=200;

%% Step 4. Calculate the correlation in the signal with a window of winlen and a shift of winshift
wincen=winlen/2;
contour=[];
while (wincen+winlen/2 <= N)
    fra=Yf(max(wincen-winlen/2,1):wincen+winlen/2); % x[n]
    frb=Yf(wincen+winlen/2+1:min(wincen+winlen/2+pitchp_high,N));
    frame2=[fra; frb]; %x_hat[n]
            
% Center Clipping at clip = % of Amax
    Amax=max(frame2);
    clip=Amax*(cl/100);
    fra(find(abs(fra) < clip))=0;
    frame2(find(abs(frame2) < clip))=0;
    % Modification for 3 level center clipper
    fra(find(fra>=clip))=1;
    fra(find(fra<=-clip))=-1;
    frame2(find(frame2>=clip))=1;
    frame2(find(frame2<=-clip))=-1;
    
% compute modified autocorrelation of fra and frame2
    lfr2=length(frame2);
    c=xcorr(fra,frame2);
    cmax=max(c); % we find the largest peak of the autocorr function
    pmax=max(c(lfr2+pitchp_low:lfr2+pitchp_high)); % fixed threshold
    if (cmax > 0)
        ploc=find(c(lfr2+pitchp_low:lfr2+pitchp_high) == pmax);
        contour=[contour ploc(1)+pitchp_low-1];
    else
        pmax=0;
        contour=[contour 0];
    end
    wincen=wincen+winshift;
end 

%% Step 5. Smooth Filter
pm=medsmth(1./contour,smooth);
pm1=medsmth(1./pm,smooth);
pm1=pm1';pm1 = circshift(pm1,100);pm1=pm1'; %Compensate the shift from autocorr
contour=contour';contour = circshift(contour,100);contour=contour'; %Compensate the shift from autocorr
subplot(211), plot(contour,'*'); title('Raw Pitch Contour');
subplot(212), plot(pm1,'*r'); title('Smooth Pitch Contour');
