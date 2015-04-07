%% Complete Project - David Lavy : Silenced, Voiced and Unvoiced (SUV) in a Speech Signal 
%% Part 1. Start and Endpoint
clc, clear all, close all

[yn1,Fs,Nbits] = wavread('H.1.wav');
yn1 = yn1';

%%%%%%%%%% Specify Gaussian noise SNR coefficient %%%%%%%%%%%%
wgn = 10; % Specify Gaussian noise SNR (0,10,20), leave [] if not noise is added
if ~isempty(wgn)   
    yn = awgn(yn1,wgn,'measured');
    % Denoising the signal
    scal = 'one'; % Use model assuming standard Gaussian white noise.
    y = wden(yn,'heursure','s',scal,3,'sym8');
    figure, plot(yn1), hold on
    figure, plot(y,'r','Linewidth',2)
    sound(y,Fs);
else y=yn1; 
end
N = length(y); % Sample length of the signal
n = 0:N-1;
ts = n*(1/Fs); % Time length for signal

% Specifications of the window
wintype = 'hamming';
winlen = 400;
winshift = 100;

%% Step 1. Sampling Rate Conversion
Fres = 10000; % 10kHz
m = lcm(Fs,Fres); %Least common multiple of Fsin, Fsout
fs = m/Fs;      % inputs for the function 
fres = m/Fres;  % resample

Y = resample(y,fs,fres); % Resample the signal to 10kHz

%% Step 2. Highpass filter the signal to eliminate DC offset and hum
% Creating a high pass filter with low freq cutoff of 200Hz
d = fdesign.highpass('N,Fc,Ast,Ap',101,.04,60,1);
designmethods(d)
Hd = design(d,'equiripple');
% Applying the filter 
Yhp = filter(Hd,Y);

%% Step 3. Calculating log energy and Zerocross rates
N = length(Yhp) %Update the number of samples since resampling
i=1;
logEn=[]; %Vector to store log energy
ZCR=[]; %Vector to store zero cross rates
while (i+winlen-1 <= N)
    frame=Yhp(i:i+winlen-1).*hamming(winlen)';
    logEn=[logEn 10*log10(sum(frame.^2))]; % Calculate log Short-time Energy
    fr1=[0 Yhp(i:i+winlen-1)];
    frn=abs(diff(sign(fr1))).*hamming(winlen)'; % Calculate Zero Crossing Rate
    ZCR=[ZCR sum(frn)];
    i=i+winshift;
end
nfrm=length(logEn);
ZCR=ZCR*winshift/(2*winlen);

% Normalize log Short-time Energy with peak at 0dB
logEm=max(logEn);
logEn(find(logEn < logEm - 60))=logEm-60;
logen=logEn-logEm;

subplot(211), plot(logen)
title('Log Short-time Energy winlen=40ms, winshift=10ms');
subplot(212), plot(ZCR)
title('Short-time Zero Crossing Rate winlen=40ms, winshift=10ms');

%% Step 4. Start and Endpoint detection
% Define threshold parameters
IZCT = mean(ZCR(1:10))+0*std(ZCR(1:10)); 
ITU = -15; % in dB
ITR = -25 % in dB

peak = find(logen==0);
B1=[]; E1=[];
i=1;
while logen(i)<ITU
    i=i+1;
end
B1=i-1;
i=1;
while logen(end-i)<ITU
    i=i+1;
end
E1=length(logen)-i+1;

ss1 = find(ZCR(1:B1)>IZCT);
if length(ss1)>4
    B2 = ss1(1);
else B2=B1;
end

ss2 = find(ZCR(E1:end)>IZCT);
if ~isempty(ss2)
    E2 = E1+ss2(1);
else E2=E1;
end
% We have to make sure that the local region around B2 and E2
% doesn't exceed the log energy threshold ITR
flag1 = find(logen(1:B2)>ITR);
if ~isempty(flag1)
    B2=flag1(1)-1;
    if B2==0 B2=1; end % Assuring we don't have 0 as index
end

flag2 = find(logen(E2:end)>ITR);
if ~isempty(flag2)
    E2_=E2+flag2(end);
    ss4 = find(ZCR(E2_:end)>IZCT);
    if length(ss4)>4
        E2 = E2_+ss4(1);
    else E2 = E2_;
    end
    if E2>length(logen) E2=length(logen)-1; end
end

UE = mean(logen)-1.5*std(logen); %Unvoiced Energy
VE = UE; %Voiced Energy
VZ = mean(ZCR)-0*std(ZCR); %Voiced ZCR
UZ = VZ; %Unvoiced ZCR
%% Step 5. Segmenting Speech
% Voiced Speech if (Energy > VE) and (Zero-cross < VZ)
% Unvoiced Speech if (Energy < UE) and (Zero-cross > UZ)
% We need to use the vectors energy and cross
% We will denote 1 - Voiced Speech, 0 - No Speech, -1 - Unvoiced Speech
voiced=zeros(1,nfrm);
unvoiced=zeros(1,nfrm);
voiced(B2:E2) = (logen(B2:E2) > VE) .* (ZCR(B2:E2) < VZ);
unvoiced(B2:E2) = -1 * ((logen(B2:E2) < UE) .* (ZCR(B2:E2) > UZ));
total = voiced + unvoiced;
total(B2-1+find(total(B2:E2)<1))=-1;

% Result of adding noise is that a missclasification of silence as unvoiced
% can be done at the beginning or the end. We do a search around to correct
% it.
total(find(total(1:25)<0))=0;
total(length(total)-20+find(total(end-20:end)<0))=0;
figure,plot(total,'*'), title('Regions of voiced and unvoiced in speech')

%% Step 6. Segmenting the audio file from Voiced and Unvoiced
% Extending the total vector to the length of the audio
vec=find(total~=0);
AV = [];
AU = [];
for i=1:vec(end)
    if total(i)==0
        AV=[AV 0*ones(1,100)];
        AU=[AU 0*ones(1,100)];
    elseif total(i)==1
        AV=[AV 1*ones(1,100)];
        AU=[AU 0*ones(1,100)];
    elseif total(i)==-1
        AV=[AV 0*ones(1,100)];
        AU=[AU 1*ones(1,100)];
    end
end
rest=N-100*length(total);
if total(vec(end))==1
    AV=[AV 1*ones(1,rest)];
    AU=[AU 0*ones(1,rest)];
else
    AV=[AV 0*ones(1,rest)];
    AU=[AU 1*ones(1,rest)];
end
AV=[AV 0*ones(1,100*(length(total)-vec(end)))];
AU=[AU 0*ones(1,100*(length(total)-vec(end)))];
audio_voiced = Yhp.*AV;
audio_unvoiced = Yhp.*AU;
figure
n = [1:length(Yhp)]; % Sample index vector
t = n./Fres; % Sample time vector
subplot(211), plot(t,audio_voiced); 
title('Voiced Speech'); xlabel('Time (sec)'); ylabel('Amplitude');
subplot(212), plot(t,audio_unvoiced); title('Unvoiced Speech');
title('Unvoiced Speech'); xlabel('Time (sec)'); ylabel('Amplitude');
