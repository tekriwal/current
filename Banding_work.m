%Pt9 frequency banding work
%Goal is to take some stretches of a given sleep stage and compare them 
%Literature supports using 4-10, 13-30, and 60-80 Hz bands as these have
%greatest differences in ON vs OFF DA-drug 

Time1 = lfp12((397*1024*30):(561*1024*30)); % Wake, epochs 397:561
Time2 = lfp12((288*1024*30):(350*1024*30)); % N2, epochs 288:350

Fs = 1024;
fs = 1024;
x = Time1;
y = Time2;

%% 4-10 Hz Band
f = [4 10];
%ff = 4; %figure out how to get variables in text
%ll = 10;

[Pxx,F] = periodogram(x,rectwin(length(x)),length(x),Fs);
pband1 = bandpower(Pxx,F,f,'psd');
ptot1 = bandpower(Pxx,F,'psd');
per_power1 = 100*(pband1/ptot1);

[Pxx,F] = periodogram(y,rectwin(length(y)),length(y),Fs);
pband2 = bandpower(Pxx,F,f,'psd');
ptot2 = bandpower(Pxx,F,'psd');
per_power2 = 100*(pband2/ptot2);

if ptot1>ptot2
    disp('ptot1 larger')
else
    disp('ptot2 larger')
end

window = [];
noverlap = [];

figure(1)
subplot(2,1,1)
pwelch(x,window,noverlap,[],fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', f);
xlabel('Frequency (Hz)')
ylabel('Power');
title('PT9 LFP12 Wake')

subplot(2,1,2)
pwelch(y,window,noverlap,[],fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', f);
xlabel('Frequency (Hz)')
ylabel('Power');
title('PT9 LFP12 N2')
%% 13- 30 Hz Band
f = [13 30];
%ff = 13; %figure out how to get variables in text
%ll = 30;

[Pxx,F] = periodogram(x,rectwin(length(x)),length(x),Fs);
pband3 = bandpower(Pxx,F,f,'psd');
ptot3 = bandpower(Pxx,F,'psd');
per_power3 = 100*(pband3/ptot3);

[Pxx,F] = periodogram(y,rectwin(length(y)),length(y),Fs);
pband4 = bandpower(Pxx,F,f,'psd');
ptot4 = bandpower(Pxx,F,'psd');
per_power4 = 100*(pband4/ptot4);

if ptot3>ptot4
    disp('ptot3 larger')
else
    disp('ptot4 larger')
end

window = [];
noverlap = [];

figure(2)
subplot(2,1,1)
pwelch(x,window,noverlap,[],fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', f);
xlabel('Frequency (Hz)')
ylabel('Power');
title('PT9 LFP12 Wake')

subplot(2,1,2)
pwelch(y,window,noverlap,[],fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', f);
xlabel('Frequency (Hz)')
ylabel('Power');
title('PT9 LFP12 N2')
%% 60-80 Hz Band
f = [60 80];
%ff = 60; %figure out how to get variables in text
%ll = 80;

[Pxx,F] = periodogram(x,rectwin(length(x)),length(x),Fs);
pband5 = bandpower(Pxx,F,f,'psd');
ptot5 = bandpower(Pxx,F,'psd');
per_power5 = 100*(pband5/ptot5);

[Pxx,F] = periodogram(y,rectwin(length(y)),length(y),Fs);
pband6 = bandpower(Pxx,F,f,'psd');
ptot6 = bandpower(Pxx,F,'psd');
per_power6 = 100*(pband6/ptot6);

if ptot5>ptot6
    disp('ptot5 larger')
else
    disp('ptot6 larger')
end

window = [];
noverlap = [];

figure(3)
subplot(2,1,1)
pwelch(x,window,noverlap,[],fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', f);
xlabel('Frequency (Hz)')
ylabel('Power');
title('PT9 LFP12 Wake')

subplot(2,1,2)
pwelch(y,window,noverlap,[],fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', f);
xlabel('Frequency (Hz)')
ylabel('Power');
title('PT9 LFP12 N2')

%% Big comparison
figure(4)
window = 1024;
overlap = 512;

subplot(3,2,1) 
spectrogram(Time1,window, overlap, 4:10, fs, 'yaxis')
title('Time1 LFP01 4-10Hz PT 9');
colorbar off

subplot(3,2,2) 
spectrogram(Time2,window, overlap, 4:10, fs, 'yaxis')
title('Time2 LFP01 4-10Hz PT 9');
colorbar off

subplot(3,2,3) 
spectrogram(Time1,window, overlap, 13:30, fs, 'yaxis')
title('Time1 LFP01 13-30 Hz PT 9');
colorbar off

subplot(3,2,4) 
spectrogram(Time2,window, overlap, 13:30, fs, 'yaxis')
title('Time2 LFP01 13-30 Hz PT 9');
colorbar off

subplot(3,2,5) 
spectrogram(Time1,window, overlap, 60:80, fs, 'yaxis')
title('Time1 LFP01 60-80 Hz PT 9');
colorbar off

subplot(3,2,6) 
spectrogram(Time2,window, overlap, 60:80, fs, 'yaxis')
title('Time2 LFP01 60-80 Hz PT 9');
colorbar off

