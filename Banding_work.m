%Goal is to take some stretches of a given sleep stage and compare them 
%Literature supports using 4-10, 13-30, and 60-80 Hz bands as these have
%greatest differences in ON vs OFF DA-drug 

PT = 9;
%% Use for comparison of specific bouts of PSG stages
%Time1_start_epoch = 397; %Wake
%Time1_end_epoch = 561;
%Time2_start_epoch = 397; %REM - wake
%Time2_end_epoch = 561;
%Time3_start_epoch = 397; %N1 - wake
%Time3_end_epoch = 561;
%Time4_start_epoch = 288; %N2
%Time4_end_epoch = 350;
%Time5_start_epoch = 397; %N3 - wake
%Time5_end_epoch = 561;

%Time1 = lfp12((Time1_start_epoch*1024*30):(Time1_end_epoch*1024*30)); % Wake, epochs 397:561
%Time2 = lfp12((Time2_start_epoch*1024*30):(Time2_end_epoch*1024*30)); % REM, 
%Time3 = lfp12((Time3_start_epoch*1024*30):(Time3_end_epoch*1024*30)); % N1, 
%Time4 = lfp12((Time4_start_epoch*1024*30):(Time4_end_epoch*1024*30)); % N2, epochs 397:561
%Time5 = lfp12((Time5_start_epoch*1024*30):(Time5_end_epoch*1024*30)); % N3, 

%Wake = Time1;
%REM = Time2;
%N1 = Time3;
%N2 = Time4;
%N3 = Time5;
%% Use for comparison of ALL PSG stages
load('IndexPT9.mat');
IndexPT9 = IndexPT9(1:epochs,:); 
IndexPT9_score = IndexPT9(:,3);
% use to change values
IndexPT9_score( IndexPT9_score==0 )=90; % Wake
IndexPT9_score( IndexPT9_score==1 )=50;  % N1
IndexPT9_score( IndexPT9_score==2 )=40;  % N2
IndexPT9_score( IndexPT9_score==3 )=35;  % N3
IndexPT9_score( IndexPT9_score==5 )=80;  % REM

Signal0_block = zeros(epochs,30720);  %post processing from main script for LFP0

for k = 1:epochs; 
    Signal0_block(k,:) = Signal0( 1+(k-1)*30*1024 : k*30*1024);
end

if signal_3_lead0(30721,1) == Signal0_block(2,1)
    disp('GOOD')
else signal_3_lead0(30721,1) ~= Signal0_block(2,1);
    disp('ERROR');
end

if signal_3_lead0(614401,1) == Signal0_block(21,1)
    disp('GOOD')
else signal_3_lead0(614401,1) ~= Signal0_block(21,1);
    disp('ERROR');
end

Signal0_block_index = [Signal0_block, IndexPT9_score]; %epochs with indexing at end

%% Grouping by sleep score 

Signal0_block_1 = zeros(epochs, 30720); %Wake

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 90
     Signal0_block_1(k,:) = Signal0_block_index( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
PT9_Wake = reshape(Signal0_block_2,[needed_rows,1]); %back in vector form
%^ is it okay to lump everything together? cause pwelch's window might
%overlap on data from separate bouts of a given PSG stage

Signal0_block_1 = zeros(epochs, 30720); %REM

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 80
     Signal0_block_1(k,:) = Signal0_block_index( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
PT9_REM = reshape(Signal0_block_2,[needed_rows,1]);

Signal0_block_1 = zeros(epochs, 30720); %N1

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 50
     Signal0_block_1(k,:) = Signal0_block_index( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
PT9_N1 = reshape(Signal0_block_2,[needed_rows,1]);

Signal0_block_1 = zeros(epochs, 30720); %N2

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 40
     Signal0_block_1(k,:) = Signal0_block_index( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
PT9_N2 = reshape(Signal0_block_2,[needed_rows,1]);

Signal0_block_1 = zeros(epochs, 30720); %N3

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 35
     Signal0_block_1(k,:) = Signal0_block_index( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
PT9_N3 = reshape(Signal0_block_2,[needed_rows,1]);

fs = 1024;
Fs = fs;

Wake = PT9_Wake;
REM = PT9_Wake;
N1 = PT9_N1;
N2 = PT9_N2;
N3 = PT9_Wake;

%% 4-10 Hz Band
band_start = 4;
band_end = 10;
band = [band_start band_end];

%Each PSG stage powers'
[Pxx,F] = periodogram(Wake,rectwin(length(Wake)),length(Wake),Fs);
pband1 = bandpower(Pxx,F,band,'psd');
ptot1 = bandpower(Pxx,F,'psd');
per_power1 = 100*(pband1/ptot1);

[Pxx,F] = periodogram(REM,rectwin(length(REM)),length(REM),Fs);
pband2 = bandpower(Pxx,F,band,'psd');
ptot2 = bandpower(Pxx,F,'psd');
per_power2 = 100*(pband2/ptot2);

[Pxx,F] = periodogram(N1,rectwin(length(N1)),length(N1),Fs);
pband3 = bandpower(Pxx,F,band,'psd');
ptot3 = bandpower(Pxx,F,'psd');
per_power3 = 100*(pband3/ptot3);

[Pxx,F] = periodogram(N2,rectwin(length(N2)),length(N2),Fs);
pband4 = bandpower(Pxx,F,band,'psd');
ptot4 = bandpower(Pxx,F,'psd');
per_power4 = 100*(pband4/ptot4);

[Pxx,F] = periodogram(N3,rectwin(length(N3)),length(N3),Fs);
pband5 = bandpower(Pxx,F,band,'psd');
ptot5 = bandpower(Pxx,F,'psd');
per_power5 = 100*(pband5/ptot5);

Per_Power = per_power1+per_power2+per_power3+per_power4+per_power5;
disp(Per_Power);

if ptot1>ptot4
    disp('ptot1 larger')
else
    disp('ptot4 larger')
end

window = [];
noverlap = [];
f = window; %mathworks tutorial said that its safest to set to 'window'

figure(1) %subplots for each PSG stage for 'band'
subplot(5,1,1)
pwelch(Wake,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for Wake ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,2)
pwelch(REM,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power')
title(['PT' num2str(PT) ', power spectrum for REM ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,3)
pwelch(N1,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N1 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,4)
pwelch(N2,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N2 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,5)
pwelch(N3,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N3 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

%pwelch overlay chaos
%window = (band_end - band_start)*50;
%noverlap = (band_end - band_start)*25;

figure(101) 
[pxx_1, f] = pwelch(Wake,window,noverlap,window,fs);
plot(f, 10*log10(pxx_1))
hold on
[pxx_2, f] = pwelch(REM,window,noverlap,window,fs);
plot(f, 10*log10(pxx_2))
hold on
[pxx_3, f] = pwelch(N1,window,noverlap,window,fs);
plot(f, 10*log10(pxx_3))
hold on
[pxx_4, f] = pwelch(N2,window,noverlap,window,fs);
plot(f, 10*log10(pxx_4))
hold on
[pxx_5, f] = pwelch(N3,window,noverlap,window,fs);
plot(f, 10*log10(pxx_5))
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for freqencies ' num2str(band_start) ' - ' num2str(band_end) 'Hz'])
legend('Wake','REM','N1', 'N2', 'N3')
legend('boxoff')


%% 13- 30 Hz Band
band_start = 13;
band_end = 30;
band = [band_start band_end];


%Each PSG stage powers'
[Pxx,F] = periodogram(Wake,rectwin(length(Wake)),length(Wake),Fs);
pband1 = bandpower(Pxx,F,band,'psd');
ptot1 = bandpower(Pxx,F,'psd');
per_power1 = 100*(pband1/ptot1);

[Pxx,F] = periodogram(REM,rectwin(length(REM)),length(REM),Fs);
pband2 = bandpower(Pxx,F,band,'psd');
ptot2 = bandpower(Pxx,F,'psd');
per_power2 = 100*(pband2/ptot2);

[Pxx,F] = periodogram(N1,rectwin(length(N1)),length(N1),Fs);
pband3 = bandpower(Pxx,F,band,'psd');
ptot3 = bandpower(Pxx,F,'psd');
per_power3 = 100*(pband3/ptot3);

[Pxx,F] = periodogram(N2,rectwin(length(N2)),length(N2),Fs);
pband4 = bandpower(Pxx,F,band,'psd');
ptot4 = bandpower(Pxx,F,'psd');
per_power4 = 100*(pband4/ptot4);

[Pxx,F] = periodogram(N3,rectwin(length(N3)),length(N3),Fs);
pband5 = bandpower(Pxx,F,band,'psd');
ptot5 = bandpower(Pxx,F,'psd');
per_power5 = 100*(pband5/ptot5);

Per_Power = per_power1+per_power2+per_power3+per_power4+per_power5;
disp(Per_Power);

if ptot1>ptot4
    disp('ptot1 larger')
else
    disp('ptot4 larger')
end

window = [];
noverlap = [];
f = window; %mathworks tutorial said that its safest to set to 'window'

figure(2) %subplots for each PSG stage for 'band'
subplot(5,1,1)
pwelch(Wake,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for Wake ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,2)
pwelch(REM,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power')
title(['PT' num2str(PT) ', power spectrum for REM ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,3)
pwelch(N1,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N1 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,4)
pwelch(N2,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N2 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,5)
pwelch(N3,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N3 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

%pwelch overlay chaos
%window = (band_end - band_start)*50;
%noverlap = (band_end - band_start)*25;

figure(201) 
[pxx_1, f] = pwelch(Wake,window,noverlap,window,fs);
plot(f, 10*log10(pxx_1))
hold on
[pxx_2, f] = pwelch(REM,window,noverlap,window,fs);
plot(f, 10*log10(pxx_2))
hold on
[pxx_3, f] = pwelch(N1,window,noverlap,window,fs);
plot(f, 10*log10(pxx_3))
hold on
[pxx_4, f] = pwelch(N2,window,noverlap,window,fs);
plot(f, 10*log10(pxx_4))
hold on
[pxx_5, f] = pwelch(N3,window,noverlap,window,fs);
plot(f, 10*log10(pxx_5))
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for freqencies ' num2str(band_start) ' - ' num2str(band_end) 'Hz'])
legend('Wake','REM','N1', 'N2', 'N3')
legend('boxoff')

%% 60-80 Hz Band
band_start = 60;
band_end = 80;
band = [band_start band_end];


%Each PSG stage powers'
[Pxx,F] = periodogram(Wake,rectwin(length(Wake)),length(Wake),Fs);
pband1 = bandpower(Pxx,F,band,'psd');
ptot1 = bandpower(Pxx,F,'psd');
per_power1 = 100*(pband1/ptot1);

[Pxx,F] = periodogram(REM,rectwin(length(REM)),length(REM),Fs);
pband2 = bandpower(Pxx,F,band,'psd');
ptot2 = bandpower(Pxx,F,'psd');
per_power2 = 100*(pband2/ptot2);

[Pxx,F] = periodogram(N1,rectwin(length(N1)),length(N1),Fs);
pband3 = bandpower(Pxx,F,band,'psd');
ptot3 = bandpower(Pxx,F,'psd');
per_power3 = 100*(pband3/ptot3);

[Pxx,F] = periodogram(N2,rectwin(length(N2)),length(N2),Fs);
pband4 = bandpower(Pxx,F,band,'psd');
ptot4 = bandpower(Pxx,F,'psd');
per_power4 = 100*(pband4/ptot4);

[Pxx,F] = periodogram(N3,rectwin(length(N3)),length(N3),Fs);
pband5 = bandpower(Pxx,F,band,'psd');
ptot5 = bandpower(Pxx,F,'psd');
per_power5 = 100*(pband5/ptot5);

Per_Power = per_power1+per_power2+per_power3+per_power4+per_power5;
disp(Per_Power);

if ptot1>ptot4
    disp('ptot1 larger')
else
    disp('ptot4 larger')
end

window = [];
noverlap = [];
f = window; %mathworks tutorial said that its safest to set to 'window'

figure(3) %subplots for each PSG stage for 'band'
subplot(5,1,1)
pwelch(Wake,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for Wake ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,2)
pwelch(REM,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power')
title(['PT' num2str(PT) ', power spectrum for REM ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,3)
pwelch(N1,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N1 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,4)
pwelch(N2,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N2 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

subplot(5,1,5)
pwelch(N3,window,noverlap,f,fs);
%set(gca, 'YLim', []);
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for N3 ' num2str(band_start) '-' num2str(band_end) ' Hz'])

%pwelch overlay chaos
%window = (band_end - band_start)*50;
%noverlap = (band_end - band_start)*25;

figure(301) 
[pxx_1, f] = pwelch(Wake,window,noverlap,window,fs);
plot(f, 10*log10(pxx_1))
hold on
[pxx_2, f] = pwelch(REM,window,noverlap,window,fs);
plot(f, 10*log10(pxx_2))
hold on
[pxx_3, f] = pwelch(N1,window,noverlap,window,fs);
plot(f, 10*log10(pxx_3))
hold on
[pxx_4, f] = pwelch(N2,window,noverlap,window,fs);
plot(f, 10*log10(pxx_4))
hold on
[pxx_5, f] = pwelch(N3,window,noverlap,window,fs);
plot(f, 10*log10(pxx_5))
set(gca, 'XLim', band);
xlabel('Frequency (Hz)')
ylabel('Power');
title(['PT' num2str(PT) ', power spectrum for freqencies ' num2str(band_start) ' - ' num2str(band_end) 'Hz'])
legend('Wake','REM','N1', 'N2', 'N3')
legend('boxoff')


%% Big comparison
figure(10)
window = 1024;
overlap = 512;

subplot(3,2,1) 
spectrogram(Time1,window, overlap, 4:10, fs, 'yaxis')
title('Time1 LFP01 4-10Hz PT 9');
colorbar off

subplot(3,2,2) 
spectrogram(Time3,window, overlap, 4:10, fs, 'yaxis')
title('Time2 LFP01 4-10Hz PT 9');
colorbar off

subplot(3,2,3) 
spectrogram(Time1,window, overlap, 13:30, fs, 'yaxis')
title('Time1 LFP01 13-30 Hz PT 9');
colorbar off

subplot(3,2,4) 
spectrogram(Time3,window, overlap, 13:30, fs, 'yaxis')
title('Time2 LFP01 13-30 Hz PT 9');
colorbar off

subplot(3,2,5) 
spectrogram(Time1,window, overlap, 60:80, fs, 'yaxis')
title('Time1 LFP01 60-80 Hz PT 9');
colorbar off

subplot(3,2,6) 
spectrogram(Time3,window, overlap, 60:80, fs, 'yaxis')
title('Time2 LFP01 60-80 Hz PT 9');
colorbar off

