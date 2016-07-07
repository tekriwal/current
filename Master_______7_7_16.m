%To-Do
%for lfp comparison set axis to total length of lfp01/whatever
%need to change IndexPTX_EMG according to PT
%need to find a better way to remove 60Hz noise
%% code for one segment of data, as opposed to two
%% if two segments, then add: below where you load PTX.matfile (just copy Pt10)
%lfpData1 = data(1).data(:,31:34);
%lfpData2 = data(2).data(:,31:34); %this is right
%allLFPs = [lfpData1 ; lfpData2];
%clear lfpData1
%clear lfpData2
%lfplead0 = allLFPs(:,1);
%lfplead1 = allLFPs(:,2);
%lfplead2 = allLFPs(:,3);
%lfplead3 = allLFPs(:,4);
%clear allLfPs

%also for Chin_reference, need to change, as well as last plot with FDC etc
%%
%Chin_referenced_1 = data(1).data(:,25)-data(1).data(:,26);
%Chin_referenced_2 = data(2).data(:,25)-data(2).data(:,26); %check to be sure
%Chin_referenced = [Chin_referenced_1 ; Chin_referenced_2];
%clear Chin_referenced_1
%clear Chin_referenced_2

%Key for figures
%if starts with '1__' then it is for reference/checking
%if starts with '2__' then it is prelim analysis (spectrogram, stairs)
%if starts with '3__' then it is more advanced (quantifying PSD, newer stuff)
clf
tic
%% change next four variables
epochs = 1149; %import Index_PTX up to this #, rest we don't have LFPs for, try and automate epoch detection - refer to the evalepochs function for this line 21-23
PT = 2;
load('IndexPT2'); %temporary
IndexPT2 = IndexPT2(1:epochs,:); 
%line 68 change too
%line 495 is where you input lfp01,12,23 for downstream calculations,
%graphs towards end as well
%MUST CHECK YOUR LFP CHANNELS! & Chin_reference etc too
load(['Patient_',num2str(PT),'_Sleep_LFP.mat'])
allLFPs = data.data(1:(epochs*30720),31:34); %need to restrict LFP's up front
lfplead0 = allLFPs(:,1);
lfplead1 = allLFPs(:,2);
lfplead2 = allLFPs(:,3);
lfplead3 = allLFPs(:,4);
clear allLfPs

figure(100) %plot to visualize raw LFP's
subplot(4,1,1)
plot(lfplead3)
title('Raw LFPs - lfplead3')
subplot(4,1,2)
plot(lfplead2)
title('Raw LFPs - lfplead2')
subplot(4,1,3)
plot(lfplead1)
title('Raw LFPs - lfplead1')
subplot(4,1,4)
plot(lfplead0)
title('Raw LFPs - lfplead0')
%% Stair plots
%load(['IndexPT',num2str(PT),'.mat']); % still need to work on, unclear how to set equal to something
%['IndexPT',num2str(PT)] = IndexPTX
IndexPTX = IndexPT2; %temporary
% for spectrogram stairs, x axis is per hour (time, fs depedendent)
IndexPTX_score = IndexPT2(:,3);
IndexPTX_score( IndexPTX_score==0 )=90; % Wake
IndexPTX_score( IndexPTX_score==1 )=50;  % N1
IndexPTX_score( IndexPTX_score==2 )=40;  % N2
IndexPTX_score( IndexPTX_score==3 )=35;  % N3
IndexPTX_score( IndexPTX_score==5 )=80;  % REM
IndexPTX_1 = IndexPTX(:,1)/(1024*60*60); % adjust scale to match that of 'spectrogram' output
IndexPTX_2 = IndexPTX(:,2)/(1024*60*60);
IndexPTX = [IndexPTX_1,IndexPTX_2,IndexPTX_score];

clear IndexPTX_1 
clear IndexPTX_2 
clear IndexPTX_score

% for EMG stairs, x axis is per sample
IndexPTX_score_sm = IndexPT2(:,3); %score_sm keeps sleep scoring small numbers
IndexPTX_score_sm( IndexPTX_score_sm==0 )=4.99; % Wake
IndexPTX_score_sm( IndexPTX_score_sm==1 )=2.99;  % N1
IndexPTX_score_sm( IndexPTX_score_sm==2 )=1.99;  % N2
IndexPTX_score_sm( IndexPTX_score_sm==3 )=0.99;  % N3
IndexPTX_score_sm( IndexPTX_score_sm==5 )=3.99;  % REM
IndexPT2_EMG = [IndexPT2(:,1),IndexPT2(:,2),IndexPTX_score_sm]; %xaxis is by sample and y is small numbers

% for sample number on X and big stairs
IndexPTX_persample = [IndexPT2_EMG(:,1),IndexPT2_EMG(:,2),IndexPTX(:,3)];

clear IndexPTX_score_sm

% plot check to make sure only difference between Index's are scales
figure(101)
subplot(2,1,1)
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
subplot(2,1,2)
stairs(IndexPT2_EMG(:,1), IndexPT2_EMG(:,3),'LineWidth',2, 'Color','b' );

%% lfp0
signal = lfplead0;

outData = evalGoodEpochsLFP(signal, 1024, 30);
eventDet = outData.events; %boolean vector, use for logical indexing
Boolean0 = eventDet;
signal_block = zeros(epochs,30720); 
for k = 1:epochs; 
    signal_block(k,:) = signal( 1+(k-1)*30*1024 : k*30*1024);
end

eventDet_block = repmat(eventDet,1,30720); %extends boolean array to match dimensions
signal_block(~eventDet_block) = 0; %logical index
signal_blockT = signal_block';
signal_2 = reshape(signal_blockT,[(epochs*30720),1]);

tubs = 10000; % notch filter for 60 Hz
f1 = 0.1165;
f2 = 0.1175;
notch = fir1(tubs,[f1 f2], 'stop'); % needs adjusting
signal_3 = filter(notch, 1, signal_2); % signal_3 is post notch

[b,a] = butter(2,1/512,'high');
signal_3B=filter(b,a,signal_3); % signal_3B is post notch&butter
signal_3_lead0 = signal_3B;

%spectrogram(signal_3B,700000,10000 , 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%spectrogram(signal_3B,700000,10000 , 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%colorbar off
%hold on
%plot(IndexPTX,'LineWidth',2, 'Color','r' );

clear signal
clear signal_block
clear signal_blockT
clear signal_2
clear tubs
clear f1
clear f2
clear a
clear b
clear signal_3
clear signal_block_3
clear signal_3B
clear outData
%% lfp1

signal = lfplead1;

outData = evalGoodEpochsLFP(signal, 1024, 30);
eventDet = outData.events; %boolean vector, use for logical indexing
Boolean1 = eventDet;
signal_block = zeros(epochs,30720); 
for k = 1:epochs; 
    signal_block(k,:) = signal( 1+(k-1)*30*1024 : k*30*1024);
end

eventDet_block = repmat(eventDet,1,30720); %extends boolean array to match dimensions
signal_block(~eventDet_block) = 0; %logical index
signal_blockT = signal_block';
signal_2 = reshape(signal_blockT,[(epochs*30720),1]);

tubs = 10000; % notch filter for 60 Hz
f1 = 0.1165;
f2 = 0.1175;
notch = fir1(tubs,[f1 f2], 'stop'); % needs adjusting
signal_3 = filter(notch, 1, signal_2); % signal_3 is post notch

[b,a] = butter(2,1/512,'high');
signal_3B = filter(b,a,signal_3); % signal_3B is post notch&butter
signal_3_lead1 = signal_3B;

%spectrogram(signal_3B,700000,10000 , 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%spectrogram(signal_3_lead1,100000, 60000, 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%colorbar off
%hold on
%stairs(IndexPTX,'LineWidth',2, 'Color','r' );

clear signal
clear signal_block
clear signal_blockT
clear signal_2
clear tubs
clear f1
clear f2
clear a
clear b
clear signal_3
clear signal_block_3
clear signal_3B
clear outData
%% lfp2

signal = lfplead2;

outData = evalGoodEpochsLFP(signal, 1024, 30);
eventDet = outData.events; %boolean vector, use for logical indexing
Boolean2 = eventDet;
signal_block = zeros(epochs,30720); 
for k = 1:epochs; 
    signal_block(k,:) = signal( 1+(k-1)*30*1024 : k*30*1024);
end

eventDet_block = repmat(eventDet,1,30720); %extends boolean array to match dimensions
signal_block(~eventDet_block) = 0; %logical index
signal_blockT = signal_block';
signal_2 = reshape(signal_blockT,[(epochs*30720),1]);

tubs = 10000; % notch filter for 60 Hz
f1 = 0.1165;
f2 = 0.1175;
notch = fir1(tubs,[f1 f2], 'stop'); % needs adjusting
signal_3 = filter(notch, 1, signal_2); % signal_3 is post notch

[b,a] = butter(2,1/512,'high');
signal_3B = filter(b,a,signal_3); % signal_3B is post notch&butter
signal_3_lead2 = signal_3B;

%spectrogram(signal_3B,700000,10000 , 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%spectrogram(signal_3_lead2,100000, 60000, 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%colorbar off
%hold on
%stairs(IndexPTX,'LineWidth',2, 'Color','r' );

clear signal
clear signal_block
clear signal_blockT
clear signal_2
clear tubs
clear f1
clear f2
clear a
clear b
clear signal_3
clear signal_block_3
clear signal_3B
clear outData
%% lfp3

signal = lfplead3;

outData = evalGoodEpochsLFP(signal, 1024, 30);
eventDet = outData.events; %boolean vector, use for logical indexing
Boolean3 = eventDet;
signal_block = zeros(epochs,30720); 
for k = 1:epochs; 
    signal_block(k,:) = signal( 1+(k-1)*30*1024 : k*30*1024);
end

eventDet_block = repmat(eventDet,1,30720); %extends boolean array to match dimensions
signal_block(~eventDet_block) = 0; %logical index
signal_blockT = signal_block';
signal_2 = reshape(signal_blockT,[(epochs*30720),1]);

tubs = 10000; % notch filter for 60 Hz
f1 = 0.1165;
f2 = 0.1175;
notch = fir1(tubs,[f1 f2], 'stop'); % needs adjusting
signal_3 = filter(notch, 1, signal_2); % signal_3 is post notch

fs = 1024;
[b,a] = butter(2,1/512,'high');
signal_3B = filter(b,a,signal_3); % signal_3B is post notch&butter
signal_3_lead3 = signal_3B;

%spectrogram(signal_3B,700000,10000 , 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%spectrogram(signal_3_lead3,100000, 60000, 0:100, fs, 'yaxis')
%title(['Spectrogram; Window=700000; Overlap=10000; LFP0 PT',num2str(PT)]);
%colorbar off
%hold on
%stairs(IndexPTX,'LineWidth',2, 'Color','r' );

clear signal
clear signal_block
clear signal_blockT
clear signal_2
clear tubs
clear f1
clear f2
clear a
clear b
clear signal_3
clear signal_block_3
clear signal_3B
clear outData
%%
fB0 = find(Boolean0==0);
fB1 = find(Boolean1==0);
fB2 = find(Boolean2==0);
fB3 = find(Boolean3==0);
fB_all = [fB0; fB1; fB2; fB3];
ufB_all = unique(fB_all); % unique epoches to index out
Boolean_master = ones(epochs,1);
Boolean_master(ufB_all) = 0;
FbT = find(Boolean_master==0); % so use BooleanT to index out
if ufB_all ~= FbT;
    disp('ERROR');
end

clear Boolean0
clear Boolean1
clear Boolean2
clear Boolean3
clear Fb0
clear Fb1
clear Fb2
clear Fb3
clear FbC
clear FbT
clear uFbC

signalblock_3_lead0 = zeros(epochs,30720); 
for k = 1:epochs;
    signalblock_3_lead0(k,:) = signal_3_lead0( 1+(k-1)*30*1024 : k*30*1024);
end
signalblock_3_lead1 = zeros(epochs,30720); 
for k = 1:epochs; 
    signalblock_3_lead1(k,:) = signal_3_lead1( 1+(k-1)*30*1024 : k*30*1024);
end
signalblock_3_lead2 = zeros(epochs,30720);
for k = 1:epochs; 
    signalblock_3_lead2(k,:) = signal_3_lead2( 1+(k-1)*30*1024 : k*30*1024);
end
signalblock_3_lead3 = zeros(epochs,30720);
for k = 1:epochs; 
    signalblock_3_lead3(k,:) = signal_3_lead3( 1+(k-1)*30*1024 : k*30*1024);
end

%Need below to match Boolean output from eval.epoch across each LFP
Boolean_master_block = repmat(eventDet,1,30720);
signalblock_3_lead0(~Boolean_master_block) = 0;
signalblock_3_lead1(~Boolean_master_block) = 0;
signalblock_3_lead2(~Boolean_master_block) = 0;
signalblock_3_lead3(~Boolean_master_block) = 0;

fBsb1 = find(signalblock_3_lead1==0);
fBsb3 = find(signalblock_3_lead0==0); 
if fBsb1 ~= fBsb3;
    disp('ERROR');
end
clear fBsb1
clear fBsb3
%%
BlockT0 = signalblock_3_lead0';
Signal0 = reshape(BlockT0,[(epochs*30720),1]);
BlockT1 = signalblock_3_lead1';
Signal1 = reshape(BlockT1,[(epochs*30720),1]);
BlockT2 = signalblock_3_lead2';
Signal2 = reshape(BlockT2,[(epochs*30720),1]);
BlockT3 = signalblock_3_lead3';
Signal3 = reshape(BlockT3,[(epochs*30720),1]);

figure(103) %plot to visualize raw LFP's vs pre-processed
subplot(4,1,1)
plot(lfplead3)
hold on
plot(Signal3)
title('Raw vs. Pre-processed LFPs - lfplead3')
subplot(4,1,2)
plot(lfplead2)
hold on
plot(Signal2)
title('Raw vs. Pre-processed LFPs -  lfplead2')
subplot(4,1,3)
plot(lfplead1)
hold on
plot(Signal1)
title('Raw vs. Pre-processed LFPs - lfplead1')
subplot(4,1,4)
plot(lfplead0)
hold on
plot(Signal0)
title('Raw vs. Pre-processed LFPs - lfplead0')

lfp01 = Signal0 - Signal1;
lfp12 = Signal1 - Signal2;
lfp23 = Signal2 - Signal3;
%%
figure(200)
subplot(4,1,4) 
spectrogram(signal_3_lead0,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP0 PT' num2str(PT)]);
colorbar off
subplot(4,1,3) 
spectrogram(signal_3_lead1,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP1 PT' num2str(PT)]);
colorbar off
subplot(4,1,2) 
spectrogram(signal_3_lead2,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP2 PT' num2str(PT)]);
colorbar off
subplot(4,1,1) 
spectrogram(signal_3_lead3,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP3 PT' num2str(PT)]);
colorbar off
%%
figure(201)
subplot(3,1,3) 
spectrogram(lfp01,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP0-LFP1 PT' num2str(PT)]);
colorbar off
subplot(3,1,2) 
spectrogram(lfp12,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP1-LFP2 PT' num2str(PT)]);
colorbar off
subplot(3,1,1) 
spectrogram(lfp23,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP2-LFP3 PT' num2str(PT)]);
colorbar off

% can use to optimize window length
%N = length(lfp01);
%figure(115)
%for len = 100000:100000:N; %starts at first numb, increased by next till N
%    spectrogram(lfp01,len,0, 0:100, fs, 'yaxis')
%    hold on
%    stairs(IndexPTX,'LineWidth',2, 'Color','r' );    
%    title(sprintf('Hamming Window Size :: %d', len))
%    pause();
%end

%% Goal is to take some stretches of a given sleep stage and compare them 
%Literature supports using 4-10, 13-30, and 60-80 Hz bands as these have
%greatest differences in ON vs OFF DA-drug 

%% Use for comparison of specific bouts of PSG stages
%Time1_start_epoch = ; %Wake
%Time1_end_epoch = ;
%Time2_start_epoch = ; %REM 
%Time2_end_epoch = ;
%Time3_start_epoch = ; %N1 
%Time3_end_epoch = ;
%Time4_start_epoch = ; %N2
%Time4_end_epoch = ;
%Time5_start_epoch = ; %N3 
%Time5_end_epoch = ;

%signal_ID = lfp12;
%Time1 = signal_ID((Time1_start_epoch*1024*30):(Time1_end_epoch*1024*30)); % Wake, epochs 397:561
%Time2 = signal_ID((Time2_start_epoch*1024*30):(Time2_end_epoch*1024*30)); % REM, 
%Time3 = signal_ID((Time3_start_epoch*1024*30):(Time3_end_epoch*1024*30)); % N1, 
%Time4 = signal_ID((Time4_start_epoch*1024*30):(Time4_end_epoch*1024*30)); % N2, epochs 397:561
%Time5 = signal_ID((Time5_start_epoch*1024*30):(Time5_end_epoch*1024*30)); % N3, 

%Wake = Time1;
%REM = Time2;
%N1 = Time3;
%N2 = Time4;
%N3 = Time5;
%% Use for comparison of ALL PSG stages

% input_lfp below is key, where data is coming from, must replicate for
% other referemced lfps

input_lfp = lfp01;
Signal0_block = zeros(epochs,30720);  

for k = 1:epochs; 
    Signal0_block(k,:) = input_lfp( 1+(k-1)*30*1024 : k*30*1024);
end

if input_lfp(30721,1) == Signal0_block(2,1)
    disp('GOOD')
else input_lfp(30721,1) ~= Signal0_block(2,1);
    disp('ERROR');
end

if input_lfp(614401,1) == Signal0_block(21,1)
    disp('GOOD')
else input_lfp(614401,1) ~= Signal0_block(21,1);
    disp('ERROR');
end

Signal0_block_index = [Signal0_block, IndexPTX(:,3)]; %epochs with indexing at end

%% Grouping by sleep score 

Signal0_block_1 = zeros(epochs, 30720); %Wake
%this is probably just a bad way to index
for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 90
     Signal0_block_1(k,:) = input_lfp( 1+(k-1)*30*1024 : k*30*1024);
    end
end
%
Signal0_block_2 = Signal0_block_1; %
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros


if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
Signal0_block_2T = Signal0_block_2';
PTX_Wake = reshape(Signal0_block_2T,[needed_rows,1]); %back in vector form
%^ is it okay to lump everything together? cause pwelch's window might
%overlap on data from separate bouts of a given PSG stage

%% REM, N1, N2, N3
%need to make some kind of if/then statement to determine whether analysis
%proceeds, if no values for a sleep stage need to not compute it
Signal0_block_1 = zeros(epochs, 30720); %REM

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 80
     Signal0_block_1(k,:) = input_lfp( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
Signal0_block_2T = Signal0_block_2';
PTX_REM = reshape(Signal0_block_2T,[needed_rows,1]);

Signal0_block_1 = zeros(epochs, 30720); %N1

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 50
     Signal0_block_1(k,:) = input_lfp( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
Signal0_block_2T = Signal0_block_2';
PTX_N1 = reshape(Signal0_block_2T,[needed_rows,1]);

Signal0_block_1 = zeros(epochs, 30720); %N2

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 40
     Signal0_block_1(k,:) = input_lfp( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
Signal0_block_2T = Signal0_block_2';
PTX_N2 = reshape(Signal0_block_2T,[needed_rows,1]);

Signal0_block_1 = zeros(epochs, 30720); %N3

for k = 1:epochs; 
    if Signal0_block_index(k,30721) == 35
     Signal0_block_1(k,:) = input_lfp( 1+(k-1)*30*1024 : k*30*1024);
    end
end

Signal0_block_2 = Signal0_block_1;
Signal0_block_2( ~any(Signal0_block_2,2), : ) = [];  % delete rows made of zeros

if sum(Signal0_block_2) ~= sum(Signal0_block_1)
    disp('ERROR')
end

[m,n]  = size(Signal0_block_2);
needed_rows = m*n;
Signal0_block_2T = Signal0_block_2';
PTX_N3 = reshape(Signal0_block_2T,[needed_rows,1]);

fs = 1024;
Fs = fs;

Wake = PTX_Wake;
REM = PTX_REM; %No REM Pt9
N1 = PTX_N1;
N2 = PTX_N2;
N3 = PTX_N3; %No N3 Pt9
%% plot total LFP vs LFP indexed by sleep stage
Ry1 = -0.01; %y axis range start
Ry2 = 0.01;
figure(104)
subplot(6,1,1)
plot(input_lfp)
set(gca, 'YLim', [Ry1 Ry2]);
title('LFP - referenced')
subplot(6,1,2)
plot(Wake)
set(gca, 'YLim', [Ry1 Ry2]);
title('Wake')
subplot(6,1,3)
plot(REM)
set(gca, 'YLim', [Ry1 Ry2]);
title('REM')
subplot(6,1,4)
plot(N1)
set(gca, 'YLim', [Ry1 Ry2]);
title('N1')
subplot(6,1,5)
plot(N2)
set(gca, 'YLim', [Ry1 Ry2]);
title('N2')
subplot(6,1,6)
plot(N3)
set(gca, 'YLim', [Ry1 Ry2]);
title('N3')

clear Ry1
clear Ry2

%% plot total LFP pwelch vs LFP indexed by sleep stage pwelch
window = [];
noverlap = [];
f = window; %mathworks tutorial said that its safest to set to 'window'

Ry1 = -100; %y axis range start
Ry2 = 0;
Rx1 = 0; %x axis range start
Rx2 = 30;

figure(105)
subplot(6,1,1)
pwelch(input_lfp,window,noverlap,f,fs)
set(gca, 'XLim', [Rx1 Rx2]);
%set(gca, 'YLim', [Ry1 Ry2]);
title('LFP - referenced')
subplot(6,1,2)
pwelch(Wake,window,noverlap,f,fs)
set(gca, 'XLim', [Rx1 Rx2]);
%set(gca, 'YLim', [Ry1 Ry2]);
title('Wake')
subplot(6,1,3)
pwelch(REM,window,noverlap,f,fs)
set(gca, 'XLim', [Rx1 Rx2]);
%set(gca, 'YLim', [Ry1 Ry2]);
title('REM')
subplot(6,1,4)
pwelch(N1,window,noverlap,f,fs)
set(gca, 'XLim', [Rx1 Rx2]);
%set(gca, 'YLim', [Ry1 Ry2]);
title('N1')
subplot(6,1,5)
pwelch(N2,window,noverlap,f,fs)
set(gca, 'XLim', [Rx1 Rx2]);
%set(gca, 'YLim', [Ry1 Ry2]);
title('N2')
subplot(6,1,6)
pwelch(N3,window,noverlap,f,fs)
set(gca, 'XLim', [Rx1 Rx2]);
%set(gca, 'YLim', [Ry1 Ry2]);
title('N3')

%% 4-10 Hz Band
band_start = 4;
band_end = 10;
band = [band_start band_end];

%Each PSG stage powers
% need to plot something meaningful here
% come back and automate the 4to10 with band_start/end
[Pxx,F] = periodogram(input_lfp,rectwin(length(input_lfp)),length(input_lfp),Fs);
pband0_4to10 = bandpower(Pxx,F,band,'psd');
ptot0_4to10 = bandpower(Pxx,F,'psd');
per_power0_4to10 = 100*(pband0_4to10/ptot0_4to10);

[Pxx,F] = periodogram(Wake,rectwin(length(Wake)),length(Wake),Fs);
pband1_4to10 = bandpower(Pxx,F,band,'psd');
ptot1_4to10 = bandpower(Pxx,F,'psd');
per_power1_4to10 = 100*(pband1_4to10/ptot1_4to10);

[Pxx,F] = periodogram(REM,rectwin(length(REM)),length(REM),Fs);
pband2_4to10 = bandpower(Pxx,F,band,'psd');
ptot2_4to10 = bandpower(Pxx,F,'psd');
per_power2_4to10 = 100*(pband2_4to10/ptot2_4to10);

[Pxx,F] = periodogram(N1,rectwin(length(N1)),length(N1),Fs);
pband3_4to10 = bandpower(Pxx,F,band,'psd');
ptot3_4to10 = bandpower(Pxx,F,'psd');
per_power3_4to10 = 100*(pband3_4to10/ptot3_4to10);

[Pxx,F] = periodogram(N2,rectwin(length(N2)),length(N2),Fs);
pband4_4to10 = bandpower(Pxx,F,band,'psd');
ptot4_4to10 = bandpower(Pxx,F,'psd');
per_power4_4to10 = 100*(pband4_4to10/ptot4_4to10);

[Pxx,F] = periodogram(N3,rectwin(length(N3)),length(N3),Fs);
pband5_4to10 = bandpower(Pxx,F,band,'psd');
ptot5_4to10 = bandpower(Pxx,F,'psd');
per_power5_4to10 = 100*(pband5_4to10/ptot5_4to10);

Per_Power_4to10 = per_power1_4to10+per_power2_4to10+per_power3_4to10+per_power4_4to10+per_power5_4to10;
disp(Per_Power_4to10);

if ptot1_4to10>ptot4_4to10
    disp('ptot1_4to10 larger')
else
    disp('ptot4_4to10 larger')
end

figure(30000)
y_pband_4to10 = [pband0_4to10,pband1_4to10,pband2_4to10,pband3_4to10,pband4_4to10,pband5_4to10];
bar(y_pband_4to10);
grid on
l = cell(1,6);
l{1}='All'; l{2}='Wake'; l{3}='REM'; l{4}='N1'; l{5}='N2'; l{6}='N3';    %funky because of N2
set(gca,'xticklabel', l) 
title(['Power bands for',num2str(band_start),'-',num2str(band_end),'Hz'])

%
window = 30720;
noverlap = 0;
f = window; %mathworks tutorial said that its safest to set to 'window'

figure(210) %subplots for each PSG stage for 'band'
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

figure(211) 
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
%legend('Wake','N1', 'N2')
legend('boxoff')


%% 13- 30 Hz Band
band_start = 13;
band_end = 30;
band = [band_start band_end];

[Pxx,F] = periodogram(input_lfp,rectwin(length(input_lfp)),length(input_lfp),Fs);
pband0_13to30 = bandpower(Pxx,F,band,'psd');
ptot0_13to30 = bandpower(Pxx,F,'psd');
per_power0_13to30 = 100*(pband0_13to30/ptot0_13to30);

[Pxx,F] = periodogram(Wake,rectwin(length(Wake)),length(Wake),Fs);
pband1_13to30 = bandpower(Pxx,F,band,'psd');
ptot1_13to30 = bandpower(Pxx,F,'psd');
per_power1_13to30 = 100*(pband1_13to30/ptot1_13to30);

[Pxx,F] = periodogram(REM,rectwin(length(REM)),length(REM),Fs);
pband2_13to30= bandpower(Pxx,F,band,'psd');
ptot2_13to30 = bandpower(Pxx,F,'psd');
per_power2_13to30 = 100*(pband2_13to30/ptot2_13to30);

[Pxx,F] = periodogram(N1,rectwin(length(N1)),length(N1),Fs);
pband3_13to30 = bandpower(Pxx,F,band,'psd');
ptot3_13to30 = bandpower(Pxx,F,'psd');
per_power3_13to30 = 100*(pband3_13to30/ptot3_13to30);

[Pxx,F] = periodogram(N2,rectwin(length(N2)),length(N2),Fs);
pband4_13to30 = bandpower(Pxx,F,band,'psd');
ptot4_13to30 = bandpower(Pxx,F,'psd');
per_power4_13to30 = 100*(pband4_13to30/ptot4_13to30);

[Pxx,F] = periodogram(N3,rectwin(length(N3)),length(N3),Fs);
pband5_13to30 = bandpower(Pxx,F,band,'psd');
ptot5_13to30 = bandpower(Pxx,F,'psd');
per_power5_13to30 = 100*(pband5_13to30/ptot5_13to30);

Per_Power_13to30 = per_power1_13to30+per_power2_13to30+per_power3_13to30+per_power4_13to30+per_power5_13to30;
disp(Per_Power_13to30);

if ptot1_13to30>ptot4_13to30
    disp('ptot1_13to30 larger')
else
    disp('ptot4_13to30 larger')
end

figure(30001)
y_pband_13to30 = [pband0_13to30,pband1_13to30,pband2_13to30,pband3_13to30,pband4_13to30,pband5_13to30];
bar(y_pband_13to30);
grid on
l = cell(1,6);
l{1}='All'; l{2}='Wake'; l{3}='REM'; l{4}='N1'; l{5}='N2'; l{6}='N3';    %funky because of N2
set(gca,'xticklabel', l) 
title(['Power bands for',num2str(band_start),'-',num2str(band_end),'Hz'])
%
window = 30720;
noverlap = 0;
f = window; %mathworks tutorial said that its safest to set to 'window'

figure(220) %subplots for each PSG stage for 'band'
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

figure(221) 
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
%legend('Wake','N1', 'N2')
legend('boxoff')

%% 60-80 Hz Band
band_start = 60;
band_end = 80;
band = [band_start band_end];

%Each PSG stage powers'
[Pxx,F] = periodogram(input_lfp,rectwin(length(input_lfp)),length(input_lfp),Fs);
pband0_60to80 = bandpower(Pxx,F,band,'psd');
ptot0_60to80 = bandpower(Pxx,F,'psd');
per_power0_60to80 = 100*(pband0_60to80/ptot0_60to80);

[Pxx,F] = periodogram(Wake,rectwin(length(Wake)),length(Wake),Fs);
pband1_60to80= bandpower(Pxx,F,band,'psd');
ptot1_60to80 = bandpower(Pxx,F,'psd');
per_power1_60to80= 100*(pband1_60to80/ptot1_60to80);

[Pxx,F] = periodogram(REM,rectwin(length(REM)),length(REM),Fs);
pband2_60to80 = bandpower(Pxx,F,band,'psd');
ptot2_60to80= bandpower(Pxx,F,'psd');
per_power2_60to80 = 100*(pband2_60to80/ptot2_60to80);

[Pxx,F] = periodogram(N1,rectwin(length(N1)),length(N1),Fs);
pband3_60to80 = bandpower(Pxx,F,band,'psd');
ptot3_60to80 = bandpower(Pxx,F,'psd');
per_power3_60to80 = 100*(pband3_60to80/ptot3_60to80);

[Pxx,F] = periodogram(N2,rectwin(length(N2)),length(N2),Fs);
pband4_60to80= bandpower(Pxx,F,band,'psd');
ptot4_60to80 = bandpower(Pxx,F,'psd');
per_power4_60to80 = 100*(pband4_60to80/ptot4_60to80);

[Pxx,F] = periodogram(N3,rectwin(length(N3)),length(N3),Fs);
pband5_60to80 = bandpower(Pxx,F,band,'psd');
ptot5_60to80 = bandpower(Pxx,F,'psd');
per_power5_60to80 = 100*(pband5_60to80/ptot5_60to80);

Per_Power_60to80 = per_power1_60to80+per_power2_60to80+per_power3_60to80+per_power4_60to80+per_power5_60to80;
disp(Per_Power_60to80);

if ptot1_60to80>ptot4_60to80
    disp('ptot1_60to80 larger')
else
    disp('ptot4_60to80larger')
end

figure(30003)
y_pband_60to80 = [pband0_60to80,pband1_60to80,pband2_60to80,pband3_60to80,pband4_60to80,pband5_60to80];
bar(y_pband_60to80);
grid on
l = cell(1,6);
l{1}='All'; l{2}='Wake'; l{3}='REM'; l{4}='N1'; l{5}='N2'; l{6}='N3';    %funky because of N2
set(gca,'xticklabel', l) 
title(['Power bands for',num2str(band_start),'-',num2str(band_end),'Hz'])

%
window = 30720;
noverlap = 0;
f = window; %mathworks tutorial said that its safest to set to 'window'

figure(230) %subplots for each PSG stage for 'band'
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

figure(231) 
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
%legend('Wake','N1', 'N2')
legend('boxoff')

%% Comparison across two stages
figure(305)
window = 1024;
overlap = 512;
subplot(3,2,1) 
spectrogram(Wake,window, overlap, 4:10, fs, 'yaxis')
title(['Wake LFP01 4-10Hz PT 9' num2str(PT)]);
colorbar off
subplot(3,2,2) 
spectrogram(N3,window, overlap, 4:10, fs, 'yaxis')
title(['N3 LFP01 4-10Hz PT' num2str(PT)]);
colorbar off
subplot(3,2,3) 
spectrogram(Wake,window, overlap, 13:30, fs, 'yaxis')
title(['Wake LFP01 13-30 Hz PT' num2str(PT)]);
colorbar off
subplot(3,2,4) 
spectrogram(N3,window, overlap, 13:30, fs, 'yaxis')
title(['N3 LFP01 13-30 Hz PT' num2str(PT)]);
colorbar off
subplot(3,2,5) 
spectrogram(Wake,window, overlap, 60:80, fs, 'yaxis')
title(['Wake LFP01 60-80 Hz PT' num2str(PT)]);
colorbar off
subplot(3,2,6) 
spectrogram(N3,window, overlap, 60:80, fs, 'yaxis')
title(['N3 LFP01 60-80 Hz PT' num2str(PT)]);
colorbar off
%%
Chin_referenced = (data.data(:,25))-(data.data(:,26));

%%
figure(3970)
subplot(10,1,3) 
spectrogram(lfp01,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP0-LFP1 PT' num2str(PT)]);
colorbar off
subplot(10,1,2) 
spectrogram(lfp12,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP1-LFP2 PT' num2str(PT)]);
colorbar off
subplot(10,1,1) 
spectrogram(lfp23,100000, 60000, 0:100, fs, 'yaxis')
hold on
stairs(IndexPTX(:,1), IndexPTX(:,3),'LineWidth',2, 'Color','r' );
title(['LFP2-LFP3 PT' num2str(PT)]);
colorbar off
% adding in subplots to determine optimal comparison
subplot(10,1,4) 
x = 1;
plot((data.data(:,x))-(data.data(:,x+11)))
hold on
stairs(IndexPTX_persample(:,1), IndexPTX_persample(:,3),'LineWidth',2, 'Color','r' );
title(['PT' num2str(PT) ' : channel number' num2str(x) '- referenced FDC ']);
subplot(10,1,5) 
x = 2;
plot((data.data(:,x))-(data.data(:,x+11)))
hold on
stairs(IndexPTX_persample(:,1), IndexPTX_persample(:,3),'LineWidth',2, 'Color','r' );
title(['PT' num2str(PT) ' : channel number' num2str(x) '- referenced bicep ']);
subplot(10,1,6) 
x = 3;
plot((data.data(:,x))-(data.data(:,x+11)))
hold on
stairs(IndexPTX_persample(:,1), IndexPTX_persample(:,3),'LineWidth',2, 'Color','r' );
title(['PT' num2str(PT) ' : channel number' num2str(x) '- referenced tricep ']);
subplot(10,1,7) 
x = 4;
plot((data.data(:,x))-(data.data(:,x+11)))
hold on
stairs(IndexPTX_persample(:,1), IndexPTX_persample(:,3),'LineWidth',2, 'Color','r' );
title(['PT' num2str(PT) ' : channel number' num2str(x) '- referenced EDC ']);
subplot(10,1,8) 
x = 25;
plot(data.data(:,x)-data.data(:,(x+1)))
hold on
stairs(IndexPTX_persample(:,1), IndexPTX_persample(:,3),'LineWidth',2, 'Color','r' );
title(['PT' num2str(PT) ' : channel number' num2str(x) ' - referenced chin ']);
subplot(10,1,9) 
x = 27;
plot(data.data(:,x)-data.data(:,(x+1)))
hold on
stairs(IndexPT2_EMG(:,1), IndexPT2_EMG(:,3),'LineWidth',2, 'Color','r' );
set(gca, 'YLim', [0 5]);
title(['PT' num2str(PT) ' : channel number' num2str(x) ' - referenced anterior tibia']);
subplot(10,1,10) 
x = 35;
plot(data.data(:,x))
hold on
stairs(IndexPT2_EMG(:,1), IndexPT2_EMG(:,3),'LineWidth',2, 'Color','r' );
set(gca, 'YLim', [0 5]);
title(['PT' num2str(PT) ' : channel number' num2str(x)]);

%% here add something using the 'pband_band' values to compare across groups

toc