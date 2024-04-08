%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source localization in EEG recordings using frequency domain and time 
% domain. Input data are general visual response at harmonics f0 = 6, 12, 
% 18, 24 Hz and face individuation response at harmonics f1 = 1.2, 2.4, 3.6, 
% 4.8 Hz. 
% Frequency domain: FFT EEG transformation using peaks taken at the harmonic 
% frequencies of general visual response (f0) and face individuation response 
% (f1). 
% Time domain: Band pass filtering at the harmonic frequencies of general 
% visual (f0) and face individuation (f1) responses   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%% Forward model 
load LF3D_Colin27_9mm vol full_dip LF_full_scalp sensBS128 scalp_sens_idxBS128 
sens = sensBS128;
scalp = vol.bnd(1);
brain = vol.bnd(3);
LF = LF_full_scalp(scalp_sens_idxBS128,:);
source_space = size(LF,2)/3;

%% EEG recordings 
% check if the path is correct!
modeldir = 'C:\Users\admin-sm\Documents\MATLAB\data\Data-review-FPVS-oddball\*.mat';
matFiles = dir(modeldir);
load eeg_info % sr and freq for each EEG
load eeg_header
nchan = 128;
elnames = {eeg_header.chanlocs(1:nchan).labels}.';
velecf1_avgview=zeros(128,4);
velecf0_avgview=zeros(128,4);

%% Normalization of lead-field (forward model)
%LF_Norm = LF_normalization(LF);

%% Localization to each EEG file
for i = 1:length(matFiles)

baseFileName = matFiles(i).name;
load(baseFileName)   
data = squeeze(data);
nsess = size(data,1);
n = size(data,3);
fs = [eeg_info(i).fs].';
sessvelec = [];

%% Organize the sessions recorded 
for isess=1:nsess
    sessvelec(:,:,isess) = double(squeeze(data(isess,1:nchan,:)));
end
%% Average the sessions (optional) or process sessions separately
velec = mean(sessvelec,3);

%% Normalization on the averaged sessions
%velec = velec/(norm(velec,'fro'));

%% Frequency transformation
velecf = conj(fft(velec')');
Rf = fs/n; %freq resolution

%% Signal frequency domain visualization
%EEGviewer(abs(velecf), 1/Rf, elnames)

%% Electrodes selection with resolute peaks at frequencies of interests
%[velecf1,velecf0,ixf1,ixf0,LFself1,LFself0] = electrodes_selection(LF,velecf,indf1,indf0);

%% Frequency peaks selection
odd = [eeg_info(i).oddball_freq].';
f1 = [odd odd*2 odd*3 odd*4]; %freq peaks at face individuation response
base = [eeg_info(i).base_freq].';
f0 = [base base*2 base*3 base*4];%freq peaks at general visual response
indf1 = round(f1/Rf+1);
indf0 = round(f0/Rf+1);   

%% Sum up the power spectrum of the frequencies in f1 and f0 for each EEG 
%% file for scalp map view purpose 
velecf1_avgview=velecf1_avgview+velecf(:,indf1);
velecf0_avgview=velecf0_avgview+velecf(:,indf0);

%% Data in frequency
% Harmonics of face individuation f1(1.2 2.4 3.6 4,8 Hz)
Vreh = real(velecf(:,indf1));
Vimh = imag(velecf(:,indf1));
Yf1h = [Vreh Vimh];

% Harmonics of general visual response f0(6 12 18 24 Hz)
Vre0h = real(velecf(:,indf0));
Vim0h = imag(velecf(:,indf0));
Yf0h = [Vre0h Vim0h];

%% Source localization 

%% in freq, face individuation response f1 = 1.2Hz, 2.4Hz, 3.6Hz, 4,8Hz
%[freq.inddipsbrf1h,freq.orrsbrf1h,freq.ampsbrf1h,freq.GOFsbrf1h] = SBRfreqloc9_2_fin_corr(Yf1h,LF,0.01);
[freq.inddipsbrf1h,freq.orrsbrf1h,freq.ampsbrf1h,freq.GOFsbrf1h] = mySBRR1_AltOpt(Yf1h,LF,0.01);

%% in freq, general visual response f0 = 6Hz, 12Hz, 18Hz, 24Hz
%[freq.inddipsbrf0h,freq.orrsbrf0h,freq.ampsbrf0h,freq.GOFsbrf0h,Xhat0] = SBRfreqloc9_2_fin_corr(Yf0h,LF,0.01);
[freq.inddipsbrf0h,freq.orrsbrf0h,freq.ampsbrf0h,freq.GOFsbrf0h] = mySBRR1_AltOpt(Yf0h,LF,0.01);

% Set the number of sources (SBR reference) to be localized by (t)rap algs in time domain 
nsources_f1 = numel(freq.inddipsbrf1h);
nsources_f0 = numel(freq.inddipsbrf0h);

%% Data in time

%% x axis adjustment to integer n of periods,
nw = fix(n/fs);
velec_adjf1 = velec(:,1:nw*fs);
velec_adjf0 = velec(:,1:nw*fs);

%% Band pass filtering
% f1 [1.2 2.4 3.6 4.8] Hz
[b,a] = butter(1,[1 5]/(fs/2));
velecavgf1 = mean(reshape(filtfilt(b,a,velec_adjf1')',nchan,fs,length(velec_adjf1)/fs),3);
% f0 [6 12 18 24] Hz
[b,a] = butter(1,[5 25]/(fs/2));
velecavgf0 = mean(reshape(filtfilt(b,a,velec_adjf0')',nchan,fs,length(velec_adjf0)/fs),3);

%% in time, face individuation response f1 = 1.2Hz, 2.4Hz, 3.6Hz, 4,8Hz
[time.inddiprapf1th,time.orrrapf1th,time.amprapf1th,time.GOFrapf1th] = allmusic(velecavgf1,LF,1,'rap',1,nsources_f1);
[time.inddiptrapf1th,time.orrtrapf1th,time.amptrapf1th,time.GOFtrapf1th] = allmusic(velecavgf1,LF,1,'trap',1,nsources_f1);

%% in time, general visual response f0 = 6Hz, 12Hz, 18Hz, 24Hz
[time.inddiprapf0th,time.orrrapf0th,time.amprapf0th,time.GOFrapf0th] = allmusic(velecavgf0,LF,1,'rap',1,nsources_f0);
[time.inddiptrapf0th,time.orrtrapf0th,time.amptrapf0th,time.GOFtrapf0th] = allmusic(velecavgf0,LF,1,'trap',1,nsources_f0);

%% Save localization for each EEG file
save(['EEGloc_fvps_SBR_freq_T_RAP_time_f1_f0_',num2str(nchan),'el_',num2str(source_space),'s_','file_',num2str(i)],'baseFileName','freq','time')

end

%% Scalp map view representing the average from EEG files, summing up the 
%% power spectrum of the frequencies in f1 and f0 across EEG files
figure,
scalpmap_mesh(double(abs(sum(velecf1_avgview,2))),sens.elecpos,scalp.pnt,scalp.tri)
axis off, view(-90,0), colorbar off

figure,
scalpmap_mesh(double(abs(sum(velecf0_avgview,2))),sens.elecpos,scalp.pnt,scalp.tri)
axis off, view(-90,0), colorbar off