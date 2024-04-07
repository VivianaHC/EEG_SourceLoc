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
modeldir = 'C:\Users\admin-sm\Documents\MATLAB\EEG-Source-localization\data\Data-review-FPVS-oddball\*.mat';
matFiles = dir(modeldir);
load eeg_info % sr and freq for each EEG
load eeg_header
nchan = 128;
elnames = {header.chanlocs(1:nchan).labels}.';

%% Normalization of lead-field (forward model)
%LF_Norm = LF_normalization(LF);

%% Localization to each EEG file
for i = 1:length(matFiles)

baseFileName = matFiles(i).name;
load(baseFileName)   
data = squeeze(data);
nsess = size(data,1);
n = size(data,3);
fs = [eeg_info(i).fs].';%256/250
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

%% frequency peaks selection
odd = [eeg_info(i).oddball_freq].';
f1 = [odd odd*2 odd*3 odd*4]; %frequencies at face individuation response
base = [eeg_info(i).base_freq].';
f0 = [base base*2 base*3 base*4];%frequencies at general visual response
indf1 = round(f1/Rf+1);
indf0 = round(f0/Rf+1);   

% %% Electrodes selection with resolute peaks at oddball harmonics 1.2,2.4,3.6,4.8 Hz
% q=1;
% velecf1=[];
% ixf1=[];
% LFself1=[];
% for j=1:size(velecf,1)
%     if  (abs(velecf(j,indf1(1)-1)) < abs(velecf(j,indf1(1)))) && (abs(velecf(j,indf1(1))) > abs(velecf(j,indf1(1)+1))) ...
%         || (abs(velecf(j,indf1(2)-1)) < abs(velecf(j,indf1(2)))) && (abs(velecf(j,indf1(2))) > abs(velecf(j,indf1(2)+1))) ...
%         || (abs(velecf(j,indf1(3)-1)) < abs(velecf(j,indf1(3)))) && (abs(velecf(j,indf1(3))) > abs(velecf(j,indf1(3)+1)))
%         velecf1(q,:)=velecf(j,:);
%         ixf1(q,1)=j;
%         q=q+1;
%     end
% end
% freq.ixf1=ixf1;
% LFself1=LF(ixf1,:);
% %LFself1=LF_Norm(ixf1,:);

% %% Electrodes selection with resolute peaks at visual response harmonics 6,12,18,24 Hz
% q=1;
% velecf0=[];
% ixf0=[];
% LFself0=[];
% for j=1:size(velecf,1)
%     if  (abs(velecf(j,indf0(1)-1)) < abs(velecf(j,indf0(1)))) && (abs(velecf(j,indf0(1))) > abs(velecf(j,indf0(1)+1))) ...
%         || (abs(velecf(j,indf0(2)-1)) < abs(velecf(j,indf0(2)))) && (abs(velecf(j,indf0(2))) > abs(velecf(j,indf0(2)+1))) ...
%         || (abs(velecf(j,indf0(3)-1)) < abs(velecf(j,indf0(3)))) && (abs(velecf(j,indf0(3))) > abs(velecf(j,indf0(3)+1)))
%         velecf0(q,:)=velecf(j,:);
%         ixf0(q,1)=j;
%         q=q+1;
%     end
% end
% freq.ixf0=ixf0;
% LFself0=LF(ixf0,:);
% %LFself0=LF_Norm(ixf0,:);

%% GOF 
interva=indf1(1)-4:indf1(end)+4;
expGOFf1=1-(norm(mean(abs(velecf(:,setdiff(interva,indf1))),2),'fro')/norm(mean(abs(velecf(:,indf1)),2),'fro'))^2;
expGOFf1=min(expGOFf1,0.99);

interva=indf0(1)-4:indf0(end)+4;
expGOFf0=1-(norm(mean(abs(velecf(:,setdiff(interva,indf0))),2),'fro')/norm(mean(abs(velecf(:,indf0)),2),'fro'))^2;
expGOFf0=min(expGOFf0,0.99);

%% Data in frequency
% Harmonics of face individuation f1(1.2 2.4 3.6 4,8 Hz)
Vreh=real(velecf(:,indf1));
Vimh=imag(velecf(:,indf1));
Yf1h=[Vreh Vimh];

% Harmonics of general visual response f0(6 12 18 24 Hz)
Vre0h=real(velecf(:,indf0));
Vim0h=imag(velecf(:,indf0));
Yf0h=[Vre0h Vim0h];

%% Scalp map view, sum of power spectrum over the frequencies in f1 and f0
velecf1_rep=zeros(128,4);
velecf0_rep=zeros(128,4);
velecf1_rep=velecf1_rep+velecf(:,indf1);
velecf0_rep=velecf0_rep+velecf(:,indf0);
figure,
scalpmap_mesh(double(abs(sum(velecf1_rep,2))),sens.elecpos,scalp.pnt,scalp.tri)
axis off
view(-90,0)
colorbar off
% saveas(1,['scalpmap_avg_f1h.fig'])
% saveas(1,['ScalpMap_File' num2str(i) '_' num2str(baseFileName) '.png'])
%close
figure,
scalpmap_mesh(double(abs(sum(velecf0_rep,2))),sens.elecpos,scalp.pnt,scalp.tri)
axis off
view(-90,0)
colorbar off
%saveas(1,['scalpmap_avg_f0h.fig'])

%% Source localization 

%% in freq, oddball harms 1.2Hz, 2.4Hz,3.6Hz, 4,8Hz
%[freq.inddiprapf1h,freq.orrrapf1h,freq.amprapf1h,freq.GOFrapf1h] = allmusic(Yfh,LF,1,'rap',1,'mdl');
%[freq.inddiptrapf1h,freq.orrtrapf1h,freq.amptrapf1h,freq.GOFtrapf1h] = allmusic(Yfh,LF,1,'trap',1,'mdl');
%[freq.inddipolsf1h,freq.orrolsf1h,freq.ampolsf1h,freq.GOFolsf1h] = OLSfreqloc10_fin(Yf1h,LF_Norm,expGOFf1);
%[freq.inddipsbrf1h,freq.orrsbrf1h,freq.ampsbrf1h,freq.GOFsbrf1h] = SBRfreqloc9_2_fin_corr(Yf1h,LF,0.01);
[freq.inddipsbrf1h,freq.orrsbrf1h,freq.ampsbrf1h,freq.GOFsbrf1h] = mySBRR1_AltOpt(Yf1h,LF,0.01);

%% in freq, harms of visual base response 6Hz, 12Hz, 18Hz, 24Hz
%[freq.inddiprapf0h,freq.orrrapf0h,freq.amprapf0h,freq.GOFrapf0h] = allmusic(Yf0h,LF,1,'rap',1,'mdl');
%[freq.inddiptrapf0h,freq.orrtrapf0h,freq.amptrapf0h,freq.GOFtrapf0h] = allmusic(Yf0h,LF,1,'trap',1,'mdl');
%[freq.inddipolsf0h,freq.orrolsf0h,freq.ampolsf0h,freq.GOFolsf0h] = OLSfreqloc10_fin(Yf0h,LF_Norm,expGOFf0);
%[freq.inddipsbrf0h,freq.orrsbrf0h,freq.ampsbrf0h,freq.GOFsbrf0h,Xhat0] = SBRfreqloc9_2_fin_corr(Yf0h,LF,0.01);
[freq.inddipsbrf0h,freq.orrsbrf0h,freq.ampsbrf0h,freq.GOFsbrf0h] = mySBRR1_AltOpt(Yf0h,LF,0.01);

nsources_f1=numel(freq.inddipsbrf1h);
nsources_f0=numel(freq.inddipsbrf0h);

%% Data in time
% velecf1=velec(ixf1,:);
% f1=length(ixf1);
% 
% velecf0=velec(ixf0,:);
% mf0=length(ixf0);

%% x axis adjustment to integer n of periods,
nw=fix(n/fs);
velec_adjf1=velec(:,1:nw*fs);
velec_adjf0=velec(:,1:nw*fs);

%% band pass filtering
% band pass harms oddball(f1) [1.2 2.4 3.6 4.8]Hz
[b,a]=butter(1,[1 5]/(fs/2));
velecavgf1=mean(reshape(filtfilt(b,a,velec_adjf1')',nchan,fs,length(velec_adjf1)/fs),3);
% band pass harms visual response(f0) [6 12 18 24]Hz
[b,a]=butter(1,[5 25]/(fs/2));
velecavgf0=mean(reshape(filtfilt(b,a,velec_adjf0')',nchan,fs,length(velec_adjf0)/fs),3);

%% in time, oddball harms 1.2Hz, 2.4Hz, 3.6Hz, 4.8Hz
%[time.inddiprapf1th,time.orrrapf1th,time.amprapf1th,time.GOFrapf1th] = allmusic(velecavgf1,LF,1,'rap',1,nsources_f1);
%[time.inddiptrapf1th,time.orrtrapf1th,time.amptrapf1th,time.GOFtrapf1th] = allmusic(velecavgf1,LF,1,'trap',1,nsources_f1);

%[time.inddiprapf1th,time.orrrapf1th,time.amprapf1th,time.GOFrapf1th] = allmusic(velecavgf1,LF_Norm,1,'rap',1,'mdl');
%[time.inddiptrapf1th,time.orrtrapf1th,time.amptrapf1th,time.GOFtrapf1th] = allmusic(velecavgf1,LF_Norm,1,'trap',1,'mdl');
%[time.Iolsf1th,time.orrolsf1th,time.ampolsf1th,time.GOFolsf1th] = OLSfreqloc10_fin(velecavgf1,LF,expGOFf1); 
%[time.Isbrf1h,time.orrsbrf1th,time.ampsbrf1th,time.GOFsbrf1th] = SBRfreqloc9_2_fin_corr(velecavgf1,LF,0.01);

%% in time, harms of visual base response 6Hh, 12Hz, 18Hz, 24Hz
%[time.inddiprapf0th,time.orrrapf0th,time.amprapf0th,time.GOFrapf0th] = allmusic(velecavgf0,LF,1,'rap',1,nsources_f0);
%[time.inddiptrapf0th,time.orrtrapf0th,time.amptrapf0th,time.GOFtrapf0th] = allmusic(velecavgf0,LF,1,'trap',1,nsources_f0);

%[time.inddiprapf0th,time.orrrapf0th,time.amprapf0th,time.GOFrapf0th] = allmusic(velecavgf0,LF_Norm,1,'rap',1,'mdl');
%[time.inddiptrapf0th,time.orrtrapf0th,time.amptrapf0th,time.GOFtrapf0th] = allmusic(velecavgf0,LF_Norm,1,'trap',1,'mdl');
%[time.Iolsf0th,time.orrolsf0th,time.ampolsf0th,time.GOFolsf0th] = OLSfreqloc10_fin(velecavgf0,LF,expGOFf0); 
%[time.Isbrf0th,time.orrsbrf0th,time.ampsbrf0th,time.GOFsbrf0th] = SBRfreqloc9_2_fin_corr(velecavgf0,LF,0.01);

save(['EEGloc_fvps_128el_2634s_newSBR_freq_f1_f0_128el_2634s_file_',num2str(nchan),'el_',num2str(source_space),'s_','file_',num2str(i)],'baseFileName','freq')

end
%% Histogram
%% freq
source_space=2634;
histf1_sbr=zeros(1,source_space);
histf0_sbr=zeros(1,source_space);
ampf1_sbr=zeros(source_space,4);
ampf0_sbr=zeros(source_space,4);

histf1_ols=zeros(1,source_space);
histf0_ols=zeros(1,source_space);
ampf1_ols=zeros(source_space,4);
ampf0_ols=zeros(source_space,4);

for ifiles=1:129
    load(['EEGloc_fvps_128el_2634s_newSBR_freq_f1_f0_128el_2634s_file_128el_2634s_file_',num2str(ifiles),'.mat'])
    histf1_sbr(freq.inddipsbrf1h)=histf1_sbr(freq.inddipsbrf1h)+1;
    histf0_sbr(freq.inddipsbrf0h)=histf0_sbr(freq.inddipsbrf0h)+1;
    ampf1_sbr(freq.inddipsbrf1h,:)=ampf1_sbr(freq.inddipsbrf1h,:)+((freq.ampsbrf1h(:,1:4).^2)+(freq.ampsbrf1h(:,5:8).^2)).^(1/2);
    ampf0_sbr(freq.inddipsbrf0h,:)=ampf0_sbr(freq.inddipsbrf0h,:)+((freq.ampsbrf0h(:,1:4).^2)+(freq.ampsbrf0h(:,5:8).^2)).^(1/2);
    
    histf1_ols(freq.inddipolsf1h)=histf1_ols(freq.inddipolsf1h)+1;
    histf0_ols(freq.inddipolsf0h)=histf0_ols(freq.inddipolsf0h)+1;
    ampf1_ols(freq.inddipolsf1h,:)=ampf1_ols(freq.inddipolsf1h,:)+((freq.ampolsf1h(:,1:4).^2)+(freq.ampolsf1h(:,5:8).^2)).^(1/2);
    ampf0_ols(freq.inddipolsf0h,:)=ampf0_ols(freq.inddipolsf0h,:)+((freq.ampolsf0h(:,1:4).^2)+(freq.ampolsf0h(:,5:8).^2)).^(1/2);   
end
% powerf1=sum((ampf1.*ampf1)');
% figure('Name','Hist SBR f1f'),bar(histf1);
% figure('Name','Power SBR f1f'),bar(powerf1);

%% time
histf1t_rap=zeros(1,source_space);
histf0t_rap=zeros(1,source_space);
ampf1t_rap=zeros(source_space,256);
ampf0t_rap=zeros(source_space,256);

histf1t_trap=zeros(1,source_space);
histf0t_trap=zeros(1,source_space);
ampf1t_trap=zeros(source_space,256);
ampf0t_trap=zeros(source_space,256);

for ifiles=1:129
load(['EEGloc_fvps_128el_2634s_T_RAP_MUSIC_time_f1_f0_128el_2634s_file_128el_2634s_file_',num2str(ifiles),'.mat'])
    histf1t_rap(time.inddiprapf1th)=histf1t_rap(time.inddiprapf1th)+1; 
    histf0t_rap(time.inddiprapf0th)=histf0t_rap(time.inddiprapf0th)+1; 
    histf1t_trap(time.inddiptrapf1th)=histf1t_trap(time.inddiptrapf1th)+1; 
    histf0t_trap(time.inddiptrapf0th)=histf0t_trap(time.inddiptrapf0th)+1;
    ampf1t_r=[]; ampf0t_r=[]; ampf1t_t=[]; ampf0t_t=[];
    if size(time.amprapf1th,2)<=250
        for l=1:size(time.amprapf1th,1)
            ampf1t_r(l,:)=interp1(linspace(0,1,length(time.amprapf1th)),time.amprapf1th(l,:),linspace(0,1,size(ampf1t_rap,2)));
        end
        ampf1t_rap(time.inddiprapf1th,:)=ampf1t_rap(time.inddiprapf1th,:)+((ampf1t_r.^2).^(1/2));
    else
        ampf1t_rap(time.inddiprapf1th,:)=ampf1t_rap(time.inddiprapf1th,:)+((time.amprapf1th.^2).^(1/2));
    end 
    if size(time.amprapf0th,2)<=250
        for l=1:size(time.amprapf0th,1)
            ampf0t_r(l,:) = interp1(linspace(0,1,length(time.amprapf0th)),time.amprapf0th(l,:),linspace(0,1,size(ampf0t_rap,2)));
        end
        ampf0t_rap(time.inddiprapf0th,:)=ampf0t_rap(time.inddiprapf0th,:)+((ampf0t_r.^2).^(1/2));
    else
        ampf0t_rap(time.inddiprapf0th,:)=ampf0t_rap(time.inddiprapf0th,:)+((time.amprapf0th.^2).^(1/2));
    end
    if size(time.amptrapf1th,2)<=250
        for l=1:size(time.amptrapf1th,1)
            ampf1t_t(l,:)=interp1(linspace(0,1,length(time.amptrapf1th)),time.amptrapf1th(l,:),linspace(0,1,size(ampf1t_trap,2)));
        end
            ampf1t_trap(time.inddiptrapf1th,:)=ampf1t_trap(time.inddiptrapf1th,:)+((ampf1t_t.^2).^(1/2));
    else
            ampf1t_trap(time.inddiptrapf1th,:)=ampf1t_trap(time.inddiptrapf1th,:)+((time.amptrapf1th.^2).^(1/2));
    end 
    if size(time.amptrapf0th,2)<=250
       for l=1:size(time.amptrapf0th,1)
           ampf0t_t(l,:)=interp1(linspace(0,1,length(time.amptrapf0th)),time.amptrapf0th(l,:),linspace(0,1,size(ampf0t_trap,2)));
       end
       ampf0t_trap(time.inddiptrapf0th,:)=ampf0t_trap(time.inddiptrapf0th,:)+((ampf0t_t.^2).^(1/2));
    else
       ampf0t_trap(time.inddiptrapf0th,:)=ampf0t_trap(time.inddiptrapf0th,:)+((time.amptrapf0th.^2).^(1/2));
    end
end

%%
thrnbs=6;

%%
figure('Name','SBR1 f1f_dips')
ft_plot_mesh(vol.bnd(:),'facealpha',0.1, 'edgealpha', 0.1);hold on
scatter3(full_dip(histf1_sbr>thrnbs,1),full_dip(histf1_sbr>thrnbs,2),full_dip(histf1_sbr>thrnbs,3),100,histf1_sbr(histf1_sbr>thrnbs)','filled')
colorbar

figure('Name','SBR f0f_dips')
ft_plot_mesh(vol.bnd(:),'facealpha',0.1, 'edgealpha', 0.1);hold on
scatter3(full_dip(histf0_sbr>thrnbs,1),full_dip(histf0_sbr>thrnbs,2),full_dip(histf0_sbr>thrnbs,3),100,histf0_sbr(histf0_sbr>thrnbs)','filled')
colorbar

figure('Name','f1t')
ft_plot_mesh(vol.bnd(:),'facealpha',0.1, 'edgealpha', 0.1);hold on
scatter3(full_dip(histf1t>thrnbs,1),full_dip(histf1t>thrnbs,2),full_dip(histf1t>thrnbs,3),100,histf1t(histf1t>thrnbs)','filled')
colorbar

figure('Name','f0t')
ft_plot_mesh(vol.bnd(:),'facealpha',0.1, 'edgealpha', 0.1);hold on
scatter3(full_dip(histf0t>thrnbs,1),full_dip(histf0t>thrnbs,2),full_dip(histf0t>thrnbs,3),100,histf0t(histf0t>thrnbs)','filled')
colorbar

%% convolution
load LF3D_Colin27_9mm vol full_dip LF_full_scalp
for ip=1:length(full_dip)
    ind=find(sqrt(sum((full_dip-full_dip(ip,:)).^2,2))<13);
    hist2f1_sbr(ip)=mean(histf1_sbr(ind));
    hist2f0_sbr(ip)=mean(histf0_sbr(ind));
    amp2f1_sbr(ip)=sum(sum(ampf1_sbr(ind,:),2));
    amp2f0_sbr(ip)=sum(sum(ampf0_sbr(ind,:),2));  

    hist2f1_ols(ip)=mean(histf1_ols(ind));
    hist2f0_ols(ip)=mean(histf0_ols(ind));
    amp2f1_ols(ip)=sum(sum(ampf1_ols(ind,:),2));
    amp2f0_ols(ip)=sum(sum(ampf0_ols(ind,:),2));  

    hist2f1t_rap(ip)=mean(histf1t_rap(ind));
    hist2f0t_rap(ip)=mean(histf0t_rap(ind));
    amp2f1t_rap(ip)=sum(sum(ampf1t_rap(ind,:),2));
    amp2f0t_rap(ip)=sum(sum(ampf0t_rap(ind,:),2));  
    
    hist2f1t_trap(ip)=mean(histf1t_trap(ind));
    hist2f0t_trap(ip)=mean(histf0t_trap(ind));
    amp2f1t_trap(ip)=sum(sum(ampf1t_trap(ind,:),2));
    amp2f0t_trap(ip)=sum(sum(ampf0t_trap(ind,:),2));  
end

%% interpolation
load('brain.mat')

x=[min(brain.pnt(:,1)):max(brain.pnt(:,1))];
y=[min(brain.pnt(:,2)):max(brain.pnt(:,2))];
z=[min(brain.pnt(:,3)):max(brain.pnt(:,3))];
[X,Y,Z] = meshgrid(x,y,z);

brain2.faces=brain.tri;
brain2.vertices=brain.pnt;

IN=inpolyhedron(brain2,x,y,z);
IN=double(IN);
IN(IN==0)=NaN;

%% SBR
F1=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f1_sbr','natural','linear');
vq1=F1(X,Y,Z);

figure('Name','SBR f1f')
s=slice(X,Y,Z,vq1.*IN,[-100:100],[],[]); axis equal; set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

F1a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f1_sbr','natural','linear');
vq1a=F1a(X,Y,Z);

figure('Name','SBR f1f amplitude')
s=slice(X,Y,Z,vq1a.*IN,[-100:100],[],[]); axis equal; set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 10000000000000])
%view(-40,20)

F0=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f0_sbr','natural','linear');
vq0=F0(X,Y,Z);

figure('Name','SBR f0f')
s=slice(X,Y,Z,vq0.*IN,[-100:100],[],[]);axis equal;set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)
%print -dpng f1fa

F0a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f0_sbr','natural','linear');
vq0a=F0a(X,Y,Z);

figure('Name','SBR f0f amplitude')
s=slice(X,Y,Z,vq0a.*IN,[-100:100],[],[]);axis equal;set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

%% RAP
F1=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f1t_rap','natural','linear');
vq1=F1(X,Y,Z);

figure('Name','RAP f1t')
s=slice(X,Y,Z,vq1.*IN,[-100:100],[],[]); axis equal; set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

F1a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f1t_rap','natural','linear');
vq1a=F1a(X,Y,Z);

figure('Name','RAP f1t amplitude')
s=slice(X,Y,Z,vq1a.*IN,[-100:100],[],[]); axis equal; set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 3])
%view(-40,20)

F0=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f0t_rap','natural','linear');
vq0=F0(X,Y,Z);

figure('Name','RAP f0t')
s=slice(X,Y,Z,vq0.*IN,[-100:100],[],[]);axis equal;set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)
%print -dpng f1fa

F0a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f0t_rap','natural','linear');
vq0a=F0a(X,Y,Z);

figure('Name','RAP f0t amplitude')
s=slice(X,Y,Z,vq0a.*IN,[-100:100],[],[]);axis equal;set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 1500000000000000000])
%view(-40,20)

%% TRAP
F1=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f1t_trap','natural','linear');
vq1=F1(X,Y,Z);

figure('Name','TRAP f1t')
s=slice(X,Y,Z,vq1.*IN,[-100:100],[],[]); axis equal; set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

F1a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f1t_trap','natural','linear');
vq1a=F1a(X,Y,Z);

figure('Name','TRAP f1t amplitude')
s=slice(X,Y,Z,vq1a.*IN,[-100:100],[],[]); axis equal; set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

F0=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f0t_trap','natural','linear');
vq0=F0(X,Y,Z);

figure('Name','TRAP f0t')
s=slice(X,Y,Z,vq0.*IN,[-100:100],[],[]);axis equal;set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

F0a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f0t_trap','natural','linear');
vq0a=F0a(X,Y,Z);

figure('Name','TRAP f0t amplitude')
s=slice(X,Y,Z,vq0a.*IN,[-100:100],[],[]);axis equal;set(s,'Edgecolor','None')
axis off
%set(gca,'clim',[0 30])
%view(-40,20)

%% OLS
% F1=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f1_ols','natural','linear');
% vq1=F1(X,Y,Z);
% 
% figure('Name','OLS f1f')
% s=slice(X,Y,Z,vq1.*IN,[-90:90],[],[]); axis equal; set(s,'Edgecolor','None')
% %set(gca,'clim',[0 30])
% axis off
% view(-40,20)
% 
% F1a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f1_ols','natural','linear');
% vq1a=F1a(X,Y,Z);
% 
% figure('Name','OLS f1f amplitude')
% s=slice(X,Y,Z,vq1a.*IN,[-50:50],[],[]); axis equal; set(s,'Edgecolor','None')
% set(gca,'clim',[0 30])
% axis off
% view(-40,20)
% 
% F0=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),hist2f0_ols','natural','linear');
% vq0=F0(X,Y,Z);
% 
% figure('Name','OLS f0f')
% s=slice(X,Y,Z,vq0.*IN,[-90:90],[],[]);axis equal;set(s,'Edgecolor','None')
% %set(gca,'clim',[0 30])
% axis off
% view(-40,20)
% %print -dpng f1fa
% 
% F0a=scatteredInterpolant(full_dip(:,1),full_dip(:,2),full_dip(:,3),amp2f0_ols','natural','linear');
% vq0a=F0a(X,Y,Z);
% 
% figure('Name','OLS f0f amplitude')
% s=slice(X,Y,Z,vq0a.*IN,[-90:90],[],[]);axis equal;set(s,'Edgecolor','None')
% set(gca,'clim',[0 30])
% axis off
% view(-40,20)

%% Loop for brain slices representation
figure
l=10; m=10;
for i=1:10
    s=slice(X,Y,Z,vq1a.*IN,[l:m],[],[]); axis equal; set(s,'Edgecolor','None')
    l=l+(-10);
    %hold on
    %s=slice(X,Y,Z,vq0a,[m],[],[]);axis equal;set(s,'Edgecolor','None')
    m=m+(10);
end
axis off

%% Results quantification : n dips and GOF 
modeldir='C:\Users\admin-sm\Documents\MATLAB\article\EEGLoc_fpvs_128el_2634_All_sess_avg_MUSICalgs/*.mat';
matFiles=dir(modeldir);

for i=1:length(matFiles)
%for i=1:129   
    %load(['EEGloc_fvps_128el_2634s_SBR_and_OLSfreq_(T)RAPtime_f1_f0_128el_2634s_file_',num2str(i)])
    baseFileName=matFiles(i).name;
    load(baseFileName)
    
    freq_SBR_f1(i)=numel(freq.inddipsbrf1h);
    freq_SBR_f0(i)=numel(freq.inddipsbrf0h);
    freq_OLS_f1(i)=numel(freq.inddipolsf1h);
    freq_OLS_f0(i)=numel(freq.inddipolsf0h);
    
    time_RAP_f1(i)=numel(time.inddiprapf1th);
    time_RAP_f0(i)=numel(time.inddiprapf0th);
    time_TRAP_f1(i)=numel(time.inddiptrapf1th);
    time_TRAP_f0(i)=numel(time.inddiptrapf0th);
     
    freq_SBR_GOFf1(i)=freq.GOFsbrf1h;    
    freq_SBR_GOFf0(i)=freq.GOFsbrf0h;
    freq_OLS_GOFf1(i)=freq.GOFolsf1h;    
    freq_OLS_GOFf0(i)=freq.GOFolsf0h;
    
    time_RAP_GOFf1(i)=time.GOFrapf1th;    
    time_RAP_GOFf0(i)=time.GOFrapf0th;
    time_TRAP_GOFf1(i)=time.GOFtrapf1th;    
    time_TRAP_GOFf0(i)=time.GOFtrapf0th;
end

%% freq
mean_freq_dips_SBR_f1=mean(freq_SBR_f1);
std_freq_dips_SBR_f1=std(freq_SBR_f1);
mean_freq_dips_SBR_f0=mean(freq_SBR_f0);
std_freq_dips_SBR_f0=std(freq_SBR_f0);
mean_freq_SBR_GOFf1=mean(freq_SBR_GOFf1);
std_freq_SBR_GOFf1=std(freq_SBR_GOFf1);
mean_freq_SBR_GOFf0=mean(freq_SBR_GOFf0);
std_freq_SBR_GOFf0=std(freq_SBR_GOFf0);

mean_freq_dips_OLS_f1=mean(freq_OLS_f1);
std_freq_dips_OLS_f1=std(freq_OLS_f1);
mean_freq_dips_OLS_f0=mean(freq_OLS_f0);
std_freq_dips_OLS_f0=std(freq_OLS_f0);
mean_freq_OLS_GOFf1=mean(freq_OLS_GOFf1);
std_freq_OLS_GOFf1=std(freq_OLS_GOFf1);
mean_freq_OLS_GOFf0=mean(freq_OLS_GOFf0);
std_freq_OLS_GOFf0=std(freq_OLS_GOFf0);

%% time   
mean_time_dips_RAP_f1=mean(time_RAP_f1);
std_time_dips_RAP_f1=std(time_RAP_f1);
mean_time_dips_RAP_f0=mean(time_RAP_f0);
std_time_dips_RAP_f0=std(time_RAP_f0);
mean_time_RAP_GOFf1=mean(time_RAP_GOFf1);
std_time_RAP_GOFf1=std(time_RAP_GOFf1);
mean_time_RAP_GOFf0=mean(time_RAP_GOFf0);
std_time_RAP_GOFf0=std(time_RAP_GOFf0);

mean_time_dips_TRAP_f1=mean(time_TRAP_f1);
std_time_dips_TRAP_f1=std(time_TRAP_f1);
mean_time_dips_TRAP_f0=mean(time_TRAP_f0);
std_time_dips_TRAP_f0=std(time_TRAP_f0);
mean_time_TRAP_GOFf1=mean(time_TRAP_GOFf1);
std_time_TRAP_GOFf1=std(time_TRAP_GOFf1);
mean_time_TRAP_GOFf0=mean(time_TRAP_GOFf0);
std_time_TRAP_GOFf0=std(time_TRAP_GOFf0);
