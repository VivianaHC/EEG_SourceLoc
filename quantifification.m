%% Results quantification: n dips and GOF 
modeldir = 'C:\Users\admin-sm\Documents\MATLAB\data\EEGloc_fvps_SBR_freq_T_RAP_time_f1_f0_/*.mat';
matFiles = dir(modeldir);

for i = 1:length(matFiles)
    baseFileName = matFiles(i).name;
    load(baseFileName)    
    freq_SBR_f1(i) = numel(freq.inddipsbrf1h);
    freq_SBR_f0(i) = numel(freq.inddipsbrf0h);
    time_RAP_f1(i) = numel(time.inddiprapf1th);
    time_RAP_f0(i) = numel(time.inddiprapf0th);
    time_TRAP_f1(i) = numel(time.inddiptrapf1th);
    time_TRAP_f0(i) = numel(time.inddiptrapf0th);
     
    freq_SBR_GOFf1(i) = freq.GOFsbrf1h;    
    freq_SBR_GOFf0(i) = freq.GOFsbrf0h;
    time_RAP_GOFf1(i) = time.GOFrapf1th;    
    time_RAP_GOFf0(i) = time.GOFrapf0th;
    time_TRAP_GOFf1(i) = time.GOFtrapf1th;    
    time_TRAP_GOFf0(i) = time.GOFtrapf0th;
end

%% freq
mean_freq_dips_SBR_f1 = mean(freq_SBR_f1);
std_freq_dips_SBR_f1 = std(freq_SBR_f1);
mean_freq_dips_SBR_f0 = mean(freq_SBR_f0);
std_freq_dips_SBR_f0 = std(freq_SBR_f0);
mean_freq_SBR_GOFf1 = mean(freq_SBR_GOFf1);
std_freq_SBR_GOFf1 = std(freq_SBR_GOFf1);
mean_freq_SBR_GOFf0 = mean(freq_SBR_GOFf0);
std_freq_SBR_GOFf0 = std(freq_SBR_GOFf0);

%% time   
mean_time_dips_RAP_f1 = mean(time_RAP_f1);
std_time_dips_RAP_f1 = std(time_RAP_f1);
mean_time_dips_RAP_f0 = mean(time_RAP_f0);
std_time_dips_RAP_f0 = std(time_RAP_f0);
mean_time_RAP_GOFf1 = mean(time_RAP_GOFf1);
std_time_RAP_GOFf1 = std(time_RAP_GOFf1);
mean_time_RAP_GOFf0 = mean(time_RAP_GOFf0);
std_time_RAP_GOFf0 = std(time_RAP_GOFf0);

mean_time_dips_TRAP_f1 = mean(time_TRAP_f1);
std_time_dips_TRAP_f1 = std(time_TRAP_f1);
mean_time_dips_TRAP_f0 = mean(time_TRAP_f0);
std_time_dips_TRAP_f0 = std(time_TRAP_f0);
mean_time_TRAP_GOFf1 = mean(time_TRAP_GOFf1);
std_time_TRAP_GOFf1 = std(time_TRAP_GOFf1);
mean_time_TRAP_GOFf0 = mean(time_TRAP_GOFf0);
std_time_TRAP_GOFf0 = std(time_TRAP_GOFf0);
