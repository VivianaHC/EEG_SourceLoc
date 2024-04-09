# EEG_SourceLoc

Matlab implementation for performing source localization in EEG recordings using the algorithms Single Best Replacement (SBR) and recursive version of MUSIC (RAP-MUSIC and TRAP-MUSIC) in the frequency domain and the time domain. 

# Preprocessing
Preprocessing was applied to 129 EEG files previous localization procedure. The preprocessing on each EEG file consisted of averaging across sessions, as a way of noise reduction. For frequency domain localization, fast Fourier transform (FFT) was applied to EEG averaged data and frequencies of interest were selected from this transformation. Real and imaginary data were selected at the frequency and harmonics of face individuation response (f1 = 1.2 Hz, 2.4 Hz, 3.6 Hz, 4.8 Hz), as well at the frequency ahd harmonics of general visual response (f0 = 6 Hz, 12 Hz, 18 Hz, 24 Hz). For time domain localization, band pass filtering was applied to EEG averaged data at the range of the frequencies of interest and their harmonics for face individuation response (f1 = from 1 Hz to 5 Hz) and general visual response (f0 = from 5 Hz to 25 Hz).  

# Localization
Localization was applied to the preprocessed data in frequency and time domain using SBR and TRAP and RAP MUSIC algorithms. For algorithms description, please refer to the sources (Castanon, 2023) and (). Review in the corresponding matlab file the description and parameters used for each algorithm. 

