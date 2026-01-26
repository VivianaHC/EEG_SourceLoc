# EEG Data
Matlab implementation for performing source localization in EEG recordings using the algorithms Single Best Replacement (SBR-R1) and TRAP MUSIC in the frequency domain and the time domain respectively. The EEG database is presented in Rossion, B., Retter, T.L., & Liu-Shuang, J. (2020). "Understanding human individuation of unfamiliar faces with oddball fast periodic visual stimulation and electroencephalography", and the access to the database is from: https://data.4tu.nl/articles/_/12707438/1
For details about methodology for data collection and data specifications, please refer to those sources.

# Preprocessing 
Preprocessing and localization were applied to all EEG recordings found in the data folder by pairing each .lw5 header file with its corresponding .mat data file (one per subject). For each file, EEG epochs were first preprocessed in the time domain by (i) removing the DC offset independently for each channel and (ii) re-referencing to a common average reference (CAR) at each time simple, which reduces slow drifts and global noise shared across electrodes.

# Frequency Domain Localization
Each preprocessed epoch was transformed with an FFT, then average of complex values was applied across epochs to preserve phase-locked FPVS responses. The analyzed time window was adjusted to contain an integer number of oddball cycles to improve frequency-bin alignment. Harmonics were defined for the base responses (6 Hz, 12 Hz, 18 Hz, 24 Hz) and oddball responses (1.2 Hz, 2.4 Hz, 3.6 Hz, 4.8 Hz). For each harmonic, a complex harmonic estimate was obtained by averaging a few adjacent FFT bins to reduce spectral leakage, and a local baseline correction was applied by substracting the mean magnitude from adjacent frequency bins. These complex harmonic stacks were then used as input to SBR-R1 algorithm for source localization.

# Time Domain Localization
Epochs were first averaged in time and then band-pass filtered with a zero-phase Butterworth filter in the frequency ranges targeting FPVS components: oddball response ~0.85-5 Hz and base response ~5-25 Hz. The filtered signals were localized using (T)RAP MUSIC algorithm, with the number of sources guided by the number of dipoles returned by the frequency domain solution.  


