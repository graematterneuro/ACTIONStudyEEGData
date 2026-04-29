%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTION_Preprocessing_Pipeline_ICA.m
%
% This script is designed to run in EEGLAB and provide semi-automated     
% preprocessing of EEG data from the ACTION study and partially fulfills  
% the requirements for the author's Doctor of Philosophy program at the   
% University of Sydney, Australia.    
% This script was coded without the use of generative AI.
%
% Author: Grace Harvie, WIMR/USyd, 2026
%
% Copyright (C) 2026 Grace Harvie, Westmead Institute for Medical Research
% and The University of Sydney, grace.harvie@sydney.edu.au
%
% License: PolyForm Noncommercial License 1.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% OPEN EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Start parallel pool
parpool;

% Create a table of Participant IDs from a prepared .csv  
PIDs = readtable('C:\\Users\\Grae\\OneDrive - Westmead Institute for Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\ParticipantIDsTest.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                            VARIABLES                            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of rows in PID table so that the for loop will iterate
NoRows = height(PIDs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                BEGIN PREPROCESSING FOR LOOP                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:NoRows
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  VARIABLES FOR PIPING INTO CODE  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Current participant ID
    CurrPID = char(string(PIDs{row, 'ParticipantID'}));

    % Write directory file path
    mkdir(append('C:\\Users\\Grae\\OneDrive - Westmead Institute for Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\', CurrPID, '\\'));
    WriteDir = append('C:\\Users\\Grae\\OneDrive - Westmead Institute for Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\', CurrPID, '\\');

    % Read directory file path
    ReadFile = append('C:\Users\Grae\OneDrive - Westmead Institute for Medical Research\Documents\EEGLAB_MyFiles\ACTION\', CurrPID, '.EO.edf');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  IMPORT CurrPID.edf INTO EEGLAB  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    EEG = pop_biosig(ReadFile);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, ...
        'setname', append(CurrPID, '.EO'), ... % Name the imported file
        'gui', 'off');
    eeglab redraw;

    % Check that the file imported correctly
    userInput  = input(append("Was there a warning about the size of the channel location structure not matching with number of channels during the import process?", ...
        newline, "[Y = 1/N = 0]", newline));

    if userInput == 1

        % If the file didn't import correctly, we need to add the channel
        % label data to EEG.chanlocs manually
        % Create a structure array of the channels
        EEG.chanlocs = [struct('labels','Fp1'); struct('labels','Fp2'); ...
	    struct('labels','F7'); struct('labels','F3'); struct('labels','Fz'); ...
        struct('labels','F4'); struct('labels','F8'); struct('labels','FC3'); ...
	    struct('labels','FCz'); struct('labels','FC4');struct('labels','T3'); ...
        struct('labels','C3'); struct('labels','Cz'); struct('labels','C4'); ...
	    struct('labels','T4'); struct('labels','CP3'); struct('labels','CPz'); ...
        struct('labels','CP4'); struct('labels','T5'); struct('labels','P3'); ...
        struct('labels','Pz'); struct('labels','P4'); struct('labels','T6'); ...
        struct('labels','O1'); struct('labels','Oz'); struct('labels','O2'); ...
        struct('labels','VPVA'); struct('labels','VNVB'); ...
        struct('labels','HPHL'); struct('labels','HNHR'); ...
        struct('labels','Erbs'); struct('labels','OrbOcc'); ...
        struct('labels','Mass'); struct('labels','EDA'); ...
        struct('labels','Resp'); struct('labels','ECG'); ...
        struct('labels','Events'); struct('labels','A1A2'); ...
        struct('labels','A2'); struct('labels','Cer7')];

        % Add the other fields to the first channel in the structure array
        % which will add these fields empty to the other channels.
        EEG.chanlocs(1).ref = [];
        EEG.chanlocs(1).theta = [];
        EEG.chanlocs(1).radius = [];
        EEG.chanlocs(1).X = [];
        EEG.chanlocs(1).Y = [];
        EEG.chanlocs(1).Z = [];
        EEG.chanlocs(1).sph_theta = [];
        EEG.chanlocs(1).sph_phi = [];
        EEG.chanlocs(1).sph_radius = [];
        EEG.chanlocs(1).type = [];

    else
    end
    
    % Add channel location data to the dataset and save the dataset
    EEG = pop_chanedit(EEG, {'lookup',['C:\Users\Grae\OneDrive - Westmead ' ...
        'Institute for Medical Research\Documents\EEGLAB_MyFiles\Neuroscan ' ...
        'Chan Locs\THIS-ONE-NuAmps40-forACTION-Centred.ced']});
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    % Save the imported data as an EEGLAB data set (.set)
    FileName = append(CurrPID, '.EO.set');
    EEG = pop_saveset(EEG, ...
        'filename', FileName, ...
        'filepath', WriteDir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    eeglab redraw;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  START OF PREPROCESSING PIPELINE  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%  SELECT THE CEPHALIC CHANNELS  %%%%%
    EEG = pop_select(EEG, ...
        'channel',{'Fp1','Fp2','F7','F3','Fz','F4','F8','FC3','FCz','FC4', ...
        'T3','C3','Cz','C4','T4','CP3','CPz','CP4','T5','P3','Pz','P4', ...
        'T6','O1','Oz','O2'});

    % Save as new dataset
    SetName = append(CurrPID, '.Ceph.set');
    SaveNew = append(WriteDir, CurrPID, '.Ceph.set');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew); 

    eeglab redraw;

    %%%%%  HIGH-PASS FILTER  %%%%%
    EEG = pop_eegfiltnew(EEG, ...
        'locutoff',1);

    % Save as new dataset
    SetName = append(CurrPID, '.HiPass.set');
    SaveNew = append(WriteDir, CurrPID, '.HiPass.set');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew);

    eeglab redraw;

    %%%%%  CHECK FOR BAD CHANNELS  %%%%%
    figure; pop_spectopo(EEG, 1, ...
        [0  119998], ...            % Epoch time range to analyze [min_ms max_ms]
        'EEG' , ...
        'percent', 100, ...         % Percent data to sample (1 to 100)
        'freq', [6 10 22 60], ...   % Frequencies to plot as scalp maps (Hz)
        'freqrange',[2 70], ...     % Plotting frequency range [lo_Hz hi_Hz]
        'electrodes','off' ...      % Spectral and scalp map options (see topoplot)
        );

    userInput  = input(append("Are there bad channels which requires removing?", ...
        newline, "[Y = 1/N = 0]", newline));

    if userInput == 1
         % Remove bad channels and save as new dataset
         % EEG = pop_cleanline(EEG, ...
         %    'bandwidth',2, ...      % Bandwidth
         %    'chanlist', [1 26], ... % Indices to clean (channel/components)
         %    'computepower',1, ...   % Visualise original and cleaned spectra (?) (T/F)
         %    'linefreqs',50, ...     % Line noise freqs to remove
         %    'newversion',1, ...     % Use new implementation (T/F)
         %    'normSpectrum',0, ...   % Normalise log spectrum by detrending (T/F)
         %    'p',0.01, ...           % P-val for detection of sinusoid
         %    'pad',2, ...            % FFT padding factor
         %    'plotfigures',0, ...    % Plot individual figures (T/F)
         %    'scanforlines',0, ...   % Scan for line noise (T/F)
         %    'sigtype', ...
         %    'Channels', ...         % Type of signal to clean (channels/components)
         %    'taperbandwidth',2, ... % Taper bandwidth
         %    'tau',100, ...          % Window overlap smoothing factor
         %    'verb',0, ...           % Produce verbose output
         %    'winsize',4, ...        % Sliding window length (s)
         %    'winstep',1 ...         % Sliding window step size (s)
         %    );
    
         % Save as new dataset
         SetName = append(CurrPID, '.BadChanRej.set');
         SaveNew = append(WriteDir, CurrPID, '.BadChanRej.set');
         [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
             'setname', SetName, ...
             'savenew', SaveNew, ...
             'gui','off'); 

         eeglab redraw;
    else
    end

    % Plot continuous data for comparison with next step
    % pop_eegplot(EEG, 1, 1, 1);

    %%%%%  LINE NOISE  %%%%%
    figure; pop_spectopo(EEG, 1, ...
        [0  119998], ...            % Epoch time range to analyze [min_ms max_ms]
        'EEG' , ...
        'percent', 100, ...         % Percent data to sample (1 to 100)
        'freq', [6 10 22 60], ...   % Frequencies to plot as scalp maps (Hz)
        'freqrange',[2 70], ...     % Plotting frequency range [lo_Hz hi_Hz]
        'electrodes','off' ...      % Spectral and scalp map options (see topoplot)
        );

    userInput  = input(append("Is there line noise which requires cleaning?", ...
        newline, "[Y = 1/N = 0]", newline));

    if userInput == 1

        % Plot channel spectra close up to determine line noise frequencies
        figure; pop_spectopo(EEG, 1, ...
            [0  119998], ...            % Epoch time range to analyze [min_ms max_ms]
            'EEG' , ...
            'percent', 100, ...         % Percent data to sample (1 to 100)
            'freq', [50 65], ...   % Frequencies to plot as scalp maps (Hz)
            'freqrange',[45 70], ...     % Plotting frequency range [lo_Hz hi_Hz]
            'electrodes','off' ...      % Spectral and scalp map options (see topoplot)
            );

        userInput = input(append("Is there more than one frequency to clean?", ...
        newline, "[Y = 1/N = 0]", newline));

        if userInput == 1

            lnFreq1 = input("Frequency 1: ");
            lnFreq2 = input ("Frequency 2: ");
            
            % Clean line noise from cephalic channel and save as new dataset
            EEG = pop_cleanline(EEG, ...
               'bandwidth',2, ...      % Bandwidth
               'chanlist', [1:26], ... % Indices to clean (channel/components)
               'computepower',1, ...   % Visualise original and cleaned spectra (?) (T/F)
               'linefreqs', [lnFreq1 lnFreq2], ...     % Line noise freqs to remove
               'newversion',1, ...     % Use new implementation (T/F)
               'normSpectrum',0, ...   % Normalise log spectrum by detrending (T/F)
               'p',0.01, ...           % P-val for detection of sinusoid
               'pad',2, ...            % FFT padding factor
               'plotfigures',0, ...    % Plot individual figures (T/F)
               'scanforlines',0, ...   % Scan for line noise (T/F)
               'sigtype', ...
               'Channels', ...         % Type of signal to clean (channels/components)
               'taperbandwidth',2, ... % Taper bandwidth
               'tau',100, ...          % Window overlap smoothing factor
               'verb',0, ...           % Produce verbose output
               'winsize',4, ...        % Sliding window length (s)
               'winstep',1 ...         % Sliding window step size (s)
               );

        else

            lnFreq = input("Frequency: ");

            % Clean line noise from cephalic channel and save as new dataset
            EEG = pop_cleanline(EEG, ...
               'bandwidth',2, ...      % Bandwidth
               'chanlist', [1:26], ... % Indices to clean (channel/components)
               'computepower',1, ...   % Visualise original and cleaned spectra (?) (T/F)
               'linefreqs', lnFreq, ...     % Line noise freqs to remove
               'newversion',1, ...     % Use new implementation (T/F)
               'normSpectrum',0, ...   % Normalise log spectrum by detrending (T/F)
               'p',0.01, ...           % P-val for detection of sinusoid
               'pad',2, ...            % FFT padding factor
               'plotfigures',0, ...    % Plot individual figures (T/F)
               'scanforlines',0, ...   % Scan for line noise (T/F)
               'sigtype', ...
               'Channels', ...         % Type of signal to clean (channels/components)
               'taperbandwidth',2, ... % Taper bandwidth
               'tau',100, ...          % Window overlap smoothing factor
               'verb',0, ...           % Produce verbose output
               'winsize',4, ...        % Sliding window length (s)
               'winstep',1 ...         % Sliding window step size (s)
               );
        end
    
        % Save as new dataset
        SetName = append(CurrPID, '.LineNoise.set');
        SaveNew = append(WriteDir, CurrPID, '.LineNoise.set');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew, ...
            'gui','off'); 
        
        eeglab redraw;
        else
    end

    %%%%%  FULL-RANK AVERAGE REFERENCE  %%%%%
    EEG = fullRankAveRef(EEG);

    % Save as new dataset
    SetName = append(CurrPID, '.Ref.set');
    SaveNew = append(WriteDir, CurrPID, '.Ref.set');
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
             'setname', SetName, ...
             'savenew', SaveNew, ...
             'gui','off');

    eeglab redraw;

    % Plot continuous data for comparison with previous step
    % pop_eegplot(EEG, 1, 1, 1);

    %%%%%  BAD DATA SEGMENT REMOVAL  %%%%%
    % EEG = pop_clean_rawdata(EEG, ...
    %    'FlatlineCriterion', 'off', ...
    %    'ChannelCriterion', 'off', ...
    %    'LineNoiseCriterion', 'off', ...
    %    'Highpass', 'off', ...
    %    'BurstCriterion', 'off', ...
    %    'WindowCriterion', 0.25, ...
    %    'BurstRejection', 'off', ...
    %    'Distance', 'Euclidian', ...
    %    'WindowCriterionTolerances',[-Inf 7] ...
    %    );
    % 
    % % Save as new dataset
    % SetName = append(CurrPID, '.ArtRej.set');
    % SaveNew = append(WriteDir, CurrPID, '.ArtRej.set');
    % [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
    %          'setname', SetName, ...
    %          'savenew', SaveNew, ...
    %          'gui','off'); 
    % 
    % eeglab redraw;
    
    %%%%  ICA DECOMPOSITION  %%%%%
    EEG = pop_runica(EEG, ...
       'icatype', 'runica', ...
       'extended', 1, ...
       'lrate', 1e-05, ...
       'maxsteps', 2000, ...
       'interrupt','off');
    
    EEG = pop_saveset(EEG, ...
        'filename', SetName, ...
        'filepath', WriteDir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    eeglab redraw;

   STUDY = []; 
   CURRENTSTUDY = 0; 
   ALLEEG = []; 
   EEG=[]; 
   CURRENTSET=[];

   eeglab redraw;
end
