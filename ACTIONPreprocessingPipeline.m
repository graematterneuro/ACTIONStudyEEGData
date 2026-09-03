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
% Other contributions noted in comments where appropriate.
%
% Copyright (C) 2026 Grace Harvie, Westmead Institute for Medical Research
% and The University of Sydney, grace.harvie@sydney.edu.au
%
% License: PolyForm Noncommercial License 1.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                        OPEN EEGLAB                              %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Create a table of Participant IDs from a prepared .csv
PIDs = readtable(['C:\\Users\\Grae\\OneDrive - Westmead Institute for Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\ParticipantIDs' ...
    '.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                            VARIABLES                            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the number of rows in PID table so that the for loop will iterate
NoRows = height(PIDs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                BEGIN PREPROCESSING FOR LOOP                     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for row = 1:NoRows

    % Check if there is currently a parallel pool, if one does not exist,
    % create one
    pool = gcp("nocreate");
    if isempty(pool)
        parpool;
    end

    passCounter = 1;
    attempt = num2str(passCounter);
    repeatWReject = 1;

    while repeatWReject == 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  VARIABLES FOR PIPING INTO CODE  %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Current participant ID
        CurrPID = char(string(PIDs{row, 'ParticipantID'}));

        % Write directory file path
        mkdir(append('C:\\Users\\Grae\\OneDrive - Westmead Institute for ', ...
            'Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\', ...
            CurrPID, '\\Attempt', attempt, '\\'));
        WriteDir = append('C:\\Users\\Grae\\OneDrive - Westmead Institute ', ...
            'for Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\', ...
            CurrPID, '\\Attempt', attempt, '\\');

        % Read directory file path
        ReadFile = append('C:\Users\Grae\OneDrive - Westmead Institute for ', ...
            'Medical Research\Documents\EEGLAB_MyFiles\ACTION\', CurrPID, ...
            '.EO.edf');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  IMPORT CurrPID.edf INTO EEGLAB  %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        EEG = pop_biosig(ReadFile);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, ...
            'setname', append(CurrPID, '.EO'), ... % Name the imported file
            'gui', 'off');
        eeglab redraw;

        % Check that the file imported correctly
        userInput  = input(append("Is the EEG file missing channel labels?", ...
            newline, "[Y = 1/N = anything else]", newline));

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

        %%%%%  HIGH-PASS FILTER AT 1Hz  %%%%%
        EEG = pop_eegfiltnew(EEG, ...
            'locutoff', 1);

        % Save as new dataset
        SetName = append(CurrPID, '.HiPass.set');
        SaveNew = append(WriteDir, CurrPID, '.HiPass.set');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew);

        eeglab redraw;


        %%%%%  CLEAN LINE NOISE USING CLEANLINE  %%%%%
        figure; pop_spectopo(EEG, 1, ...
            [0  119998], ...            % Epoch time range to analyze [min_ms max_ms]
            'EEG' , ...
            'percent', 100, ...         % Percent data to sample (1 to 100)
            'freq', [50 60 65], ...   % Frequencies to plot as scalp maps (Hz)
            'freqrange',[45 70], ...     % Plotting frequency range [lo_Hz hi_Hz]z]
            'electrodes','off' ...      % Spectral and scalp map options (see topoplot)
            );

        userInput  = input(append("Is there line noise which requires cleaning?", ...
            newline, "[Y = 1/N = anything else]", newline));

        if userInput == 1

            % Plot channel spectra close up to determine line noise frequencies
            userInput = input(append("Is there more than one frequency to clean?", ...
                newline, "[Y = 1/N = anything else]", newline));

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
        end

        fprintf("Press any key to continue to ASR.");

        pause();

        %%%%%  CLEAN BAD DATA SEGMENTS USING ASR  %%%%%
        cleanEEG = clean_asr(EEG, 10);
        vis_artifacts(cleanEEG, EEG);
        EEG = cleanEEG;
        SetName = append(CurrPID, '.ASR.set');
        SaveNew = append(WriteDir, CurrPID, '.ASR.set');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew, ...
            'gui','off');

        eeglab redraw;


        %%%%%  REJECT BAD CHANNELS  %%%%%
        if passCounter > 1

            fprintf("Press any key to continue to bad channel rejection.");

            pause();

            figure; pop_spectopo(EEG, 1, ...
                [0  119998], ...            % Epoch time range to analyze [min_ms max_ms]
                'EEG' , ...
                'percent', 100, ...         % Percent data to sample (1 to 100)
                'freq', [6 10 22 60], ...   % Frequencies to plot as scalp maps (Hz)
                'freqrange',[2 70], ...     % Plotting frequency range [lo_Hz hi_Hz]
                'electrodes','off' ...      % Spectral and scalp map options (see topoplot)
                );
            pop_eegplot(EEG, 1, 1, 1);

            % Get the original channel numbers and a copy of the full 
            % channel structure - Needed for eeg_icainterp() on line 348
            chanlocs = EEG.chanlocs;

            BadChan = input("Using numeric notation, list the channels which require interpolation in square brackets:");

            EEG = pop_select(EEG, 'rmchannel', BadChan);

            SetName = append(CurrPID, '.ChanRej.set');
            SaveNew = append(WriteDir, CurrPID, '.ChanRej.set');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
                'setname', SetName, ...
                'savenew', SaveNew, ...
                'gui','off');

            eeglab redraw;

        end

        fprintf("Press any key to continue to full rank average referencing.");

        pause();

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

        fprintf("Press any key to continue to ICA decomposition and IC labelling.");

        pause();

        %%%%  ICA DECOMPOSITION AND LABELLING  %%%%%
        EEG = pop_runica(EEG, ...
            'icatype', 'runica', ...
            'extended', 1, ...
            'lrate', 1e-05, ...
            'maxsteps', 2000, ...
            'interrupt','off');

        EEG = pop_iclabel(EEG, 'default');

        SetName = append(CurrPID, '.ICALabelled.set');
        SaveNew = append(WriteDir, CurrPID, '.ICALabelled.set');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew, ...
            'gui','off');

        eeglab redraw;

        %%%%%  ICA INTERPOLATION  %%%%%

        if passCounter > 1
            fprintf("Press any key to interpolate icawinv.");

            pause();

            addpath('C:\Users\Grae\OneDrive - Westmead Institute for Medical Research\Documents\EEGLAB_MyFiles\GitHubRepo');

            EEG = eeg_icainterp(EEG, BadChan, chanlocs);

            SetName = append(CurrPID, '.ICAInterp.set');
            SaveNew = append(WriteDir, CurrPID, '.ICAInterp.set');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
                'setname', SetName, ...
                'savenew', SaveNew, ...
                'gui','off');

            eeglab redraw;

        end

        fprintf("Press any key to label independent components.");

        pause();

        %%%%%  IC REJECTION  %%%%%%
        fprintf("Select ICs for rejection then press any key to continue.");
        addpath('C:\Users\Grae\OneDrive - Westmead Institute for Medical Research\Documents\eeglab_current\eeglab2026.0.0\plugins\ICLabel\viewprops');
        pop_viewprops(EEG, 0);
        pop_selectcomps(EEG, [1:size(EEG.icawinv,2)]);

        pause();

        % Plot single trials before and after IC rejection.
        % Code on lines 345 to 356 (inclusive) merged from pop_subcomp.m
        % (https://github.com/sccn/eeglab/blob/develop/functions/popfunc/pop_subcomp.m)
        components = [];
        components = find(EEG.reject.gcompreject == 1);
        components = components(:)';
        component_keep = setdiff_bc(1:size(EEG.icaweights,1), components);
        compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
        compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);

        eegplot(EEG.data(EEG.icachansind,:,:), ...
            'srate', EEG.srate, ...
            'title', 'Black = channel before rejection; red = after rejection -- eegplot()', ...
            'limits', [EEG.xmin EEG.xmax]*1000, ...
            'data2', compproj);

        fprintf("Press any key to remove the rejected components.");

        pause();

        EEG = pop_subcomp( EEG, [], 0);

        SetName = append(CurrPID, '.ICACleaned.set');
        SaveNew = append(WriteDir, CurrPID, '.ICACleaned.set');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
            'setname', SetName, ...
            'savenew', SaveNew, ...
            'gui','off');

        eeglab redraw;

        repeatWReject = input(append("Are there channels which are not adequately cleaned by ICA?", ...
        newline, "[Y = 1/N = anything else]", newline));

        if any(repeatWReject) == 0
            repeatWReject = 0;
            fprintf("Press any key to continue to epoching.");

            pause();

            %%%%%  EPOCH DATA  %%%%%
            EEG = eeg_regepochs(EEG, ...
                'recurrence', 2);

            % Save as new dataset
            SetName = append(CurrPID, '.Epoched.set');
            SaveNew = append(WriteDir, CurrPID, '.Epoched.set');
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET+1, ...
                'setname', SetName, ...
                'savenew', SaveNew);

            eeglab redraw;

            fprintf("Press any key to export cleaned data to text.");

            pause();


            %%%%% EXPORT AS TEXT FILE %%%%%
            % Write directory file path
            WriteDir = append('C:\\Users\\Grae\\OneDrive - Westmead Institute for ', ...
                'Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\', ...
                'Preprocessing\\ExportedEEGs');

            fileName = append(WriteDir, CurrPID, '.Cleaned.txt');
            pop_export(EEG, fileName, ...
                'transpose', 'on', ...
                'precision', 4);

            fprintf("Press any key to clear all and load next participant.");

            pause();

            STUDY = [];
            CURRENTSTUDY = 0;
            ALLEEG = [];
            EEG=[];
            CURRENTSET=[];

            eeglab redraw;

            repeatWReject = 0;

        elseif repeatWReject == 1;
            fprintf("Press any key to clear all and try again.");

            pause();

            STUDY = [];
            CURRENTSTUDY = 0;
            ALLEEG = [];
            EEG=[];
            CURRENTSET=[];

            eeglab redraw;

            passCounter = passCounter+1;

            attempt = num2str(passCounter);

            mkdir(append('C:\\Users\\Grae\\OneDrive - Westmead Institute for Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\', CurrPID, '\\Attempt'));

        end

    end
    
end
