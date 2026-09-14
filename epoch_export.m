% function export = epoch_export(EEG,currPID,writeDir,)
%
% Exports good epochs over multiple files.
%
% Inputs:
%
%     EEG         -     a dataset with rejected epochs
%     currPID     -     Current participant ID
%
% Outputs:
%
%     export      -     EEG exported as .txt
%
% Copyright (C) 2026 Grace Harvie, Westmead Institute for Medical Research
% and The University of Sydney, grace.harvie@sydney.edu.au

function export = epoch_export(EEG,CurrPID)

writeDir = ['C:\\Users\\Grae\\OneDrive - Westmead Institute for Medical ' ...
    'Research\\Documents\\EEGLAB_MyFiles\\ACTION\\Preprocessing\\ExportedEEGs\\'];

% Exports one good epoch per data file
for i = 1:size(EEG.reject.rejmanual, 2)
    if EEG.reject.rejmanual(i) == 1
        % Skip the rejected epoch
    else
        temp = pop_select(EEG, 'trial', [i:i]);
        epoch = num2str(i);
        filename = append(writeDir, 'Epoch', epoch, '//', CurrPID, '.Epoch', epoch, '.txt');
        pop_export(temp, filename,'transpose','on','precision',4);
    end
end


mkdir(append('C:\\Users\\Grae\\OneDrive - Westmead Institute for ', ...
            'Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\', ...
            'Preprocessing\\ExportedEEGs\\', CurrPID, '\\'));

writeDir = append('C:\\Users\\Grae\\OneDrive - Westmead Institute for ', ...
            'Medical Research\\Documents\\EEGLAB_MyFiles\\ACTION\\', ...
            'Preprocessing\\ExportedEEGs\\', CurrPID, '\\');

% Exports the runs of good epochs as discrete chunks

 rejEpochs = find(EEG.reject.rejmanual);

if length(rejEpochs) == 1 & rejEpochs(1) == 1
    fprintf("Exporting epochs 2-60")
    temp = pop_select(EEG, 'trial', [2:60]);
    filename = append(writeDir, CurrPID, '.Epoch2-60.txt');
    pop_export(temp, filename,'transpose','on','precision',4);

elseif length(rejEpochs) == 1 & rejEpochs(1) == 60
    fprintf("Exporting epochs 1-59")
    temp = pop_select(EEG, 'trial', [1:59]);
    filename = append(writeDir, CurrPID, '.Epoch1-59.txt');
    pop_export(temp, filename,'transpose','on','precision',4);

elseif length(rejEpochs) == 1 & rejEpochs(1) ~= 1 & rejEpochs(1) ~= 6
    stopEpoch = rejEpochs(1)-1;
    temp = pop_select(EEG, 'trial', [1:stopEpoch]);
    stopEpoch = num2str(stopEpoch);
    filename = append(writeDir, CurrPID, '.1-', stopEpoch, '.txt');
    pop_export(temp, filename,'transpose','on','precision',4);

    startEpoch = rejEpochs(1)+1;
    temp = pop_select(EEG, 'trial', [startEpoch:60])
    startEpoch = num2str(startEpoch);
    filename = append(writeDir, CurrPID, '.', startEpoch, '-60.txt');
    pop_export(temp, filename,'transpose','on','precision',4);

    fprintf(append("Exporting epochs 1 to ", stopEpoch, " and epochs ", startEpoch , " to 60"));


elseif length(rejEpochs) > 1 & rejEpochs(1) == 1 & rejEpochs(length(rejEpochs)) == 60
    for i=1:length(rejEpochs)-1 % Requires one less iteration

        if i >= 1 & i < length(rejEpochs)+1 & rejEpochs(i) == rejEpochs(i+1)-1
            % Skip epochs which are consecutive.

        elseif i == 1 %& rejEpochs(i) == 1 
            startEpoch = rejEpochs(i)+1;
            stopEpoch = rejEpochs(2)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        else
            startEpoch = rejEpochs(i)+1;
            stopEpoch = rejEpochs(i+1)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));
        end
    end

elseif length(rejEpochs) > 1 & rejEpochs(1) ~= 1 & rejEpochs(length(rejEpochs)) ~= 60
    for i = 1:length(rejEpochs)+1 % Requires one extra iteration

        if i > 1 & i < length(rejEpochs)+1 & rejEpochs(i) == rejEpochs(i-1)+1
            % Skip epochs which are consecutive.

        elseif i == 1
            startEpoch = 1;
            stopEpoch = rejEpochs(i)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        elseif i == length(rejEpochs)+1
            startEpoch = rejEpochs(i-1)+1;
            stopEpoch = 60;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        else
            startEpoch = rejEpochs(i-1)+1;
            stopEpoch = rejEpochs(i)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));
        end
    end

elseif length(rejEpochs) > 1 & rejEpochs(1) == 1 & rejEpochs(length(rejEpochs)) ~= 60
    for i=1:length(rejEpochs)

        if i >= 1 & i < length(rejEpochs) & rejEpochs(i) == rejEpochs(i+1)-1
            % Skip epochs which are consecutive.

        elseif i == 1
            startEpoch = rejEpochs(1)+1;
            stopEpoch = rejEpochs(i+1)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs 2 to ", stopEpoch, newline));

        elseif i == length(rejEpochs)
            startEpoch = rejEpochs(i)+1;
            stopEpoch = 60;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        else
            startEpoch = rejEpochs(i)+1;
            stopEpoch = rejEpochs(i+1)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        end
    end

elseif length(rejEpochs) > 1 & rejEpochs(1) ~= 1 & rejEpochs(length(rejEpochs)) == 60
    for i=1:length(rejEpochs)
        if i > 1 & i < length(rejEpochs) & rejEpochs(i) == rejEpochs(i+1)-1
            % Skip epochs which are consecutive.

        elseif i == length(rejEpochs) & rejEpochs(i) == rejEpochs(i-1)+1
            % Skip epochs which are consecutive.

        elseif i == 1
            startEpoch = 1;
            stopEpoch = rejEpochs(i)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        elseif i == length(rejEpochs)
            startEpoch = rejEpochs(i-1)+1;
            stopEpoch = rejEpochs(i)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));

        else
            startEpoch = rejEpochs(i)+1;
            stopEpoch = rejEpochs(i+1)-1;
            temp = pop_select(EEG, 'trial', [startEpoch:stopEpoch]);

            startEpoch = num2str(startEpoch);
            stopEpoch = num2str(stopEpoch);
            filename = append(writeDir, CurrPID, '.', startEpoch, '-', stopEpoch, '.txt');
            pop_export(temp, filename,'transpose','on','precision',4);

            fprintf(append("Exporting epochs ", startEpoch, " to ", stopEpoch, newline));
        end
    end
end
