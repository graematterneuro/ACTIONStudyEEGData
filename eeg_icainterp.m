% function EEGout = eeg_icainterp(EEG,rejchan,fullchanlocs)
%
% Interpolate missing channels of EEG.icawinv (ICA maps) and EEG data
%
% Inputs:
%
%     EEG         -     a dataset with ICA from data with rejected / missing channels
%     rejchan     -     the channel/s which has/have beem rejected
%     fullchanlocs  -     complete channel structure which will be used to
%                       interpolate missing data
%
% Outputs:
%
%     EEGout        -    dataset with ICA maps and data with rejected /
%                        missing channels interpolated
%
% Author: Jason Palmer, West Virginia University, 2026
% https://sccn.ucsd.edu/pipermail/eeglablist/2026/019041.html
%
% Edited by Grace Harvie to work with the ACTION EEG preprocessing pipeline
% and use the Spherical Kang method, Sept 2026

function EEGout = eeg_icainterp(EEG,rejchan,fullchanlocs)

% get the original channel numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE: I couldn't get this to work, I think because my dataset has data 
%%%       labels. Instead I created a copy of the complete data structure 
%%%       in my preprocessing pipeline and get the chanlist from that - GH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k = 1:EEG.nbchan
%     v(k) = EEG.urchanlocs(k).urchan;
% end

fullChanList = 1:size(fullchanlocs,2);

% create fake dataset to interpolate ica maps
EEGtmp = EEG;
EEGtmp.data = EEG.icawinv;
EEGtmp.icaact = [];
EEGtmp.trials = 1;
EEGtmp.pnts = size(EEG.icawinv,2);[]
EEGtmp2 = eeg_interp(EEGtmp,fullchanlocs,'sphericalKang'); % GH edit here

% copy the icaact and interpolated icawinv into output
EEGout = EEG;
EEGout.nbchan = length(fullChanList);
EEGout.icawinv = EEGtmp2.data;
EEGout.icaact = EEG.icaact;

% add zero columns to icasphere where (new) interpolated data channels are
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE: This wouldn't work for me either, so I wrote a new method to a
%%%       zero-filled columns to icasphere                            - GH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EEGout.icasphere(:,fullChanList) = EEG.icasphere;
% EEGout.icasphere(:,setdiff(1:EEGout.nbchan,fullChanList)) = 0;

% New method to add zero columns to icasphere where (new) interpolated data channels are
add = length(EEG.icasphere)+1;
EEGout.icasphere(:,add:EEGout.nbchan) = 0;
EEG.icasphere(:,add:EEGout.nbchan) = 0;
readCount = 1;
for writeCount = 1:EEGout.nbchan
    if ismember(writeCount,rejchan) == 1
        readCount = readCount-1;
    else
        EEGout.icasphere(:,writeCount) = EEG.icasphere(:,readCount);
    end
    readCount = readCount+1;
end

for zeroCount = 1:length(rejchan)
    EEGout.icasphere(:,rejchan(zeroCount)) = 0;
end

% add data with interpolated channels
EEGout.data = reshape(EEGout.icawinv * double(EEGout.icaact(:,:)),EEGout.nbchan,EEGout.pnts,EEGout.trials);

% copy full chanlocs and icachansind
EEGout.chanlocs = fullchanlocs;
EEGout.icachansind = 1:EEGout.nbchan;