function [onsets,peaks] = find_OURAPPGonsets(ppg,sf,ploton)

%% Function to analyze PPG signal and find pulse onsets
% Previous version was findPPGpeaks.m which was used in  PPG-features
% analysis in 2022. Here, movement is excluded from input variables.


% INPUTS:
%             ppg: 1st ppg signal,
%             sf: sampling frequency of ppg signal
%             ploton: 1 for plot on , 0 for plot off

% DEPENDENCIES:
% 1)Physionet Cardiovascular toolbox's modified scripts (qppg)
%

%``````````````````````````````````````````````````````````````
% Last updated by Gizem Y.D. @ 2023.08.10


%%%%--------------------------------------------------------------------------------------------------------------
if nargin < 3
    ploton=0;

elseif nargin < 2
    error('Wrong number of input arguments');
end

onsets=[];
peaks=[];


ppg1=ppg;
SF=sf;


% Define PPG
%     ppgZ=((ppg1-min(ppg1))/(max(ppg1)-min(ppg1))); %min max normalisation
ppgZ=ppg1;
ppg_t=(0:1/SF:(length(ppgZ)-1)/SF)';

% Set windowing parameters & Add 0 to artefactual regions
win_len=30;  %sec
overlap=10;    %sec
increment=win_len-overlap;
nx=numel(ppgZ)/SF;
Nwinds = fix((nx-overlap)/(win_len-overlap));    % number of sliding windows
win_start = (0:(Nwinds-1))*(win_len-overlap);  % starting index of each windows

% Start peak finding loop for PPG1
onsets_ann=[];
peaks_ann=[];

for ii=1:numel(win_start)

    [sig,tseg,ppgseg,ind_wind,ppgOnsets,ppgPeaks]=deal([]);

    %         if ii~=numel(win_start)

    ind_wind=(ppg_t>=win_start(ii) & ppg_t< win_start(ii)+win_len);
    ppgseg=ppgZ(ind_wind);
    tseg= ppg_t(ind_wind);

    % Adaptive threshols based on envelope and signal amplitude
    % %             [UP,LO]=envelope((ppgseg),200);
    %              if numel(find(UP>3.5*std(ppgseg)))<3
    %                 if numel(find(isnan(ppgseg)))/numel(ppgseg)< 0.2 % ratio of NANs should be < 0.2

    % PPG Detection - qppg - inside the window
    ppgOnsets = qppg(ppgseg,SF);


    if ~isempty(ppgOnsets)
        ppgOnsets= win_start(ii)*SF +ppgOnsets;
        onsets_ann=[onsets_ann;ppgOnsets'];
    end

    %         end
end




onsets_ann=unique(onsets_ann);


% rr=diff(peaks_ann)./SF;
% t_rr=peaks_ann(2:end)./SF;


%% Clean RR intervals (Here RR means IBI)
% [cleantOO,cleanOO] = cleanPPGpeaks(ppgZ,SF,onsets_ann);
[clean_onsets,clean_peaks,cleantNN,cleanNN]= cleanPPGpeaks(ppgZ,SF,onsets_ann);

if ~isempty(clean_onsets)


% Remove unrealistic cyclelengths (if there are any missed by cleaning algo):
clean_onsets(diff([clean_onsets(1);clean_onsets])<0.5*SF | diff([clean_onsets(1);clean_onsets])> 2.5*SF)=nan;

clean_onsets(isnan(clean_onsets))=[]; %this step is to remova NANs from annotations.

% Find 0.5 distance from each onset:
cycles=clean_onsets+ [floor(diff([clean_onsets]).*0.5); SF/2] ; %25 is approx duration of a hald cycle, added for the last pulse


newpeaks=nan(size(clean_onsets));
for u=1:numel(clean_onsets);
    if ~isnan(clean_onsets(u)) && ~isnan(cycles(u))
        [pulse,peakval,imax1]=deal([]);
        pulse=ppgZ(clean_onsets(u):cycles(u));
        [peakval,imax1]=max(pulse);
        newpeaks(u)=clean_onsets(u)+imax1;
    end
end


% Plot PPG segment and extracted peaks/onsets together, if selected
if ploton

    figure;
    s0=subplot(2,1,1);
    plot(ppg_t,ppgZ,'k');hold on;
    legend(s0,'OuraPPG');
    ylabel(s0,'PPG');
    s1=subplot(2,1,2);
    plot(ppg_t,ppgZ,'k');hold on;
    %     plot(onsets_ann./SF,ppgZ(onsets_ann),'bo');hold on
    plot(clean_onsets./SF,ppgZ(clean_onsets),'bo');hold on;
    plot(newpeaks./SF,ppgZ(newpeaks),'r*');hold on;

    legend(s1,'PPG','onsets');
    ylabel('PPG');
    linkaxes([s0,s1],'x');

end


% Output:
onsets=clean_onsets;
peaks=newpeaks;

% disp(['Cleaned PPG onsets and peaks are ready for window <3 '])



end






% Prepare for extraction
% PPG.signal=ppgZ;
% PPG.signal_t=ppg_t;
% PPG.acc=accdiff;

% PPG.type=type;
% PPG.SF=SF;
% % PPG.peaks=peaks_ann;
% PPG.onsets=onsets_ann;
% PPG.tNN=cleantNN;
% % PPG.NN=cleanNN;
% PPG.cleanPeaks=newpeaks;
% PPG.cleanOnsets=cleantOO;



end


