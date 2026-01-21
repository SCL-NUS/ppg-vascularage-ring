function [cleanOnsets,cleanPeaks,cleantNN,cleanNN] = cleanPPGpeaks(ppg,sf,onsets)
% This function cleans PPG peaks/ onsets according to the established
% physiological ranges, so that, values that are outside of
% physiological/normal ranges are removed.
% Cleaning process is similar to what we do in cleaning RR time series, but
% with different parameters. 


% INPUTS:
%             ppg: 1st ppg signal,
%             sf: sampling frequency of ppg signal
%             type: string, 'Oura' or 'Pleth'
%             peaks: peaks/onsets in ppg signal (as output of findPPGpeaks.m)
%            

% DEPENDENCIES:
% 1)Physionet Cardiovascular toolbox's modified scripts (FindSpikesInRRseries.m)
%

% Created by Gizem Yilmaz Durmusoglu @ 09.06.22

%%%%--------------------------------------------------------------------------------------------------------------
% Display current file:

if nargin < 3
    error('Wrong number of input arguments');
end

ppg1=ppg;
SF=sf;
onsets_ann=onsets;

% Define PPG
ppgZ=ppg1;
ppg_t=(0:1/SF:(length(ppgZ)-1)/SF)';

%% clean using Onsets
rr=diff(onsets_ann) ./SF; %
t_rr=onsets_ann(2:end)./SF; %

%--------------Preprocessing, Artefact Removal, Peak Detection settings fro
%              Physionet scripts

HRVparams.windowlength = 300;	      % Default: 300, seconds
HRVparams.increment = 30;             % Default: 30, seconds increment
HRVparams.numsegs = 5;                % Default: 5, number of segments to collect with lowest HR
HRVparams.RejectionThreshold = .20;   % Default: 0.2, amount (%) of data that can be rejected before a
% window is considered too low quality for analysis
HRVparams.MissingDataThreshold = .15; % Default: 0.15, maximum percentage of data allowable to be missing
% from a window .15 = 15%

HRVparams.preprocess.figures = 0;                   % Figures on = 1, Figures off = 0
HRVparams.preprocess.gaplimit = 2;                  % Default: 2, seconds; maximum believable gap in rr intervals
HRVparams.preprocess.per_limit = 0.4;               % Default: 0.2, Percent limit of change from one interval to the next
HRVparams.preprocess.forward_gap = 3;	            % Default: 3, Maximum tolerable gap at beginning of timeseries in seconds
HRVparams.preprocess.method_outliers = 'rem';       % Default: 'rem', Method of dealing with outliers
% 'cub' = replace outlier points with cubic spline method
% 'rem' = remove outlier points
% 'pchip' = replace with pchip method
HRVparams.preprocess.lowerphysiolim = 60/160;       % Default: 60/160
HRVparams.preprocess.upperphysiolim = 60/30;        % Default: 60/30
HRVparams.preprocess.method_unphysio = 'rem';       % Default: 'rem', Method of dealing with unphysiologically low beats
% 'cub' = replace outlier points with cubic spline method
% 'rem' = remove outlier points
% 'pchip' = replace with pchip method

% The following settings do not yet have any functional effect on
% the output of preprocess.m:
HRVparams.preprocess.threshold1 = 0.90 ;	        % Default: 0.9, Threshold for which SQI represents good data
HRVparams.preprocess.minlength = 30;            % Default: 30, The minimum length of a good data segment in seconds



%% Clean RR intervals (Here RR means IBI)

% Remove Large intervals Caused by Gaps
% These are not counted towards the total signal removed
preprocess_upperphysiolim=2 ;%*SF; % maximum believable gap in rr intervals
preprocess_lowerphysiolim=0.5 ;% *SF;
idx_remove = find(rr >= preprocess_upperphysiolim);
rr(idx_remove) = [];
t_rr(idx_remove) = [];
clear idx_remove;

% Find RR Over Given Percentage Change : TRESHOLD HERE. necessary???
perLimit=0.4;   % Percent limit of change from one interval to the next
idxRRtoBeRemoved = FindSpikesInRRseries(rr, perLimit);

% Combine Annotations and Percentage Outliers
% Combination of methods
outliers_combined = idxRRtoBeRemoved(:); % + ann_outliers(:);
outliers = logical(outliers_combined);

% Extrapolation of outliers
preprocess_method_outliers='rem'; %'cub';
idx_outliers = find(outliers == 1);
% Keep count of outliers
numOutliers = length(idx_outliers);
%
rr_original = rr;
rr(idx_outliers) = NaN;


switch preprocess_method_outliers
    case 'cub'
        NN_Outliers = interp1(t_rr,rr,t_rr,'spline','extrap');
        t_Outliers = t_rr;
    case 'pchip'
        NN_Outliers = interp1(t_rr,rr,t_rr,'pchip');
        t_Outliers = t_rr;
    case 'lin'
        NN_Outliers = interp1(t_rr,rr,t_rr,'linear','extrap');
        t_Outliers = t_rr;
    case 'rem'
        NN_Outliers = rr;
        NN_Outliers(idx_outliers) = [];
        t_Outliers = t_rr;
        t_Outliers(idx_outliers) = [];
    otherwise % By default remove outliers
        NN_Outliers = rr;
        NN_Outliers(idx_outliers) = [];
        t_Outliers = t_rr;
        t_Outliers(idx_outliers) = [];
end


% Identify Non-physiologic Beats and Remove
preprocess_method_unphysio='rem';

toolow = NN_Outliers < preprocess_lowerphysiolim;
idx_toolow = find(toolow == 1);
NN_NonPhysBeats = NN_Outliers;

%number
numOutliers = numOutliers + length(idx_toolow);

switch preprocess_method_unphysio

    case 'cub'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'spline','extrap');
        t_NonPhysBeats = t_Outliers;
        flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
    case 'pchip'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'pchip');
        t_NonPhysBeats = t_Outliers;
        flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
    case 'lin'
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'linear','extrap');
        t_NonPhysBeats = t_Outliers;
        flagged_beats = logical(outliers(:) + toohigh(:)+ toolow(:));
    case 'rem'
        NN_NonPhysBeats(idx_toolow) = [];
        t_NonPhysBeats = t_Outliers;
        t_NonPhysBeats(idx_toolow) = []; % Review this line of code for improvement
    otherwise % use cubic spline interpoletion as default
        NN_NonPhysBeats = interp1(t_Outliers,NN_NonPhysBeats,t_Outliers,'pchip');
        t_NonPhysBeats = t_Outliers;
end

%  Interpolate Through Beats that are Too Fast (if we need to interpolate)
toohigh = NN_NonPhysBeats > preprocess_upperphysiolim;
idx_outliers_2ndPass = find(toohigh==1);
NN_TooFastBeats = NN_NonPhysBeats;
NN_TooFastBeats(idx_outliers_2ndPass) = NaN;

% numOutliers = numOutliers + length(idx_outliers_2ndPass);
% if strcmp(preprocess_method_unphysio,'rem')
%     flagged_beats = numOutliers;
% end

switch preprocess_method_outliers     % switch HRVparams.preprocess.method_outliers
    case 'cub'
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'spline','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
    case 'pchip'
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'pchip');
        t_TooFasyBeats = t_NonPhysBeats;
    case 'lin'
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'linear','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
    case 'rem'
        NN_TooFastBeats(idx_outliers_2ndPass) = [];
        t_TooFasyBeats = t_NonPhysBeats;
        t_TooFasyBeats(idx_outliers_2ndPass) = []; % Review this line of code for improvement
    otherwise % USe cubic spline interpoletion as default
        NN_TooFastBeats = interp1(t_NonPhysBeats,NN_TooFastBeats,t_NonPhysBeats,'spline','extrap');
        t_TooFasyBeats = t_NonPhysBeats;
end

% Remove erroneous data at the end of a record
%       (i.e. a un-physiologic point caused by removing data at the end of
%       a record)

while NN_TooFastBeats(end) > preprocess_upperphysiolim	% equivalent to RR = 2
    NN_TooFastBeats(end) = [];
    t_TooFasyBeats(end) = [];
end

cleanNN = NN_TooFastBeats;
cleantNN = t_TooFasyBeats; % time of clean NN intervarls, in seconds, from 0 to 30 sec range 


%% Detect pulse peaks based on distance from onset
  
[cleantOO,newpeaks,cycles]=deal([]);

cleantOO=floor(cleantNN*SF);
newpeaks=nan(size(cleantOO));

for u=1:numel(cleantOO)-1;
    [pulse,peakval,imax1,cyclelen]=deal([]);
    cyclelen=((cleantOO(u+1)-cleantOO(u))/SF);
    if cyclelen> 0.5 && cyclelen < 2
        pulse=ppgZ([cleantOO(u):cleantOO(u)+floor(cyclelen*0.5*SF)-1]);%within 0.5 of the cycle length
        [pkt,lct]=findpeaks(pulse,SF);
        if ~isempty(lct)
        newpeaks(u)=cleantOO(u)+((lct(1)*SF)-1);
        end
    end
end

cleantOO(isnan(newpeaks))=nan;
cleantOO(isnan(cleantOO))=[];
newpeaks(isnan(newpeaks))=[];

% Output:
cleanOnsets=cleantOO; % annotations (sample) for cleaned onsets
cleanPeaks=newpeaks;  % annotations (sample) for cleaned peaks


end



