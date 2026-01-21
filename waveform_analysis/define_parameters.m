% How to define parameters for analysis:

%% Define parameters to be used in analysis
% Parameters for main analysis
parameters.dataset="FINGERTIP_PPG"; % 
parameters.template_method='allSTAGES-NsecOverlap' ; % Method to get template: 'allSTAGES-5minEPOCHs','fromDEEP','fromBASELINE' or 'allSTAGES'
parameters.SFdevice=256; % if device is Oura, use Oura's sampling rate (e.g. 50 Hz). If the device is Finger-tip pleth, use PSG-sampling rate (e.g. 256 Hz)
parameters.SFoura=50; % if device is Oura, use Oura's sampling rate (e.g. 50 Hz). If the device is Finger-tip pleth, use PSG-sampling rate (e.g. 256 Hz)
parameters.getFEATURESforBP=1; 
parameters.data_ratio=0.5;  % minimum ratio of available data expected to be in the window to be considered for further analysis
% % parameters.consecNmin=nan; % Window length for consecutive stage search (minutes)
% % parameters.minN_highQ=0.3; %Window acceptace criteria: minimum ratio of high quality pulses within a window
parameters.quality_threshold=0.95; % Normalized cross correlation threshold to accept pulse as good quality. This value is arbitrary, can be increased or decreased!
parameters.data_ratio=0.5;  % minimum ratio of available data expected to be in the window to be considered for further analysis
% % parameters.consecNmin=nan; % Window length for consecutive stage search (minutes)
% % parameters.minN_highQ=0.3; %Window acceptace criteria: minimum ratio of high quality pulses within a window
parameters.extract_features.do_extract=1;
parameters.extract_features.SFdevice=256;
parameters.extract_features.pulse_width_max=1.7; %Max posssible Pulse Width (sec)
parameters.extract_features.normalized_width=200; % samples, fixed width
parameters.extract_features.upsample=0; % upsample to 256 hz if 1, leave in original SF if 0
parameters.extract_features.SFupsample=256; %upsampling frequency if upsampling paramater is selected
parameters.extract_features.ploton=0;
parameters.extract_features.save_pulses=1; % save individual pulse waveforms
parameters.extract_features.normalize_amplitude=0;  %if 1, normalize amplitude to 1. If 0, leave as z-scored original amplitude
parameters.extract_features.normalize_amplitude_range=[0,1];  %
parameters.extract_features.win_len=30; %sec
parameters.extract_features.overlap=0; %sec
parameters.extract_features.envelope_factor=10;
parameters.main_figureon=0; % Plot the timeseries figure with different channels and stages?
parameters.saveresults=1;