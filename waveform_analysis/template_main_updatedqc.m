function [Results_Nmin] = template_main_updatedqc(ppg_clean,SCORES2use,parameters,SubjInfo)
% This function finds template and matching pulses

% This function generates the pulse template based on the input parameters:
%% Dependencies:
% find_OURAPPGonsets.m
% create_pulseTEMPLATE_ver1_updatedqc.m
% extract_pulseFEATURES_updatedqc.m

%% Inputs:
% ppg_clean: filtered ppg signal
% SCORES2use: sleep staging file to get timings of 30 sec windows
% parameters: pre-defined analysis parameters (see define_parameters.m)
% SubjInfo: a struct containing participant name, age, bmi, height and weight



%% Created at: 10 Sept by GYD
%% Last modified: @ 10 sept 2024 by Gizem:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    error('Wrong number of input arguments');
end

Results_Nmin=[];
SFdevice=parameters.SFdevice; % sampling rate of the device
SFupsample=parameters.extract_features.SFupsample;

ppg_t=(0:1/SFdevice:(length(ppg_clean)-1)/SFdevice)';

switch parameters.template_method

    case "allSTAGES-NsecOverlap"


        % Important points of pulse: % Different PPG waveforn features;
        
        pointsPPG_varnames={'onset','max_slope_point','systolic_peak','dicrotic_notch','diastolic_peak','offset'};
        feats_varnames={'PW','CT','dT',...
            'sysAmp','dicAmp','diasAmp',...
            'dT_norm', 'CT_norm',...
            'AS','AS_norm',... % Ascending slope
            'RefInd','RefInd_norm',...      % Reflection Index,Formerly used as:'AIx2dias_zsc','AIx2dia_zsc100' in ppg manuscript
            'SI','SI_norm',...
            'AIx', 'AIx_norm'};


        %% Here, we will go trough overlapping windows to extract PPG features:
        % Create 30-sec array of scores
        [starting_ind,win_start,scores_t]=deal([]);
        scores_t=seconds(SCORES2use.timestamp-SCORES2use.timestamp(1));

        % Setup overlapping windows
        [nx,Nwinds,tWin_start,tWin]=deal([]);
        %         Nmin=1 ; %length of windows, in minutes
        %         win_len=60*Nmin; % Nmin converted to seconds.
        %         overlap=win_len-30; %sec
        %         increment=win_len-overlap;

        win_len=parameters.extract_features.win_len; % Nmin converted to seconds.
        overlap=parameters.extract_features.overlap; %sec
        increment=win_len-overlap;

        nx = numel(ppg_clean)/SFdevice;               % length of sequence
        Nwinds = fix((nx-overlap)/(win_len-overlap));    % number of sliding windows
        tWin_start= [(0:(Nwinds-1))*(win_len-overlap)]';  % starting time (sec) of each window

        tWin= cell2mat(arrayfun(@(a) (linspace(a,a+(overlap),(win_len/increment))), tWin_start,'uni',false)); % within an Nmin window, these are starting times of 30 sec windows. This might be necessary while storing data from 30 sec windows

        % Preallocate nan/empty arrays:
        STAGE_features=[];
        STAGE_NCCvalues=[];

        % Store for N-min only;
        pulse_templates= num2cell(nan(size(tWin,1),1));
        NCC_fortemplate=(nan(size(tWin,1),1));
        Window_ratio_miss=(nan(size(tWin,1),1));

        % Store for every 'increment- sec within Nmin (crucial especially if we define overlapping windows):
        stages=num2cell(nan(size(tWin)));stages(:)={"NaN"};
        onsets_realtime=num2cell(nan(size(tWin)));
        % onsets_extracted_final=num2cell(nan(size(tWin)));

        pointsPPG_extracted_final=num2cell(nan(size(tWin)));
        pointsPPG_extracted_final(:)={nan(1,size(pointsPPG_varnames,2))};

        features_extracted_final=num2cell(nan(size(tWin)));
        features_extracted_final(:)={nan(1,size(feats_varnames,2))};

        NCCmax_mean_final=num2cell(nan(size(tWin)));
        NCCmax_fitted_final=num2cell(nan(size(tWin)));
        dtw_norm_final=num2cell(nan(size(tWin)));
        sim_score_final=num2cell(nan(size(tWin)));
        raw_pulse_length_final=num2cell(nan(size(tWin)));

        pulse_waveforms_extracted_final=num2cell(nan(size(tWin)));pulse_waveforms_extracted_final(:)={{nan(1,parameters.extract_features.normalized_width)}};


        for d=1:size(tWin,1)

            display(['... current window: ' num2str(d) ' / ' num2str(size(tWin,1))  ' ...'])

            clear win_start STAGES;
            [pulse_waveforms,pointsPPG,feats_derived, NCC_max,TEMPLATE,gaussfit]=deal([]);

            % Get window-timing info:
            win_start=tWin(d,1);% starttime of  window, in seconds

            % Store window's stage information
            stages{d,1}=string(unique(SCORES2use.categorical(scores_t >= win_start & scores_t < win_start + increment)));

            % if overlap ~=0
            % stages{d,2}=string(unique(SCORES2use.categorical(scores_t >= (win_start+increment) & scores_t < win_start + win_len)));
            % end

            % Isolate PPG seg in the Nmin window first:
            clear ppg_d;
            ppg_d=ppg_clean(ppg_t>=win_start & ppg_t<win_start+win_len);

            ppg_d_t=(0:1/SFdevice:(length(ppg_d)-1)/SFdevice)';

            % Calculate the ratio of missing data within the 5 min window.
            % Continue if there is enough data
            ratio_miss=[];
            ratio_miss= numel(find(isnan(ppg_d))) / numel(ppg_d);
            Window_ratio_miss(d)=ratio_miss;

            %% Extract the pulse template and SQI:
            if 1 - ratio_miss >= parameters.data_ratio % continue only if the amount of data > data ratio defined in the parameters

                % Extract onsets for this ppg segment :
                [onsets_d,peaks_d]=deal([]);

                if range(ppg_d)>= 1
                    try
                        [onsets_d,peaks_d] = find_OURAPPGonsets(ppg_d,SFdevice,0);
                    end
                end

                % % Check ppg segment and extracted onsents:
                % figure;plot(ppg_d_t,ppg_d,'k');hold on;plot(ppg_d_t(onsets_d),ppg_d(onsets_d),'ro');


                if numel(onsets_d)>=5  %  at least 5 pulses required within a window

                    % Detrend // Remove baseline drift usign cubic spline (Yang et al 2014; Gonzales et al 2013);
                    [onsets_up_t,onsets_amp,onsets_amp_up,baselinedrift,ppg_i_win_nodrift]=deal([]);
                    onsets_amp=ppg_d(onsets_d);
                    onsets_up_t = (onsets_d(1)/SFdevice):1/SFdevice:(onsets_d(end)/SFdevice);
                    onsets_amp_up = interp1((onsets_d/SFdevice),onsets_amp, onsets_up_t,'spline')'; % cubic spline interpolation
                    %
                    % figure;plot(ppg_d_t,ppg_d,'k');hold on;plot(onsets_up_t,onsets_amp_up,'m'); hold on;
                    % xline((tWin(d,:)-tWin(d,1)),'b--')

                    baselinedrift=zeros(size(ppg_d));
                    baselinedrift(onsets_d(1):onsets_d(end))=onsets_amp_up;
                    ppg_d_nodrift=ppg_d-baselinedrift;
                    %

                    %% Find the template pulse for this ppg segment:
                    TEMPLATE=nan(parameters.extract_features.normalized_width,3);

                    try
                        [TEMPLATE,gaussfit, ncc_meanvsfitted] = create_pulseTEMPLATE_ver1_updatedqc(ppg_d_nodrift,onsets_d,parameters);
                    end
                    %
                    % figure;plot([1:parameters.extract_features.normalized_width],TEMPLATE(:,2),'k');
                    % hold on;
                    % plot([1:parameters.extract_features.normalized_width],TEMPLATE(:,3),'m');
                    % ylabel('Amplitude');xlabel('Pulse Width');
                    % legend({'Mean Pulse','Gaussian-fitted template'});
                    % title(['Pulse template from a ',num2str(SubjInfo.subjAge),' yo, CC=',num2str(round(gaussfit.GOF(2),3)) ]);
                    % set(gca,'FontSize',14);


                    %% Extract pulse waveform features:
                    if all(~isnan(TEMPLATE))
                        % Change any parameters?

                        % Extract pulses and waveform features:
                        % [feats_derived,pointsPPG,pulse_waveforms,NCC_max,pointsPPG_varnames,feats_varnames,~]=...
                        %     extract_pulseFEATURES(ppg_d_nodrift,onsets_refined,TEMPLATE,gaussfit,SubjInfo.subjH,parameters); % we use the z-scored ppg segment here!!

                        parameters.extract_features.ploton=0;
                        % stages{d,1}

                        [feats_derived,pointsPPG,pulse_waveforms,NCCmax_mean, NCCmax_fitted,dtw_norm,sim_score,pointsPPG_varnames,feats_varnames,raw_pulse_length] = ...
                            extract_pulseFEATURES_updatedqc(ppg_d_nodrift,onsets_d, SubjInfo.subjH,parameters,TEMPLATE,gaussfit);



                        % Store results:
                        w=[];
                        for w=1:size(tWin,2)

                            % ind_ons=[];
                            % ind_ons=ppg_d_t(onsets_refined(1:end-1)) >= tWin(1,w) & ppg_d_t(onsets_refined(1:end-1))< tWin(1,w)+increment;
                            %
                            onsets_realtime{d,w}=onsets_d+ ((win_start*SFdevice));
                            % onsets_extracted_final{d,w}=onsets_refined;
                            pointsPPG_extracted_final{d,w}=pointsPPG; % actually this is not exactly correct!
                            features_extracted_final{d,w}=feats_derived;
                            pulse_waveforms_extracted_final{d,w}=pulse_waveforms;

                            NCCmax_mean_final{d,w}=NCCmax_mean;
                            NCCmax_fitted_final{d,w}=NCCmax_fitted;
                            dtw_norm_final{d,w}=dtw_norm;
                            sim_score_final{d,w}=sim_score;
                            raw_pulse_length_final{d,w}=raw_pulse_length;

                        end

                        %                         % Save results:
                        %                         STAGE_features=[STAGE_features;feats_derived];
                        %                         STAGE_NCCvalues=[STAGE_NCCvalues;NCC_max];

                    end

                else

                    display(['.. no onsets, skipping to the next window']);
                end
            end


            Results_Nmin.window_templates{d}= TEMPLATE;
            Results_Nmin.window_gaussfit{d}= gaussfit;


        end

        %
        % % Following part is to correct for an error in saving results
        % % after feature extraction. Basically, here I am keeping values if they
        % % have respective PPG-points for the same window
        % for u=1:numel(tWin);
        %
        %     clear curr_onsets curr_points curr_feats curr_NCC curr_waveforms;
        %     curr_onsets_realtime= onsets_realtime{u};
        %     curr_onsets= onsets_extracted_final{u};
        %     curr_points= pointsPPG_extracted_final{u};
        %     curr_feats=features_extracted_final{u};
        %     curr_NCC=NCC_max_extracted_final{u};
        %     curr_waveforms=pulse_waveforms_extracted_final{u};
        %
        %     if any(~isnan(curr_onsets)) && any(~isnan(curr_points(:,1)))
        %
        %         onsets_realtime{u}=curr_onsets_realtime(~isnan(curr_points(:,1)));
        %         onsets_extracted_final{u}=curr_onsets(~isnan(curr_points(:,1)));
        %         pointsPPG_extracted_final{u}=curr_points(~isnan(curr_points(:,1)),:);
        %         features_extracted_final{u}=curr_feats(~isnan(curr_points(:,1)),:);
        %         NCC_max_extracted_final{u}=curr_NCC(~isnan(curr_points(:,1)),:);
        %         pulse_waveforms_extracted_final{u}=cell2mat(curr_waveforms(~isnan(curr_points(:,1)))')';
        %
        %     else
        %
        %         onsets_extracted_final{u}=nan;
        %         pointsPPG_extracted_final{u}=nan(1,size(pointsPPG_varnames,2));
        %         features_extracted_final{u}=nan(1,size(feats_varnames,2));
        %         NCC_max_extracted_final{u}=nan;
        %         pulse_waveforms_extracted_final{u}=nan(1,200);
        %     end
        %
        % end
        %
        %

        %% Save as: 'Nmin windows' format (as is): ========================
        Results_Nmin.Stages=stages;
        Results_Nmin.Onsets_realtime=onsets_realtime;
        Results_Nmin.PointsPPG_out=pointsPPG_extracted_final;
        Results_Nmin.Features_out=features_extracted_final;

        Results_Nmin.NCCmax_mean_out= NCCmax_mean_final;
        Results_Nmin.NCCmax_fitted_out= NCCmax_fitted_final;
        Results_Nmin.dtw_norm_out=dtw_norm_final;
        Results_Nmin.sim_score_out=sim_score_final;
        Results_Nmin.raw_pulse_length_out=raw_pulse_length_final;

        Results_Nmin.PulseWaveforms_out=pulse_waveforms_extracted_final;


        %         if overlap ~= 0
        %
        %         %% Save as: 'overlapping windows are combined' format =============
        %         clear reshaped_onsets reshaped_pointsPPG reshaped_features reshaped_NCC_max ...
        %             reshaped_pulse_waveforms tWin_reshaped stages_reshaped;
        %         clear reshaped_onsets_cellarray reshaped_pointsPPG_cellarray reshaped_features_cellarray reshaped_NCC_max_cellarray...
        %             reshaped_pulse_waveforms_cellarray stages_reshaped_cellarray
        %
        %         % Windows used in the analysis
        %         tWin_reshaped=[tWin(1):increment:tWin(end)]';
        %
        %         % Sleep stages for every 30 sec window:
        %         stages_reshaped_cellarray=combine_overlaps(stages,tWin,"stages");
        %         stages_reshaped=table2array(cell2table(cellfun(@(x) (x(1)),stages_reshaped_cellarray(:),'UniformOutput',false)));
        %         stages_reshaped(find(stages_reshaped=='WakeAfter')) = categorical({'Wake'});
        %         stages_reshaped=removecats(stages_reshaped,'WakeAfter');
        %
        %         % Some definitions:
        %         % Individual-pulse level onsets (annotations) which were used to get features:   onsets_extracted_final
        %         % Individual-pulse level PPG fidducial points:                                   pointsPPG_extracted_final
        %         % Individual feature values :                                                    features_extracted_final
        %         % Individual-pulse level NCC values:                                             NCC_max_extracted_final
        %         % Individual-pulse waveforms which were used to get features:                    pulse_waveforms_extracted_final
        %
        %         % Reshape onsets_extracted_final
        %         reshaped_onsets_cellarray=combine_overlaps(onsets_extracted_final,tWin,"onsets");
        %         reshaped_onsets= cell2mat(cellfun(@(x) (x),reshaped_onsets_cellarray(:),'UniformOutput',false));
        %         reshaped_onsets=reshaped_onsets(~isnan(reshaped_onsets));
        %         % Reshape pointsPPG_extracted_final
        %         reshaped_pointsPPG_cellarray=combine_overlaps(pointsPPG_extracted_final,tWin,"pointsPPG");
        %         reshaped_pointsPPG= cell2mat(cellfun(@(x) (x), reshaped_pointsPPG_cellarray(:),'UniformOutput',false));
        %         reshaped_pointsPPG=reshaped_pointsPPG(~isnan(reshaped_pointsPPG(:,1)),:);
        %         % Reshape features_extracted_final
        %         reshaped_features_cellarray=combine_overlaps(features_extracted_final,tWin,"features");
        %         % % First get median of each 30 sec window:
        %         reshaped_features_30secMEDIAN= cell2mat(cellfun(@(x) median(x,1,'omitnan'), reshaped_features_cellarray(:),'UniformOutput',false));
        %         % % Second, get individual values:
        %         reshaped_features= cell2mat(cellfun(@(x) (x),   reshaped_features_cellarray(:),'UniformOutput',false));
        %         reshaped_features=  reshaped_features(~isnan(  reshaped_features(:,1)),:);
        %         % Reshape NCC_max_extracted_final
        %         reshaped_NCC_max_cellarray=combine_overlaps(NCC_max_extracted_final,tWin,"NCC_max");
        %         reshaped_NCC_max= cell2mat(cellfun(@(x) (x), reshaped_NCC_max_cellarray(:),'UniformOutput',false));
        %         reshaped_NCC_max= reshaped_NCC_max(~isnan( reshaped_NCC_max));
        % %         % Reshape pulse_waveforms_extracted_final
        % %         reshaped_pulse_waveforms_cellarray=combine_overlaps(pulse_waveforms_extracted_final,tWin,"pulse_waveforms");
        % %         reshaped_pulse_waveforms= cell2mat(cellfun(@(x) (x), reshaped_pulse_waveforms_cellarray(:),'UniformOutput',false));
        % %         reshaped_pulse_waveforms=reshaped_pulse_waveforms(~isnan(reshaped_pulse_waveforms(:,1)),:);
        %
        %         % Store inside the main struct:
        %         % % As a double-array (matrix):
        %         OUT.Results_withOverlap.Stages=stages_reshaped;
        %         OUT.Results_withOverlap.tWin=tWin_reshaped;
        %         OUT.Results_withOverlap.asMTX.Onsets_out=reshaped_onsets;
        %         OUT.Results_withOverlap.asMTX.PointsPPG_out=reshaped_pointsPPG;
        %         OUT.Results_withOverlap.asMTX.Features_out=reshaped_features;
        %         OUT.Results_withOverlap.asMTX.NCC_pulses_out=reshaped_NCC_max;
        %         % % Also in the cell array format (in case we need to trace back to stages)
        %         OUT.Results_withOverlap.asCELL.Onsets_out=reshaped_onsets_cellarray;
        %         OUT.Results_withOverlap.asCELL.PointsPPG_out=reshaped_pointsPPG_cellarray;
        %         OUT.Results_withOverlap.asCELL.Features_out=reshaped_features_cellarray;
        %         OUT.Results_withOverlap.asCELL.NCC_pulses_out=reshaped_NCC_max_cellarray;
        %
        % % %         % Save : 'all extracted features are concatenated in an array'
        % % %         OUT.Features_array=STAGE_features; % all PPG features concatenated for current stage, as a matrix
        % % %         OUT.NCCvalues_array= STAGE_NCCvalues; % all corresponding PPG features concatenated for current stage, as a matrix
        %         end


        %% Save other info ===============================================

        Results_Nmin.scores_categorical=SCORES2use;
        Results_Nmin.window_times=tWin;
        Results_Nmin.used_feature_variablenames=feats_varnames;
        Results_Nmin.used_pointPPG_variablenames=pointsPPG_varnames;
        Results_Nmin.used_parameters=parameters;
        Results_Nmin.extractdate=datetime('today');



end



if ~isempty(Results_Nmin)
    display(['PPG Feature extraction for ',SubjInfo.SUBID,'  is COMPLETE']);
    display(['***********************************************************************************************']);

else

    display(['NO PPG Feature extraction for ',SubjInfo.SUBID,'  !']);
    display(['x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x ']);
end

end

