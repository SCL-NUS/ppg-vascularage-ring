function [feats_derived,pointsPPG,pulse_waveforms, NCCmax_mean, NCCmax_fitted,dtw_norm,sim_score,pointsPPG_varnames,feats_varnames,raw_pulse_length] = ...
    extract_pulseFEATURES_updatedqc(ppg_seg,onsets,pt_height,parameters,TEMPLATE,gaussfit)


%% INPUTS:
%        pulses_array: array of pulses from oura
%        pt_height: participant's height (cm), needed to calculate stiffness index
%        paramaters:
%        TEMPLATE:

%% OUTPUTS



%``````````````````````````````````````````````````````````````
% Last updated by Gizem Y.D. @ 2024.10.10


%%%%--------------------------------------------------------------------------------------------------------------
if nargin < 5
    error('Wrong number of input arguments');
end

warning('off','all');


%% Define and set parameters:
SFdevice= parameters.extract_features.SFdevice;
upsample= parameters.extract_features.upsample; % upsample to SFupsample Hz if 1. Do nothing if 0.
SFupsample= parameters.extract_features.SFupsample;
ploton= parameters.extract_features.ploton;
save_pulses= parameters.extract_features.save_pulses;
pulse_width_max=parameters.extract_features.pulse_width_max;
normalized_width=parameters.extract_features.normalized_width;

% Important points of pleth
pointsPPG_varnames={'onset','max_slope_point','systolic_peak','dicrotic_notch','diastolic_peak','offset'};
% pointsPPG=num2cell(nan(numel(onsets)-1,1));
pointsPPG=deal(nan(numel(onsets),size(pointsPPG_varnames,2)));
error=num2cell(nan(2,numel(onsets)));


feats_varnames={'PW','CT','dT',...
    'sysAmp','dicAmp','diasAmp',...
    'dT_norm', 'CT_norm',...
    'AS','AS_norm',... % Ascending slope
    'RefInd','RefInd_norm',...      % Reflection Index,Formerly used as:'AIx2dias_zsc','AIx2dia_zsc100' in ppg manuscript
    'SI','SI_norm',...
    'AIx', 'AIx_norm'};

feats_derived=deal(nan(numel(onsets),size(feats_varnames,2)));

NCCmax_mean=deal(nan(numel(onsets),1));
NCCmax_fitted=deal(nan(numel(onsets),1));
dtw_norm=deal(nan(numel(onsets),1));
sim_score=deal(nan(numel(onsets),1));

raw_pulse_length=deal(nan(numel(onsets),1));

pulse_waveforms=cell(numel(onsets),1);pulse_waveforms(:)={NaN(1,normalized_width)};



if isempty(find(isnan(TEMPLATE)))
    % Use TEMPLATE to select which pulses to include:
    template_t=TEMPLATE(:,1)';
    template_original=TEMPLATE(:,2)';
    template_fitted=TEMPLATE(:,3)';

    % get approximate peak locations:
    % [FitResults,GOF,baseline,coeff,residual,xi,yi,BootstrapErrors]=peakfit([[1:numel(TEMPLATE(:,3))]',TEMPLATE(:,3)],0,0,2,1,0, 0, 0, 0, 0,0);

    FitResults=gaussfit.FitResults;
    [~,i1]=min(abs(template_t-FitResults(1,2)));
    [~,i2]=min(abs(template_t-FitResults(2,2)));

    FitResults(1,2)=i1;
    FitResults(2,2)=i2;


    PDs=[];
    PDs=diff(onsets)/SFdevice; %pulse durations

    medianPD=median(PDs(PDs < pulse_width_max));
    PD_upper=prctile(PDs(PDs < pulse_width_max),[97]);


    %% Extract  waveform features for each individual PPG-pulse:
    for p= 1:numel(onsets)-1;

        try

            % [pulsePL,pulsePL_t]=deal([]);
            % [pulsePL,pulsePL_t,pulsePL_1dv,pulsePL_2dv]=deal([]);
            %
            % pulsePL=pulses_array(p,:);
            % pulsePL_t=[1:numel(pulsePL)];

            [currpulse,pulsePL,pulsePL_t]=deal([]);
            [pulsePL,pulsePL_t,pulsePL_1dv,pulsePL_2dv]=deal([]);


            if (onsets(p+1)-onsets(p))/SFdevice < PD_upper % does curren pulse have the offset point? No means, next pulse was discarded during SQI thresholding.

                currpulse=ppg_seg(onsets(p):onsets(p+1));

            else

                % if this pulse is a lonely pulse without any high quality
                % neighbouring pulses, we can follow the onset using simple find
                % peaks function
                pulse_template=[]; onset_defined=[];
                delay=PD_upper/3; %sec
                pulse_template=ppg_seg(onsets(p)+delay*SFdevice:onsets(p)+(1+delay)*SFdevice); % assign a fixed width

                loc=[];
                [~,loc]=findpeaks(-pulse_template,'SortStr','descend');

                onset_defined=  loc(1)+ (delay*SFdevice)-1;

                currpulse=ppg_seg(onsets(p):onsets(p)+onset_defined);

            end

            % Check if onset points is TRUE onset:
            [val, ind]=min(currpulse(1:10));
            if ind~=1 % it means onset points doesn't have the lowest amplitude -> so FALSE onset
                currpulse=currpulse(ind:end);
            end


            currpulse_t=(0:1/SFdevice:(length(currpulse)-1)/SFdevice)'; %time array for the pulse
            
            raw_pulse_length(p)=length(currpulse);

            %% Resample currpulse to normalized_width samples (i.e. normalize duration)

            [pulseZrawdur,tpulse,tpulse_new]= deal([]);
            [F,pulse200,sampratio]=deal([]);

            pulseZrawdur=currpulse;
            tpulse = currpulse_t; %original time arrat of the pulse
            tpulse_new = (linspace(0,max(tpulse),normalized_width))'; % to normalize to 200 equally sampled intervals
            F = griddedInterpolant(tpulse,pulseZrawdur,'spline','spline'); %interpolation and extrapolation grid
            pulse200 = (F(tpulse_new));
            sampratio=normalized_width/(length(pulseZrawdur)-1);

  
            %% Normalize pulse amplitude to 0 -100 range (if needed)
            [pulsePL,pulsePL_ampnormed]=deal([]);
            pulsePL=pulse200';
            pulsePL_ampnormed=normalize(pulsePL,"range",[0 100]);pulsePL=pulsePL-pulsePL(1);
            %pulsePL_ampnormed=normalize(pulsePL,"range",[0 1]);pulsePL=pulsePL-pulsePL(1);

            if parameters.extract_features.normalize_amplitude
            pulsePL=pulsePL_ampnormed;
            end

            pulsePL_t=[1:normalized_width]; % unit in samples

            %% Check cross-corr between current pulse and template:

            cut1=1;
            cut2=120;

            % Nrmalized cross-corr between 2 duration normalized signals
            [Cmean,Cfitted,max_NCC_temp,max_NCC_temp2]=deal([]);
            
            % [Cmean]=normxcorr2(template_original(cut1:cut2),pulsePL_ampnormed);
            % max_NCC_temp=round(max(Cmean),2);
            % NCCmax_mean(p)=max_NCC_temp;

            [Cmean, lags_mean] = xcorr(template_original(cut1:cut2),pulsePL_ampnormed(cut1:cut2), 'coeff');
            max_NCC_temp=round(max(Cmean),2);
            NCCmax_mean(p)=max_NCC_temp;

            % [Cfitted]=normxcorr2(template_fitted(cut1:cut2),pulsePL_ampnormed);
            % max_NCC_temp2=round(max(Cfitted),2);
            % NCCmax_fitted(p)=max_NCC_temp2;

            [Cfitted, lags_fitted] = xcorr(template_fitted(cut1:cut2),pulsePL_ampnormed(cut1:cut2), 'coeff');
            max_NCC_temp2=round(max(Cfitted),2);
            NCCmax_fitted(p)=max_NCC_temp2;






            % Compute the DTW distance and alignment
            [dist,ix,iy,dist_raw,path_length]=deal([]);
            [dist,ix,iy] = dtw(template_original(cut1:cut2),pulsePL_ampnormed(cut1:cut2),'euclidean');
            dist_raw=dist;
            path_length=length(ix);
            dtw_norm(p)= round(dist_raw/ path_length,2);

            % [dist,ix,iy,dist_raw,path_length]=deal([]);
            % [dist,ix,iy] = dtw(template_fitted(cut1:cut2),pulsePL(cut1:cut2));
            % dist_raw=dist;
            % path_length=length(ix);
            % dtw_norm(p)= round(dist_raw/ path_length,2);


            % Compute similarity score:
            % anything below 0 is unacceptable
            clear curr_similarity_score;
            [curr_similarity_score] = SQrank(NCCmax_mean(p),dtw_norm(p));

            sim_score(p)=curr_similarity_score;


            % figure('Position',[400 400 800 800]);
            % % subplot(1,2,1)
            % plot(pulsePL_t,pulsePL_ampnormed,'k','LineWidth',1);hold on;
            % plot(pulsePL_t(cut1:cut2),pulsePL_ampnormed(cut1:cut2),'k.','MarkerSize',8);hold on;
            % plot(pulsePL_t,template_original,'r','LineWidth',1);hold on;
            % plot(pulsePL_t(cut1:cut2),template_original(cut1:cut2),'r.','MarkerSize',8);
            % plot(pulsePL_t,template_fitted,'m','LineWidth',0.5);hold on;
            % % plot(template_t(cut1:cut2),template_fitted(cut1:cut2),'m.','MarkerSize',8); ylim([0 105]);
            % ylim([0 105]);
            % xlim([0 normalized_width]);
            % xlabel('Duration (samples)');ylabel('Normalized amplitude');
            % legend({'pulse',"","mean","","fitted"});
            % % title(['pulse',num2str(p), ' NCCmean ',num2str(NCCmax_mean(p)), ' NCCfit ',num2str(NCCmax_fitted(p)),' DTW ',num2str(distance)]);
            % title(['p',num2str(p),' ncc: ',num2str(  NCCmax_mean(p)),' dtw:',num2str(dtw_norm(p)),...
            %     ' score:',num2str(curr_similarity_score)]);


            if NCCmax_fitted(p) >= 0.98 && ...
                    sim_score(p)>=0 % && ... !! UPDATE according to needs of project!!!
                    %abs(pulsePL(end)-pulsePL(1)) <=10 && ... % %10 difference is arbitrary
                    % max(pulsePL)>=90
                % % %
              

                    %% Prepare Pulse for feature extraction
                    [pulse200]=deal([]);
                    pulse200=pulsePL;

                    % First & sevond derivatives of the 'pulse'
                    pulse_1dv =smooth(diff(pulse200),0.1,'sgolay',3);
                    pulse_2dv = smooth(diff(pulse_1dv),0.1,'sgolay',3);
                    pulse_3dv = smooth(diff(pulse_2dv),0.1,'sgolay',3);
                    pulse_4dv = smooth(diff(pulse_3dv),0.1,'sgolay',3);

                    pulse_1dv=[nan;pulse_1dv];
                    pulse_2dv=[nan;nan;pulse_2dv];
                    pulse_3dv=[nan;nan;nan;pulse_3dv];
                    pulse_4dv=[nan;nan;nan;nan;pulse_4dv];

                    [indons,indmsl,indsys,inddic,inddias,indend,inddic_alt]=deal(nan);

                    %%  Extract main points on PPG Waveform
                    %% Find Onset and offset-------------------------
                    indons=1;
                    indend=numel(pulse200);

                    %% Find max slope point------------------------------------
                    % Point where 1st derivative has 1st peak
                    [pk,loc,indmsl]=deal([]);
                    [pk,loc]=findpeaks(pulse_1dv,'NPeaks',1);
                    indmsl=loc;

                    %% Find Systolic Peak-----------------------
                    % Point where 1st derivative crosses zero for the 1st time (from  max slope point until end)
                    [indzc,indsys]=deal([]);
                    indzc=find(movprod(pulse_1dv(indmsl:end),2,'Endpoints','fill')<=0,1,'first');
                    indsys=(indzc-1)+indmsl-1;

                    % use the peak location of fitted template to check for
                    % the correctness of indsys:
                    % Using peak location of 1st gaussian,:
                    
                    [sys_loc,sys_amp]=deal([]);
                    [~,sys_loc]=min(abs(pulsePL_t-FitResults(1,2)));
                    sys_amp=pulsePL(sys_loc);

                    % % % we don't want indsys and sys_loc  to be too further apart (approx 0.2 sec
                    % % if (abs(indsys-sys_loc)>0.2*SFdevice) %added on sept 2024
                    % %     indsys=nan;
                    % % end
                    % % 
                    % % [syst_loc,indpeak]=deal([]);
                    % % [valx,syst_loc]=max(template_fitted);
                    % % [valpeak,indpeak]=max(pulsePL); % to check if real peak of pulsePL is close to estimated peak
                    % % % 
                    % 
                    % figure;
                    % sb1=subplot(3,1,1);plot(pulsePL_t,pulse200);
                    % sb2=subplot(3,1,2);plot([pulsePL_t],pulse_1dv);yline(0);
                    % sb3=subplot(3,1,3);plot([pulsePL_t],pulse_2dv);yline(0);




                    %% Find Diastolic peak and dicrotic notch ---------------------
                    % Force diastolic peak to be close to template's diastolic peak

                    [indw,indw2,inddic_ori,pk2,loc2,loc5,loc6]=deal([]);

                    % use the location of second gaussian peak diastolc
                    % peak: dia_loc (with 10% margin) is the estimated location of diastolic peak
                    % based on current window's template
                    [loc3,flip_pulse_3dv,dia_loc]=deal([]);
                    [~,dia_loc]=min(abs(pulsePL_t-FitResults(2,2))); 
                    dia_loc=dia_loc+round(dia_loc*0.1);


                    try

                        %                 % Exclude last 1/4 of the pulse from search
                        %                 indw2=round(numel(pulse200)*3/4);

                        % Limit dicrotic-notch search window to until estimated
                        % location of diastolic peak:
                        indw2=dia_loc;
                        % Find the 1st zero-crossing after b point/dip in the 2nd-der
                        indw=find(movprod(pulse_2dv((indmsl+2):indw2),2,'Endpoints','fill')<=0,1,'First')-1; %indmss+2 is just to make sure we skip first crossing
                        %% indw currently is not in use!!
                        
                        %% Alternative 1 for  dicrotic nothc detection:

                        % Find the last zero-crossing on 3rd derivative signal just
                        % before the predicted diastolic peak time (peak location of
                        % gauss2 )
                        flip_pulse_3dv=flipud(pulse_3dv);
                        flip_pulse_3dv_seg=flip_pulse_3dv(numel(pulse_3dv)-dia_loc+1:end);
                        loc3=find(movprod(flip_pulse_3dv_seg,2,'Endpoints','fill')<=0); %here we have
                        % locations where zero crossing occurs but we need the crossing to be from negative to positive.
                        % hence, the Nth value (indexes in loc3) itself must be positive and (N-1)th value must be negative:
                        loc3a=loc3(find(flip_pulse_3dv_seg(loc3)>=0,1,'First'));

                        % figure;subplot(2,1,1);plot(pulse200,'.');yline(0);xline(dia_loc);xline(indsys);xline(loc3a,'r');
                        % subplot(2,1,2);plot(pulse_3dv,'.');yline(0);xline(dia_loc);xline(indsys);xline(loc3a,'r');
                        % 
                        if ~isempty(loc3) && (dia_loc-loc3a+1) < dia_loc && indsys < (dia_loc-loc3a+1) %check if dia_loc is located between predicted syspeak and predicted diaspeak
                            inddic_alt=dia_loc-loc3a+1;
                        end

                        inddic=inddic_alt;

                        %% Diastolic peak:
                        % From dicrotic nothch to dicrotic notch+ 1/5 of descending wave,
                        [indzc,indneg,seg,pk,loc]=deal([]);
                        seg=inddic+round((indsys+indend)/5);

                        % Find zero crossings on the first derivative on
                        % selected segment
                        indzc=find(movprod(pulse_1dv(inddic:seg),2,'Endpoints','fill')<=0);
                        % On second derivative signal, find zero crossing
                        % where 
                        indneg=find(pulse_2dv((indzc+inddic-1)-1)<0);

                        if numel(indneg)>1
                            [~,inddias]=min(pulse_2dv(indzc+inddic-1-1));
                            inddias=indzc(inddias)+(inddic-1);
                        elseif isempty(indneg)
                            [pk,loc]=findpeaks(-pulse_2dv(inddic:seg),'NPeaks',1);
                            if isempty(loc)
                                inddias=round((inddic+indend)/2);
                            else
                                inddias=inddic+loc -1;
                            end
                        else
                            inddias=indzc(indneg)+(inddic-1);
                        end



                    end


                    if isempty(find(isnan([indons,indmsl,indsys,indend]))) %&& ...  % MODIFY this section according to the needs of data/project!!!
                                                        %abs(indsys - syst_loc) < SFdevice*0.2  %0.2 sec !! check if it applies to pleth-ppg

                            % pulse200(indsys) > pulse200(inddias) && pulse200(indsys) > pulse200(inddic) 
                            %%... && abs(indpeak-indsys) < SFdevice*0.1 % 0.1 sec !! check if it applies to pleth-ppg
                      

                        if ploton

                         fig=figure('Name',['pulse',num2str(p)],'Position',[200 400 1200 600],'WindowStyle','normal');
                         tiledlayout(fig,1,2,'TileSpacing','none');
                         t1=nexttile(1);
                        template_t=[1:normalized_width];
                       
                        plot(template_t,pulsePL_ampnormed,'k','LineWidth',1);hold on;
                        plot(template_t(cut1:cut2),pulsePL_ampnormed(cut1:cut2),'k.','MarkerSize',8);hold on;
                        plot(template_t,template_original,'r','LineWidth',1);hold on;
                        plot(template_t(cut1:cut2),template_original(cut1:cut2),'r.','MarkerSize',8);
                        plot(template_t,template_fitted,'m','LineWidth',0.5);hold on;
                        % plot(template_t(cut1:cut2),template_fitted(cut1:cut2),'m.','MarkerSize',8); ylim([0 105]);
                        ylim([0 105]);xlim([0 normalized_width]);
                        legend({'pulse',"","mean","","fitted"});
                        % title(['pulse',num2str(p), ' NCCmean ',num2str(NCCmax_mean(p)), ' NCCfit ',num2str(NCCmax_fitted(p)),' DTW ',num2str(distance)]);
                        title(t1,['p',num2str(p),' ncc: ',num2str(  NCCmax_mean(p)),' dtw:',num2str(dtw_norm(p)),...
                            ' score:',num2str(curr_similarity_score)]);

                        %Pulse Shape Plot-2: For displayin fiducial points!
                        t2=nexttile(2);
                        % fig=figure('WindowStyle','normal','visible','on');
                        plot(pulsePL_t,pulse200,'k.');hold on;
                        xlim([0 normalized_width]); hold on; %ylim([0 4]);
                        xline([round(FitResults(1,2)) ,round(FitResults(2,2))]);
                        title(t2,['PPG Pulse and Fiducial Points p=',num2str(p)]);
                        % plot(pulse_gauss_t,pulse_gauss,'g');hold on;
                        plot(pulsePL_t([indons,indend]),pulse200([indons,indend]),'ko','MarkerFaceColor','k','MarkerSize',8);
                        plot(pulsePL_t([indmsl]),pulse200([indmsl]),'ko');
                        plot(pulsePL_t([indsys]),pulse200([indsys]),'mo','MarkerFaceColor','m','MarkerSize',8);
                        % xline(pulsePL_t([indsys]),'m-.');
                        plot(pulsePL_t([inddic]),pulse200([inddic]),'ro','MarkerFaceColor','r','MarkerSize',8);
                        % xline(pulsePL_t([inddic]),'r-.');
                        plot(pulsePL_t([inddias]),pulse200([inddias]),'bo','MarkerFaceColor','b','MarkerSize',8);
                        set(gca,'Fontsize',14);

    %
                % % Pulse Shape and derivatives Plot: diagnostic & for displayin fiducial points!
                % figder=figure('position',[800,150,560,2000],'WindowStyle','normal','visible','on');
                % tiledlayout(figder,4,1,'TileSpacing','tight');
                % s1=nexttile(1);
                % plot(pulsePL_t,pulse200,'k','LineWidth',1.5);hold on; ylim([0 1.1]);
                % % plot(tpulse_new,TEMPLATE(:,3),'g');
                % legend('Pulse','Box', 'off');
                % subtitle(s1,['p:' ,num2str(p),' Pulse']);
                % s2=nexttile(2);
                % plot(pulsePL_t,[pulse_1dv],'b','LineWidth',1.5);yline(0);
                % legend(s2,'1st der.','Box', 'off');
                % s3=nexttile(3);
                % plot(pulsePL_t,[pulse_2dv],'m','LineWidth',1.5);yline(0);
                % hold on;
                % plot(pulsePL_t([indw]),pulse_2dv([indw]),'rd');
                % legend(s3,{'2nd der.','',''},'Box', 'off');
                % s4=nexttile(4);
                % plot(pulsePL_t,[pulse_3dv],'g','LineWidth',1.5);yline(0);
                % xline([s4],pulsePL_t(dia_loc),'b--');
                % %xline([s4],pulsePL_t(inddic),'b--');
                % legend(s4,{'3rd der.','',''},'Box', 'off');
                % linkaxes([s1,s2,s3,s4],'x');
                % grid([s1,s2,s3,s4], 'on');
                % grid([s1,s2,s3,s4], 'minor');
                % xlabel(s4,'Pulse width (samples)');
                % set([s1,s2,s3,s4],'Fontsize',16);



                %             %             xline([s1],loc5(1)/SF);
                %             %             xline([s2],loc5(1)/SF);
                %             %             xline([s3],loc5(1)/SF);
                %             %             xline([s1],loc6/SF);
                %             %             xline([s2],loc6/SF);
                %             %             xline([s3],loc6/SF);
                %             plot(s1,tpulse_new([indons,indend]),pulse200([indons,indend]),'ko','MarkerSize',12);
                %             plot(s1,tpulse_new([indmsl]),pulse200([indmsl]),'ko');
                %             plot(s1,tpulse_new([indsys]),pulse200([indsys]),'m|','MarkerSize',26);
                %             plot(s1,tpulse_new([inddic]),pulse200([inddic]),'r|','MarkerSize',26);
                %
                % %             %             plot(s1,tpulse_new([inddic_alt]),pulse200([inddic_alt]),'rd','MarkerFaceColor','r','MarkerSize',8);
                % %
           












                        end

                        %Store main fiducial points:
                        [main_points] =[];
                        main_points=[indons,indmsl,indsys,inddic,inddias,indend];


                        % save pulse waveforms and fiducial points for plotting later:
                        if save_pulses
                            pulse_waveforms{p}=pulse200;
                            pointsPPG(p,:)=main_points;

                        end



                        %% Extract Features:

                        if parameters.extract_features.do_extract

                        % Calculate next set of waveform features using the waveform and fiducial points:
                                feats_varnames={'PW','CT','dT',...
                                            'sysAmp','dicAmp','diasAmp',...
                                            'dT_norm', 'CT_norm',...
                                            'AS','AS_norm',... % Ascending slope
                                            'RefInd','RefInd_norm',...      % Reflection Index,Formerly used as:'AIx2dia_zsc','AIx2dia_zsc200' in ppg manuscript
                                            'SI','SI_norm',...
                                            'AIx', 'AIx_norm'};

                                % 1) Pulse Width (or Pulse Duration)
                                PW=nan;
                                PW=round(ceil((indend-indons)*1/sampratio)/SFdevice,3);
                                feats_derived(p,find(strcmp(feats_varnames,'PW')))=PW;

                                % 2) Crest Time (or Rise Time)
                                CT=nan;
                                CT=round(ceil((indsys-indons)*1/sampratio)/SFdevice,3);
                                feats_derived(p,find(strcmp(feats_varnames,'CT')))=CT;

                                % 3) Delta T (time difference between sys and dias peaks)
                                dT=nan;
                                try
                                    dT=round(ceil((inddias- indsys)*1/sampratio)/SFdevice,3);
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'dT')))=dT;

                                % 4) Systolic Amplitude
                                sysAmp=nan;
                                sysAmp=pulse200(indsys)-pulse200(indons);
                                feats_derived(p,find(strcmp(feats_varnames,'sysAmp')))=sysAmp;

                                % 5) Dicrotic Amplitude
                                dicAmp=nan;
                                try
                                    dicAmp=pulse200(inddic)-pulse200(indons);
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'dicAmp')))=dicAmp;

                                % 6) Diastolic Amplitude
                                diasAmp=nan;
                                try
                                    diasAmp=pulse200(inddias)-pulse200(indons);
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'diasAmp')))=diasAmp;

                                % 7) dT_norm: deltaT in duration-normalized pulse
                                dT_norm=nan;
                                try
                                    dT_norm=(inddias- indsys);
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'dT_norm')))=dT_norm;

                                % 8) CT_norm: crest time in duration-normalized pulse
                                CT_norm=nan;
                                CT_norm=(indsys-indons);
                                feats_derived(p,find(strcmp(feats_varnames,'CT_norm')))=CT_norm;

                                % 9) AS: Ascending Slope
                                AS=nan;
                                AS=[pulse200(indsys)-pulse200(indons)]/CT;
                                feats_derived(p,find(strcmp(feats_varnames,'AS')))=AS;

                                % 10) AS_norm: Ascending Slope in duration normalized pulse
                                AS_norm=nan;
                                AS_norm= [pulse200(indsys)-pulse200(indons)]/[(indsys-indons)];
                                feats_derived(p,find(strcmp(feats_varnames,'AS_norm')))=AS_norm;

                                % 11) RefInd: Reflection Index, dias amp/ sys amp (formerly used as: AIx2dia_zsc):
                                RefInd=nan;
                                try
                                    RefInd=(pulsePL(ceil(inddias*1/sampratio))-pulsePL(ceil(indons*1/sampratio)))/...
                                        (pulsePL(ceil(indsys*1/sampratio))-pulsePL(ceil(indons*1/sampratio)));
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'RefInd')))=RefInd;

                                % 12) RefInd_norm: Reflection Index in duration normalized pulse (formerly used as:AIx2dia_zsc200)
                                RefInd_norm=nan;
                                try
                                    RefInd_norm=(pulse200(inddias)-pulse200(indons))/(pulse200(indsys)-pulse200(indons));
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'RefInd_norm')))=RefInd_norm;

                                % 13) SI: Stiffness Index (SI)
                                SI=nan;
                                try
                                    SI=(pt_height/100)/(dT);
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'SI')))=SI;

                                % 14) SI_norm: Stiffness Index (SI)on duration normalized pulse
                                SI_norm=nan;
                                try
                                    SI_norm=(pt_height/100)/(dT_norm);
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'SI_norm')))=SI_norm;

                                % 15) AIx: Augmentation Index, as ratio between dicrotic
                                % notch amplitude and systolic peak amplitude on original
                                % pulse
                                AIx=nan;
                                try
                                    AIx=(pulsePL(ceil(inddic*1/sampratio))-pulsePL(ceil(indons*1/sampratio)))/...
                                        (pulsePL(ceil(indsys*1/sampratio))-pulsePL(ceil(indons*1/sampratio)));
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'AIx')))=AIx;

                                % 16) AIx_norm: Augmentation Index, as ratio between dicrotic
                                % notch amplitude and systolic peak amplitude on duration
                                % normalized pulse
                                AIx_norm=nan;
                                try
                                    AIx_norm=(pulse200(inddic)-pulse200(indons))/(pulse200(indsys)-pulse200(indons));
                                end
                                feats_derived(p,find(strcmp(feats_varnames,'AIx_norm')))=AIx_norm;

                        end
                        %              end
                    end

            end

        end

    end


    % if ploton
    %     title(gca,['PPG Pulses and Fiducial Points (n=',num2str(numel(~isnan(feats_derived(:,1)))),' pulses)']);
    %     set(gca,'Fontsize',14);
    %
    % end
    %     catch ME
    %
    %         error{1,p}=[p];
    %         error{2,p}=[ME.message];
    %
    % %         display('error:')
    % %         display(ME.message)
    %
    %
    %      end


    % display(['.....Feature extraction for the current segment is COMPLETE']);
    % display(['.......................................................... ']);





end

end







