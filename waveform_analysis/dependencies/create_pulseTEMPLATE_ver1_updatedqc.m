function [TEMPLATE,gaussfit, ncc_meanvsfitted] = create_pulseTEMPLATE_ver1_updatedqc(ppg_seg,onsets,parameters)

% Function to create template within selected window using the updated QC
% appoach which I developed during sg70ppg analysis (in Aug 2024) . Main
% difference from previous versions (which wre used in OV feature
% extractions) is inclusion of dtw -based similarity scores.


% INPUTS:
%             ppg_seg: ppg segment signal,
%             SF: sampling frequency of ppg signal
%             onsets: annotations of pulse onsets
%             peaks: anotations of pulse peaks
%             ploton: 1 for plot on , 0 for plot off

%``````````````````````````````````````````````````````````````
% Last updated by Gizem Y.D. @ 2024.10.10


%%%%--------------------------------------------------------------------------------------------------------------
if nargin < 3
    error('Wrong number of input arguments');
end

warning('off','all');

% disp(['..... Generating a TEMPLATE ....'])



%% Define and set parameters:
pulse_width_max=parameters.extract_features.pulse_width_max;
SFdevice=parameters.extract_features.SFdevice;
NCC_threshold=parameters.quality_threshold; %arbitrary, can be increased or decreased!
upsample=parameters.extract_features.upsample; % upsample to SFupsample Hz if 1. Do nothing if 0.
SFupsample=parameters.extract_features.SFupsample;
normalized_width=parameters.extract_features.normalized_width;


template_pulse=[]; 
clean_pulses=[];
clean_onsets=[];
clean_NCC=[];
ratio_of_dataloss=[];

ncc_meanvsfitted=[];
ncc_pulsevsmean=[];
ncc_pulsevsfitted=[];

gaussfit=[];
TEMPLATE=nan(normalized_width,3);
onsets_refined=[];




SF=SFdevice;
ppg_t=(0:1/SF:(length(ppg_seg)-1)/SF)';

if ~isempty(onsets) & ~isempty(ppg_seg) % Continue only if these are non-empty!!


        % Find individual pulse waveforms:
        reshape_factor=pulse_width_max ;  % to equalize length of pulses to calculate median. unit is sec
        
        cac = cell( 1, length(onsets)-1 );
        
        PDs=[];
        PDs=diff(onsets)/SF; %pulse durations
        
        medianPD=median(PDs(PDs < pulse_width_max));
        PD_upper=prctile(PDs(PDs < pulse_width_max),[95]);
        
        for jj = 1 :numel(onsets)-1;
          
            clear pulsePL ;
            pulsewave_temp=zeros(1,floor(SF*reshape_factor)); %zeros array to equalize pulse widths - to calculate mean/median

             % Does curren pulse have the offset point? No means, next pulse was discarded during SQI thresholding.
            if (onsets(jj+1)-onsets(jj))/SF <= PD_upper 
                pulsePL=ppg_seg(onsets(jj):onsets(jj+1),:)';

                if numel(pulsePL)<=numel(pulsewave_temp)
                pulsewave_temp(1:numel(pulsePL))=pulsePL;
               
                else
                pulsewave_temp=nan(1,floor(SF*reshape_factor)); %nan array
                PDs(jj)=nan;
                end

            else

                pulsewave_temp=nan(1,floor(SF*reshape_factor)); %nan array
                PDs(jj)=nan;
            end

            cac{jj} = pulsewave_temp;
            % cac_amp1{jj}=(pulsewave_temp-min(pulsewave_temp))/...
            %     (max(pulsewave_temp)-min(pulsewave_temp));

            cac_amp1{jj}=normalize(pulsewave_temp,"range",[0 100]);

            % figure;plot( cac_amp1{jj});title('pulseno',num2str(jj))
        end

       %% Option-1: Template PULSE (Calculate median of the pulses)
       [mean_pulse,allpulses,temp_pulse]=deal([]);
       
%         allpulses=reshape(cell2mat(cac(:)),size(cac,2),SF*reshape_factor);
        allpulses=reshape(cell2mat(cac_amp1(:)),size(cac_amp1,2),floor(SF*reshape_factor));
%         for i=1:size(allpulses,1);plot(allpulses(i,:));pause;end
       
       sk=[];
       for i=1:size(allpulses,1)
           currpulse=[];
           currpulse=allpulses(i,:);
           sk(i) = skewness(currpulse);
%              plot([1:1:75],currpulse,'k'); 
%              title(['Pulse no= ',num2str(i),' sk= ',num2str(sk(i))]); pause;
%              
       end

[sk_refined,allpulses_refined,onsets_refined, mean_pulse,temp_pulse]=deal([]);
skewness_threshold=-0.1;
sk_refined=sk(sk>=skewness_threshold );
allpulses_refined=allpulses(find(sk>= skewness_threshold),:);

onsets_refined=onsets(1:numel(onsets)-1);
onsets_refined=onsets(find(sk>= skewness_threshold ));

mean_pulse=mean(allpulses_refined,1,'omitnan'); % 
median_pulse=median(allpulses_refined,1,'omitnan'); % 

temp_pulse=[];

if ~isempty(onsets_refined)

% Amplitude-normalize temp_pulse;
temp_pulse=reshape((mean_pulse-min(mean_pulse))/(max(mean_pulse)-min(mean_pulse)),[],1); % min max normalized signal
temp_pulse=normalize(temp_pulse,"range",[0 100]); % normalize scale to max 100 (ro be compatible with Oura)
pulsePL_t=(0:1/SFdevice:(size(allpulses_refined,2)-1)/SFdevice)';


%% Should I check for presence of 2 peaks by fitting gaussians?
pulsePL=[];
pulsePL_t=[];
temp2up=[];
temp2up=temp_pulse;
temp2up(medianPD*SF+1:end)=[];

pulsePL_t=(0:1/SFdevice:(numel(temp2up)-1)/SFdevice)';

% pulsePL=temp2up;
if upsample==1

    [pulsePL_t_up,template_up]=deal([]);
    pulsePL_t_up = pulsePL_t(1):1/SFupsample:pulsePL_t(end);
    template_up = interp1(pulsePL_t,temp2up, pulsePL_t_up,'spline')'; % cubic spline interpolation

    pulsePL=template_up;
    pulsePL_t=pulsePL_t_up';
    SF=SFupsample;

else
   SF=SFdevice;
   pulsePL=temp_pulse;
   pulsePL(PD_upper*SF+1:end)=[];
   pulsePL_t=(0:1/SF:(numel(pulsePL)-1)/SF);

end

pulsePL_t=reshape(pulsePL_t,[],1);
%% Use peakfit.m to fit 2 gaussians:
% https://terpconnect.umd.edu/~toh/spectrum/functions.html
% Copyright (c) 2019, Thomas C. O'Haver

% Define result arrays:
[FitResults, GOF, baseline, coeff, residual, xi, yi, BootstrapErrors, gaussfit]=deal([]);
[template_fitted,template_original,template_t]=deal([]);
[gauss1,gauss2,pulse_gauss,pulse_refit_t]=deal([]);


% Define the signal for peak detection:
signal=[];

if ceil(SF*PD_upper) <= numel(pulsePL) %if template pulse has a long zeros tail;shorten it to PD-upper
signal=[pulsePL_t(1:ceil(SF*PD_upper)),...
    pulsePL(1:ceil(SF*PD_upper))];

else %use as it is

signal=[pulsePL_t,pulsePL];

end

% % Define approximate center and width parameters for gaussians (inital guess)
% center1=0.2*SF; 
% center2=0.2*SF+0.4*SF;
% window1= 0; %(max(pulsePL_t)-min(pulsePL_t))/3;
% window2= 0 ; %(max(pulsePL_t)-min(pulsePL_t))/3;

%Fit 2 gaussians:
[FitResults,GOF,baseline,coeff,residual,xi,yi,BootstrapErrors]=peakfit(signal,0,0,2,1,0, 0, 0, 0, 0,0);
% figure;hold on;
% [FitResults,GOF,baseline,coeff,residual,xi,yi,BootstrapErrors]=peakfit(signal,0,0,3,1,0, 0, 0, 0, 0,1);

% Using peak location of 1st gaussian, estimate sys peak for the
% mean-pulse-template
[sys_loc,sys_amp]=deal([]);
[~,sys_loc]=min(abs(pulsePL_t-FitResults(1,2))); 
sys_amp=pulsePL(sys_loc);

% Using peak location of 2nd gaussian, estimate dias peak for the
% mean-pulse-template
[dia_loc,dias_amp]=deal([]);
[~,dia_loc]=min(abs(pulsePL_t-FitResults(2,2))); 
dias_amp=pulsePL(dia_loc);

 % figure;plot(pulsePL_t,pulsePL,'k');hold on;xline(pulsePL_t(sys_loc),'r');xline(pulsePL_t(dia_loc),'b')
%% Continue creating template if conditions below are met:

if ~isempty(FitResults)  &&... % we must have non-empty FitResults
        sys_amp >= dias_amp % systolic amplitude must be larger than diastolic amplitude for MEAN pulse from this segment

gauss1=yi(1,:);
gauss2=yi(2,:);
pulse_gauss=gauss1+gauss2;
pulse_gauss_t=xi;

gaussfit.FitResults=FitResults; %peak number, peak position, height, width, and area 
gaussfit.GOF=GOF;
gaussfit.coeff=coeff;
gaussfit.xi=xi;
gaussfit.yi=yi;

%% Duration normalization
%  Interpolte upsampled-pulsePL and pulse_refit both into equal sizes
% (width normalization!)

% Width-normalize pulsePL
[pulseZrawdur,tpulse,tpulse_new]= deal([]);
[F,pulsePL_200,sampratio]=deal([]);
pulseZrawdur=pulsePL;
tpulse = pulsePL_t; %original time array of the pulse
tpulse_new = linspace(0,max(tpulse),normalized_width); % to normalize to 200 equally sampled intervals
F = griddedInterpolant(tpulse,pulseZrawdur,'spline','spline'); %interpolation and extrapolation grid
pulsePL_200 = [F(tpulse_new)]';
pulsePL_200=normalize(pulsePL_200,"range",[0 100]);

% Width-normalize gaussian-fitted pulse (pulse_gauss)
[pulseZrawdur,tpulse,tpulse_new]= deal([]);
[F,pulse_gauss_200,sampratio]=deal([]);
pulseZrawdur=pulse_gauss;
tpulse = pulse_gauss_t; %original time arrat of the pulse
tpulse_new = linspace(0,max(tpulse),normalized_width); % to normalize to 2000 equally sampled intervals
F = griddedInterpolant(tpulse,pulseZrawdur,'spline','spline'); %interpolation and extrapolation grid
pulse_gauss_200 = [F(tpulse_new)]';
pulse_gauss_200=normalize(pulse_gauss_200,"range",[0 100]);


% % Do you want to plot the results?
% fig=figure('Position',[600 400 1000 600]);
% % plot(pulsePL_t,pulsePL,'k');hold on;
% plot(tpulse_new,pulsePL_200,'k','LineWidth',1.5);hold on;
% ylim([0 102]);
% plot(tpulse_new,pulse_gauss_200,'r','LineWidth',1.5);
% plot(pulse_gauss_t,gauss1,'r-.','LineWidth',1.5);hold on;
% plot(pulse_gauss_t,gauss2,'r--','LineWidth',1.5);hold on;
% xline(pulse_gauss_t(find(gauss1==max(gauss1))),'b--','LineWidth',0.7);
% xline(pulse_gauss_t(find(gauss2==max(gauss2))),'b--','LineWidth',0.7);
% legend({'Average pulse','Sum of gaussians','Gaussian 1','Gaussian 2','',''});
% % title(['RMS fitting error: ',num2str(GOF(1)),'... R-squared: ',num2str(GOF(2)), '...NCC with meanpulse: ',num2str(max_NCC_refit)]);
% fontsize(gcf,18,"points")

% 
% for i=1:size(allpulses_refined,1)
% 
%     currpulse=[];
%     currpulse=allpulses_refined(i,:);
%     curr_t=(0:1/SFdevice:(size(allpulses_refined,2)-1)/SFdevice)';
%   
%     close(gcf);
%     plot(curr_t,currpulse,'k');hold on;
%     plot(curr_t,temp_pulse,'b');ylim([0 1.1]); hold on;
%     plot(tpulse_new,pulse_gauss_200,'r');
% 
%     legend({'current pulse','mean template','gauss template'});
% pause;
% 
% end


% Use pulse_gauss_2000 as a template to check for quality!
template_original=pulsePL_200;
template_fitted=pulse_gauss_200;
template_t=tpulse_new;

 
cut1=1;
cut2=120;

% % Calculate NCC between mean-pulse-template and gauss-based-template
% [NCC_refit,lags]=normxcorr2_general(template_original(cut1:cut2),template_fitted(cut1:cut2));
% ncc_meanvsfitted=round(max(NCC_refit),3);

[NCC_refit, lags] = xcorr(template_original(cut1:cut2),template_fitted(cut1:cut2), 'coeff');
ncc_meanvsfitted=round(max(NCC_refit),3);


 % Compute the DTW distance and alignment

 [dist,ix,iy,dist_raw,path_length]=deal([]);
 [dist,ix,iy] = dtw(template_original(cut1:cut2),template_fitted(cut1:cut2),'euclidean');
 dist_raw=dist;
 path_length=length(ix);
 dtw_norm= round(dist_raw/ path_length,2);


% figure;
% plot(template_t,template_original,'k');hold on;
% plot(template_t(cut1:cut2),template_original(cut1:cut2),'k.','MarkerSize',8);hold on;
% plot(template_t,template_fitted,'m');hold on;
% plot(template_t(cut1:cut2),template_fitted(cut1:cut2),'m.','MarkerSize',8);
% legend("mean","","fitted","");
% title(['NCC: ',num2str(ncc_meanvsfitted)," dtw: ",num2str(dtw_norm)]);
% 

if round(ncc_meanvsfitted,2)>=0.98 %&& GOF(2)>= 0.98
    TEMPLATE=[template_t',template_original,template_fitted];
else
    TEMPLATE=nan(normalized_width,3);

end
%% Running through each pulse and check against template?

% 
% 
% for i=1:size(allpulses_refined,1)
% 
%     currpulse=[];
%     currpulse=allpulses_refined(i,:);
%     currpulse=currpulse(1:(medianPD*SFdevice));
% 
% % Width-normalize current pulse: ---
% [pulseZrawdur,tpulse,tpulse_new,F,currpulse_200,sampratio]= deal([]);
% pulseZrawdur=currpulse;
% tpulse = (0:1/SFdevice:(numel(currpulse)-1)/SFdevice)'; %original time array of the pulse
% tpulse_new = linspace(0,max(tpulse),200); % to normalize to equally sampled intervals
% F = griddedInterpolant(tpulse,pulseZrawdur,'spline','spline'); %interpolation and extrapolation grid
% currpulse_200 = [F(tpulse_new)]';
% currpulse_200=normalize(currpulse_200,"range",[0 100]);
% 
%  % Normalized cross correlation as a measure of SQI --
%     [C_ncc,lags_ncc,Cpeaks]=deal([]);
%     [C_ncc,lags_ncc] = normxcorr2_general(template_fitted(cut1:cut2),currpulse_200(cut1:cut2),100);
%        % figure;plot(C_ncc)
%     NCC_withtemplate(i)=round(max(C_ncc),2);
% 
% 
%     [C_ncc2,lags_ncc2,Cpeaks2]=deal([]);
%     [C_ncc2,lags_ncc2] = normxcorr2_general(template_original(cut1:cut2),currpulse_200(cut1:cut2),100);
%        % figure;plot(C_ncc)
%     NCC_withtemplate2(i)=round(max(C_ncc2),2);
% 
% 
%     figure('Position',[400 400 1500 500]);    subplot(1,2,1)
%     plot(template_t,currpulse_200,'k');hold on;
%     plot(template_t(cut1:cut2),currpulse_200(cut1:cut2),'k.','MarkerSize',8);hold on;
%     plot(template_t,template_fitted,'m');hold on;
%     plot(template_t(cut1:cut2),template_fitted(cut1:cut2),'m.','MarkerSize',8); ylim([0 105]);
%     legend({'pulse',"",'template-gauss',""});
%     title(['template and pulse-',num2str(i), ' -->  NCC is ',num2str(NCC_withtemplate(i))]);
%     subplot(1,2,2)
%     plot(template_t,currpulse_200,'k');hold on;
%     plot(template_t(cut1:cut2),currpulse_200(cut1:cut2),'k.','MarkerSize',8);hold on;
%     plot(template_t,template_original,'b');hold on;
%     plot(template_t(cut1:cut2),template_original(cut1:cut2),'b.','MarkerSize',8); ylim([0 105]);
%     legend({'pulse',"",'template-mean',""});
%     title(['template and pulse-',num2str(i), ' -->  NCC is ',num2str(NCC_withtemplate2(i))]);
% % 
% end

% 


end

end


end


end

















