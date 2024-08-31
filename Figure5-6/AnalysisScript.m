miceopt = {'Galilei','Haydn','Jemison','Faraday','Lully'};%%options for mice
MouseID = [1,2,3,5,7]; % ID in paper
GenOpt = {'Drd2','Drd1','Drd2','Drd2','Drd1'};%
LocalFolder = '\Path\To\Data\'
ServerDir = '\Path\To\Server'
traininglogdir = '\Path\To\LogFiles\'
Stim2Check = 'SpatialWM' %Name of the stimulus as written in the LOG-file

AnalysisParameters.motioncorrection=0;
AnalysisParameters.NotMotivatedDate = {'20200428','20200921','20200428','20200921','20201028'};
AnalysisParameters.WindowSize = 25; %Nr. trials around window
AnalysisParameters.baselinemethod = 1; %1: trial by trial baseline, 2: filtered baseline (F0), 3: weighted average baseline after detrending (over all trials)
AnalysisParameters.smoothfact = 2; %Smoothing with neighbours (Gaussian); fill in standard deviation
AnalysisParameters.ScaleFact = 0.5; %between 0-1 is scale down
AnalysisParameters.AREASGrouped = {{'V1'},{'Vl','Vli','Val','Vpor'},{'Vrl','Va'},{'Vam','Vpm'},{'RSP'},{'M1'},{'M2'}};%note: l includes li and l{'V1','HVLat_Vl_Val','HVAnterior_Vrl_Va','HVMedial_Vam_Vpm','RSP','M1','M2'};%% %
AREAGROUPNAMES = {'Primary Visual Area','Lateral Visual Areas','Anterior Visual Areas',...
    'Medial Visual Areas','Retrosplenial Area','Primary Motor Cortex','Secondary Motor Cortex'};
AnalysisParameters.AREAS = {'V1','Vl','Val','Vrl','Va','Vam','Vpm','RSP','M1','M2'};%note: l includes li and l{'V1','HVLat_Vl_Val','HVAnterior_Vrl_Va','HVMedial_Vam_Vpm','RSP','M1','M2'};%% %
AnalysisParameters.ReactionOpt = {'Hit','Error','Miss'} %Which response types to include
AnalysisParameters.OptoOpt = {'Opto Off','Opto Baseline','Opto Visual','Opto Delay','Opto Response','Opto Post-trial'};
AnalysisParameters.SideOpt = {'left','right'}
AnalysisParameters.trialtypes = {'1500'}
AnalysisParameters.TimeWindowNames = {'Baseline','BaselineOpt','StimOn','DelayEarly','DelayOpto','DelayLate','Response','Post-Trial','ITI'};
AnalysisParameters.TW = {[-850 -500],[-500 0],[0 500],[650 1000],[1000 1500],[1500 1850],[1950 2500],[3500 4000],[4500 7500]};
AnalysisParameters.templatemouse = 'Jemison'


HemOpt = {'RightHem','LeftHem'};
tasknames = {'StimulusSide','ResponseSide'};
variablenames = {'Trial_ID','Session_ID','Reaction','Side','Delay','Opto','GavePassive','DifMotion','TimePoint','Pixel_Idx','dFF'};
OptoOpt = AnalysisParameters.OptoOpt;
OptoOptGrouped = {'Off','Behavioural Effect','Opto Late'};
GenoTypeOpt = unique(GenOpt);
optogroupedid={1,[2:4],[5:6]};
perfthresh=0.55;
bootstrapMVPA = 0;
searchwin=15;
SideOpt = AnalysisParameters.SideOpt; %Always same order!
ReactionOpt = AnalysisParameters.ReactionOpt; %Always same order!
DelayOpt = AnalysisParameters.trialtypes;
AREAS = AnalysisParameters.AREASGrouped;
% motioncorrection = AnalysisParameters.motioncorrection;
workwithallcollapsed = 0; %Work with 'collapsed' data (projected onto one hemisphere) or both hemispheres?
cols = cat(2,[108; 99; 91], [140; 20; 0], [110;81;76], [140;20;160], [189;105;41], [229;176;46], [154;192;209], [73;121;159],[111;60;151], [46;141;67], [209;211;212], [37;34;25])';
cols = cols./255;
NewDim = 40;
nback = 50;
nboot = 500;
if length(AnalysisParameters.OptoOpt)==1
    traininglogdir = '\\vs02\VandC\WorkingMemory_EB\WMTask_Widefield\TrainingLogs\'
    TW = {[-300 0],[0 500],[700 1200],[1200 1900],[1950 2500],[3000 4000]};
    OptoTW = {[nan nan]};
    TimeWindowNames = {'Baseline','StimOn','DelayEarly','DelayLate','Response', 'Post-Reward'};
elseif length(AnalysisParameters.OptoOpt)>1
    traininglogdir = '\\vs02\VandC\WorkingMemory_EB\WMTask_Opto\LogFiles\';
    TW = {[-850 -500],[-500 0],[0 500],[700 1200],[1200 1900],[1950 2500],[3000 4000],[4000 5500],[5500 7000]};
    TimeWindowNames = {'Baseline','BaselineOpt','StimOn','DelayEarly','DelayLate','Response','Post-Reward','Motor Down','ITI'};
    OptoTW = {[nan nan],[-500,0],[0 500],[1000 1500],[2000 2500],[4000 4500]};
    % 0 = off
    % 1 = baseline (-500 to -0ms)
    % 2 = visual (0 to 500ms)
    % 3 = delay (1000 to 1500ms)
    % 4 = response (2000 to 2500ms)
    % 5 = post-trial (3500 to 4000)
    TWRewardBiasAna = {[500 2000],[3000 4500]};
    TWRewardBiasNames = {'Delay','ITI'};
end

% Get timeline
tmp = matfile(fullfile(LocalFolder,'\Faraday\TimeCourses.mat'));
timeline = tmp.timeline;

%% Linear Model
% Build Design Matrix with variables Stimulus Left, Stimulus Right, Stimulus Offset, Left
% Choice, Right Choice, Reward. These are shifted in time as well.
timevec=downsample(timeline,2); % downsample
SR = unique(diff(timevec)); if length(SR)>1; keyboard;end
newpixsz = 20; %Downsample to this
ReloadDatPerPix = 0;
Redo = 0


%% %Define kernel indices
DummyIdx = 1;
VisKernelIdx = 2:length(timevec)*2+1;

RewardKernelIdx = length(timevec)*2+2:length(timevec)*4+1;
ChoiceKernelIdx = length(timevec)*2+2:length(timevec)*6+1;
WholeBrainIdx = length(timevec)*6+2;

LeftLicksKernelIdx = length(timevec)*6+3;
RightLicksKernelIdx = length(timevec)*6+4;
EyeMovKernelIdx = length(timevec)*6+5;
BodyMovKernelIdx = length(timevec)*6+6;
AllKernels = {DummyIdx,VisKernelIdx,RewardKernelIdx,ChoiceKernelIdx,WholeBrainIdx,...
    LeftLicksKernelIdx,RightLicksKernelIdx,EyeMovKernelIdx,BodyMovKernelIdx};

AllKernelNames = {'Dummy','Stimulus','Reward','Choice','WholeBrain','LeftLicks','RightLicks','EyeMovement','BodyMovement'};
vecvector=2:5; % Which to shuffle analysis?

for midx=1:length(miceopt)
    for optid=1
        if exist(fullfile(LocalFolder,miceopt{midx},['LinearModels_Opto' num2str(optid) '.mat'])) && ~Redo
            disp('existing data, continue...')
            continue
        end
        % Load in Processed Data
        % Copy to local, if not yet done
        if ~exist(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'))
            if ~isdir(fullfile(LocalFolder,miceopt{midx}))
                mkdir(fullfile(LocalFolder,miceopt{midx}))
            end
            disp('Copying data to local')
            copyfile(fullfile(ServerDir,miceopt{midx},'ProcessedData.mat'),fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'))
            copyfile(fullfile(ServerDir,miceopt{midx},[miceopt{midx} 'BrainModel.mat']),fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']))
        end
        TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
        model = load(fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));
        model = model.BrainModel;
        % Get dFF
        disp('Loading in data, which will take a bit of time..')
        tmptrials = TMPDat.Trials2Keep;

        LM = matfile(fullfile(LocalFolder,miceopt{midx},['LinearModels_Opto' num2str(optid) '.mat']));
        % we can load dat per pix
        if exist('LM') && ~ReloadDatPerPix
            DatPerPix = LM.DatPerPix;%nan(newpixsz,newpixsz,length(timevec),0);
            ReExtractDatPix = 0;
            DatPerPix = reshape(DatPerPix,newpixsz,newpixsz,length(timevec),[]);
        else
            AllDat = TMPDat.AllDat; %get dFF %Pre/Post Training, RewardOptions (1/0) , Stimulus Options (HRP/Catch/LRP)
            DatPerPix = nan(newpixsz,newpixsz,length(timevec),0);
            ReExtractDatPix = 1;
        end
        SideVec = nan(1,0);
        RewardVec = nan(1,0);
        ChoiceVec = nan(1,0);
        LeftActionVec = nan(length(timevec),0);
        RightActionVec = LeftActionVec;
        DifMotionVec = LeftActionVec;
        EyeDataVec = LeftActionVec;
        TrialVec = nan(1,0);
        SesVec = nan(1,0);
        for ridx=1:3
            for sidx=1:length(AnalysisParameters.SideOpt)
                disp(['Downsampling data for ridx = ' num2str(ridx) ', sidx=' num2str(sidx) '...'])

                if ReExtractDatPix

                    tmp = nan(newpixsz,newpixsz,size(AllDat{optid,ridx,sidx},3),size(AllDat{optid,ridx,sidx},4));
                    % Downsample in space
                    for trid = 1:size(tmp,4)
                        tmp(:,:,:,trid) = Nan_imresize(AllDat{optid,ridx,sidx}(:,:,:,trid),newpixsz./size(AllDat{optid,ridx,sidx},1));
                    end

                    if isempty(tmp)
                        continue
                    end
                    % Get trials      - downsample in time: add low-pass 9Hz
                    % filter? (anti-aliasing) (resample includes this already!)
                    tmp = reshape(permute(tmp,[3,1,2,4]),length(timeline),[]);
                    tmp = single(resample(double(tmp),length(timevec),length(timeline)));
                    tmp = permute(reshape(tmp,length(timevec),newpixsz,newpixsz,[]),[2,3,1,4]);

                    DatPerPix = cat(4,DatPerPix,tmp);
                end


                % Extract action
                tmpdatFiles = dir(fullfile(LocalFolder,miceopt{midx},'**',['Log.mat']));
                tmpses = unique(tmptrials(2,:));
                if length(tmpdatFiles)<length(tmpses)
                    keyboard
                end
                for sesid = 1:length(tmpses)
                    % Get correct log file
                    Log = matfile(fullfile(tmpdatFiles(sesid).folder,tmpdatFiles(sesid).name));
                    TimeLineHere = Log.TimeLineHere;
                    Log = Log.Log;
                    takethesetrials = tmptrials(1,tmptrials(2,:)==tmpses(sesid));
                    side = Log.Side(takethesetrials);
                    reaction = Log.Reaction(takethesetrials);
                    opto = Log.Opto(takethesetrials);
                    gavepassive = Log.Gavepassive(takethesetrials);
                    delay = Log.currentdelay(takethesetrials);
                    %Extract trials included in dataset
                    takethesetrials = takethesetrials(strcmp(side,SideOpt{sidx})&strcmp(reaction,ReactionOpt{ridx})&opto==optid-1&gavepassive==0&delay==1500);
                    if isempty(takethesetrials)
                        continue
                    end
                    % Extract Motion and downsample in time
                    if isfield(Log,'DifMotion')
                        DifMotionVec = cat(2,DifMotionVec,resample(Log.DifMotion(TimeLineHere>=timevec(1)&TimeLineHere<=timevec(end),takethesetrials),length(timevec),length(TimeLineHere)));
                    else
                        DifMotionVec = cat(2,DifMotionVec,nan(length(timevec),length(takethesetrials)));
                    end
                    % Extract eye movements
                    if isfield(Log,'EyeMovement')
                        tmpeye = nanmean(Log.EyeMovement(TimeLineHere>=timevec(1)&TimeLineHere<=timevec(end),takethesetrials,:),3); %take average eyemotion across the 4 coordinates
                        EyeDataVec = cat(2,EyeDataVec,resample(tmpeye,length(timevec),length(TimeLineHere)));
                    else
                        EyeDataVec = cat(2,EyeDataVec,nan(length(timevec),length(takethesetrials)));
                    end

                    % Extract licks
                    tmplicksRight = arrayfun(@(X) Log.RTrightVec{X},takethesetrials,'UniformOutput',0);
                    tmplicksLeft = arrayfun(@(X) Log.RTleftVec{X},takethesetrials,'UniformOutput',0);

                    tmpactionL = zeros(length(timevec),size(tmplicksLeft,2)); % Left first
                    tmpactionR = zeros(length(timevec),size(tmplicksRight,2)); % Now right

                    % Set to 'event' at correct time
                    for trid=1:size(tmplicksLeft,2)
                        tmp=arrayfun(@(Y) abs(timevec-Y)<=SR/2,tmplicksLeft{trid},'UniformOutput',0);
                        if ~isempty(tmp)
                            tmpactionL(:,trid)=nansum(cat(1,tmp{:}),1);
                        end
                        tmp=arrayfun(@(Y) abs(timevec-Y)<=SR/2,tmplicksRight{trid},'UniformOutput',0);
                        if ~isempty(tmp)
                            tmpactionR(:,trid)=nansum(cat(1,tmp{:}),1);
                        end
                    end
                    LeftActionVec = cat(2,LeftActionVec,tmpactionL);
                    RightActionVec = cat(2,RightActionVec,tmpactionR);
                    ntrialstmp = length(takethesetrials);
                    % Assign condition
                    SideVec = cat(2,SideVec,repmat(sidx,1,ntrialstmp));
                    RewardVec = cat(2,RewardVec,repmat(ridx,1,ntrialstmp));
                    if ridx==1
                        ChoiceVec = cat(2,ChoiceVec,repmat(sidx,1,ntrialstmp));
                    elseif ridx==2
                        if sidx==1
                            ChoiceVec = cat(2,ChoiceVec,repmat(2,1,ntrialstmp));
                        elseif sidx==2
                            ChoiceVec = cat(2,ChoiceVec,repmat(1,1,ntrialstmp));
                        end
                    elseif ridx==3 % No choice made
                        ChoiceVec = cat(2,ChoiceVec,repmat(0,1,ntrialstmp));
                    end
                    TrialVec = cat(2,TrialVec,takethesetrials);
                    SesVec = cat(2,SesVec,repmat(tmpses(sesid),1,length(takethesetrials)));

                end
            end
        end
        clear AllDat
        %

        %Create Design Matrix for events and their shifts: HRP, LRP,
        %Reward
        disp('Creating Design Matrix')
        ntrials = size(DatPerPix,4);
        DM = zeros(4,length(timevec),ntrials);
        %Extract left stimulus trials
        tmpDM = zeros(length(timevec),ntrials);
        tmpDM(1,SideVec==1)=1;
        DM(1,:,:) = tmpDM;

        %Extract right stimulus trials
        tmpDM = zeros(length(timevec),ntrials);
        tmpDM(1,SideVec==2)=1;
        DM(2,:,:) = tmpDM;

        %Extract Reward event trials
        tmpDM = zeros(length(timevec),ntrials);
        tmpDM(1,RewardVec==1)=1;
        DM(3,:,:) = tmpDM;

        %Extract Error event trials
        tmpDM = zeros(length(timevec),ntrials);
        tmpDM(1,RewardVec==2)=1;
        DM(4,:,:) = tmpDM;

        % Extract choice event
        tmpDM = zeros(length(timevec),ntrials);
        tmpDM(1,ChoiceVec==1)=1;
        DM(5,:,:) = tmpDM;

        % Extract choice event
        tmpDM = zeros(length(timevec),ntrials);
        tmpDM(1,ChoiceVec==2)=1;
        DM(6,:,:) = tmpDM;

        % Create shifted copies
        fullDM = nan(length(timevec),size(DM,1),size(DM,2),size(DM,3));
        for tp=1:length(timevec)
            fullDM(tp,:,:,:) = circshift(DM,tp-1,2);
        end
        DM = reshape(fullDM,length(timevec)*size(DM,1),size(DM,2),size(DM,3));
        clear fullDM

        % Add the left licks
        %Normalize Action Vec (0-1)
        LeftActionVec = (LeftActionVec-nanmin(LeftActionVec(:)))./(nanmax(LeftActionVec(:))-nanmin(LeftActionVec(:)));
        LeftActionVec(isnan(LeftActionVec))=0;
        DM = cat(1,DM,reshape(LeftActionVec,[1,size(LeftActionVec,1),size(LeftActionVec,2)]));

        % Add the right licks
        RightActionVec = (RightActionVec-nanmin(RightActionVec(:)))./(nanmax(RightActionVec(:))-nanmin(RightActionVec(:)));
        RightActionVec(isnan(RightActionVec))=0;
        DM = cat(1,DM,reshape(RightActionVec,[1,size(RightActionVec,1),size(RightActionVec,2)]));

        % Add eye movements
        % Remove outliers
        EyeDataVecZ = (EyeDataVec-nanmean(EyeDataVec(:)))./nanstd(EyeDataVec(:));
        EyeDataVec(abs(EyeDataVecZ)>5)=nan;
        EyeDataVec = (EyeDataVec-nanmin(EyeDataVec(:)))./(nanmax(EyeDataVec(:))-nanmin(EyeDataVec(:)));
        EyeDataVec = fillmissing(EyeDataVec,'linear',1,'EndValues','nearest');
        DM = cat(1,DM,reshape(EyeDataVec,[1,size(EyeDataVec,1),size(EyeDataVec,2)]));

        % Add body movements
        % Remove outliers
        DifMotionVecZ = (DifMotionVec-nanmean(DifMotionVec(:)))./nanstd(DifMotionVec(:));
        DifMotionVec(abs(DifMotionVecZ)>5)=nan;
        DifMotionVec = (DifMotionVec-nanmin(DifMotionVec(:)))./(nanmax(DifMotionVec(:))-nanmin(DifMotionVec(:)));
        DifMotionVec = fillmissing(DifMotionVec,'linear',1,'EndValues','nearest');
        DM = cat(1,DM,reshape(DifMotionVec,[1,size(DifMotionVec,1),size(DifMotionVec,2)]));

        % Reshape
        DM = reshape(DM,size(DM,1),size(DM,2)*size(DM,3));

        % Add Dummy variable
        DM = cat(1,ones(1,size(DM,2)),DM);

        %Reshape data
        TotalTrialVec = TrialVec + SesVec*1000; %add surreal number to do ordering
        [~,sortidx] = sort(TotalTrialVec);

        tmpfileinfo = dir(fullfile(LocalFolder,miceopt{midx},['LinearModels_Opto' num2str(optid) '.mat']));
        if datetime(tmpfileinfo.date)<datetime('27-04-2023','InputFormat','dd-MM-yyyy') %otherwise already sorted
            DatPerPix = DatPerPix(:,:,:,sortidx);
        end
        % Check for logical reasons
        figure;
        subplot(1,3,1)
        h=imagesc(nanmean(nanmean(DatPerPix(:,:,timevec>200&timevec<700,SideVec(sortidx)==1),3),4) - nanmean(nanmean(DatPerPix(:,:,timevec>200&timevec<700,SideVec(sortidx)==2),3),4),[-0.03 0.03]);
        colormap redblue
        title('Stimulus')
        axis square
        subplot(1,3,2)
        h=imagesc(nanmean(nanmean(DatPerPix(:,:,timevec>2000&timevec<3000,ChoiceVec(sortidx)==1),3),4) - nanmean(nanmean(DatPerPix(:,:,timevec>2000&timevec<3000,ChoiceVec(sortidx)==2),3),4),[-0.03 0.03]);
        colormap redblue
        title('Choice')
        axis square

        subplot(1,3,3)
        h=imagesc(nanmean(nanmean(DatPerPix(:,:,timevec>2000&timevec<4000,RewardVec(sortidx)==1),3),4) - nanmean(nanmean(DatPerPix(:,:,timevec>2000&timevec<4000,RewardVec(sortidx)==2),3),4),[-0.03 0.03]);
        colormap redblue
        title('Reward')
        axis square



        % Check if brain top view is aligned
        CheckBrainAlignment

        DM = reshape(DM,size(DM,1),length(timevec),ntrials);
        DM = DM(:,:,sortidx);
        DM = reshape(DM,size(DM,1),length(timevec)*ntrials);
        %Reshape data
        DatPerPix = reshape(DatPerPix,newpixsz*newpixsz,[]);

        % Add in whole brain
        WholeBrainResp = nanmean(DatPerPix,1);
        % Normalize between 0 and 1 like the other DMs
        WholeBrainResp = (WholeBrainResp-nanmin(WholeBrainResp))./(nanmax(WholeBrainResp)-nanmin(WholeBrainResp));
        DM(WholeBrainIdx+1:size(DM,1)+1,:) = DM(WholeBrainIdx:size(DM,1),:); % Move over the rest
        DM(WholeBrainIdx,:) = WholeBrainResp;

        figure('name','DesignMatrix')
        imagesc(DM)
        colormap(flipud(gray))
        xlabel('TimexTrial')
        ylabel('Predictor')
        FigureDefault
        drawnow
        %%
        RunLinearModels

        %Save stuff here
        disp('Saving...')
        try
            save(fullfile(LocalFolder,miceopt{midx},['LinearModels_Opto' num2str(optid) '.mat']),'Rsq_PerVar','RExcl_PerVar','TotalR','ChosenK','PredicteddFF','RExcl_PerVar_Shuffled','PredicteddFF_LeaveOut','DatPerPix','ntrials')
        catch ME
            disp(ME)
            keyboard
        end

        %%
        % Visualize Results
        figure('name',[miceopt{midx} ' ' OptoOpt{optid} ' TotalR'])
        lims = [quantile(TotalR(:),0.05) quantile(TotalR(:),0.95)];
        for cv = 1:size(TotalR,2)
            subplot(1,size(TotalR,2),cv)
            h=imagesc(reshape(TotalR(:,cv),newpixsz,newpixsz),lims);
            colorbar
            axis square
            colormap(flipud(gray))
            FigureDefault
            axis off
        end
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR_' OptoOpt{optid} '.fig']))
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR_' OptoOpt{optid} '.bmp']))

        % Visualize Results
        figure('name',[miceopt{midx} ' ' OptoOpt{optid} ' TotalR'])
        h=imagesc(reshape(nanmean(TotalR,2),newpixsz,newpixsz),lims);
        colorbar
        axis square
        colormap(flipud(gray))
        FigureDefault
        axis off
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR_Avg' OptoOpt{optid} '.fig']))
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR_Avg' OptoOpt{optid} '.bmp']))

        % Visualize Results
        DeltaVar = (repmat(TotalR,[1,1,size(RExcl_PerVar,2)])-permute(RExcl_PerVar,[1,3,2]))./repmat(TotalR,[1,1,size(RExcl_PerVar,2)]); % Calculate difference in explained variance (normalized by total variance explained0
        DeltaVar = permute(DeltaVar,[1,3,2]);


        for kernelid=1:length(AllKernels)
            try
                figure('name',[miceopt{midx} ' ' OptoOpt{optid} ' ' AllKernelNames{kernelid}])
                HrpRsq = squeeze(Rsq_PerVar(:,kernelid,:)); % Take the maximum RsqPerVar across the HRP kernel
                HRPDelta = squeeze(DeltaVar(:,kernelid,:)); % Take the maximum difference in explained variance
                lims1 = [quantile(HrpRsq(:),0.05) quantile(HrpRsq(:),0.95)];
                lims2 = [quantile(HRPDelta(:),0.05) quantile(HRPDelta(:),0.95)];
                for cv = 1:size(HrpRsq,2)
                    subplot(2,size(HrpRsq,2),cv)
                    h=imagesc(reshape(HrpRsq(:,cv),newpixsz,newpixsz),lims1);
                    colorbar
                    axis square

                    subplot(2,size(HRPDelta,2),cv+size(HRPDelta,2))
                    h=imagesc(reshape(HRPDelta(:,cv),newpixsz,newpixsz),lims2);
                    colorbar
                    axis square
                end
                colormap hot

            catch ME
                disp(ME)
            end

        end

        figure('name',[miceopt{midx} ' ' OptoOpt{optid} 'KernelRepresentations'])
        for kernelid=1:length(AllKernels)
            colormap(flipud(gray))
            HRPDelta = squeeze(nanmean(DeltaVar(:,kernelid,:),3)); % Take the maximum difference in explained variance
            lims2 = [quantile(HRPDelta(:),0.05) quantile(HRPDelta(:),0.95)];
            if all(isnan(lims2))
                continue
            end
            subplot(ceil(sqrt(length(AllKernels))),round(sqrt(length(AllKernels))),kernelid)
            h=imagesc(reshape(HRPDelta,newpixsz,newpixsz),lims2);
            colorbar
            axis square

            axis off
            title(AllKernelNames{kernelid})
        end


        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' OptoOpt{optid} '.fig']))
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' OptoOpt{optid} '.bmp']))


        %%
        cols = [0 0 1; 0 0.6 0.8; 1 0 0; 0.8 0.6 0];
        linstls = {'-','--'};
        figure('name',['Explained variance across time Full model']);
        for areaid=1:length(AREAS)

            for hemid=1:2
                % Find area pixels
                idx = find(ismember(model.Rnames,AREAS{areaid}));
                if isempty(idx)
                    idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),model.Rnames,'UniformOutput',0)));
                    areamask = false(400,400);
                    for ridx = 1:length(idx)
                        areamask(logical(round(model.Regions{idx(ridx)}))')=true;
                    end
                else
                    areamask = logical(any(cat(3,model.Regions{idx}),3))';
                end
                if hemid==2 %left hemisphere
                    areamask(1:round(model.Lambda(2)),:)=0;
                else
                    areamask(round(model.Lambda(2)):end,:)=0;
                end
                % Downsample areamask
                areamask = Nan_imresize(areamask,newpixsz./400);
                areamask(areamask>0)=1;
                pixidx = find(areamask);
                %Extract prediction for every timepoint and trial
                tmpPred = squeeze(nanmean(reshape(PredicteddFF(pixidx,:),length(pixidx),length(timevec),ntrials),1));
                tmpActual = squeeze(nanmean(reshape(DatPerPix(pixidx,:),length(pixidx),length(timevec),ntrials),1));

                % Calculate UniqueExplainedVariancePerTP
                TotalR = 1 - (nansum((tmpActual - tmpPred).^2,2)./nansum((tmpActual - nanmean(tmpActual,2)).^2,2));

                % Draw Left Stimulus Rewarded
                subplot(length(AREAS),2,(areaid-1)*2+hemid)
                hold on

                h(1)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==1&RewardVec==1),2)),'color',cols(1,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==1&RewardVec==1),2)),'color',cols(1,:),'LineStyle',linstls{2})

                % Right stimulus Rewarded
                h(2)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==2&RewardVec==1),2)),'color',cols(2,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==2&RewardVec==1),2)),'color',cols(2,:),'LineStyle',linstls{2})

                % Left not rewarded
                h(3)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==1&RewardVec==2),2)),'color',cols(3,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==1&RewardVec==2),2)),'color',cols(3,:),'LineStyle',linstls{2})

                % Right stimulus not Rewarded
                h(4)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==2&RewardVec==2),2)),'color',cols(4,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==2&RewardVec==2),2)),'color',cols(4,:),'LineStyle',linstls{2})


                %             title(AREAGROUPNAMES{areaid})
                if hemid==1
                    ylabel([AREAGROUPNAMES{areaid}])
                end
                if areaid==1
                    title(HemOpt{hemid})
                end
                FigureDefault
                xlim([-100 timevec(end)])
            end
        end
        legend([h(:)],{'Left Rewarded','Right Rewarded','Left Error','Right Error'})
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' OptoOpt{optid} '_ModelPerformance.fig']))
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' OptoOpt{optid} '_ModelPerformance.bmp']))


        % Extract prediction for different timings and areas

        MSE = nan(length(timevec),length(AllKernels),length(AREAS),2);
        UniqEV = MSE; % Difference between totalR and leaveoutR
        clear h
        for kernelid=1:length(AllKernels)

            figure('name',['Explained variance across time ' AllKernelNames{kernelid}]);
            for areaid=1:length(AREAS)

                for hemid=1:2
                    % Find area pixels
                    idx = find(ismember(model.Rnames,AREAS{areaid}));
                    if isempty(idx)
                        idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),model.Rnames,'UniformOutput',0)));
                        areamask = false(400,400);
                        for ridx = 1:length(idx)
                            areamask(logical(round(model.Regions{idx(ridx)}))')=true;
                        end
                    else
                        areamask = logical(any(cat(3,model.Regions{idx}),3))';
                    end
                    if hemid==2 %left hemisphere
                        areamask(1:round(model.Lambda(2)),:)=0;
                    else
                        areamask(round(model.Lambda(2)):end,:)=0;
                    end
                    % Downsample areamask
                    areamask = Nan_imresize(areamask,newpixsz./400);
                    areamask(areamask>0)=1;
                    pixidx = find(areamask);
                    %Extract prediction for every timepoint and trial
                    tmpPred = squeeze(nanmean(reshape(PredicteddFF_LeaveOut(pixidx,:,kernelid),length(pixidx),length(timevec),ntrials),1));
                    tmpActual = squeeze(nanmean(reshape(DatPerPix(pixidx,:),length(pixidx),length(timevec),ntrials),1));


                    % Calculate UniqueExplainedVariancePerTP
                    MSE(:,kernelid,areaid,hemid) = nanmean(squeeze(nanmean(tmpPred,1)-nanmean(tmpActual,1)).^2,2);
                    tmpPrdFF = squeeze(nanmean(reshape(PredicteddFF(pixidx,:),length(pixidx),length(timevec),ntrials),1));

                    TotalR = 1 - (nansum((tmpActual - tmpPrdFF).^2,2)./nansum((tmpActual - nanmean(tmpActual,2)).^2,2));
                    TotalRLeaveOut = 1 - (nansum((tmpActual - tmpPred).^2,2)./nansum((tmpActual - nanmean(tmpActual,2)).^2,2));
                    UniqEV(:,kernelid,areaid,hemid) = TotalR-TotalRLeaveOut;

                    % Draw Left Stimulus Rewarded
                    subplot(length(AREAS),2,(areaid-1)*2+hemid)
                    hold on
                    h(1)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==1&RewardVec==1),2)),'color',cols(1,:),'LineStyle',linstls{1});
                    plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==1&RewardVec==1),2)),'color',cols(1,:),'LineStyle',linstls{2})

                    % Right stimulus Rewarded
                    h(2)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==2&RewardVec==1),2)),'color',cols(2,:),'LineStyle',linstls{1});
                    plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==2&RewardVec==1),2)),'color',cols(2,:),'LineStyle',linstls{2})

                    % Left not rewarded
                    h(3)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==1&RewardVec==2),2)),'color',cols(3,:),'LineStyle',linstls{1});
                    plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==1&RewardVec==2),2)),'color',cols(3,:),'LineStyle',linstls{2})

                    % Right stimulus not Rewarded
                    h(4)=plot(timevec,squeeze(nanmean(tmpActual(:,SideVec==2&RewardVec==2),2)),'color',cols(4,:),'LineStyle',linstls{1});
                    plot(timevec,squeeze(nanmean(tmpPred(:,SideVec==2&RewardVec==2),2)),'color',cols(4,:),'LineStyle',linstls{2})


                    %             title(AREAGROUPNAMES{areaid})
                    if hemid==1
                        ylabel([AREAGROUPNAMES{areaid}])
                    end
                    if areaid==1
                        title(HemOpt{hemid})
                    end
                    FigureDefault

                    xlim([-100 timevec(end)])
                end
            end
            legend([h(:)],{'Left Rewarded','Right Rewarded','Left Error','Right Error'})
            saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' OptoOpt{optid} '_' AllKernelNames{kernelid} '.fig']))
            saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' OptoOpt{optid} '_' AllKernelNames{kernelid} '.bmp']))


        end


        cols = distinguishable_colors(length(AllKernels));
        clear h
        for areaid=1:length(AREAS)
            figure('name',[AREAGROUPNAMES{areaid} ' VarianceExplainedDifferentKernels'])

            for hemid=1:2
                subplot(1,2,hemid)
                hold on
                for kernelid=1:length(AllKernels)
                    h(kernelid) = plot(timevec./1000,UniqEV(:,kernelid,areaid,hemid).*100,'-','color',cols(kernelid,:));
                end
                FigureDefault
                xlabel('Time (s)')
                ylabel('\DeltaR^2')
                if areaid==1
                    title(HemOpt{hemid})
                end
            end
            legend(h(:),AllKernelNames)
            saveas(gcf,fullfile(LocalFolder,miceopt{midx},['ExplainedKernelVariance_LinearModels_' AREAGROUPNAMES{areaid} '_' OptoOpt{optid} '.fig']))
            saveas(gcf,fullfile(LocalFolder,miceopt{midx},['ExplainedKernelVariance_LinearModels_' AREAGROUPNAMES{areaid} '_' OptoOpt{optid} '.bmp']))

        end
    end
end

TW = {[-500 0],[0 1000],[1000 2000],[2000 4000],[4000 6000]};
TimeWindowNames = {'Baseline','Stimulus','Delay','Response','Post-trial'};
%% Summarize
nshuffle = 500
vecvector=2:4; % Which to shuffle analysis?
SwitchNames = {'Switch','Same'};
for optid=1
    p_NextTrialPerf = nan(length(AllKernels),length(miceopt));

    MSEAllMice = nan(length(timevec),length(AllKernels),length(AREAS),2,length(miceopt));
    UniqEVAllMice = MSEAllMice; % Difference between totalR and leaveoutR
    %     RExcl_PerVar_Shuffled
    RtotalAcrossMice = nan(length(AREAS),2,length(miceopt));
    RShuffleAcrossMice = nan(length(AllKernels),nshuffle,length(AREAS),2,length(miceopt));
    for midx = 1:length(miceopt)
        % Load LM Data
        tmp = matfile(fullfile(LocalFolder,miceopt{midx},['LinearModels_Opto' num2str(optid) '.mat']));

        % Extract all the fields
        disp('Extracting data for this mouse')
        s=fieldnames(tmp);
        for id = 1:length(s)
            eval([s{id} '=tmp.' s{id} ';'])
        end
        model = load(fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));
        model = model.BrainModel;
        MSE = nan(length(timevec),length(AllKernels),length(AREAS),2);
        UniqEV = MSE; % Difference between totalR and leaveoutR
        for kernelid=vecvector
            for areaid=1:length(AREAS)
                for hemid=1:2
                    % Find area pixels
                    idx = find(ismember(model.Rnames,AREAS{areaid}));
                    if isempty(idx)
                        idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),model.Rnames,'UniformOutput',0)));
                        areamask = false(400,400);
                        for ridx = 1:length(idx)
                            areamask(logical(round(model.Regions{idx(ridx)}))')=true;
                        end
                    else
                        areamask = logical(any(cat(3,model.Regions{idx}),3))';
                    end
                    if hemid==2 %left hemisphere
                        areamask(1:round(model.Lambda(2)),:)=0;
                    else
                        areamask(round(model.Lambda(2)):end,:)=0;
                    end
                    % Downsample areamask
                    areamask = Nan_imresize(areamask,newpixsz./400);
                    areamask(areamask>0)=1;
                    pixidx = find(areamask);
                    %Extract prediction for every timepoint and trial
                    tmpPred = reshape(PredicteddFF_LeaveOut(pixidx,:,kernelid),length(pixidx),length(timevec),ntrials); %Predicted leave out this kernel
                    tmpActual = reshape(DatPerPix(pixidx,:),length(pixidx),length(timevec),ntrials); %Actual dFF
                    tmpPrdFF = reshape(PredicteddFF(pixidx,:),length(pixidx),length(timevec),ntrials); % predicted total model

                    % Calculate UniqueExplainedVariancePerTP
                    MSEAllMice(:,kernelid,areaid,hemid,midx) = nanmean(squeeze(nanmean(tmpPred,1)-nanmean(tmpActual,1)).^2,2);
                    % R2 for leave out
                    % R2 total model
                    tmpR = 1 - nansum((squeeze(nanmean(tmpActual,1)) - squeeze(nanmean(tmpPred,1))).^2,2)./nansum((squeeze(nanmean(tmpActual,1)) - nanmean(squeeze(nanmean(tmpActual,1)),2)).^2,2);
                    TotalRtmp = 1 - nansum((squeeze(nanmean(tmpActual,1)) - squeeze(nanmean(tmpPrdFF,1))).^2,2)./nansum((squeeze(nanmean(tmpActual,1)) - nanmean(squeeze(nanmean(tmpActual,1)),2)).^2,2);
                    UniqEVAllMice(:,kernelid,areaid,hemid,midx) = TotalRtmp-tmpR;

                    % Save totalR
                    RtotalAcrossMice(areaid,hemid,midx) = nanmean(nanmean(TotalR(pixidx,:),2),1);
                    % Save shuffled R
                    RShuffleAcrossMice(kernelid,:,areaid,hemid,midx) = squeeze(nanmean(nanmean(RExcl_PerVar_Shuffled(pixidx,kernelid,:,:),3),1));

                end
            end
        end

        % Behavior
        TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
        tmptrials = TMPDat.Trials2Keep;

        SideVec = nan(1,0);
        RewardVec = nan(1,0);
        RewardNextTrial = nan(1,0);
        ChoiceVec = nan(1,0);
        LeftActionVec = nan(length(timevec),0);
        RightActionVec = LeftActionVec;
        DifMotionVec = LeftActionVec;
        EyeDataVec = LeftActionVec;
        TrialVec = nan(1,0);
        SesVec = nan(1,0);
        for ridx=1:3
            for sidx=1:length(AnalysisParameters.SideOpt)

                % Extract action
                tmpdatFiles = dir(fullfile(LocalFolder,miceopt{midx},'**',['Log.mat']));
                tmpses = unique(tmptrials(2,:));
                if length(tmpdatFiles)<length(tmpses)
                    keyboard
                end
                for sesid = 1:length(tmpses)
                    % Get correct log file
                    Log = matfile(fullfile(tmpdatFiles(sesid).folder,tmpdatFiles(sesid).name));
                    TimeLineHere = Log.TimeLineHere;
                    Log = Log.Log;
                    takethesetrials = tmptrials(1,tmptrials(2,:)==tmpses(sesid));
                    side = Log.Side(takethesetrials);
                    reaction = Log.Reaction(takethesetrials);
                    opto = Log.Opto(takethesetrials);
                    gavepassive = Log.Gavepassive(takethesetrials);
                    delay = Log.currentdelay(takethesetrials);
                    %Extract trials included in dataset
                    takethesetrials = takethesetrials(strcmp(side,SideOpt{sidx})&strcmp(reaction,ReactionOpt{ridx})&opto==optid-1&gavepassive==0&delay==1500);


                    if isempty(takethesetrials)
                        continue
                    end
                    % Extract Motion and downsample in time
                    if isfield(Log,'DifMotion')
                        DifMotionVec = cat(2,DifMotionVec,resample(Log.DifMotion(TimeLineHere>=timevec(1)&TimeLineHere<=timevec(end),takethesetrials),length(timevec),length(TimeLineHere)));
                    else
                        DifMotionVec = cat(2,DifMotionVec,nan(length(timevec),length(takethesetrials)));
                    end
                    % Extract eye movements
                    if isfield(Log,'EyeMovement')
                        tmpeye = nanmean(Log.EyeMovement(TimeLineHere>=timevec(1)&TimeLineHere<=timevec(end),takethesetrials,:),3); %take average eyemotion across the 4 coordinates
                        EyeDataVec = cat(2,EyeDataVec,resample(tmpeye,length(timevec),length(TimeLineHere)));
                    else
                        EyeDataVec = cat(2,EyeDataVec,nan(length(timevec),length(takethesetrials)));
                    end

                    % Extract licks
                    tmplicksRight = arrayfun(@(X) Log.RTrightVec{X},takethesetrials,'UniformOutput',0);
                    tmplicksLeft = arrayfun(@(X) Log.RTleftVec{X},takethesetrials,'UniformOutput',0);

                    tmpactionL = zeros(length(timevec),size(tmplicksLeft,2)); % Left first
                    tmpactionR = zeros(length(timevec),size(tmplicksRight,2)); % Now right

                    % Set to 'event' at correct time
                    for trid=1:size(tmplicksLeft,2)
                        tmp=arrayfun(@(Y) abs(timevec-Y)<=SR/2,tmplicksLeft{trid},'UniformOutput',0);
                        if ~isempty(tmp)
                            tmpactionL(:,trid)=nansum(cat(1,tmp{:}),1);
                        end
                        tmp=arrayfun(@(Y) abs(timevec-Y)<=SR/2,tmplicksRight{trid},'UniformOutput',0);
                        if ~isempty(tmp)
                            tmpactionR(:,trid)=nansum(cat(1,tmp{:}),1);
                        end
                    end
                    LeftActionVec = cat(2,LeftActionVec,tmpactionL);
                    RightActionVec = cat(2,RightActionVec,tmpactionR);
                    ntrialstmp = length(takethesetrials);
                    % Assign condition
                    SideVec = cat(2,SideVec,repmat(sidx,1,ntrialstmp));
                    RewardVec = cat(2,RewardVec,repmat(ridx,1,ntrialstmp));
                    if ridx==1
                        ChoiceVec = cat(2,ChoiceVec,repmat(sidx,1,ntrialstmp));
                    elseif ridx==2
                        if sidx==1
                            ChoiceVec = cat(2,ChoiceVec,repmat(2,1,ntrialstmp));
                        elseif sidx==2
                            ChoiceVec = cat(2,ChoiceVec,repmat(1,1,ntrialstmp));
                        end
                    elseif ridx==3 % No choice made
                        ChoiceVec = cat(2,ChoiceVec,repmat(0,1,ntrialstmp));
                    end
                    Log.Reaction(cellfun(@isempty,Log.Reaction)) = {'none'};

                    tmpchoice = zeros(1,length(Log.Reaction));
                    tmpchoice(ismember(Log.Reaction,'Hit')&ismember(Log.Side,'right')) = 1; %right
                    tmpchoice(ismember(Log.Reaction,'Error')&ismember(Log.Side,'left')) = 1; %right
                    tmpchoice(ismember(Log.Reaction,'Hit')&ismember(Log.Side,'left')) = 2; %left
                    tmpchoice(ismember(Log.Reaction,'Error')&ismember(Log.Side,'right')) = 2; %left
                    takenexttrial = takethesetrials+1;
                    if max(takenexttrial)>length(Log.Reaction)
                        addone = 1;
                        takenexttrial(takenexttrial>length(Log.Reaction)) = [];
                    else
                        addone = 0;
                    end
                    same = arrayfun(@(X) ismember(tmpchoice(takenexttrial(X)),tmpchoice(takethesetrials(X))),1:length(takenexttrial));
                    RewardNextTrial = cat(2,RewardNextTrial,SwitchNames(same+1));
                    if addone
                        RewardNextTrial = cat(2,RewardNextTrial,{'Unknown'});
                    end


                    TrialVec = cat(2,TrialVec,takethesetrials);
                    SesVec = cat(2,SesVec,repmat(tmpses(sesid),1,length(takethesetrials)));

                end
            end
        end

        %Reshape data
        TotalTrialVec = TrialVec + SesVec*1000; %add surreal number to do ordering
        [~,sortidx] = sort(TotalTrialVec);

        % Performance ordered by trial
        orderedreward = RewardNextTrial(sortidx);
        %Name to number
        orderedrewardNr = nan(size(orderedreward));
        orderedrewardNr(ismember(orderedreward,'Hit')) = 1;
        orderedrewardNr(ismember(orderedreward,'Error')) = 2;
        orderedrewardNr(ismember(orderedreward,'Miss')) = 3;

        tmpActual = squeeze(nanmean(reshape(DatPerPix,[],length(timevec),ntrials),1)); %Actual dFF
        tmpPrdFF = squeeze(nanmean(reshape(PredicteddFF,[],length(timevec),ntrials),1)); % predicted total model
        twid = 5;

        GeneralPerf = 1-nansum((tmpActual(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:)-tmpPrdFF(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:)).^2,1)./nansum((tmpActual(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:)-nanmean(tmpActual(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:),1)).^2,1);

        figure('name','Performance per trial')
        %Extract prediction for every timepoint and trial
        for kernelid = vecvector
            tmpPred = squeeze(nanmean(reshape(PredicteddFF_LeaveOut(:,:,kernelid),[],length(timevec),ntrials),1)); %Predicted leave out this kernel
            % Difference in activity
            KernelPerf = 1-nansum((tmpActual(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:)-tmpPred(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:)).^2,1)./nansum((tmpActual(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:)-nanmean(tmpActual(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:),1)).^2,1);

            deltaR = GeneralPerf - KernelPerf;

            takethese = ismember(orderedreward,{'Switch','Same'}) & ~isnan(nanmean(tmpPred,1));

            subplot(2,2,kernelid-1)
            violinplot(deltaR(takethese),orderedreward(takethese))
            xlabel('Condition'); ylabel('\DeltaR')

            title(AllKernelNames{kernelid})


            p_NextTrialPerf(kernelid,midx) = anovan(deltaR(takethese),{orderedreward(takethese)},'display','off');

        end

        drawnow
    end

    %% Shuffle versus total R
    figure('name','ShuffleVersusRDistribution')
    for kernelid=1:length(vecvector)
        subplot(1,length(vecvector),kernelid)
        %                         RtotalAcrossMice = nan(length(AREAS),2,length(miceopt));
        %     RShuffleAcrossMice = nan(length(AllKernels),nshuffle,length(AREAS),2,length(miceopt));fleAcrossMice = nan(length(timevec),length(AllKernels),nshuffle,length(AREAS),2,length(miceopt));
        tmpRtotal = squeeze(nanmean(RtotalAcrossMice,2));
        tmpShuff = permute(reshape(nanmean(RShuffleAcrossMice(vecvector(kernelid),:,:,:,:),4),nshuffle,length(AREAS),[]),[2,1,3]);
        tmpShuff = repmat(tmpRtotal,1,1,nshuffle)-permute(tmpShuff,[1,3,2]);
        h=Groupedviolinplot(permute(tmpShuff,[3,1,2]).*100,string(AREAGROUPNAMES),string(miceopt));
        legend off
        hold on
        for hid=1:length(h)
            h(hid).ViolinAlpha={[1]};
            h(hid).ShowMedian = 0;
        end
        set(gca,'XTickLabelRotation',90)
        if kernelid>1
            set(gca,'XTickLabel','')
        end
        line(get(gca,'xlim'),[0 0],'color',[0 0 0],'LineStyle','--')
        ylim([-1 6])
        for areaid = 1:length(AREAS)
            for midx = 1:length(miceopt)
                if quantile(tmpShuff(areaid,midx,:),0.001)>0
                    text((areaid-1)*10+midx,5.9-midx*0.075,'***')
                elseif quantile(tmpShuff(areaid,midx,:),0.01)>0
                    text((areaid-1)*10+midx,5.9-midx*0.075,'**')
                elseif quantile(tmpShuff(areaid,midx,:),0.05)>0
                    text((areaid-1)*10+midx,5.9-midx*0.075,'*')
                end
            end
        end

        FigureDefault
        ylabel('\DeltaR^2')
        title(AllKernelNames{vecvector(kernelid)})


    end
    %% Time courses (zoomed in)
    UniqEVAllMice = reshape(UniqEVAllMice,length(timevec),length(AllKernelNames),length(AREAS),[]); % Take area as N

    cols = [1 0 0; 0 0 1; 0 0.5 0; 0 0 0; 0.5 0.5 0.5];
    clear h
    figure('name',['VarianceExplainedDifferentKernels'])
    for areaid=1:length(AREAS)
        for twid=1:length(TW)
            subplot(length(AREAGROUPNAMES),length(TW),(areaid-1)*length(TW)+twid)
            hold on
            for kernelid=vecvector
                h(kernelid) = shadedErrorBar(timevec(timevec>=TW{twid}(1)&timevec<=TW{twid}(2))./1000,nanmean(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),kernelid,areaid,:).*100,4),nanstd(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),kernelid,areaid,:).*100,[],4)./sqrt(length(miceopt)*2),{'-','color',cols(kernelid-min(vecvector)+1,:)},0);
            end
            FigureDefault
            if areaid==length(AREAGROUPNAMES)
                xlabel('Time (s)')
            else
                set(gca,'XTickLabel','')
            end
            if twid==1
                ylabel(AREAGROUPNAMES{areaid})
            else
                %                 ylabel('\DeltaR^2 (%)')
            end
            if areaid==1
                title(TimeWindowNames{twid})
            end
            if ismember(twid,[1,3])
                ylim([-1 5])
            elseif ismember(twid,[2,5])
                ylim([-1 10])
            else
                ylim([-1 20])
            end
        end
    end
    legend([h(:).mainLine],AllKernelNames(vecvector))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels_' OptoOpt{optid} '.fig']))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels_' OptoOpt{optid} '.bmp']))

    %%
    UniqEVAllMice = reshape(UniqEVAllMice,length(timevec),length(AllKernels),length(AREAS),2,length(miceopt));
    SumperTW = nan(length(TW),length(AllKernels),length(AREAS),2,length(miceopt));
    for twid=1:length(TW)
        SumperTW(twid,:,:,:,:) = squeeze(nanmean(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:,:,:,:),1));
    end

    %%
    kernelidx = vecvector;
    figure('name',['barplots'])

    %     kernelidx = [2,3];
    for twid=1:length(TW)
        subplot(1,length(TW),twid)
        tmp = reshape(SumperTW(twid,kernelidx,:,:,:).*100,length(kernelidx),length(AREAS),[]);
        h=bar(nanmean(tmp,3));
        hold on
        for areaid=1:length(AREAS)
            h(areaid).EdgeColor = 'none';
            h(areaid).BarWidth = 1;
            for kernelid=1:length(kernelidx)
                scatter(repmat(h(areaid).XEndPoints(kernelid),1,length(miceopt)*2),reshape(SumperTW(twid,kernelidx(kernelid),areaid,:,:).*100,1,[]),10,[0.2 0.2 0.2],'filled')
            end
        end
        set(gca,'XTick',1:length(kernelidx),'XTickLabel',AllKernelNames(kernelidx),'YAxisLocation','right')
        ylabel(['\DeltaR^2 %'])
        ylim([-2 22])
        FigureDefault

        % Statistics across hemispheres
        Ytmp = permute(tmp,[3,2,1]);
        tmpstr = arrayfun(@(X) arrayfun(@(Y) [AllKernelNames{X} '_AREA' num2str(Y)],1:length(AREAS),'Uni',0),kernelidx,'UniformOutput',0);
        tmpstr = [tmpstr{:}];


        kernelvec =repmat([1:length(kernelidx)]',[1,length(AREAGROUPNAMES)])';
        kernelvec = kernelvec(:);

        areavec = repmat([1:length(AREAGROUPNAMES)]',[1,length(kernelidx)]);
        areavec = areavec(:);


        ttmp = array2table(reshape(Ytmp,size(Ytmp,1),[]),'VariableNames',tmpstr);
        withinDesign = table(kernelvec,areavec,'VariableNames',{'Kernel','Area'});
        withinDesign.Kernel=categorical(withinDesign.Kernel);
        withinDesign.Area=categorical(withinDesign.Area);

        try
            Mdl = fitrm(ttmp,[tmpstr{1} '-' tmpstr{end} '~1'],'WithinDesign',withinDesign);
            ranovatbl = ranova(Mdl,'WithinModel','Kernel*Area');
            title([TimeWindowNames{twid}])
        catch
        end


        writetable(ttmp,fullfile(LocalFolder,['TABLE_' TimeWindowNames{twid} '_' OptoOpt{optid} '.csv']))
    end
    legend(h,AREAGROUPNAMES)
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BarPlot_' OptoOpt{optid} '.fig']))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BarPlot_' OptoOpt{optid} '.bmp']))

    %% Test whether these are significantly more than 0 across mice
    p = nan(length(TW),length(kernelidx),length(AREAS));
    for twid=1:length(TW)
        for kernelid=1:length(kernelidx)
            for areaid=1:length(AREAS)
                tmp = squeeze(SumperTW(twid,kernelidx(kernelid),areaid,:,:));
                h = kstest(tmp(:));
                if ~h
                    [~,p(twid,kernelid,areaid)] = ttest(tmp(:),0,'tail','right');
                else
                    p(twid,kernelid,areaid) = signrank(tmp(:),0,'tail','right');
                end
                %                end
            end
            p(twid,kernelid,:) = bonf_holm(p(twid,kernelid,:));

        end

    end
    % Information for figure 1B
    for areaid = 1:length(AREAS)
        for kernelid = 1:length(kernelidx)
            if any(p(1,kernelid,areaid)<0.05)
                disp([TimeWindowNames{1}  ' ' AREAGROUPNAMES{areaid} ' ' AllKernelNames{kernelidx(kernelid)} ' significant'])
            end
        end
    end

    % Information for figure 1B
    for areaid = 1:length(AREAS)
        for kernelid = 1:length(kernelidx)
            if any(p(1,kernelid,areaid)<0.05)
                disp(['Baseline Encoding: ' AREAGROUPNAMES{areaid} ' ' AllKernelNames{kernelidx(kernelid)} ' significant'])
            end
        end
    end
    % Information for figure 1B
    for areaid = 1:length(AREAS)
        for kernelid = 1:length(kernelidx)
            if any(p(2:4,kernelid,areaid)<0.05)
                disp(['Event Encoding: ' AREAGROUPNAMES{areaid} ' ' AllKernelNames{kernelidx(kernelid)} ' significant'])
            end
        end
    end
    % Information for figure 1B
    for areaid = 1:length(AREAS)
        for kernelid = 1:length(kernelidx)

            if any(p(5,kernelid,areaid)<0.05)
                disp(['Post-trial: ' AREAGROUPNAMES{areaid} ' ' AllKernelNames{kernelidx(kernelid)} ' significant'])
            end
        end
    end
end

%% Behavior Plots (Fig 5A)
if 1
    for midx = 1:length(miceopt)
        TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
        tmptrials = TMPDat.Trials2Keep;
        for optid = 1:length(OptoOpt)
            for sidx=1:length(AnalysisParameters.SideOpt)

                SideVec = nan(1,0);
                RewardVec = nan(1,0);
                TrialVec = nan(1,0);
                for ridx=1:3

                    % Extract action
                    tmpdatFiles = dir(fullfile(ServerDir,miceopt{midx},'*','*',[miceopt{midx} '*RawData_All.mat']));
                    tmpses = unique(tmptrials(2,:));
                    if length(tmpdatFiles)<length(tmpses)
                        keyboard
                    end
                    for sesid = 1:length(tmpses)
                        % Get correct log file
                        Log = matfile(fullfile(tmpdatFiles(tmpses(sesid)).folder,tmpdatFiles(tmpses(sesid)).name));
                        TimeLineHere = Log.timeline;
                        Log = Log.Log;
                        takethesetrials = tmptrials(1,tmptrials(2,:)==tmpses(sesid));
                        side = Log.Side(takethesetrials);
                        reaction = Log.Reaction(takethesetrials);
                        opto = Log.Opto(takethesetrials);
                        gavepassive = Log.Gavepassive(takethesetrials);
                        delay = Log.currentdelay(takethesetrials);
                        %Extract trials included in dataset
                        takethesetrials = takethesetrials(strcmp(side,SideOpt{sidx})&strcmp(reaction,ReactionOpt{ridx})&opto==optid-1&gavepassive==0&delay==1500);

                        SideVec = cat(2,SideVec,repmat(sidx,1,length(takethesetrials)));
                        RewardVec = cat(2,RewardVec,repmat(ridx,1,length(takethesetrials)));
                        TrialVec = cat(2,TrialVec,takethesetrials);

                    end
                end
                PerformanceAcrossmice(sidx,optid,midx) = sum(RewardVec==1)./sum(RewardVec~=3);
                OmissionsAcrossmice(sidx,optid,midx)  = sum(RewardVec==3)./length(RewardVec);
                figure; scatter(TrialVec,RewardVec==3,10,[0 0 0],'filled');
                title([AnalysisParameters.SideOpt{sidx} ' ' OptoOpt{optid} ' ' miceopt{midx}])
                drawnow
            end
        end
    end


    % Statistics
    PerformanceAcrossmice
    g1 = repmat(SideOpt',[1,length(OptoOpt),length(miceopt)]);
    g2= repmat(OptoOpt',[1,length(SideOpt),length(miceopt)]);
    g2 = permute(g2,[2,1,3]);
    g3 = repmat(GenOpt',[1,length(SideOpt),length(OptoOpt)]);
    g3 = permute(g3,[2,3,1]);
    g4 = repmat(miceopt',[1,length(SideOpt),length(OptoOpt)]);
    g4 = permute(g4,[2,3,1]);

    % First across all mice
    tbl = table(PerformanceAcrossmice(:),g1(:),g2(:),g4(:),'VariableNames',{'Performance','Side','OptoEpoch','MOUSEID'});
    glme = fitglme(tbl,'Performance ~ 1 + OptoEpoch + Side + OptoEpoch:Side + (1|MOUSEID)','Distribution','Normal');
    p=anova(glme)

    %
    figure('name','Performance no opto')
    h = bar(squeeze(nanmean(PerformanceAcrossmice(:,1,:),3)));
    hold on
    line(get(gca,'xlim'),[0.5 0.5],'LineStyle','--','color',[0 0 0])
    h.EdgeColor = 'none';
    h.FaceColor = [0.5 0.5 0.5];
    cols = lines(length(miceopt));
    for midx=1:length(miceopt)
        plot([1,2],squeeze(PerformanceAcrossmice(:,1,midx)),'color',cols(midx,:))
    end
    set(gca,'XTickLabel',SideOpt)
    ylabel('Behavioral accuracy')

    p=nan(1,2);
    for sidx = 1:2

        [~,p(sidx)] = ttest(squeeze(PerformanceAcrossmice(sidx,1,:))-0.5);
    end
    bonf_holm(p)
    FigureDefault



    % Now take your pick
    idx = ismember(g3(:),'Drd1')

    tbl = table(OmissionsAcrossmice(idx),g1(idx),g2(idx),g4(idx),'VariableNames',{'Performance','Side','OptoEpoch','MOUSEID'});
    glme = fitglme(tbl,'Performance ~ 1 + OptoEpoch + Side + OptoEpoch:Side + (1|MOUSEID)','Distribution','Normal');
    p=anova(glme)


    figure('name','OptoEffect')
    linestls = {'-','--'};
    tmplabel = {OptoOpt{:} 'ITI Opto'};
    optorder = [1,7,2:6];
    for gidx = 1:length(GenoTypeOpt)
        subplot(2,1,gidx)
        h=plot(1:length(OptoOpt)+1,repmat(nanmean(nanmean(PerformanceAcrossmice(:,1,ismember(genotype,GenoTypeOpt{gidx})),1),3),1,length(OptoOpt)+1),...
            'Color',[0 0 0],'LineStyle',linestls{2});
        hold on
        for sidx=1:2
            h=shadedErrorBar(1:length(OptoOpt)+1,nanmean(PerformanceAcrossmice(sidx,optorder,ismember(genotype,GenoTypeOpt{gidx})),3),...
                nanstd(PerformanceAcrossmice(sidx,optorder,ismember(genotype,GenoTypeOpt{gidx})),[],3)./sqrt(sum(ismember(genotype,GenoTypeOpt{gidx}))-1),...
                {'Color',cols(sidx,:),'LineStyle',linestls{1}},1);

            h=shadedErrorBar(1:length(OptoOpt)+1,nanmean(OmissionsAcrossmice(sidx,optorder,ismember(genotype,GenoTypeOpt{gidx})),3),...
                nanstd(OmissionsAcrossmice(sidx,optorder,ismember(genotype,GenoTypeOpt{gidx})),[],3)./sqrt(sum(ismember(genotype,GenoTypeOpt{gidx}))-1),...
                {'Color',cols(sidx,:),'LineStyle',linestls{2}},1);
        end
        title(GenoTypeOpt{gidx})
        box off
        set(gca,'XTick',[1:length(OptoOpt)+1],'XTickLabel',tmplabel(optorder))

    end
    saveas(gcf,fullfile(storepath,'Figures','Behavior',['Behavior_OptoWM.fig']))
end

%% Across mice - save out time course per area, keep trial information
if exist(fullfile(LocalFolder,'AllMiceTC.mat'))
    load(fullfile(LocalFolder,'AllMiceTC.mat'))
else
    for midx=1:length(miceopt)

        TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
        Trials2Keep = TMPDat.Trials2Keep;
        AllDat = TMPDat.AllDat;

        % get model
        model = load(fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));


        % AllTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(ReactionOpt),length(SideOpt),length(HemOpt),length(OptoOpt));
        if midx==1
            timeline = TMPTC.timeline;
            AllMiceTC = cell(length(OptoOpt),2,2,length(AREAS),2,length(miceopt));
        end

        % Average Across Grouped AreasTake
        for areaid=1:length(AnalysisParameters.AREASGrouped)
            for hemid=1:2
                % Find area pixels
                idx = find(ismember(model.BrainModel.Rnames,AREAS{areaid}));
                if isempty(idx)
                    idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),model.BrainModel.Rnames,'UniformOutput',0)));
                    areamask = false(400,400);
                    for i = 1:length(idx)
                        areamask(logical(round(model.BrainModel.Regions{idx(i)}))')=true;
                    end
                else
                    areamask = logical(any(cat(3,model.BrainModel.Regions{idx}),3))';
                end
                if hemid==2 %left hemisphere
                    areamask(1:round(model.BrainModel.Lambda(2)),:)=0;
                else
                    areamask(round(model.BrainModel.Lambda(2)):end,:)=0;
                end

                pixidx = find(areamask);
                %Extract data
                tmp = squeeze(AllDat(:,ismember(ReactionOpt,{'Hit','Error'}),:));
                %reshape to index with pixidx
                tmp = cellfun(@(X) reshape(X,400*400,size(X,3),size(X,4)),tmp,'UniformOutput',0);
                %Average over the pixels
                tmp = cellfun(@(X) squeeze(nanmean(X(pixidx,:,:),1)),tmp,'UniformOutput',0);

                % Save out
                AllMiceTC(:,:,:,areaid,hemid,midx) = tmp;
            end
        end
        save(fullfile(LocalFolder,'AllMiceTC.mat'),'AllMiceTC')
    end
end
% Time courses across hemispheres
ColOpt = [0 0 0; 1 0 0; 0.5 0.5 0.5; 0 0 1];
LineOpt = {'-','-- '}; % ContraOpto / IpsiOpto

% % For WM paper
ylims = [-.02 0.02];
xlims = [-600,2000];
optiidvec = {1, 2};


clear h
maxnrt = 700;
for areaid=1:length(AREAGROUPNAMES)
    figure('name',['AllMiceTCEarlyvsOff ' [AREAS{areaid}{:}]])
    for genidx=1:2
        for hemid=1:2
            subplot(2,2,(genidx-1)*2+hemid)
            tmp4ana = nan(length(timeline),2,2,maxnrt);
            mousentrials = nan(2,2,maxnrt);
            for optiid = 1:2
                for sidx=1:2 %First should be contra
                    tmp = reshape(cat(2,AllMiceTC(optiidvec{optiid},:,sidx,areaid,hemid,ismember(GenOpt,GenoTypeOpt{genidx}))),[],sum(ismember(GenOpt,GenoTypeOpt{genidx})));
                    if any(any(any(cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0)))))
                        tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))} = tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))}';
                    end
                    %                     % Concatenate trials of different conditions
                    %                     tmp = cellfun(@(X) smooth(nanmean(X,2),3),tmp,'UniformOutput',0); %Average across trials

                    % Get number of trials per mouse
                    mousen = cell2mat(cellfun(@(X) size(X,2),arrayfun(@(X) cat(2,tmp{:,X}),1:sum(ismember(GenOpt,GenoTypeOpt{genidx})),'UniformOutput',0),'UniformOutput',0));
                    mousen = arrayfun(@(X) repmat(X,1,mousen(X)),1:length(mousen),'UniformOutput',0);
                    mousentrials(optiid,sidx,1:length(cat(2,mousen{:}))) = cat(2,mousen{:});
                    tmpall = cat(2,tmp{:}); %Concatenate different mice/Trials
                    tmp = reshape(tmpall,length(timeline),[]);
                    tmp4ana(:,optiid,sidx,1:size(tmp,2)) = tmp;
                    if size(tmp,2)>maxnrt
                        disp('Choose bigger maxnrt')
                        keyboard
                    end
                    h(optiid,sidx) = shadedErrorBar(timeline,(squeeze(nanmean(tmp,2))),(squeeze(nanstd(tmp,[],2)))./sqrt(size(tmp,2)-1),...
                        {'color',ColOpt((genidx-1)*2+optiid,:),'LineStyle',LineOpt{sidx},'LineWidth',2},1)
                    hold on
                end
            end

            box off
            xlim(xlims)
            title([HemOpt{hemid} ' ' AREAGROUPNAMES{areaid}])
            ylim(ylims)
            line([0 0],ylims,'LineStyle','--','color', [0 0 0])
            line([500 500],ylims,'LineStyle','--','color', [0 0 0])
            line([2000 2000],ylims,'LineStyle','--','color', [0 0 0])
            line([4000 4000],ylims,'LineStyle','--','color', [0 0 0])
            line([OptoTW{optiidvec{2}}(1) OptoTW{optiidvec{2}}(1)],ylims,'LineStyle','--','color',[1 0 0])
            line([OptoTW{optiidvec{2}}(2) OptoTW{optiidvec{2}}(2)],ylims,'LineStyle','--','color',[1 0 0])


            % statistical test - Mixed Linear Model
            opto = repmat({'OFF','ON'}',[1,2,size(tmp4ana,4)]);
            side = repmat({'Left','Right'}',[1,2,size(tmp4ana,4)]);
            side = permute(side,[2,1,3]);


            p = nan(length(timeline),3);
            parfor tp = 1:length(timeline)
                tmp = tmp4ana(tp,:,:,:);
                idx = find(~isnan(tmp(:)));
                tbl = table(double(tmp(idx)),opto(idx),side(idx),mousentrials(idx),'VariableNames',{'dFF','Opto','Side','MOUSEID'});
                try

                    glme = fitglme(tbl,'dFF ~ 1 + Opto + Side + Opto:Side + (1|MOUSEID)','Distribution','Normal');
                    tmp = anova(glme);
                    p(tp,:) = tmp.pValue(2:end);
                catch ME
                    disp(ME)
                end
            end

            %Opto
            plot(timeline(find(p(:,1)<0.05)),repmat(ylims(2)*0.9,[1,sum(p(:,1)<0.05)]),'.','color',ColOpt(genidx*2,:),'MarkerSize',8)

            %Stimulus side
            plot(timeline(find(p(:,2)<0.05)),repmat(ylims(2)*0.85,[1,sum(p(:,2)<0.05)]),'.','color',[0 0 0],'MarkerSize',8)

            %Interaction side
            plot(timeline(find(p(:,3)<0.05)),repmat(ylims(2)*0.80,[1,sum(p(:,3)<0.05)]),'.','color',[0.6 0.6 0.6],'MarkerSize',8)

            drawnow
        end
    end

    saveas(gcf,fullfile(LocalFolder,['AllHemTCOpto ' OptoOpt{optiidvec{1}} 'vs'  OptoOpt{optiidvec{2}} [AREAS{areaid}{:}] '.fig']))


end

%% Check if response to stimulus increases over time
clear h
ColOpt = distinguishable_colors(length(AREAGROUPNAMES));
maxnrt = 700;
% optiidvec = {1, [2,3,4]};
figure('name',['StimulusResponseOverTrials '])
r = nan(length(miceopt),length(AREAGROUPNAMES));
p = r;
Slopes = r;
SlopesnShuffle = nan(length(miceopt),length(AREAGROUPNAMES),nboot);
quantilep = r;
for midx=1:length(miceopt)
    subplot(ceil(sqrt(length(miceopt))),round(sqrt(length(miceopt))),midx)
    hold on
    for areaid=1:length(AREAGROUPNAMES)
        tmpvecdff = [];
        tmpvecid = [];
        for hemid=1:2


            tmp = AllMiceTC{1,1,1,areaid,hemid,midx}; %only hits
            tmpvec = nanmean(tmp(timeline>0&timeline<500,:),1);
            tmpvec(isnan(tmpvec)) = [];
            tmpvecid = [tmpvecid 1:length(tmpvec)];
            tmpvecdff = [tmpvecdff tmpvec];


        end
        % Normalize
        scatter(tmpvecid,tmpvecdff,10,ColOpt(areaid,:),'filled');

        pol = polyfit(tmpvecid,tmpvecdff,1);
        Slopes(midx,areaid) = pol(1);
        for bootid = 1:nboot
            pol = polyfit(datasample(tmpvecid,numel(tmpvecid),'replace',false),tmpvecdff,1);
            SlopesnShuffle(midx,areaid,bootid) = pol(1);
        end
        quantilep(midx,areaid) = 1-(invprctile(SlopesnShuffle(midx,areaid,:),Slopes(midx,areaid))/100);

        plot(1:max(tmpvecid),polyval(pol,1:max(tmpvecid)),'-','Color',ColOpt(areaid,:))

        [rtmp,ptmp] = corrcoef(tmpvecid,tmpvecdff);

        r(midx,areaid) = rtmp(1,2);
        p(midx,areaid) = ptmp(1,2);


    end

    xlabel('Trial')
    ylabel('VisualResponse')
    makepretty
    offsetAxes
    drawnow



end

subplot(ceil(sqrt(length(miceopt))),round(sqrt(length(miceopt))),midx+1)
scatter(quantile(SlopesnShuffle,0.95,3),Slopes,15,[0 0 0],'filled')
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
lims = [min([xlims(1) ylims(1)]) max([xlims(2) ylims(2)])];
set(gca,'xlim',lims,'ylim',lims)
hold on
line(lims,lims)
axis square
xlabel('95th percentile')
ylabel('Real Slope')
axis square
makepretty
offsetAxes

saveas(gcf,fullfile(LocalFolder,['StimulusResponseOverTrials.fig']))

%% Time courses across hemispheres - Opto X Response X Stim
ColOpt = [0 0 0; 1 0 0; 0.5 0.5 0.5; 0.8 0.5 0]; %Contra Off, Contra Opto, Ipsi Off, Ipsi Opto
LineOpt = {'-','-- '}; % Rewarded / Not Rewarded



% % For WM paper
ylims = [-.02 0.02];
xlims = [-600,2000];
optiidvec = {1, 4};

cols = lines(7);
clear h
maxnrt = 700;
% optiidvec = {1, [2,3,4]};
for areaid=1:length(AREAS)
    figure('name',['AllMiceTCEarlyvsOff ' [AREAS{areaid}{:}]])
    for genidx=1:2
        for hemid=1:2
            subplot(2,2,(genidx-1)*2+hemid)
            tmp4ana = nan(length(timeline),2,2,2,maxnrt);
            mousentrials = nan(2,2,2,maxnrt);
            for optiid = 1:2
                for sidx=1:2 %First should be contra
                    for ridx = 1:2
                        tmp = squeeze(AllMiceTC(optiidvec{optiid},ridx,sidx,areaid,hemid,ismember(GenOpt,GenoTypeOpt{genidx})));
                        if any(any(any(cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0)))))
                            tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))} = tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))}';
                        end
                        %                     % Concatenate trials of different conditions
                        %                     tmp = cellfun(@(X) smooth(nanmean(X,2),3),tmp,'UniformOutput',0); %Average across trials

                        % Get number of trials per mouse
                        mousen = cell2mat(cellfun(@(X) size(X,2),tmp,'UniformOutput',0));
                        mousen = arrayfun(@(X) repmat(X,1,mousen(X)),1:length(mousen),'UniformOutput',0);

                        mousentrials(optiid,sidx,ridx,1:length(cat(2,mousen{:}))) = cat(2,mousen{:});
                        tmpall = cat(2,tmp{:}); %Concatenate different mice/Trials
                        tmp = reshape(tmpall,length(timeline),[]);
                        tmp4ana(:,optiid,sidx,ridx,1:size(tmp,2)) = tmp;
                        if size(tmp,2)>maxnrt
                            disp('Choose bigger maxnrt')
                            keyboard
                        end
                        h(optiid,sidx,ridx) = shadedErrorBar(timeline,(squeeze(nanmean(tmp,2))),(squeeze(nanstd(tmp,[],2)))./sqrt(size(tmp,2)-1),...
                            {'color',ColOpt((sidx-1)*2+optiid,:),'LineStyle',LineOpt{ridx},'LineWidth',2},1);
                        hold on
                    end
                end
            end

            box off
            xlim(xlims)
            title([HemOpt{hemid} ' ' AREAGROUPNAMES{areaid}])
            ylim(ylims)
            line([0 0],ylims,'LineStyle','--','color', [0 0 0])
            line([500 500],ylims,'LineStyle','--','color', [0 0 0])
            line([2000 2000],ylims,'LineStyle','--','color', [0 0 0])
            line([4000 4000],ylims,'LineStyle','--','color', [0 0 0])
            line([OptoTW{optiidvec{2}}(1) OptoTW{optiidvec{2}}(1)],ylims,'LineStyle','--','color',[1 0 0])
            line([OptoTW{optiidvec{2}}(2) OptoTW{optiidvec{2}}(2)],ylims,'LineStyle','--','color',[1 0 0])


            % statistical test - Mixed Linear Model
            opto = repmat({'OFF','ON'}',[1,2,2,size(tmp4ana,5)]);
            side = repmat({'Left','Right'}',[1,2,2,size(tmp4ana,5)]);
            side = permute(side,[2,1,3,4]);
            reaction = repmat({'Hit','Error'}',[1,2,2,size(tmp4ana,5)]);
            reaction = permute(reaction,[2,3,1,4]);

            p = nan(length(timeline),7);
            parfor tp = 1:length(timeline)
                tmp = tmp4ana(tp,:,:,:);
                idx = find(~isnan(tmp(:)));
                tbl = table(double(tmp(idx)),opto(idx),side(idx),reaction(idx),mousentrials(idx),'VariableNames',{'dFF','Opto','Side','Reaction','MOUSEID'});
                try

                    glme = fitglme(tbl,'dFF ~ 1 + Opto*Side*Reaction + (1|MOUSEID)','Distribution','Normal');
                    tmp = anova(glme);
                    p(tp,:) = tmp.pValue(2:end);
                catch ME
                    disp(ME)
                end
            end
            for i = 1:7
                if any(find(p(:,i)<0.05))
                    hsig(i) = plot(timeline(find(p(:,i)<0.05)),repmat(ylims(2)*0.95-(0.001*i),[1,sum(p(:,i)<0.05)]),'.','color',cols(i,:),'MarkerSize',8);
                end
            end

            drawnow
        end
    end
    legend([hsig(:)],{'Opto','Side','Reaction','Opto*Side','Opto*Reaction','Side*Reaction','Side*Reaction*Opto'})
    saveas(gcf,fullfile(LocalFolder,['AllHemTCOpto ' OptoOpt{optiidvec{1}} 'vs'  OptoOpt{optiidvec{2}} [AREAS{areaid}{:}] '.fig']))


end

%% TC for Opto OFF all mice ()
ColOpt = [0 0 0; 0.5 0.5 0.5];
LineOpt = {'-','-- '};
ylims =[-0.035 0.052];
for optid=1 %:length(OptoOpt)
    for areaid=1:length(AREAS)
        figure('name',['AllMiceTC ' OptoOpt{optid} ' ' [AREAS{areaid}{:}]])
        for hemid=1:2
            for genid = 1:2
                subplot(2,2,(genid-1)*2+hemid)
                tmp4ana = nan(length(timeline),2,2,length(miceopt));
                for ridx=1:2
                    for sidx=1:2
                        tmp = squeeze(AllMiceTC(optid,ridx,sidx,areaid,hemid,:));
                        if any(cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0)))
                            tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))} = tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))}';
                        end
                        tmp = squeeze(cellfun(@(X) nanmean(X,2),tmp,'UniformOutput',0));
                        tmp = cat(2,tmp{:});
                        h(ridx,sidx) = shadedErrorBar(timeline,squeeze(nanmean(tmp,2)),squeeze(nanstd(tmp,[],2))./sqrt(length(miceopt)-1),...
                            {'color',ColOpt(ridx,:),'LineStyle',LineOpt{sidx},'LineWidth',2},1)
                        hold on
                        tmp4ana(:,ridx,sidx,:)=tmp;
                    end
                end
                box off
                xlim([-600,6000])
                title(HemOpt{hemid})
                ylim(ylims)
                line([0 0],ylims,'LineStyle','--','color', [0 0 0])
                line([500 500],ylims,'LineStyle','--','color', [0 0 0])
                line([2000 2000],ylims,'LineStyle','--','color', [0 0 0])
                line([4000 4000],ylims,'LineStyle','--','color', [0 0 0])
                ylabel(GenoTypeOpt{genid})

                % statistical test
                response = repmat({'Hit','Error'}',[1,2,length(miceopt)]);
                side = repmat({'Left','Right'}',[1,2,length(miceopt)]);
                side = permute(side,[2,1,3]);
                choice = side;
                choice(strcmp(response,'Error')&strcmp(side,'Right')) = {'Left'};
                choice(strcmp(response,'Error')&strcmp(side,'Left'))={'Right'};

                p = nan(length(timeline),3);
                parfor tp = 1:length(timeline)
                    tmp = tmp4ana(tp,:,:,:);
                    if any(~isnan(tmp(:)))
                        p(tp,:) = anovan(tmp(:),{choice(:),side(:)},'model','interaction','display','off');
                    end
                end

                %choice
                plot(timeline(find(p(:,1)<0.05)),repmat(ylims(2)*0.9,[1,sum(p(:,1)<0.05)]),'*','color',[0 0 0])

                %Stimulus side
                plot(timeline(find(p(:,2)<0.05)),repmat(ylims(2)*0.85,[1,sum(p(:,2)<0.05)]),'*','color',[0.8 0.8 0.8])

                %Interaction side
                plot(timeline(find(p(:,3)<0.05)),repmat(ylims(2)*0.80,[1,sum(p(:,3)<0.05)]),'*','color',[0.6 0.6 0.6])


            end
        end
        saveas(gcf,fullfile(LocalFolder,['AllMiceTC' OptoOpt{optid} [AREAS{areaid}{:}] '.fig']))
    end
end

%% Timecourses opto effect
ColOpt = [0 0 0; 0.5 0.5 0.5];
LineOpt = {'-','-- '};
ylims =[-0.035 0.052];
for optid=1 %:length(OptoOpt)
    for areaid=1:length(AREAS)
        figure('name',['AllMiceTC ' OptoOpt{optid} ' ' [AREAS{areaid}{:}]])
        for hemid=1:2
            for genid = 1:2
                subplot(2,2,(genid-1)*2+hemid)
                tmp4ana = nan(length(timeline),2,2,length(miceopt));
                for ridx=1:2
                    for sidx=1:2
                        tmp = squeeze(AllMiceTC(optid,ridx,sidx,areaid,hemid,:));
                        if any(cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0)))
                            tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))} = tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))}';
                        end
                        tmp = squeeze(cellfun(@(X) nanmean(X,2),tmp,'UniformOutput',0));
                        tmp = cat(2,tmp{:});
                        h(ridx,sidx) = shadedErrorBar(timeline,squeeze(nanmean(tmp,2)),squeeze(nanstd(tmp,[],2))./sqrt(length(miceopt)-1),...
                            {'color',ColOpt(ridx,:),'LineStyle',LineOpt{sidx},'LineWidth',2},1)
                        hold on
                        tmp4ana(:,ridx,sidx,:)=tmp;
                    end
                end
                box off
                xlim([-600,6000])
                title(HemOpt{hemid})
                ylim(ylims)
                line([0 0],ylims,'LineStyle','--','color', [0 0 0])
                line([500 500],ylims,'LineStyle','--','color', [0 0 0])
                line([2000 2000],ylims,'LineStyle','--','color', [0 0 0])
                line([4000 4000],ylims,'LineStyle','--','color', [0 0 0])
                ylabel(GenoTypeOpt{genid})

                % statistical test
                response = repmat({'Hit','Error'}',[1,2,length(miceopt)]);
                side = repmat({'Left','Right'}',[1,2,length(miceopt)]);
                side = permute(side,[2,1,3]);
                choice = side;
                choice(strcmp(response,'Error')&strcmp(side,'Right')) = {'Left'};
                choice(strcmp(response,'Error')&strcmp(side,'Left'))={'Right'};

                p = nan(length(timeline),3);
                parfor tp = 1:length(timeline)
                    tmp = tmp4ana(tp,:,:,:);
                    if any(~isnan(tmp(:)))
                        p(tp,:) = anovan(tmp(:),{choice(:),side(:)},'model','interaction','display','off');
                    end
                end

                %choice
                plot(timeline(find(p(:,1)<0.05)),repmat(ylims(2)*0.9,[1,sum(p(:,1)<0.05)]),'*','color',[0 0 0])

                %Stimulus side
                plot(timeline(find(p(:,2)<0.05)),repmat(ylims(2)*0.85,[1,sum(p(:,2)<0.05)]),'*','color',[0.8 0.8 0.8])

                %Interaction side
                plot(timeline(find(p(:,3)<0.05)),repmat(ylims(2)*0.80,[1,sum(p(:,3)<0.05)]),'*','color',[0.6 0.6 0.6])


            end
        end
        saveas(gcf,fullfile(LocalFolder,['AllMiceTC' OptoOpt{optid} [AREAS{areaid}{:}] '.fig']))
    end
end

%% DPrime Opto ON vs OFF
TWRewardBiasAna = {[2750 5000]};
TWRewardBiasNames = {'ITI'};
aligncell = load(fullfile(LocalFolder,'AlignToTEmplateWithAllenFrame2.mat'));
takeopto = [1,6];
DprimeStim = nan(40,40,2,length(miceopt));
DPrimeChoice = nan(40,40,2,length(miceopt));
for midx=1:length(miceopt)
    midx
    TMP = load(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
    % get model
    model = load(fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));
    TM = aligncell.AlignCell{find(ismember(aligncell.miceopt,miceopt{midx}))};
    % Take Hit/Error
    tmp = squeeze(TMP.AllDat(takeopto,ismember(ReactionOpt,{'Hit','Error'}),:));
    tmp = cellfun(@(X) squeeze(nanmean(X(:,:,timeline>=TWRewardBiasAna{1}(1)&timeline<=TWRewardBiasAna{1}(2),:),3)),tmp,'UniformOutput',0);
    % Dprime Stimulus - concatenate over response
    for sidx = 1:2
        dtmp = (nanmean(cat(3,tmp{2,:,sidx}),3)-nanmean(cat(3,tmp{1,:,sidx}),3))./...
            (0.5*sqrt(nanvar(cat(3,tmp{2,:,sidx}),[],3)+nanvar(cat(3,tmp{1,:,sidx}),[],3)));

        %Align to template mouse
        dtmp =imwarp(dtmp,TM,'OutputView',imref2d(size(dtmp)));
        dtmp(dtmp==0)=nan;
        DprimeStim(:,:,sidx,midx) = Nan_imresize(dtmp,0.1);
    end

    % Dprime Stimulus - concatenate over response
    for cidx = 1:2
        acidx = cidx+1;
        if acidx>2
            acidx=1;
        end
        tmpoff = cat(3,tmp{1,1,cidx},tmp{1,2,acidx});
        tmpon = cat(3,tmp{2,1,cidx},tmp{2,2,acidx});

        dtmp = (nanmean(tmpon,3)-nanmean(tmpoff,3))./...
            (0.5*sqrt(nanvar(tmpoff,[],3)+nanvar(tmpon,[],3)));

        %Align to template mouse
        dtmp =imwarp(dtmp,TM,'OutputView',imref2d(size(dtmp)));
        dtmp(dtmp==0)=nan;
        DPrimeChoice(:,:,cidx,midx) = Nan_imresize(dtmp,0.1);
    end

end




lims = [-.75 0.75];
figure;
for genidx = 1:length(GenoTypeOpt)
    for sidx = 1:2
        subplot(2,4,(genidx-1)*4+sidx)
        h=imagesc(imresize(nanmean(DprimeStim(:,:,sidx,ismember(GenOpt,GenoTypeOpt{genidx})),4),10),lims);
        set(h,'AlphaData', ~isnan(imresize(nanmean(DprimeStim(:,:,sidx,ismember(GenOpt,GenoTypeOpt{genidx})),4),10)))
        if sidx==1
            ylabel(GenoTypeOpt{genidx})
        end
        if genidx==1
            title(['Stimulus '  SideOpt{sidx}])
        end
        colormap redblue
        axis off
        hold on
        plot((aligncell.TemplateModel.AllX),(aligncell.TemplateModel.AllY),'k.','MarkerSize',1)

        subplot(2,4,(genidx-1)*4+sidx+2)
        h=imagesc(imresize(nanmean(DPrimeChoice(:,:,sidx,ismember(GenOpt,GenoTypeOpt{genidx})),4),10),lims);
        set(h,'AlphaData', ~isnan(imresize(nanmean(DPrimeChoice(:,:,sidx,ismember(GenOpt,GenoTypeOpt{genidx})),4),10)))

        if genidx==1
            title(['Response '  SideOpt{sidx}])
        end
        colormap redblue
        axis off
        hold on
        plot((aligncell.TemplateModel.AllX),(aligncell.TemplateModel.AllY),'k.','MarkerSize',1)
    end
end
saveas(gcf,fullfile(LocalFolder,['AllMice_DPrimeOptoOffON.fig']))

%%
ColOpt = [0 0 0; 1 0 0];
LineOpt = {'-','--'}
load(fullfile(LocalFolder,'AllMiceTC.mat'))

figure;
clear h
for midx = 1:length(miceopt)
    for hemid=1:2
        subplot(length(miceopt),2,(midx-1)*2+hemid);

        for ridx=1:2
            for sidx=1:2
                for optid = 1
                    tmp = squeeze(AllMiceTC(optid,ridx,sidx,1,hemid,midx));
                    if any(cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0)))
                        tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))} = tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))}';
                    end
                    tmp = squeeze(cellfun(@(X) nanmean(X,2),tmp,'UniformOutput',0));
                    tmp = cat(2,tmp{:});
                    h(ridx,sidx) = shadedErrorBar(timeline,squeeze(nanmean(tmp,2)),squeeze(nanstd(tmp,[],2))./sqrt(length(miceopt)-1),...
                        {'color',ColOpt(ridx,:),'LineStyle',LineOpt{sidx},'LineWidth',2},1)
                    hold on
                end
            end
        end
        xlim([-200 2000])

    end
end

figure;
clear h
for midx = 1:length(miceopt)
    for hemid=1:2
        subplot(length(miceopt),2,(midx-1)*2+hemid);

        for ridx=1:2
            for optid = 1

                for sidx=1:2
                    tmp = squeeze(AllMiceTC(optid,ridx,sidx,1,hemid,midx));
                    if any(cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0)))
                        tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))} = tmp{cell2mat(cellfun(@(X) size(X,1)==1,tmp,'UniformOutput',0))}';
                    end
                    tmp = squeeze(cellfun(@(X) nanmean(X,2),tmp,'UniformOutput',0));
                    tmps{sidx} = cat(2,tmp{:});
                end
                tmp = tmps{2} - tmps{1};
                h(ridx) = plot(timeline,squeeze(nanmean(tmp,2)),...
                    'color',ColOpt(ridx,:),'LineStyle',LineOpt{1},'LineWidth',2);
                hold on
                line([timeline(1) timeline(end)],[0 0],'color',[0.2 0.2 0.2])
            end
        end
        xlim([-500 2000])

    end
end

%% Decoder Per Pixel during delay and ITI for opto off and opto during response
part = 2; %1 = seperate learning, 2 = simultaneous
Int = 1; % Int = 2: integration of modalities, Int = 1: concatenation of modalities
takeopto = {[1],[2:4]};
TWRewardBiasAna = {[0 500],[1550 1900],[2000 2500]} %Late delay for WM Opto
TWRewardBiasNames = {'Visual','Late Delay','Response'}
sideopt = [-1 1];
newpixsz = 40;
% % Bootstrapping
ntestset = nan(length(TWRewardBiasAna),length(takeopto),length(miceopt));
ErrorLog = cell(1,length(miceopt));
currname = 'MalsarOutput_0505_V4_Balancedstim.mat'
for midx=1:length(miceopt)
    if exist(fullfile(LocalFolder,miceopt{midx},currname))
        tmp = load(fullfile(LocalFolder,miceopt{midx},currname),'PerfCrossPredPerSide','PerfWholeBrain','PerfCrossPred','PerfWholeBrainPerSide');
        % if any(any(tmp.PerfCrossPredPerSide(:,:,1,:)<0.5))
        %     disp([miceopt{midx} ' bad performance saved, redo..'])
        % else
        disp([miceopt{midx} ' already done, continue..'])
        continue
        % end
    end
    % Save whole brain performance
    % PerfWholeBrain = nan(2,length(TWRewardBiasAna),length(takeopto),nfolds,length(miceopt));
    try

        TMP = load(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'),'AllDat');
        % get model
        model = load(fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));

        %%
        BetaPerPix = nan(newpixsz*newpixsz,2,length(TWRewardBiasAna),length(takeopto),'single');
        Output = cell(length(TWRewardBiasAna),length(takeopto));
        PerfCrossPred = nan(2,length(TWRewardBiasAna),length(takeopto));
        PerfCrossPredPerSide = nan(2,length(TWRewardBiasAna),length(takeopto),length(SideOpt));
        PerfWholeBrain = nan(2,length(TWRewardBiasAna),length(takeopto));
        PerfWholeBrainPerSide = nan(2,2,length(TWRewardBiasAna),length(takeopto));

        % Sample minntr trials from each condition
        for twid=1:length(TWRewardBiasAna)
            for optid = 1:length(takeopto)
                maxntr = max(max(nansum(cell2mat(cellfun(@(X) size(X,4),TMP.AllDat(takeopto{optid},ismember(ReactionOpt,{'Hit','Error'}),:),'UniformOutput',0)),1)));

                % Prepare Data for MALSAR
                XDat = nan(newpixsz,newpixsz,maxntr,2,2,'single'); %10 times downsampled
                reactionvec = nan(maxntr,2,2);
                sidevec = nan(maxntr,2,2);
                for ridx = 1:2
                    for sidx = 1:2
                        % tmp = [];
                        tridx=1;
                        for optid2 = 1:numel(takeopto{optid})
                            if size(TMP.AllDat{takeopto{optid}(optid2),ridx,sidx},4)==0
                                continue
                            end
                            tmp =  imresize(nanmean(TMP.AllDat{takeopto{optid}(optid2),ridx,sidx}(:,:,timeline>=TWRewardBiasAna{twid}(1)&timeline<=TWRewardBiasAna{twid}(2),:),3),newpixsz./400);
                            XDat(:,:,tridx:size(tmp,4)+tridx-1,ridx,sidx) = tmp;
                            tridx = tridx+size(tmp,4);
                        end
                        reactionvec(:,ridx,sidx) = ridx;
                        sidevec(:,ridx,sidx)=sidx;
                    end
                end
                clear tmp
                XDat = reshape(XDat,newpixsz*newpixsz,[]);
                % Remove ridiculous outliers
                if optid==1
                    meanmat = nanmean(XDat,2);
                    stdmat = nanstd(XDat,[],2);
                end
                XDat = (XDat - meanmat)./stdmat; % First pass, treat opto on same


                % XDat = (XDat - nanmean(XDat,1))./nanstd(XDat,[],1);
                % Establish which pixels should run the analysis
                pixelidx = find(~isnan(nanmean(reshape(XDat,newpixsz*newpixsz,[]),2)));
                Nanpixelidx = find(isnan(nanmean(reshape(XDat,newpixsz*newpixsz,[]),2)));
                actualtrialidx = find(~isnan(nanmean(reshape(XDat,newpixsz*newpixsz,[]),1)));

                % Fill missing in 2D
                for trid = 1:numel(actualtrialidx)
                    if any(isnan(XDat(pixelidx,actualtrialidx(trid))))
                        tmp = reshape(XDat(:,actualtrialidx(trid)),newpixsz,newpixsz);
                        % tmp = fillmissing2(tmp,'cubic');
                        tmp = fillmissing2(tmp,"movmedian",15);
                        tmp = reshape(tmp,newpixsz*newpixsz,1);
                        tmp(Nanpixelidx) = nan;
                        XDat(:,actualtrialidx(trid)) = tmp;
                    end
                end

                XDat = double(XDat(pixelidx,actualtrialidx));

                if any(isnan(XDat(:)))
                    ErrorLog{midx} = 'nan pixls';

                end
                % Also normalize over pixels
                reactionvec = reactionvec(actualtrialidx);
                reactionvec(reactionvec==2)=-1; %to make error -1
                sidevec = sidevec(actualtrialidx);
                sidevec(sidevec==2)=-1; %To make right grating -1
                choicevec = sidevec;
                choicevec(reactionvec==-1&sidevec==-1) = 1;
                choicevec(reactionvec==-1&sidevec==1)=-1;
                choiceopt = unique(choicevec);
                sideopt = unique(sidevec);

                % Make sense of life
                figure('name',[miceopt{midx} ', tw=' num2str(twid) ', opto=' num2str(optid)])
                for ridx = 1:2
                    for sidx = 1:2
                        tmp1 = nan(newpixsz,newpixsz);
                        tmp1(pixelidx) = nanmean(XDat(:,choicevec==choiceopt(ridx)&sidevec==sideopt(sidx)),2);

                        subplot(2,2,(ridx-1)*2+sidx)
                        h = imagesc(tmp1,[-.5 .5]); colormap redblue
                        if ridx ==1
                            title(['Stimulus ' num2str(sideopt(sidx))])
                        end
                        if sidx==1
                            ylabel(['Choice ' num2str(choiceopt(ridx))])
                        end
                    end
                end

                T = combvec(sideopt,choiceopt);
                % Keep out data
                if optid==1
                    minntr = min(arrayfun(@(X) sum(ismember(sidevec,T(1,X))&ismember(choicevec,T(2,X))),1:size(T,2)));

                    % need a balanced training set
                    trainsetcross = [];
                    for ridx = 1:2 % Balanced
                        for sidx = 1:2
                            trainsetcross = [trainsetcross;datasample(find(choicevec==choiceopt(ridx) & sidevec==sideopt(sidx)),minntr,'replace',false)];
                        end
                    end
                    % Hold out testset
                    testsetcross = trainsetcross(:,1:floor(0.2*size(trainsetcross,2)));
                    testsetcross = testsetcross(:);
                    trainsetcross = trainsetcross(:);
                    trainsetcross(ismember(trainsetcross,testsetcross)) = [];

                    % add extra test set trials so long as it's not biased
                    % for stimulus side (already presented, not yet
                    % decided)
                    UnUsedIdx = 1:numel(choicevec);
                    UnUsedIdx(ismember(UnUsedIdx,testsetcross)|ismember(UnUsedIdx,trainsetcross)) = [];
                    nExtra = min([sum(sidevec(UnUsedIdx)==-1) sum(sidevec(UnUsedIdx)==1)]);
                    for sidx = 1:2
                        testsetcross = [testsetcross; UnUsedIdx(datasample(find(sidevec(UnUsedIdx)==sideopt(sidx)),nExtra,'replace',false))'];
                    end
                    % testsetcross = [testsetcross; UnUsedIdx'];


                    % Define parameters with 1/2 of training data
                    ncv = floor(numel(trainsetcross)/4/2)*4;
                    paramMTL_max = DefineBestParam(double(XDat(:,trainsetcross))',cat(1,sidevec(trainsetcross),choicevec(trainsetcross)),ncv,part,Int,4,1:ncv);

                    % Whole brain performance
                    % Calculate average performance (cross-validated) across pixels
                    Output{twid,optid} = MALSAR_MTL_WF_Fast(XDat(:,trainsetcross)',cat(1,sidevec(trainsetcross),choicevec(trainsetcross)),length(trainsetcross),part,Int,4,paramMTL_max,1,1,1:numel(trainsetcross));
                    PerfWholeBrain(:,twid,optid) = nansum(Output{twid,optid}.PredictionsVsActual(:,:,3)==sign(Output{twid,optid}.PredictionsVsActual(:,:,1)),2)./sum(~isnan(Output{twid,optid}.PredictionsVsActual(:,:,3)),2);
                    for sidx = 1:2
                        id2 = Output{twid,optid}.PredictionsVsActual(1,:,3) == sideopt(sidx);
                        PerfWholeBrainPerSide(:,sidx,twid,optid) = nansum(Output{twid,optid}.PredictionsVsActual(:,id2,3)==sign(Output{twid,optid}.PredictionsVsActual(:,id2,1)),2)./sum(~isnan(Output{twid,optid}.PredictionsVsActual(:,id2,3)),2);
                    end

                    tmpbeta = squeeze(nanmean(Output{twid,optid}.Weights.Modalities,1))';
                    BetaPerPix(pixelidx,:,twid,optid) = tmpbeta(2:end,:);
                    figure('name',['Beta ' miceopt{midx} ', tw=' num2str(twid)])
                    for taskid = 1:2
                        subplot(1,2,taskid)
                        h = imagesc(reshape(BetaPerPix(:,taskid,twid,optid),newpixsz,newpixsz),[-0.1 0.1]); colormap redblue
                        if taskid ==1
                            title(['Stimulus decoding'])
                        else
                            title(['Choice decoding'])
                        end
                    end
                    drawnow

                else

                    % % % Hold out test set

                    minntr = nanmin(arrayfun(@(X) sum(ismember(sidevec,T(1,X))),1:size(T,2)));

                    testsetcross = [];
                    for sidx = 1:2 % Balanced stimulus side
                        testsetcross = [testsetcross;datasample(find(sidevec==sideopt(sidx)),minntr,'replace',false)];
                    end
                    testsetcross = testsetcross(:);
                end
                ntestset(twid,optid,midx) = length(testsetcross);

                tmpbeta = squeeze(nanmean(Output{twid,1}.Weights.Modalities,1))';

                % Also cross-test; test the model of opto off on opto on
                tmp2 = XDat(:,testsetcross)';
                tmp2 = cat(2,ones(size(tmp2,1),1),tmp2);
                % Add row of ones
                PerfCrossPred(:,twid,optid) =  nansum(sign(tmp2 * tmpbeta) == cat(1,sidevec(testsetcross),choicevec(testsetcross))',1)./size(tmp2,1);
                for sidx=1:2
                    id2 = sidevec(testsetcross)==sideopt(sidx);
                    PerfCrossPredPerSide(:,twid,optid,sidx) = nansum(sign(tmp2(id2,:)* tmpbeta)== cat(1,sidevec(testsetcross(id2)),choicevec(testsetcross(id2)))',1)./sum((id2));
                end

                if any(any(PerfCrossPredPerSide(:,twid,optid,:)<0.5)) & optid==1
                    ErrorLog{midx} = 'BAD PERFORMANCE';
                end


                disp(['Whole brain performance ' miceopt{midx} ' ' TWRewardBiasNames{twid} ' ' OptoOpt{takeopto{optid}}])
                disp(['p_side normal performance = '  num2str(PerfWholeBrain(1,twid,optid)) ', cross-performance = ' num2str( PerfCrossPred(1,twid,optid))])
                disp(['p_choice normal performance = '  num2str(PerfWholeBrain(2,twid,optid)) 'cross-performance = ' num2str(PerfCrossPred(2,twid,optid))])
            end
            save(fullfile(LocalFolder,miceopt{midx},currname),'PerfCrossPred','PerfCrossPredPerSide','PerfWholeBrain','PerfWholeBrainPerSide','BetaPerPix','Output')
        end

        %%
        clear TMP BetaPerPix Output PerfCrossPred PerfCrossPredPerSide PerfWholeBrain XDat tmp tmp2
    catch ME
        ErrorLog{midx} = ME;
        clear TMP BetaPerPix Output PerfCrossPred PerfCrossPredPerSide PerfWholeBrain XDat tmp tmp2

    end
end

%% % Load for all mice
PerfCrossPredPerSideAll = nan(2,length(TWRewardBiasAna),length(takeopto),length(SideOpt),length(miceopt));
PerfWholeBrainAll = nan(2,length(TWRewardBiasAna),length(takeopto),length(miceopt));

for midx = 1:length(miceopt)
    tmp = load(fullfile(LocalFolder,miceopt{midx},currname),'PerfCrossPredPerSide','PerfWholeBrain','PerfCrossPred','PerfWholeBrainPerSide');

    PerfCrossPredPerSideAll(:,:,:,:,midx) = tmp.PerfCrossPredPerSide;
    PerfWholeBrainAll(:,:,:,midx) = tmp.PerfCrossPred;
end


for twid=1:length(TWRewardBiasNames)
    TWRewardBiasNames{twid}
    y1 = squeeze(PerfCrossPredPerSideAll(1,twid,:,:,:));
    y2 = squeeze(PerfCrossPredPerSideAll(2,twid,:,:,:));
    % Factorial MANOVA
    g3 = repmat({'Off','On'}',[1,length(SideOpt),length(miceopt)]);
    g4 = repmat(GenOpt',[1,2,length(SideOpt)]);
    g4 = permute(g4,[2,3,1]);
    g5 = repmat(miceopt',[1,2,length(SideOpt)]);
    g5 = permute(g5,[2,3,1]);
    g2 = repmat(SideOpt',[1,2,length(miceopt)]);
    g2 = permute(g2,[2,1,3]);

    t = table(y1(:),y2(:),g2(:),g3(:),g4(:),g5(:),'VariableNames',{'AccuracyVisual','AccuracyResponse','Side','Opto','GenoType','MouseID'});

    fmfit = fitrm(t,'AccuracyVisual-AccuracyResponse ~ Opto * GenoType * Side');
    manovatbl = anova(fmfit)


    cols = [0 0 0; 1 0 0; 0.5 0.5 0.5; 0 0 1];
    %Whole brain decoding
    for taskid=1:2
        figure('name',['DecodingPerformance ' TWRewardBiasNames{twid} tasknames{taskid}])

        % Stats
        g1 = repmat({'Off','On'}',[1,2,length(GenOpt)]);
        g2 = repmat(GenOpt',[1,2,2]);
        g2 = permute(g2,[2,3,1]);
        g3 = repmat(fliplr(SideOpt)',[1,2,length(GenOpt)]);
        g3 = permute(g3,[2,1,3]);

        y1= squeeze(PerfCrossPredPerSideAll(taskid,twid,:,:,:)); % opto X stimside X mouse

        t = table(y1(:),g1(:),g2(:),g3(:),'VariableNames',{'Accuracy','Opto','GenoType','Side'});

        fmfit = fitrm(t,'Accuracy ~ Opto * GenoType * Side');
        anovatbl = anova(fmfit)

        for sidx = 1:2
            subplot(1,2,sidx)
            toplot = nan(2,length(GenoTypeOpt),2);
            for gidx = 1:length(GenoTypeOpt)
                toplot(:,gidx,1) = nanmean(squeeze(PerfCrossPredPerSideAll(taskid,twid,:,sidx,ismember(GenOpt,GenoTypeOpt{gidx}))),2);
                toplot(:,gidx,2) = nanstd(squeeze(PerfCrossPredPerSideAll(taskid,twid,:,sidx,ismember(GenOpt,GenoTypeOpt{gidx}))),[],2)./sqrt(sum(ismember(GenOpt,GenoTypeOpt{gidx}))-1);
            end
            h=barwitherr(reshape(toplot(:,:,2),2,[])',reshape(toplot(:,:,1),2,[])');
            for hid = 1:length(h)
                h(hid).FaceColor = cols(hid,:);
                h(hid).EdgeColor = 'none';
                h(hid).BarWidth=1;
            end
            hold on

            box off
            hold on
            line([0.5 size(PerfWholeBrainAll,3)+0.5],[0.5 0.5],'LineStyle','--','color',[1 0 0],'LineWidth',2)
            ylim([0 1])

            set(gca,'XTickLabel',GenoTypeOpt)
            if sideopt(sidx) == -1 % Different for output decoder
                title('right')
            else
                title('left')
            end

            xlabel('Genotype')

            % Individual mice
            for gidx=1:2
                tmp2 = squeeze(PerfCrossPredPerSideAll(taskid,twid,:,sidx,ismember(GenOpt,GenoTypeOpt{gidx})));
                plot([gidx-0.15 gidx+0.15],tmp2,'color', cols(gidx*2,:),'LineWidth',2)

            end

        end
        legend({'Off','Opto'})

        saveas(gcf,fullfile(LocalFolder,[TWRewardBiasNames{twid} tasknames{taskid} 'DecodingPerformance.fig']))
    end
end

%% Plot - Figure 5D RewardBias
load(fullfile(LocalFolder,'MalsarOutput_Reward1_v5.mat'))
g1 = repmat(tasknames',[1,length(TWRewardBiasNames),2,length(miceopt)]);
g2 = repmat(TWRewardBiasNames',[1,2,2,length(miceopt)]);
g2 = permute(g2,[2,1,3,4]);
g3 = repmat({'Off','On'}',[1,length(TWRewardBiasNames),2,length(miceopt)]);
g3 = permute(g3,[2,3,1,4]);
g4 = repmat(GenOpt',[1,length(TWRewardBiasNames),2,2]);
g4 = permute(g4,[2,3,4,1]);
g5 = repmat(miceopt',[1,length(TWRewardBiasNames),2,2]);
g5 = permute(g5,[2,3,4,1]);

% idx = true(length(g1(:)),1);
% idx = ismember(g1(:),tasknames{1});
idx = ismember(g2,'Delay')
tbl = table(PerfCrossPred(idx),g1(idx),g2(idx),g3(idx),g4(idx),g5(idx),'VariableNames', {'Performance','Task','Time','Opto','Genotype','MouseID'});
mld = fitglme(tbl,'Performance ~ Task*Opto*Genotype + (1|MouseID)')
anova(mld)


% Factorial MANOVA
g3 = repmat({'Off','On'}',[1,length(miceopt)]);
g4 = repmat(GenOpt',[1,2]);
g4 = permute(g4,[2,1]);
g5 = repmat(miceopt',[1,2]);
g5 = permute(g5,[2,1]);
y1 = PerfCrossPred(1,length(TWRewardBiasNames),:,:);
y2 = PerfCrossPred(2,length(TWRewardBiasNames),:,:);

t = table(y1(:),y2(:),g3(:),g4(:),g5(:),'VariableNames',{'AccuracyVisual','AccuracyResponse','Opto','GenoType','MouseID'});

fmfit = fitrm(t,'AccuracyResponse ~ Opto * GenoType');
manovatbl = anova(fmfit)

% Divide per genotype
for taskid = 1:2

    figure('name',['DecodingPerformance Per Mouse ' tasknames{taskid}])
    for midx = 1:length(miceopt)
        tmpboot = reshape(PerfCrossPredBoot(taskid,:,:,midx,:),[length(TWRewardBiasNames),2,nboot]); %TW, Opto, boot
        tmpori = reshape(PerfCrossPred(taskid,:,:,midx),[length(TWRewardBiasNames),2]);
        for twid=1:length(TWRewardBiasNames)
            subplot(length(miceopt),length(TWRewardBiasNames),(midx-1)*length(TWRewardBiasNames)+twid)
            hold on
            for optid = 1:2
                pval(twid,optid) = 1-invprctile(squeeze(tmpboot(twid,optid,:)),tmpori(twid,optid))/100;

                hh = histogram(tmpboot(twid,optid,:),20);
                hh.FaceColor = cols(optid+(find(ismember(GenoTypeOpt,GenOpt{midx}))-1)*2,:);
                hh.EdgeColor = 'none';
                hh.FaceAlpha = 0.1;
                hold on
                line([tmpori(twid,optid) tmpori(twid,optid)],get(gca,'ylim'),'linestyle','-','Color',cols(optid+(find(ismember(GenoTypeOpt,GenOpt{midx}))-1)*2,:),'LineWidth',2)
                xlim([0 1])

                if pval(twid,optid)<0.001
                    text(tmpori(twid,optid),max(get(gca,'ylim'))*0.9,'***')
                elseif pval(twid,optid)<0.01
                    text(tmpori(twid,optid),max(get(gca,'ylim'))*0.9,'**')
                elseif pval(twid,optid)<0.05
                    text(tmpori(twid,optid),max(get(gca,'ylim'))*0.9,'*')
                end

                newboot = tmpboot(twid,2,:)-tmpboot(twid,1,:);
                newori = tmpori(twid,2)-tmpori(twid,1)
                pvaldiff = 1-invprctile(squeeze(newboot),newori)/100;

                title([miceopt{midx} ' ' TWRewardBiasNames{twid} ', p = ' num2str(pvaldiff)])
                box off

            end
        end
    end
end


for twid=1:length(TWRewardBiasNames)
    figure('name',['DecodingPerformance ' TWRewardBiasNames{twid}])
    cols = [0 0 0; 1 0 0; 0.5 0.5 0.5; 0 0 1];
    %Whole brain decoding
    for taskid=1:2

        subplot(1,2,taskid)
        toplot = nan(2,length(GenoTypeOpt),2);
        for gidx = 1:length(GenoTypeOpt)
            toplot(:,gidx,1) = nanmean(squeeze(PerfCrossPred(taskid,twid,:,ismember(GenOpt,GenoTypeOpt{gidx}))),2);
            toplot(:,gidx,2) = nanstd(squeeze(PerfCrossPred(taskid,twid,:,ismember(GenOpt,GenoTypeOpt{gidx}))),[],2)./sqrt(sum(ismember(GenOpt,GenoTypeOpt{gidx}))-1);
        end
        h=barwitherr(reshape(toplot(:,:,2),2,[])',reshape(toplot(:,:,1),2,[])');
        for hid = 1:length(h)
            h(hid).FaceColor = cols(hid,:);
            h(hid).EdgeColor = 'none';
            h(hid).BarWidth=1;
        end
        hold on

        box off
        hold on
        line([0.5 size(PerfWholeBrain,3)+0.5],[0.5 0.5],'LineStyle','--','color',[1 0 0],'LineWidth',2)
        ylim([0 1])

        set(gca,'XTickLabel',GenoTypeOpt)
        title(tasknames{taskid})

        xlabel('Genotype')

        % Individual mice
        for gidx=1:2
            tmp2 = squeeze(PerfCrossPred(taskid,twid,:,ismember(GenOpt,GenoTypeOpt{gidx})));
            plot([gidx-0.15 gidx+0.15],tmp2,'color', cols(gidx*2,:),'LineWidth',2)

        end

    end
    legend({'Off','Opto'})

    saveas(gcf,fullfile(LocalFolder,[TWRewardBiasNames{twid} 'DecodingPerformance.fig']))
end



%% weightmaps
aligncell = load(fullfile(LocalFolder,'AlignToTEmplateWithAllenFrame2.mat'));

for midx = 1:length(miceopt)
    for twid = 1:length(TWRewardBiasNames)
        for optid = 1:2
            for taskid=1:2
                tmp = imresize(reshape(BetaPerPix(:,taskid,twid,optid,midx),40,40),10);

                TM = aligncell.AlignCell{find(ismember(aligncell.miceopt,miceopt{midx}))};
                tmp =imwarp(tmp,TM,'OutputView',imref2d(size(tmp)));
                tmp(tmp==0)=nan;
                BetaPerPix(:,taskid,twid,optid,midx) = reshape(Nan_imresize(tmp,0.1),40*40,[]);
            end
        end
    end

end

tmp = nanmean(BetaPerPix,5);
lims = [-quantile(abs(tmp(:)),0.75) quantile(abs(tmp(:)),0.75)];
figure;
twid=length(TWRewardBiasNames)
optid=1;
subplot(2,1,1)
tmp = imresize(reshape(nanmean(BetaPerPix(:,1,twid,optid,:),5),40,40),10);
h=imagesc(tmp,lims);
set(h,'AlphaData', ~isnan(tmp))

title(['Stimulus '  OptoOpt{takeopto(optid)}])
colormap redblue
axis off
hold on
plot((aligncell.TemplateModel.AllX),(aligncell.TemplateModel.AllY),'k.','MarkerSize',1)

subplot(2,1,2)
tmp = imresize(reshape(nanmean(BetaPerPix(:,2,twid,optid,:),5),40,40),10);

h=imagesc(tmp,lims);
set(h,'AlphaData', ~isnan(tmp))

title(['Response '   OptoOpt{takeopto(optid)}])
colormap redblue
axis off
hold on
plot((aligncell.TemplateModel.AllX),(aligncell.TemplateModel.AllY),'k.','MarkerSize',1)

saveas(gcf,fullfile(LocalFolder,['AllMice_Weightmap.fig']))

%% Decoder Per Pixel during delay and ITI
nfolds = 5;
part = 2; %1 = seperate learning, 2 = simultaneous
Int = 1; % Int = 2: integration of modalities, Int = 1: concatenation of modalities
pixgrid = reshape(1:40*40,40,40);

% Save whole brain performance
PerfWholeBrain = nan(2,length(TW),length(miceopt));
BetaPerPix = nan(40*40,2,length(TW),length(miceopt),'single');

PerfWholeBrainBoot = nan(2,length(TW),nboot,length(miceopt));
BetaPerPixBoot = nan(40*40,2,nboot,length(TW),length(miceopt));
for midx=1:length(miceopt)
    if exist(fullfile(savefolder,['MVPA_Output' miceopt{midx} '.mat']))
        tmpresults = load(fullfile(savefolder,['MVPA_Output' miceopt{midx} '.mat']));
        PerfWholeBrainBoot(:,:,:,midx) = tmpresults.PerfWholeBrainBoot(:,:,:,midx);
    end

    % Get timeline
    try
        tmp = matfile(fullfile(LocalFolder{folderid(midx)},miceopt{midx},'TimeCourses.mat'));
        timeline = tmp.timeline;
    catch
        tmp = dir(fullfile(LocalFolder{folderid(midx)},miceopt{midx},'**','*RawData_All.mat'))
        tmp = matfile(fullfile(tmp(1).folder,tmp(1).name));
        timeline = tmp.timeline;
    end

    disp('Loading data')
    TMP = matfile(fullfile(LocalFolder{folderid(midx)},miceopt{midx},'ProcessedData.mat'));
    % Take opto OFF data, and Reaction Hit/Error
    tmp = squeeze(TMP.AllDat(find(strcmp(OptoOpt,'Opto Off')),find(ismember(ReactionOpt,{'Hit','Error'})),find(ismember(SideOpt,{'left','right'}))));
    minntr = min(min(cell2mat(cellfun(@(X) size(X,4),tmp,'UniformOutput',0))));

    % Prepare Data for MALSAR
    XDat = nan(40,40,minntr,2,2,length(TW),'single'); %10 times downsampled
    trialidx = nan(minntr,2,2);
    reactionvec = nan(minntr,2,2);
    sidevec = nan(minntr,2,2);
    for twid = 1:length(TW)
        for ridx = 1:2
            for sidx = 1:2
                trialidx(:,ridx,sidx) = randsample(size(tmp{ridx,sidx},4),minntr,false);
                XDat(:,:,:,ridx,sidx,twid) = Nan_imresize(nanmean(tmp{ridx,sidx}(:,:,timeline>=TW{twid}(1)&timeline<=TW{twid}(2),trialidx(:,ridx,sidx)),3),0.1);
                reactionvec(:,ridx,sidx) = ridx;
                sidevec(:,ridx,sidx)=sidx;
            end
        end
    end
    XDat = reshape(XDat,40*40,[],length(TW));
    reactionvec = reactionvec(:);
    reactionvec(reactionvec==2)=-1; %to make error -1
    sidevec = sidevec(:);
    sidevec(sidevec==2)=-1; %To make right grating -1
    choicevec = sidevec;
    choicevec(reactionvec==-1&sidevec==-1) = 1;
    choicevec(reactionvec==-1&sidevec==1)=-1;

    % Establish which pixels should run the analysis
    pixelidx = find(~isnan(nanmean(reshape(XDat,40*40,[]),2)));

    % Sample minntr trials from each condition
    parfor twid=1:length(TW)
        if any(isnan(PerfWholeBrain(:,twid,midx)))
            % Define parameters
            disp('Fitting best parameters')
            tp = tic;
            paramMTL_max = DefineBestParam(nanmean(XDat(pixelidx,:,twid),3)',cat(2,sidevec,choicevec)',length(choicevec),part,Int,nfolds);
            disp(['Took '  num2str(round(toc(tp)./60)) ' minutes'])
            % Whole brain performance
            % Calculate average performance (cross-validated) across pixels
            Outputboot = MALSAR_MTL_WF_Fast(XDat(pixelidx,:,twid)',cat(2,sidevec,choicevec)',length(choicevec),part,Int,nfolds,paramMTL_max,1,1);
            ModelOutputs{twid,midx} = Outputboot;
            PerfWholeBrain(:,twid,midx) = nansum(Outputboot.PredictionsVsActual(:,:,3)==sign(Outputboot.PredictionsVsActual(:,:,1)),2)./sum(~isnan(Outputboot.PredictionsVsActual(:,:,3)),2);
            tmpweights = squeeze(nanmean(Outputboot.Weights.Modalities,1))';
            BetaPerPix(pixelidx,:,twid,midx) = tmpweights(2:end,:);
            for bootid=1:nboot
                tmpy=cat(2,sidevec,choicevec)';
                tmpy = tmpy(:,randsample(size(tmpy,2),size(tmpy,2),true));
                Outputboot = MALSAR_MTL_WF_Fast(XDat(pixelidx,:)',tmpy,length(choicevec),part,Int,nfolds,paramMTL_max,1,1);
                PerfWholeBrainBoot(:,twid,bootid,midx) = nansum(Outputboot.PredictionsVsActual(:,:,3)==sign(Outputboot.PredictionsVsActual(:,:,1)),2)./sum(~isnan(Outputboot.PredictionsVsActual(:,:,3)),2);
                tmpweights = squeeze(nanmean(Outputboot.Weights.Modalities,1))';
                BetaPerPixBoot(pixelidx,:,bootid,twid,midx) = tmpweights(2:end,:);
            end
        end
    end

    save(fullfile(savefolder,['MVPA_Output2' miceopt{midx} '.mat']),'PerfWholeBrain','BetaPerPix','PerfWholeBrainBoot','BetaPerPixBoot','ModelOutputs','TW')

    figure('name',miceopt{midx})
    for twid=1:length(TW)
        for tid = 1:size(PerfWholeBrain,1)
            pvals(tid) = invprctile(squeeze(PerfWholeBrainBoot(tid,twid,:,midx)),PerfWholeBrain(tid,twid,midx));
            pvals(tid) = (100-pvals(tid))/100;
        end
        disp(['Whole brain performance ' miceopt{midx} ' ' TimeWindowNames{twid}])
        disp(['performance_side = ' num2str(PerfWholeBrain(1,twid,midx)) ', p=' num2str(pvals(1))])
        disp(['performance_choice = ' num2str(PerfWholeBrain(2,twid,midx)) ', p=' num2str(pvals(2))])

        tmpthismouse = BetaPerPix(:,:,twid,midx);
        threshold = squeeze(abs(BetaPerPixBoot(:,:,:,twid,midx)));
        threshold = quantile(threshold,0.975,3);
        lims = quantile(abs(tmpthismouse(~isnan(tmpthismouse)&tmpthismouse~=0)),0.95);

        subplot(length(TW),2,(twid-1)*2+1)
        try
            h=imagesc(imresize(reshape(BetaPerPix(:,1,twid,midx),40,40),10),[-lims lims]);
        catch
            h=imagesc(imresize(reshape(BetaPerPix(:,1,twid,midx),40,40),10));
        end
        colormap redblue
        binimg = zeros(400,400);
        binimg(~isnan(imresize(reshape(BetaPerPix(:,1,twid,midx),40,40),10)))=0.5;
        binimg(imresize(reshape(BetaPerPix(:,1,twid,midx)>threshold(:,1),40,40),10)) = 1;
        binimg(imresize(reshape(BetaPerPix(:,1,twid,midx)<-threshold(:,1),40,40),10)) = 1;
        binimg(isnan(imresize(reshape(BetaPerPix(:,1,twid,midx),40,40),10)))=0;
        set(h,'AlphaData',binimg)
        hold on
        plot((model.BrainModel.AllX),(model.BrainModel.AllY),'k.','MarkerSize',2)
        title([TimeWindowNames{twid} '  Grating Left/Right'])
        box off
        axis off
        axis square


        subplot(length(TW),2,(twid-1)*2+2)
        try
            h=imagesc(imresize(reshape(BetaPerPix(:,2,twid,midx),40,40),10),[-lims lims]);
        catch
            h=imagesc(imresize(reshape(BetaPerPix(:,2,twid,midx),40,40),10),[-lims lims]);
        end
        colormap redblue
        binimg = zeros(400,400);
        binimg(~isnan(imresize(reshape(BetaPerPix(:,2,twid,midx),40,40),10)))=0.5;
        binimg(imresize(reshape(BetaPerPix(:,2,twid,midx)>threshold(:,2),40,40),10)) = 1;
        binimg(imresize(reshape(BetaPerPix(:,2,twid,midx)<-threshold(:,2),40,40),10)) = 1;
        binimg(isnan(imresize(reshape(BetaPerPix(:,2,twid,midx),40,40),10)))=0;
        set(h,'AlphaData',binimg)
        hold on
        plot((model.BrainModel.AllX),(model.BrainModel.AllY),'k.','MarkerSize',2)
        title([TimeWindowNames{twid} '  Choice left/right'])
        box off
        axis off
        axis square
    end
    saveas(gcf,fullfile(savefolder,['MVPA_Output ' miceopt{midx} '.fig']))
end



%% Across mice performance
figure;
for tid =1:2
    subplot(1,2,tid)
    tmp = squeeze(PerfWholeBrain(tid,:,:));
    h=barwitherr(nanstd(tmp,[],2)./sqrt(length(miceopt)-1),nanmean(tmp,2));
    h.FaceColor=[0.5 0.5 0.5];
    h.BarWidth=1;
    hold on
    line([0 length(TimeWindowNames)+1],[0.5 0.5],'LineStyle','--','color',[1 0 0])
    g1 = repmat(TimeWindowNames',[1,length(miceopt)]);
    anovan(tmp(:),{g1(:)})
    ylim([0 1])
    box off
    set(gca,'XTickLabel',TimeWindowNames)

    title(tasknames{tid})

end

%% Performance for hit and errors separately
PerfWholeBrainHE = nan(length(tasknames),length(TW),2,length(miceopt));
for midx = 1:length(miceopt)
    for twid=1:length(TW)
        correctidx = find(ModelOutputs{twid,midx}.PredictionsVsActual(1,:,3)==ModelOutputs{twid,midx}.PredictionsVsActual(2,:,3)); %visual matches response side
        PerfWholeBrainHE(:,twid,1,midx) = nansum(ModelOutputs{twid,midx}.PredictionsVsActual(:,correctidx,3)==sign(ModelOutputs{twid,midx}.PredictionsVsActual(:,correctidx,1)),2)./sum(~isnan(ModelOutputs{twid,midx}.PredictionsVsActual(:,correctidx,3)),2);
        erroridx = find(ModelOutputs{twid,midx}.PredictionsVsActual(1,:,3)~=ModelOutputs{twid,midx}.PredictionsVsActual(2,:,3)); %visual doesn't match response side
        PerfWholeBrainHE(:,twid,2,midx) = nansum(ModelOutputs{twid,midx}.PredictionsVsActual(:,erroridx,3)==sign(ModelOutputs{twid,midx}.PredictionsVsActual(:,erroridx,1)),2)./sum(~isnan(ModelOutputs{twid,midx}.PredictionsVsActual(:,erroridx,3)),2);

    end
end

figure;
pall = nan(2,length(TW),2);
for tid =1:2
    subplot(1,2,tid)
    tmp = squeeze(PerfWholeBrainHE(tid,:,:,:));
    h=barwitherr(nanstd(tmp,[],3)./sqrt(length(miceopt)-1),nanmean(tmp,3));
    h(1).FaceColor=[0 0 0];
    h(2).FaceColor=[0.5 0.5 0.5];
    h(1).BarWidth=1;
    h(2).BarWidth=1;

    hold on
    line([0 length(TimeWindowNames)+1],[0.5 0.5],'LineStyle','--','color',[1 0 0])
    g1 = repmat(TimeWindowNames',[1,2,length(miceopt)]);
    g2 = repmat({'Correct','Error'}',[1,length(TW),length(miceopt)]);
    g2 = permute(g2,[2,1,3]);
    anovan(tmp(:),{g1(:),g2(:)},'model','interaction','varnames',{'TimeWindow','Correctness'})
    ylim([0 1])
    box off
    set(gca,'XTickLabel',TimeWindowNames)

    for twid=1:length(TW)
        for cid = 1:2
            [~,p] = ttest(squeeze(PerfWholeBrainHE(tid,twid,cid,:)),0.5,'tail','right');
            pall(tid,twid,cid)=p
            if p<0.001
                text(  h(cid).XData(twid)-0.1+0.2*(cid-1),0.95,'***')
            elseif p<0.01
                text(  h(cid).XData(twid)-0.1+0.2*(cid-1),0.95,'**')
            elseif p<0.05
                text(  h(cid).XData(twid)-0.1+0.2*(cid-1),0.95,'*')
            end
        end
    end

    title(tasknames{tid})

end

%% Factorial MANOVA
g2 = repmat(TimeWindowNames',[1,2,length(miceopt)]);
g3 = repmat({'Correct','Error'}',[1,length(TW),length(miceopt)]);
g3 = permute(g3,[2,1,3]);
y1 = PerfWholeBrainHE(1,:,:,:);
y2 = PerfWholeBrainHE(2,:,:,:);
t = table(y1(:),y2(:),g2(:),g3(:),'VariableNames',{'AccuracyVisual','AccuracyResponse','TimeWindow','Correctness'});

fmfit = fitrm(t,'AccuracyVisual-AccuracyResponse ~ TimeWindow * Correctness');
manovatbl = manova(fmfit)


anovan(PerfWholeBrainHE(:),{g1(:),g2(:),g3(:)},'model',3,'varnames',{'Task','TimeWindow','Correctness'})


%% Figure 5D for WM paper
% Save whole brain performance
PerfWholeBrain = nan(3,2,length(miceopt)); %Grating, Choice, Response X OptoOff,optoOn
BetaPerPix = nan(40*40,3,length(miceopt),'single');

for midx=1:length(miceopt)
    figure('name',miceopt{midx})

    TMP = load(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
    % get model
    model = load(fullfile(LocalFolder,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));


    % Sample minntr trials from each condition
    % Take opto OFF data, and Reaction Hit/Error
    tmpoff = squeeze(TMP.AllDat(optogroupedid{1},ismember(ReactionOpt,{'Hit','Error'}),:)); %opto Off
    tmpon = squeeze(TMP.AllDat(optogroupedid{2},ismember(ReactionOpt,{'Hit','Error'}),:)); %Opto on
    minntr = round(min(min(cell2mat(cellfun(@(X) size(X,4),tmpoff,'UniformOutput',0)))).*0.75);

    % Prepare Data for MALSAR
    XDattrain = nan(40,40,minntr,2,2,'single'); %10 times downsampled
    XDattest_off = cell(2,2);
    XDattest_on = cell(2,2);
    trialidx = nan(minntr,2,2);
    reactionvec = nan(minntr,2,2);
    sidevec = nan(minntr,2,2);
    for ridx = 1:2
        for sidx = 1:2
            trialidx(:,ridx,sidx) = randsample(size(tmpoff{ridx,sidx},4),minntr,false);
            %Take late delay timewindow
            XDattrain(:,:,:,ridx,sidx) = Nan_imresize(nanmean(tmpoff{ridx,sidx}(:,:,timeline>=TW{5}(1)&timeline<=TW{5}(2),trialidx(:,ridx,sidx)),3),0.1);
            reactionvec(:,ridx,sidx) = ridx;
            sidevec(:,ridx,sidx)=sidx;

            % Test data - Opto Off
            XDattest_off{ridx,sidx} = squeeze(Nan_imresize(nanmean(tmpoff{ridx,sidx}(:,:,timeline>=TW{5}(1)&timeline<=TW{5}(2),:),3),0.1));
            XDattest_off{ridx,sidx}(:,:,trialidx(:,ridx,sidx)) = []; %Exclude trainingset

            % Test data - Opto ON
            tmp2 = cat(4,tmpon{:,ridx,sidx});
            XDattest_on{ridx,sidx} = squeeze(Nan_imresize(nanmean(tmp2(:,:,timeline>=TW{5}(1)&timeline<=TW{5}(2),:),3),0.1));

        end
    end

    %Prepare training data
    XDattrain = reshape(XDattrain,40*40,[]);
    reactionvec = reactionvec(:);
    reactionvec(reactionvec==2)=-1; %to make error -1
    sidevec = sidevec(:);
    sidevec(sidevec==2)=-1; %To make right grating -1
    choicevec = sidevec;
    choicevec(reactionvec==-1&sidevec==-1) = 1;
    choicevec(reactionvec==-1&sidevec==1)=-1;

    % Establish which pixels should run the analysis
    pixelidx = find(~isnan(nanmean(reshape(XDattrain,40*40,[]),2)));

    % Define parameters
    paramMTL_max = DefineBestParam(XDattrain(pixelidx,:)',cat(2,sidevec,choicevec,reactionvec)',length(choicevec),part,Int,nfolds);

    % use nfolds = 1 to extracts weights on full training data
    output = MALSAR_MTL_WF_Fast(XDattrain(pixelidx,:)',cat(2,sidevec,choicevec,reactionvec)',length(choicevec),part,Int,1,paramMTL_max,1,1)

    %Extract weights
    w = squeeze(output.Weights.Modalities);

    % prepare testdata
    XDattest_off = cellfun(@(X) reshape(X,40*40,[]),XDattest_off,'UniformOutput',0);
    nrtrials = cell2mat(cellfun(@(X) size(X,2),XDattest_off,'UniformOutput',0));
    reactionvec = [];
    sidevec = [];
    Xtest = nan(length(pixelidx),0);
    for ridx=1:2
        for sidx=1:2
            reactionvec = [reactionvec repmat(ridx,1,nrtrials(ridx,sidx))];
            sidevec = [sidevec repmat(sidx,1,nrtrials(ridx,sidx))];
            Xtest = cat(2,Xtest,XDattest_off{ridx,sidx}(pixelidx,:));
        end
    end
    reactionvec(reactionvec==2)=-1; %to make error -1
    sidevec(sidevec==2)=-1; %To make right grating -1
    choicevec = sidevec;
    choicevec(reactionvec==-1&sidevec==-1) = 1;
    choicevec(reactionvec==-1&sidevec==1)=-1;

    % Contralateral grating trials
    Predict = w*cat(1,ones(1,size(Xtest,2)),Xtest);
    nansum(cat(1,sidevec,choicevec,reactionvec) == sign(Predict),2)./sum(~isnan(Predict),2);


    % Prepare test data opto ON
    XDattest_on = cellfun(@(X) reshape(X,40*40,[]),XDattest_on,'UniformOutput',0);


    PerfWholeBrain(:,twid,midx) = nansum(Outputboot.PredictionsVsActual(:,:,3)==sign(Outputboot.PredictionsVsActual(:,:,1)),2)./sum(~isnan(Outputboot.PredictionsVsActual(:,:,3)),2);
    tmp = squeeze(nanmean(Outputboot.Weights.Modalities,1))';
    BetaPerPix(pixelidx,:,midx) = tmp(2:end,:);

    disp(['Whole brain performance ' miceopt{midx} ' ' TWRewardBiasNames{twid}])
    disp(['p_side = ' num2str(PerfWholeBrain(1,twid,midx))])
    disp(['p_choice = ' num2str(PerfWholeBrain(2,twid,midx))])
    disp(['p_reward = ' num2str(PerfWholeBrain(3,twid,midx))])

    disp('Now running pixel based model')
    PredPerPix = nan(length(pixelidx),3,'single');
    PredPerPixBoot = nan(length(pixelidx),3,nboot,'single');
    % Calculate average performance (cross-validated) across pixels
    for pixid=1:length(pixelidx)
        % find 7x7 spotlight patch
        [r,c] = find(pixgrid==pixelidx(pixid));
        try
            pixtmp = pixgrid(r-3:r+3,c-3:c+3);
        catch
            continue
        end
        Outputboot = MALSAR_MTL_WF_Fast(XDat(pixtmp(:),:)',cat(2,sidevec,choicevec,reactionvec)',length(choicevec),part,Int,nfolds,paramMTL_max,1,1);
        PredPerPix(pixid,:) = nansum(Outputboot.PredictionsVsActual(:,:,3)==sign(Outputboot.PredictionsVsActual(:,:,1)),2)./sum(~isnan(Outputboot.PredictionsVsActual(:,:,3)),2);
        parfor bootid=1:nboot
            tmpy=cat(2,sidevec,choicevec,reactionvec)';
            tmpy = tmpy(1,randsample(size(tmpy,2),size(tmpy,2),true));
            Outputboot = MALSAR_MTL_WF_Fast(XDat(pixtmp(:),:)',tmpy,length(choicevec),part,Int,nfolds,paramMTL_max,1,1);
            PredPerPixBoot(pixid,:,bootid) = nansum(Outputboot.PredictionsVsActual(:,:,3)==sign(Outputboot.PredictionsVsActual(:,:,1)),2)./sum(~isnan(Outputboot.PredictionsVsActual(:,:,3)),2);
        end
    end
    Gmap1 = nan(40,40);
    Cmap2 = nan(40,40);
    Rmap3 = nan(40,40);
    Gmap1(pixelidx)=PredPerPix(:,1);
    Cmap2(pixelidx)=PredPerPix(:,2);
    Rmap3(pixelidx)=PredPerPix(:,3);

    Gmaptresh = nan(40,40);
    Cmaptresh = nan(40,40);
    Rmaptresh = nan(40,40);
    Gmaptresh(pixelidx)=quantile(PredPerPixBoot(:,1,:),0.95,3);
    Cmaptresh(pixelidx)=quantile(PredPerPixBoot(:,2,:),0.95,3);
    Rmaptresh(pixelidx)=quantile(PredPerPixBoot(:,3,:),0.95,3);


    % Save maps
    GratingMap(:,:,twid,midx) = reshape(Gmap1,40,40);
    ChoiceMap(:,:,twid,midx) = reshape(Cmap2,40,40);
    RewardMap(:,:,twid,midx) = reshape(Rmap3,40,40);

    GratingTMap(:,:,twid,midx) = reshape(Gmaptresh,40,40);
    ChoiceTmap(:,:,twid,midx) = reshape(Cmaptresh,40,40);
    RewardMap(:,:,twid,midx) = reshape(Rmaptresh,40,40);

    PerMouse{midx} = PredPerPix;
    subplot(2,3,(twid-1)*3+1)
    h=imagesc(imresize(reshape(Gmap1,40,40),10),[0.25 0.75]);
    colormap redblue
    binimg = zeros(400,400);
    binimg(~isnan(imresize(reshape(Gmap1,40,40),10)))=0.5;
    binimg(imresize(reshape(Gmap1>Gmaptresh,40,40),10)) = 1;
    binimg(isnan(imresize(reshape(Gmap1,40,40),10)))=0;
    set(h,'AlphaData',binimg)
    hold on
    plot((model.BrainModel.AllX),(model.BrainModel.AllY),'k.','MarkerSize',2)
    title([TWRewardBiasNames{twid} '  Grating Left/Right'])
    box off
    axis off


    subplot(2,3,(twid-1)*3+2)
    h=imagesc(imresize(reshape(Cmap2,40,40),10),[0.25 0.75]);
    colormap redblue
    binimg = zeros(400,400);
    binimg(~isnan(imresize(reshape(Cmap2,40,40),10)))=0.5;
    binimg(imresize(reshape(Cmap2>Cmaptresh,40,40),10)) = 1;
    binimg(isnan(imresize(reshape(Cmap2,40,40),10)))=0;
    set(h,'AlphaData',binimg)
    hold on
    plot((model.BrainModel.AllX),(model.BrainModel.AllY),'k.','MarkerSize',2)
    title([TWRewardBiasNames{twid} '  Choice left/right'])
    box off
    axis off


    subplot(2,3,(twid-1)*3+3)
    h=imagesc(imresize(reshape(Rmap3,40,40),10),[0.25 0.75]);
    colormap redblue
    binimg = zeros(400,400);
    binimg(~isnan(imresize(reshape(Rmap3,40,40),10)))=0.5;
    binimg(imresize(reshape(Rmap3>Rmaptresh,40,40),10)) = 1;
    binimg(isnan(imresize(reshape(Rmap3,40,40),10)))=0;
    set(h,'AlphaData',binimg)
    hold on
    plot((model.BrainModel.AllX),(model.BrainModel.AllY),'k.','MarkerSize',2)
    title([TWRewardBiasNames{twid} '  Reward On/Off'])
    box off
    axis off

    for areaid=1:length(AnalysisParameters.AREASGrouped)
        for hemid=1:2
            % Find area pixels
            idx = find(ismember(model.BrainModel.Rnames,AREAS{areaid}));
            if isempty(idx)
                idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),model.BrainModel.Rnames,'UniformOutput',0)));
                areamask = false(400,400);
                for i = 1:length(idx)
                    areamask(logical(round(model.BrainModel.Regions{idx(i)}))')=true;
                end
            else
                areamask = logical(any(cat(3,model.BrainModel.Regions{idx}),3))';
            end
            if hemid==2 %left hemisphere
                areamask(1:round(model.BrainModel.Lambda(2)),:)=0;
            else
                areamask(round(model.BrainModel.Lambda(2)):end,:)=0;
            end
            areamask = imresize(areamask,0.1);
            pixidx = find(areamask);
            %                 gratingperformance = nan(length(AREAS),2,length(miceopt));
            gratingperformance(areaid,hemid,midx) = nanmean(Gmap1(pixidx));

            %                 choiceperformance = nan(length(AREAS),2,length(miceopt));
            choiceperformance(areaid,hemid,midx) = nanmean(Cmap2(pixidx));
            rewardperformance(areaid,hemid,midx) = nanmean(Rmap3(pixidx));
        end
    end

    save(fullfile(LocalFolder,['OptovsNoOptoPredictionPerPixel_' miceopt{midx} '.mat']),'RewardMap','ChoiceTmap','GratingTMap','gratingperformance','choiceperformance','rewardperformance','GratingMap','ChoiceMap','RewardMap')
    saveas(gcf,fullfile(LocalFolder,['OptovsNoOptoPredictionPerPixel_ ' miceopt{midx} '.fig']))

    save(fullfile(LocalFolder,'MalsarOutput_WMOpto.mat'),'gratingperformance','choiceperformance','rewardperformance',...
        'GratingMap','ChoiceMap','RewardMap','GratingTMap','ChoiceTmap','RewardMap','PerfWholeBrain','BetaPerPix')

end

%%
figure;
for did = 1:length(GenoTypeOpt)
    subplot(2,2,did)

    h=bar(squeeze(PerfWholeBrain(1,:,ismember(GenOpt,GenoTypeOpt{did}))));

    set(gca,'XTickLabel',LegendName,'XTickLabelRotation',25)
    box off
    ylim([0 100])
    hold on
    line([0.5 countid-0.5],[50 50],'LineStyle','--','color',[1 0 0])
    xlim([0.5 countid-0.5])

    countid=1;

    for sidx = 1:length(GratingOpt)
        for gidx = 1:length(GenOpt)
            for optidx=1:length(StimOpt)
                plot(h.XData(countid),Data{countid},'ko','MarkerSize',10)
                countid=countid+1;
            end
        end
    end

end
%% Build linear model - using cross validation
nfolds = 5;
part = 2; %1 = seperate learning, 2 = simultaneous
Int = 1; % Int = 2: integration of modalities, Int = 1: concatenation of modalities
pixgrid = reshape(1:40*40,40,40);


for midx=1:length(miceopt)

    TMP = load(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));


    % Take opto OFF data, and Reaction Hit/Error
    tmp = squeeze(TMP.AllDat(strcmp(OptoOpt,'Opto Off'),ismember(ReactionOpt,{'Hit','Error'}),:));
    minntr = min(min(cell2mat(cellfun(@(X) size(X,4),tmp,'UniformOutput',0))));

    % Sample minntr trials from each condition
    for twid=1:length(TWRewardBiasAna)

        % Prepare Data for MALSAR
        % Prepare Data for MALSAR
        XDat = nan(40,40,minntr,2,2,'single'); %10 times downsampled
        trialidx = nan(minntr,2,2);
        reactionvec = nan(minntr,2,2);
        sidevec = nan(minntr,2,2);
        for ridx = 1:2
            for sidx = 1:2
                trialidx(:,ridx,sidx) = randsample(size(tmp{ridx,sidx},4),minntr,false);
                XDat(:,:,:,ridx,sidx) = Nan_imresize(nanmean(tmp{ridx,sidx}(:,:,timeline>=TWRewardBiasAna{twid}(1)&timeline<=TWRewardBiasAna{twid}(2),squeeze(trialidx(:,ridx,sidx))),3),0.1);
                reactionvec(:,ridx,sidx) = ridx;
                sidevec(:,ridx,sidx)=sidx;
            end
        end
        XDat = reshape(XDat,40*40,[]);
        reactionvec = reactionvec(:);
        reactionvec(reactionvec==2)=-1; %to make error -1
        sidevec = sidevec(:);
        sidevec(sidevec==2)=-1; %To make right grating -1
        choicevec = sidevec;
        choicevec(reactionvec==-1&sidevec==-1) = 1;
        choicevec(reactionvec==-1&sidevec==1)=-1;

        % Establish which pixels should run the analysis
        pixelidx = find(~isnan(nanmean(reshape(XDat,40*40,[]),2)));

        % Create DM - Grating, choice, reward
        DM = cat(2,ones(length(sidevec),1),sidevec,choicevec,reactionvec);

        ntrials = size(XDat,2);
        randomvec = randsample(ntrials,ntrials,0);
        nrperfold = floor(ntrials/nfolds);

        %Initialize
        dFFmap = nan(40*40,nfolds);
        Imap = nan(40*40,nfolds);
        Gmap = nan(40*40,nfolds);
        Cmap = nan(40*40,nfolds);
        Rmap = nan(40*40,nfolds);
        Emap = nan(40*40,nfolds);
        parfor cvid = 1:nfolds
            testid = randomvec((cvid-1)*nrperfold+1:cvid*nrperfold);
            trainid = randomvec(~ismember(randomvec,testid));
            [betavals,stdvals]=lscov(fillmissing(XDat(pixelidx,trainid)','linear',1),DM(trainid,:));

            %
            dFFPred = betavals*DM(testid,:)';
            Emap(pixelidx,cvid) = nanmean((fillmissing(XDat(pixelidx,testid)','linear',1)-dFFPred'),1);
            dFFmap(pixelidx,cvid) = nanmean(fillmissing(XDat(pixelidx,testid)','linear',1),1);

            Imap(pixelidx,cvid) = betavals(:,1);
            Gmap(pixelidx,cvid) = betavals(:,2);
            Cmap(pixelidx,cvid) = betavals(:,3);
            Rmap(pixelidx,cvid) = betavals(:,4);
        end

        figure('name',[miceopt{midx} TWRewardBiasNames{twid}])
        tmp2 = smooth2a(fliplr(reshape(nanmean(dFFmap,2),40,40)'),2);
        lims = quantile(abs(tmp2(:)),0.99);

        subplot(2,3,1)
        h=imagesc(smooth2a(fliplr(reshape(nanmean(dFFmap,2),40,40)'),2),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',fliplr(~isnan(nanmean(reshape(XDat,40,40,[]),3))'))
        title(['dFF'])
        axis off

        tmp2 = smooth2a(fliplr(reshape(nanmean(Emap,2),40,40)'),2);
        lims = quantile(abs(tmp2(:)),0.99);

        subplot(2,3,2)
        h=imagesc(smooth2a(fliplr(reshape(nanmean(Imap,2),40,40)'),2),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',fliplr(~isnan(nanmean(reshape(XDat,40,40,[]),3))'))
        title(['Intercept'])
        axis off

        subplot(2,3,3)
        h=imagesc(smooth2a(fliplr(reshape(nanmean(Gmap,2),40,40)'),2),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',fliplr(~isnan(nanmean(reshape(XDat,40,40,[]),3))'))
        title(['Grating'])
        axis off

        subplot(2,3,4)
        h=imagesc(smooth2a(fliplr(reshape(nanmean(Cmap,2),40,40)'),2),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',fliplr(~isnan(nanmean(reshape(XDat,40,40,[]),3))'))
        title(['Choice'])
        axis off

        subplot(2,3,5)
        h=imagesc(smooth2a(fliplr(reshape(nanmean(Rmap,2),40,40)'),2),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',fliplr(~isnan(nanmean(reshape(XDat,40,40,[]),3))'))
        title(['Reward'])
        axis off

        subplot(2,3,6)
        h=imagesc(smooth2a(fliplr(reshape(nanmean(Emap,2),40,40)'),2),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',fliplr(~isnan(nanmean(reshape(XDat,40,40,[]),3))'))
        title(['Error'])
        axis off

    end
end

%% Average performance
figure('name','ModelPerformance')
subplot(1,3,1)
barwitherr(nanstd(gratingperformance,[],3)./sqrt(length(miceopt)-1),nanmean(gratingperformance,3))
set(gca,'XTickLabel',AREAGROUPNAMES)
legend(TWRewardBiasNames)
box off
ylabel('average p-value')
ylim([0 1])
title('Grating')

subplot(1,3,2)
barwitherr(nanstd(choiceperformance,[],3)./sqrt(length(miceopt)-1),nanmean(choiceperformance,3))
set(gca,'XTickLabel',AREAGROUPNAMES)
legend(TWRewardBiasNames)
box off
ylabel('average p-value')
ylim([0 1])
title('Choice')

subplot(1,3,3)
barwitherr(nanstd(rewardperformance,[],3)./sqrt(length(miceopt)-1),nanmean(rewardperformance,3))
set(gca,'XTickLabel',AREAGROUPNAMES)
legend(TWRewardBiasNames)
box off
ylabel('average p-value')
ylim([0 1])
title('Reward')

saveas(gcf,fullfile(LocalFolder,['PredictionAcrossAreas.fig']))


%% Decoder All Pixels during delay and ITI
nfolds = 10;
part = 2; %1 = seperate learning, 2 = simultaneous
Int = 1; % Int = 2: integration of modalities, Int = 1: concatenation of modalities
for midx=1:length(miceopt)
    figure('name',miceopt{midx})

    TMP = load(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));

    % Take opto OFF data, and Reaction Hit/Error
    tmp = squeeze(TMP.AllDat(strcmp(OptoOpt,'Opto Off'),ismember(ReactionOpt,{'Hit','Error'}),:));
    minntr = min(min(cell2mat(cellfun(@(X) size(X,4),tmp,'UniformOutput',0))));

    % Sample minntr trials from each condition
    for twid=1:length(TWRewardBiasAna)

        % Prepare Data for MALSAR
        XDat = nan(40,40,minntr,2,2,'single'); %10 times downsampled
        trialidx = nan(2,2,minntr);
        reactionvec = nan(2,2,minntr);
        sidevec = nan(2,2,minntr);
        for ridx = 1:2
            for sidx = 1:2
                trialidx(ridx,sidx,:) = randsample(size(tmp{ridx,sidx},4),minntr,false);
                XDat(:,:,:,ridx,sidx) = Nan_imresize(nanmean(tmp{ridx,sidx}(:,:,timeline>=TWRewardBiasAna{twid}(1)&timeline<=TWRewardBiasAna{twid}(2),trialidx(ridx,sidx,:)),3),0.1);
                reactionvec(ridx,sidx,:) = ridx;
                sidevec(ridx,sidx,:)=sidx;
            end
        end
        XDat = reshape(XDat,40*40,[]);
        reactionvec = reactionvec(:);
        reactionvec(reactionvec==2)=-1; %to make error 0
        sidevec = sidevec(:);
        sidevec(sidevec==2)=-1; %To make right grating 0

        % Establish which pixels should run the analysis
        pixelidx = find(~isnan(nanmean(reshape(XDat,40*40,[]),2)));

        % Define parameters
        paramMTL_max = DefineBestParam(XDat(pixelidx,:)',cat(2,sidevec,reactionvec)',length(reactionvec),part,Int,nfolds);
        Outputboot = MALSAR_MTL_WF_Fast(XDat(pixelidx,:)',cat(2,sidevec,reactionvec)',length(reactionvec),part,Int,nfolds,paramMTL_max,1,1);
        nansum(Outputboot.PredictionsVsActual(1,:,3)==sign(Outputboot.PredictionsVsActual(1,:,1)))./sum(~isnan(Outputboot.PredictionsVsActual(1,:,3)))

        gratingmap = nan(40,40);
        gratingmap(pixelidx) = squeeze(nanmean(Outputboot.Weights.Modalities(:,1,2:end),1));

        reactionmap = nan(40,40);
        reactionmap(pixelidx) = squeeze(nanmean(Outputboot.Weights.Modalities(:,2,2:end),1));

        lims = quantile(abs(gratingmap(:)),0.99);

        subplot(2,2,(twid-1)*2+1)
        h=imagesc(reshape(gratingmap,40,40),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',~isnan(gratingmap))
        title([TWRewardBiasNames{twid} '  Grating Side'])
        axis off

        subplot(2,2,twid*2)
        h=imagesc(reshape(reactionmap,40,40),[-lims lims]);
        colormap redblue
        set(h,'AlphaData',~isnan(reactionmap))
        title([TWRewardBiasNames{twid} '  Response Side'])
        axis off
    end
end

%% Statistical Test
TMP = [];
g1 = []; %area
g2 = []; %hemisphere
g3 = []; %Side
g4 = []; %opto
g5 = []; %Genotype
for areaid=1:length(AREAS)
    for hemid = 1:2
        for sidx=1:length(SideOpt)
            for optidx = 1:length(OptoOptGrouped)
                for dsidx=1:2
                    midx = find(ismember(GenOpt,GenoTypeOpt{dsidx}));
                    if optidx==1
                        opt2take = 1;
                    elseif optidx==2
                        opt2take = [2,3,4];
                    elseif optidx==3
                        opt2take = [5,6]
                    end
                    if strcmp(HemOpt{hemid},'RightHem')
                        if strcmp(SideOpt{sidx},'left')
                            sidx2take = 1;
                        elseif strcmp(SideOpt{sidx},'right')
                            sidx2take = 2;
                        else
                            continue
                        end
                    else
                        if strcmp(SideOpt{sidx},'left')
                            sidx2take = 2;
                        elseif strcmp(SideOpt{sidx},'right')
                            sidx2take = 1;
                        else
                            continue
                        end
                    end
                    tmp = squeeze(nanmean(AllMiceTC(areaid,timeline>=TW{5}(1)&timeline<=TW{5}(2),:,sidx2take,hemid,opt2take,midx),2));
                    TMP = [TMP reshape(tmp,1,[])];
                    g1 = [g1 repmat(AREAS(areaid),[1,prod(size(tmp))])];
                    g2 = [g2 repmat(HemOpt(hemid),[1,prod(size(tmp))])];
                    g3 = [g3 repmat(SideOpt(sidx),[1,prod(size(tmp))])];
                    g4 = [g4 repmat(OptoOptGrouped(optidx),[1,prod(size(tmp))])];
                    g5 = [g5 repmat(GenoTypeOpt(dsidx),[1,prod(size(tmp))])];
                end
            end
        end
    end
end

rmid = isnan(TMP);
TMP(rmid)=[];
g1(rmid)=[];
g2(rmid)=[];
g3(rmid)=[];
g4(rmid)=[];
g5(rmid)=[];
[p,tbl,stats] = anovan(TMP,{g4,g1,g2,g3,g5},'Model',5,'VarNames',{'Opto','Area','Hemisphere','FigSide','Genotype'});



%% Load template mouse first
load(fullfile(storepath,AnalysisParameters.templatemouse, 'brainareamodel.mat'))
mask = zeros(800*AnalysisParameters.ScaleFact,800*AnalysisParameters.ScaleFact);
xpix = 800*AnalysisParameters.ScaleFact;
ypix = 800*AnalysisParameters.ScaleFact;
flagscaler = 0;
if any(size(Model.Regions{1})>size(mask))
    Model.Regions = cellfun(@(X) imresize(X,AnalysisParameters.ScaleFact),Model.Regions,'UniformOutput',0);
    Model.AllX = round(Model.AllX.*AnalysisParameters.ScaleFact);
    Model.AllY = round(Model.AllY.*AnalysisParameters.ScaleFact);
    Model.Lambda = round(Model.Lambda.*AnalysisParameters.ScaleFact);
    flagscaler=1;
end
halfidx = round(median(Model.AllY));

%Remove areas
throwawayareas = find(cellfun(@isempty,Model.Rnames));
throwawayareas = [throwawayareas; find(cellfun(@(X) ismember(X,{'OlfactoryBulb','fibrtracts','InfCol','SupColSens'}),Model.Rnames))];
keepareas = 1:length(Model.Rnames);
keepareas(throwawayareas)=[];
removepix = true(xpix,ypix);
for areaid = 1:length(keepareas)
    bounds = Model.Boundaries{keepareas(areaid)};
    for boundid = 1:length(bounds)
        if flagscaler|| any(any(bounds{boundid}>xpix*1.5))
            bounds{boundid} = round(bounds{boundid}.*AnalysisParameters.ScaleFact);
        end
        removepix(poly2mask(bounds{boundid}(:,1),bounds{boundid}(:,2),xpix,ypix)) = 0;
    end
end

storepath = 'E:\Enny\TMPMatlab\StoreFigures'

%% Across mice
for midx=1:length(miceopt)
    load(fullfile(LocalFolder,miceopt{midx},'TimeCourses.mat'))
    mouse = miceopt{midx}
    AllPaths = dir(fullfile(LocalFolder,miceopt{midx}));
    AllPaths = AllPaths(AllPaths.isdir);
    AllPaths(cell2mat(cellfun(@isempty,AllPaths,'UniformOutput',0)))=[];
    tmppath = AllPaths{1};
    date = strsplit(tmppath,mouse);
    date = date{3}(1:end-1) %Find date
    expnr = strsplit(tmppath,mouse);
    expnr = str2num(expnr{end});%find session nr
    TMPMat = matfile(fullfile(storepath, mouse, [mouse date], [mouse num2str(expnr)],[mouse num2str(expnr) '_RawData_All.mat'])); %Make a workable matfile

    timelinehere = TMPMat.timeline;
    if midx==1
        AllMiceTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(ReactionOpt),length(SideOpt),length(HemOpt),length(OptoOpt),length(miceopt));
    end
    try
        AllMiceTC(:,ismember(timeline,timelinehere),:,:,:,:,midx)=AllTC;
    catch
        AllMiceTC(:,:,:,:,:,:,midx)=AllTC;
    end
end
%% Statistical Test
TMP = [];
g1 = []; %area
g2 = []; %hemisphere
g3 = []; %Side
g4 = []; %opto
g5 = []; %Genotype
for areaid=1:length(AREAS)
    for hemid = 1:2
        for sidx=1:length(SideOpt)
            for optidx = 1:length(OptoOptGrouped)
                for dsidx=1:2
                    midx = find(ismember(GenOpt,GenoTypeOpt{dsidx}));
                    if optidx==1
                        opt2take = 1;
                    elseif optidx==2
                        opt2take = [2,3,4];
                    elseif optidx==3
                        opt2take = [5,6]
                    end
                    if strcmp(HemOpt{hemid},'RightHem')
                        if strcmp(SideOpt{sidx},'left')
                            sidx2take = 1;
                        elseif strcmp(SideOpt{sidx},'right')
                            sidx2take = 2;
                        else
                            continue
                        end
                    else
                        if strcmp(SideOpt{sidx},'left')
                            sidx2take = 2;
                        elseif strcmp(SideOpt{sidx},'right')
                            sidx2take = 1;
                        else
                            continue
                        end
                    end
                    tmp = squeeze(nanmean(AllMiceTC(areaid,timeline>=TW{5}(1)&timeline<=TW{5}(2),:,sidx2take,hemid,opt2take,midx),2));
                    TMP = [TMP reshape(tmp,1,[])];
                    g1 = [g1 repmat(AREAS(areaid),[1,prod(size(tmp))])];
                    g2 = [g2 repmat(HemOpt(hemid),[1,prod(size(tmp))])];
                    g3 = [g3 repmat(SideOpt(sidx),[1,prod(size(tmp))])];
                    g4 = [g4 repmat(OptoOptGrouped(optidx),[1,prod(size(tmp))])];
                    g5 = [g5 repmat(GenoTypeOpt(dsidx),[1,prod(size(tmp))])];
                end
            end
        end
    end
end

rmid = isnan(TMP);
TMP(rmid)=[];
g1(rmid)=[];
g2(rmid)=[];
g3(rmid)=[];
g4(rmid)=[];
g5(rmid)=[];
[p,tbl,stats] = anovan(TMP,{g4,g1,g2,g3,g5},'Model',5,'VarNames',{'Opto','Area','Hemisphere','FigSide','Genotype'});
%%
% cols((dsid-1).*3+optidx,:)
cols = [0 0 0; 1 0 0; 0.7 0.3 0; 0.5 0.5 0.5;  0 0 1; 0 0.3 0.7];
linestyls = {'-','--'};
GenoTypeOpt = unique(GenOpt);
for areaid = 1:length(AnalysisParameters.AREAS)
    for hemid = 1:length(HemOpt)

        FF(areaid) = figure('name',['AcrossMice_' AnalysisParameters.AREAS{areaid} '_' HemOpt{hemid}]);
        clear h
        miny=[];
        maxy=[];
        legendnames = {};
        count=1;
        for dsid = 1:length(GenoTypeOpt)
            for optidx = 1:length(OptoOptGrouped)%
                midx = find(ismember(GenOpt,GenoTypeOpt{dsid}));
                for sidx = 1:length(SideOpt)
                    % Per Area
                    if strcmp(HemOpt{hemid},'RightHem')
                        if strcmp(SideOpt{sidx},'left')
                            sidx2take = 1;
                            LegendInput={'Contra'};
                        elseif strcmp(SideOpt{sidx},'right')
                            sidx2take = 2;
                            LegendInput={'Ipsi'};
                        else
                            continue
                        end
                    else
                        if strcmp(SideOpt{sidx},'left')
                            sidx2take = 2;
                            LegendInput={'Ipsi'};
                        elseif strcmp(SideOpt{sidx},'right')
                            sidx2take = 1;
                            LegendInput={'Contra'};
                        else
                            continue
                        end
                    end
                    %Dont'plot
                    if strcmp(SideOpt{sidx},'none')
                        continue
                    end
                    if optidx==1
                        opt2take = 1;
                    elseif optidx==2
                        opt2take = [2,3,4];
                    elseif optidx==3
                        opt2take = [5,6];
                    end
                    tmp = reshape(AllMiceTC(areaid,:,1:2,sidx2take,hemid,opt2take,midx),length(timeline),[]);
                    if isempty(tmp)
                        continue
                    end
                    tmpmean = nanmean(tmp,2);
                    tmpsem = nanstd(tmp,[],2)./sqrt(size(tmp,2)-1);
                    if dsid==1
                        subplot(4,4,[1:3,5:7])
                    else
                        subplot(4,4,[9:11,13:15])
                    end
                    hold on

                    h{count} = shadedErrorBar(timeline,tmpmean,tmpsem,{'LineWidth',2,'color',cols((dsid-1).*3+optidx,:),'LineStyle',linestyls{sidx}},1);

                    box off
                    drawnow
                    if optidx==1
                        legendnames = {legendnames{:} [SideOpt{sidx} 'No Opto']};
                    elseif optidx==2
                        legendnames = {legendnames{:} [SideOpt{sidx} 'Opto ' GenoTypeOpt{dsid}]};
                    end
                    count=count+1;
                end
            end


            xlim([-600 6000])
            ytmp = get(gca,'ylim');
            if isempty(miny)|| ytmp(1)<miny
                miny = ytmp(1);
            end
            if isempty(maxy) || ytmp(2)>maxy
                maxy=ytmp(2);
            end

            set(gca,'ylim',[-0.02 0.035])
            ytmp = get(gca,'ylim');
            widthtmp = 0.05*(max(ytmp)-min(ytmp));
            hp(1) = patch([500 2000 2000 500],[min(ytmp) min(ytmp) min(ytmp)+widthtmp min(ytmp)+widthtmp],[1 0 0]);
            hp(2) = patch([0 500 500 0],[min(ytmp) min(ytmp) min(ytmp)+widthtmp min(ytmp)+widthtmp],[0 0.5 0]);
            hp(3) = patch([2000 3500 3500 2000],[min(ytmp) min(ytmp) min(ytmp)+widthtmp min(ytmp)+widthtmp],[0 0 0.5]);


            ylabel(['dFF ' OptoOpt{optidx}])




            subplot(4,4,[4,8])
            idx = find(ismember(Model.Rnames,AnalysisParameters.AREAS{areaid}));
            if isempty(idx) || strcmp(AREAS{areaid},'Vl')
                idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AnalysisParameters.AREAS{areaid}),Model.Rnames,'UniformOutput',0)));
                areamask = false(xpix,ypix);
                for i = 1:length(idx)
                    areamask(logical(round(Model.Regions{idx(i)}))')=true;
                end
            else
                areamask = logical(round(Model.Regions{idx}))';
            end
            if hemid==1
                areamask(halfidx:end,:)=0;
            else
                areamask(1:halfidx,:)=0;

            end
            ha=imagesc(areamask);
            set(ha,'AlphaData',areamask==1)
            colormap(jet)
            hold on
            plot(Model.AllX,Model.AllY,'k.','MarkerSize',3)
            title(AREAS{areaid})
            Pos = get(gca,'Position');
            axis square
            axis off
        end
        hsub = subplot(4,4,[12,16]);
        Pos = get(hsub,'Position');
        delete(hsub)


        saveas(gcf,fullfile(storepath,['AcrossMice_' HemOpt{hemid} AnalysisParameters.AREAS{areaid} '.fig']))
        saveas(gcf,fullfile(storepath,['AcrossMice_' HemOpt{hemid} AnalysisParameters.AREAS{areaid} '.bmp']))
    end
end
