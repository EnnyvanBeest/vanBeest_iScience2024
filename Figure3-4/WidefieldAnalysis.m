%% User-Input
Scriptsdir = '\Path\To\GitHub\Scripts';%
storepath = '\Path\To\Storage\' %Main storage --> Will contain all processed data & model of the brain
RawDataDirectory = '\Path\To\RawImagingData\' %Raw images.
DataDirectory = RawDataDirectory; %Might be different when downsampling data first
miceopt = {'Faraday','Galilei','Haydn','Jemison','Karlsen','Lully','Mahler'}%
mousegen = {'Drd2','Drd2','Drd1','Drd2','Drd1','Drd1','Drd1'};%
MappingDir = '\Path\To\pRF_Results\' % Allen Brain Mapping
LocalFolder = '\Path\To\Local\';
ServerDir = '\Path\To\Server\'
Stim2Check = 'SimpleLick' %Name of the stimulus as written in the LOG-file

AnalysisParameters.UsepRFMask = 0; % Use only those pixels that show up in the pRF mask
AnalysisParameters.UseEvokedMask = 0; %Use only pixels taht show up in the Evoked Response (to a figure on grey BG)
AnalysisParameters.WindowSize = 15; %Nr. trials around window
AnalysisParameters.baselinemethod = 1; %1: trial by trial baseline, 2: filtered baseline (F0), 3: weighted average baseline after detrending (over all trials)
AnalysisParameters.takeequalsample = 0;  %Eq nr of trials per condition
AnalysisParameters.dontthrowtrialsgenerally = 0; %Throw out trials based on motion
AnalysisParameters.smoothfact = 2; %Smoothing with neighbours
AnalysisParameters.ScaleFact = 0.5; %between 0-1 is scale down
AnalysisParameters.normalize=0; %normalize to Figure on Grey condition
AnalysisParameters.AREAS = {'V1','Vpor','Vl','Vli','Val','Vrl','Va','Vam','Vpm','RSP','SS1barrel','nose','SSmouth','M1','M2'};%{'V1','HVLat_Vl_Val','HVAnterior_Vrl_Va','HVMedial_Vam_Vpm','RSP','M1','M2'};%% %
AnalysisParameters.AREASGrouped = {{'V1'},{'Vl','Vli','Val','Vpor'},{'Vrl','Va'},{'Vam','Vpm'},{'RSP'},{'M1'},{'M2'}};%note: l includes li and l{'V1','HVLat_Vl_Val','HVAnterior_Vrl_Va','HVMedial_Vam_Vpm','RSP','M1','M2'};%% %
AREAGROUPNAMES = {'V1','Vlat','PPC','Vmed','RS','M1','M2'};
AnalysisParameters.ReactionOpt = {'Hit','Miss'}; %Which response types to include
AnalysisParameters.SideOpt = {'left','right'};
AnalysisParameters.trialtypes = {'same'};
AnalysisParameters.templatemouse='Jemison'

AnalysisParameters.UseCluster = 0; %Make use of found clusters?
global UserQuestions %Defines whether gui's show up to ask what to do with existing datasets. If 0, existing datasets are used and not overwritten
UserQuestions = false; %put to 1 to do stuff manually
OptoOpt = {'OptoOff','OptoOnBeforeStim','OptoOn1stLick'};
SideOptNames = {'ContraOpto','IpsiOpto'};
AREAS = AnalysisParameters.AREASGrouped
ReactionOpt = AnalysisParameters.ReactionOpt
HemOpt = {'right','left'};
TW = {[-760 -500],[-500 0],[0 1500],[1500 3000]};
TimeWindowNames = {'Baseline','BaselineOpto','Visual+Opto','Post-Opto'};

%% Example no opto visual stimulus
for midx = 1:length(miceopt)
    TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
    Model = TMPDat.Model;
    SideOpt = TMPDat.SideOpt;
    timeline = TMPDat.timeline;

    figure('name',['No Opto Example ' miceopt{midx}])
    for sidx = 1:2
        for ridx = 1:2
            subplot(2,2,(sidx-1)*2+ridx)
            ReldFF = TMPDat.AlldFF(sidx,ridx,find(ismember(OptoOpt,'OptoOff'))); %BigDataTMP.AlldFF = cell(length(SideOpt),length(ReactionOpt),length(OptoOpt));
            ReldFF = ReldFF{1};
            tmp = nanmean(nanmean(ReldFF(:,:,timeline>0&timeline<1000,:),3),4)';
            h = imagesc(tmp,[-0.01 0.01]);
            colormap redblue
            set(h,'AlphaData',~isnan(tmp))
            hold on
            h2 = plot(Model.AllY,Model.AllX,'k.');
            axis off

            title([SideOpt{sidx} ' ' ReactionOpt{ridx}])
        end
    end
end

%% Linear Models
TMPDat = matfile(fullfile(ServerDir,miceopt{1},'ProcessedData.mat'));
timeline = TMPDat.timeline;
timevec=downsample(timeline,2); % downsample
SR = unique(diff(timevec)); if length(SR)>1; keyboard;end
newpixsz = 20; %Downsample to this
Redo=0;
ReloadDatPerPix = 0;

%Define kernel indices
DummyIdx = 1;
VisKernelIdx = 2:length(timevec)*2+1;
OptoKernelIdx = length(timevec)*2+2:length(timevec)*3+1;
LicksKernelIdx = length(timevec)*3+[2,3];
EyeMovKernelIdx = length(timevec)*3+4;
BodyMovKernelIdx = length(timevec)*3+5;
AllKernels = {DummyIdx,VisKernelIdx,OptoKernelIdx,...
    LicksKernelIdx,EyeMovKernelIdx,BodyMovKernelIdx};

AllKernelNames = {'Dummy','Stimulus','Opto','Licks','EyeMovement','BodyMovement'};
vecvector=2:4; % Which to shuffle analysis?
BootvecOpt = vecvector;
for midx=1:length(miceopt)
    if exist(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat'])) && ~Redo
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
    end
    close all
    TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
    % get brain model
    model = TMPDat.Model;
    disp('Loading in data, which will take a bit of time..')

    % Get dFF
    LM = matfile(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']));
    % we can load dat per pix
    if exist('LM') && ~ReloadDatPerPix
        DatPerPix = LM.DatPerPix;%nan(newpixsz,newpixsz,length(timevec),0);
        ReExtractDatPix = 0;
        DatPerPix = reshape(DatPerPix,newpixsz,newpixsz,length(timevec),[]);
    else
        AllDat = TMPDat.AlldFF; %BigDataTMP.AlldFF = cell(length(SideOpt),length(ReactionOpt),length(OptoOpt)); % Stimulus Options, ReactionOpt%OptoOpt,
        DatPerPix = nan(newpixsz,newpixsz,length(timevec),0);
        ReExtractDatPix = 1;
    end
    tmptrials = TMPDat.trialidentifier;
    tmptrialsall = [tmptrials{:}];
    tmpses = unique(tmptrialsall(2,:));

    %% Initialize
    SideVec = nan(1,0);
    OptoVec = nan(1,0);
    LeftActionVec = nan(length(timevec),0);
    RightActionVec = LeftActionVec;
    DifMotionVec = LeftActionVec;
    EyeDataVec = LeftActionVec;
    TrialVec = nan(1,0);
    SesVec = nan(1,0);
    for ridx=1
        for sidx=1:2
            for optidx=1:2
                if ReExtractDatPix
                    disp('Is this order correct?')
                    keyboard
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

                % Extract condition
                for sesid = 1:length(tmpses)
                    % Get correct log file
                    Log = TMPDat.LogsThisMouse;% matfile(fullfile(tmpdatFiles(tmpses(sesid)).folder,tmpdatFiles(tmpses(sesid)).name));
                    Log = Log{sesid};
                    %                     TimeLineHere = Log.timeline;
                    %                     Log = Log.Log;
                    takethesetrials = tmptrials{sidx,ridx,optidx}(1,tmptrials{sidx,ridx,optidx}(2,:)==tmpses(sesid));
                    if isempty(takethesetrials)
                        continue
                    end
                    Gavepassive = Log.Gavepassive(takethesetrials);

                    %                     takethesetrials(GavePassive)=[];

                    side = Log.Side(takethesetrials);
                    reaction = Log.Reaction(takethesetrials);
                    opto = Log.Opto(takethesetrials);

                    % Extract Motion and downsample in time
                    if isfield(Log,'DifMotion')
                        DifMotionVec = cat(2,DifMotionVec,resample(Log.DifMotion(timeline>=timevec(1)&timeline<=timevec(end),takethesetrials),length(timevec),length(timeline)));
                    else
                        DifMotionVec = cat(2,DifMotionVec,nan(length(timevec),length(takethesetrials)));
                    end
                    % Extract eye movements
                    if isfield(Log,'EyeMovement')
                        tmpeye = nanmean(Log.EyeMovement(timeline>=timevec(1)&timeline<=timevec(end),takethesetrials,:),3); %take average eyemotion across the 4 coordinates
                        EyeDataVec = cat(2,EyeDataVec,resample(tmpeye,length(timevec),length(timeline)));
                    elseif isfield(Log,'Pupil')
                        tmpeye = nanmean(Log.Pupil(timeline>=timevec(1)&timeline<=timevec(end),takethesetrials,:),3); %take average eyemotion across the 4 coordinates
                        EyeDataVec = cat(2,EyeDataVec,resample(tmpeye,length(timevec),length(timeline)));
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

                    % Assign condition
                    SideVec = cat(2,SideVec,repmat(sidx,1,length(takethesetrials)));
                    %                 RewardVec = cat(2,RewardVec,repmat(ridx,1,size(AllDat{sidx,ridx,optidx},4)));
                    OptoVec = [OptoVec repmat(optidx,1,length(takethesetrials))];
                    TrialVec = cat(2,TrialVec,takethesetrials);
                    SesVec = cat(2,SesVec,repmat(tmpses(sesid),1,length(takethesetrials)));


                end

            end
        end
    end
    % Create Design Matrix for events and their shifts: HRP, LRP,
    %Reward
    disp('Creating Design Matrix')
    ntrials = size(DatPerPix,4);
    DM = zeros(3,length(timevec),ntrials);
    %Extract left stimulus trials
    tmpDM = zeros(length(timevec),ntrials);
    tmpDM(1,SideVec==1)=1;
    DM(1,:,:) = tmpDM;

    %Extract right stimulus trials
    tmpDM = zeros(length(timevec),ntrials);
    tmpDM(1,SideVec==2)=1;
    DM(2,:,:) = tmpDM;

    % Extract opto event
    tmpDM = zeros(length(timevec),ntrials);
    tmpDM(1,OptoVec==2)=1;
    DM(3,:,:) = tmpDM;

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

    tmpfileinfo = dir(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']));
    if datetime(tmpfileinfo.date)<datetime('27-04-2023','InputFormat','dd-MM-yyyy') %otherwise already sorted
        DatPerPix = DatPerPix(:,:,:,sortidx);
    end

    % Check if brain top view is aligned
    CheckBrainAlignment

    DM = reshape(DM,size(DM,1),length(timevec),ntrials);
    DM = DM(:,:,sortidx);
    DM = reshape(DM,size(DM,1),length(timevec)*ntrials);

    DatPerPix = reshape(DatPerPix,newpixsz*newpixsz,[]);

    figure('name','DesignMatrix')
    imagesc(DM)
    colormap(flipud(gray))
    xlabel('TimexTrial')
    ylabel('Predictor')
    FigureDefault
    drawnow
    % Load in existing results -- saved during analysis to prevent time
    % loss
    %%
    RunLinearModels

    %Save stuff here
    disp('Saving...')
    try
        save(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']),'RExcl_PerVar_Shuffled','Rsq_PerVar','RExcl_PerVar','TotalR','ChosenK','PredicteddFF','PredicteddFF_LeaveOut','DatPerPix','ntrials')
    catch ME
        disp(ME)
        keyboard
    end

    %% Shuffled results
    figure('name',[miceopt{midx} ' Significance'])

    Pvalues = nan(newpixsz*newpixsz,nfolds,length(BootvecOpt));
    for vecid = 1:length(BootvecOpt)
        for cv=1:nfolds
            tmpboot = squeeze(RExcl_PerVar_Shuffled(pixid,BootvecOpt(vecid),cv,:));
            tmp = TotalR(pixid,cv);
            pval = invprctile(tmpboot,tmp,2);
            Pvalues(pixid,cv,vecid) = 1-diag(pval)./100;
        end
        PvalLabel = zeros(newpixsz*newpixsz,nfolds);
        PvalLabel(Pvalues(:,:,vecid)<0.1) = 1;
        PvalLabel(Pvalues(:,:,vecid)<0.05) = 2;
        PvalLabel(Pvalues(:,:,vecid)<0.01) = 3;
        PvalLabel(Pvalues(:,:,vecid)<0.001) = 4;

        % Visualize Results
        for cv = 1:nfolds
            subplot(3,nfolds,(vecid-1)*3+cv)
            h=imagesc(reshape(PvalLabel(:,cv),newpixsz,newpixsz),[0 4]);
            colormap(flipud(gray))
            axis off
            title([AllKernelNames{BootvecOpt(vecid)} ' fold ' num2str(cv)])
            FigureDefault
        end
    end
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['SignificanceLinearModels.fig']))
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['SignificanceLinearModels.bmp']))

    %%    % Visualize Results
    figure('name',[miceopt{midx} ' TotalR'])
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
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR.fig']))
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR.bmp']))

    % Visualize Results
    figure('name',[miceopt{midx} ' TotalR'])
    h=imagesc(reshape(nanmean(TotalR,2),newpixsz,newpixsz),lims);
    colorbar
    axis square
    colormap(flipud(gray))
    FigureDefault
    axis off
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR_Avg.fig']))
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['TotalR_Avg.bmp']))

    % Visualize Results
    %           Rsq_PerVar = nan(newpixsz*newpixsz,size(DM,1),nfolds);
    %             RExcl_PerVar = nan(newpixsz*newpixsz,size(DM,1),nfolds);
    %             TotalR=nan(newpixsz*newpixsz,nfolds);
    %             ChosenK = nan(newpixsz*newpixsz,1);
    DeltaVar = (repmat(TotalR,[1,1,size(RExcl_PerVar,2)])-permute(RExcl_PerVar,[1,3,2]))./repmat(TotalR,[1,1,size(RExcl_PerVar,2)]); % Calculate difference in explained variance (normalized by total variance explained0
    DeltaVar = permute(DeltaVar,[1,3,2]);

    for kernelid=1:length(AllKernels)
        try
            figure('name',[miceopt{midx} ' ' AllKernelNames{kernelid}])
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

    figure('name',[miceopt{midx} ' KernelRepresentations'])
    for kernelid=1:length(AllKernels)

        colormap(flipud(gray))
        HRPDelta = squeeze(nanmean(DeltaVar(:,kernelid,:),3)); % Take the maximum difference in explained variance
        lims2 = [quantile(HRPDelta(:),0.05) quantile(HRPDelta(:),0.95)];
        if ~any(~isnan(lims2))
            continue
        end
        subplot(ceil(sqrt(length(AllKernels))),round(sqrt(length(AllKernels))),kernelid)
        h=imagesc(reshape(HRPDelta,newpixsz,newpixsz),lims2);
        colorbar
        axis square

        axis off
        title(AllKernelNames{kernelid})
    end

    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels.fig']))
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels.bmp']))

    load(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']))
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
            tmpPred = reshape(PredicteddFF(pixidx,:),length(pixidx),length(timevec),ntrials);
            tmpActual = reshape(DatPerPix(pixidx,:),length(pixidx),length(timevec),ntrials);

            % Calculate UniqueExplainedVariancePerTP
            tmpPrdFF = reshape(PredicteddFF(pixidx,:),length(pixidx),length(timevec),ntrials);
            TotalR = cell2mat(arrayfun(@(X) corr(squeeze(nanmean(tmpPrdFF(:,X,:),1)),squeeze(nanmean(tmpActual(:,X,:),1)),'rows','complete').^2,1:length(timevec),'Uni',0));

            % Draw Left Stimulus OptoOff
            subplot(length(AREAS),2,(areaid-1)*2+hemid)
            hold on

            h(1)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==1&OptoVec==1),1),3)),'color',cols(1,:),'LineStyle',linstls{1});
            plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==1&OptoVec==1),1),3)),'color',cols(1,:),'LineStyle',linstls{2})

            % Right stimulus OptoOff
            h(2)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==2&OptoVec==1),1),3)),'color',cols(2,:),'LineStyle',linstls{1});
            plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==2&OptoVec==1),1),3)),'color',cols(2,:),'LineStyle',linstls{2})

            % Left Opto On
            h(3)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==1&OptoVec==2),1),3)),'color',cols(3,:),'LineStyle',linstls{1});
            plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==1&OptoVec==2),1),3)),'color',cols(3,:),'LineStyle',linstls{2})

            % Right Opto On
            h(4)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==2&OptoVec==2),1),3)),'color',cols(4,:),'LineStyle',linstls{1});
            plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==2&OptoVec==2),1),3)),'color',cols(4,:),'LineStyle',linstls{2})


            %             title(AREAGROUPNAMES{areaid})
            if hemid==1
                ylabel([AREAGROUPNAMES{areaid}])
            end
            if areaid==1
                title(HemOpt{hemid})
            end
            FigureDefault
            xlim([-600 timevec(end)])
        end
    end
    legend([h(:)],{'Left Opto Off','Right Opto Off','Left Opto On','Right Opto Ons'})
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_ModelPerformance.fig']))
    saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_ModelPerformance.bmp']))


    % Extract prediction for different timings and areas
    kernelidx2excl = 1;
    MSE = nan(length(timevec),length(AllKernels),length(AREAS),2);
    UniqEV = MSE; % Difference between totalR and leaveoutR
    clear h
    for kernelid=1:length(AllKernels)
        if any(ismember(AllKernels{kernelid},kernelidx2excl))
            continue
        end

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
                tmpPred = reshape(PredicteddFF_LeaveOut(pixidx,:,kernelid),length(pixidx),length(timevec),ntrials);
                tmpActual = reshape(DatPerPix(pixidx,:),length(pixidx),length(timevec),ntrials);

                % Calculate UniqueExplainedVariancePerTP
                MSE(:,kernelid,areaid,hemid) = nanmean(squeeze(nanmean(tmpPred,1)-nanmean(tmpActual,1)).^2,2);
                tmpPrdFF = reshape(PredicteddFF(pixidx,:),length(pixidx),length(timevec),ntrials);
                TotalR = cell2mat(arrayfun(@(X) corr(squeeze(nanmean(tmpPrdFF(:,X,:),1)),squeeze(nanmean(tmpActual(:,X,:),1)),'rows','complete').^2,1:length(timevec),'Uni',0));
                UniqEV(:,kernelid,areaid,hemid) = TotalR-cell2mat(arrayfun(@(X) corr(squeeze(nanmean(tmpPred(:,X,:),1)),squeeze(nanmean(tmpActual(:,X,:),1)),'rows','complete').^2,1:length(timevec),'Uni',0));

                % Draw Left Stimulus Rewarded
                subplot(length(AREAS),2,(areaid-1)*2+hemid)
                hold on
                h(1)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==1&OptoVec==1),1),3)),'color',cols(1,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==1&OptoVec==1),1),3)),'color',cols(1,:),'LineStyle',linstls{2})

                % Right stimulus Rewarded
                h(2)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==2&OptoVec==1),1),3)),'color',cols(2,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==2&OptoVec==1),1),3)),'color',cols(2,:),'LineStyle',linstls{2})

                % Left not rewarded
                h(3)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==1&OptoVec==2),1),3)),'color',cols(3,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==1&OptoVec==2),1),3)),'color',cols(3,:),'LineStyle',linstls{2})

                % Right stimulus not Rewarded
                h(4)=plot(timevec,squeeze(nanmean(nanmean(tmpActual(:,:,SideVec==2&OptoVec==2),1),3)),'color',cols(4,:),'LineStyle',linstls{1});
                plot(timevec,squeeze(nanmean(nanmean(tmpPred(:,:,SideVec==2&OptoVec==2),1),3)),'color',cols(4,:),'LineStyle',linstls{2})


                %             title(AREAGROUPNAMES{areaid})
                if hemid==1
                    ylabel([AREAGROUPNAMES{areaid}])
                end
                if areaid==1
                    title(HemOpt{hemid})
                end
                FigureDefault

                xlim([-600 timevec(end)])
            end
        end
        legend([h(:)],{'Left Opto Off','Right Opto Off','Left Opto On','Right Opto Ons'})
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' AllKernelNames{kernelid} '.fig']))
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['LinearModels_' AllKernelNames{kernelid} '.bmp']))


    end

    %     save(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']),'RShuffle_PerVar_Boot','Rsq_PerVar','RExcl_PerVar','TotalR','ChosenK','PredicteddFF','PredicteddFF_LeaveOut','DatPerPix','ntrials','MSE','UniqEV')

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
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['ExplainedKernelVariance_LinearModels_' AREAGROUPNAMES{areaid} '.fig']))
        saveas(gcf,fullfile(LocalFolder,miceopt{midx},['ExplainedKernelVariance_LinearModels_' AREAGROUPNAMES{areaid} '.bmp']))
    end
end

%% Summarize - copied from WMOpto
MSEAllMice = nan(length(timevec),length(AllKernels),length(AREAS),2,length(miceopt));
UniqEVAllMice = MSEAllMice; % Difference between totalR and leaveoutR
%     RExcl_PerVar_Shuffled
RtotalAcrossMice = nan(length(AREAS),2,length(miceopt));
RShuffleAcrossMice = nan(length(AllKernels),nshuffle,length(AREAS),2,length(miceopt));
for midx = 1:length(miceopt)
    % Load LM Data
    tmp = matfile(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']));

    % Extract all the fields
    disp('Extracting data for this mouse')
    s=fieldnames(tmp);
    for id = 1:length(s)
        eval([s{id} '=tmp.' s{id} ';'])
    end
    TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
    % get brain model
    model = TMPDat.Model;

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
end

%% Shuffle versus total R
figure('name','ShuffleVersusRDistribution')
plevels = nan(length(vecvector),length(AREAGROUPNAMES),2*length(miceopt));
for kernelid=1:length(vecvector)
    subplot(1,length(vecvector),kernelid)
    %                         RtotalAcrossMice = nan(length(AREAS),2,length(miceopt));
    %     RShuffleAcrossMice = nan(length(AllKernels),nshuffle,length(AREAS),2,length(miceopt));fleAcrossMice = nan(length(timevec),length(AllKernels),nshuffle,length(AREAS),2,length(miceopt));
    tmpRtotal = reshape(RtotalAcrossMice,length(AREAS),[]);
    tmpShuff = permute(reshape(RShuffleAcrossMice(vecvector(kernelid),:,:,:,:),nshuffle,length(AREAS),[]),[2,1,3]);

    tmpShuff = repmat(tmpRtotal,1,1,nshuffle)-permute(tmpShuff,[1,3,2]);
    idx = find(~any(isnan(squeeze(nanmean(tmpShuff,3))),1)); %find non nan
    strnames = repmat(string(miceopt),2,1);
    strnames = strnames(idx);
    h=Groupedviolinplot(permute(tmpShuff(:,idx,:),[3,1,2]).*100,string(AREAGROUPNAMES),strnames);
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
    ylim([-15 35])
    for areaid = 1:length(AREAS)
        for midx = 1:length(miceopt)*2
            plevels(kernelid,areaid,midx) = invprctile(squeeze(tmpShuff(areaid,midx,:)),0)./100;
            if quantile(tmpShuff(areaid,midx,:),0.001)>0
                text((areaid-1)*10+midx,35-midx*0.075,'***')
            elseif quantile(tmpShuff(areaid,midx,:),0.01)>0
                text((areaid-1)*10+midx,35-midx*0.075,'**')
            elseif quantile(tmpShuff(areaid,midx,:),0.05)>0
                text((areaid-1)*10+midx,35-midx*0.075,'*')
            end
        end
    end

    FigureDefault
    ylabel('\DeltaR^2')
    title(AllKernelNames{vecvector(kernelid)})

end
idx2 = find(squeeze(any(any(plevels<0.05,1),2)));
idx = idx(ismember(idx,idx2));
saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels_Shuffled.fig']))
saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels_Shuffled.bmp']))

figure;
bar(sum(plevels<0.05,3)./size(plevels,3)*100,'EdgeColor','none')
ylim([0 100])
ylabel('significant \DeltaR^2 (%)')
set(gca,'XTickLabel',AllKernelNames(vecvector))
FigureDefault
%% Time courses (zoomed in)
UniqEVAllMice = reshape(UniqEVAllMice,length(timevec),length(AllKernelNames),length(AREAS),[]); % Take area as N

cols = [1 0 0; 0 0 1; 0 1 0; 0 0 0; 0.5 0.5 0.5];
clear h
figure('name',['VarianceExplainedDifferentKernels'])
for areaid=1:length(AREAS)
    for twid=1:length(TW)
        subplot(length(AREAGROUPNAMES),length(TW),(areaid-1)*length(TW)+twid)
        hold on
        for kernelid=vecvector
            h(kernelid) = shadedErrorBar(timevec(timevec>=TW{twid}(1)&timevec<=TW{twid}(2))./1000,nanmean(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),kernelid,areaid,idx).*100,4),nanstd(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),kernelid,areaid,idx).*100,[],4)./sqrt(length(idx)-1),{'-','color',cols(kernelid-min(vecvector)+1,:)},0);
        end
        FigureDefault
        xlabel('Time (s)')
        if twid==1
            ylabel(AREAGROUPNAMES{areaid})
        else
            ylabel('\DeltaR^2 (%)')
        end
        if areaid==1
            title(TimeWindowNames{twid})
        end
        if ismember(twid,[1,2])
            ylim([-2 12])
        elseif ismember(twid,[3])
            ylim([-1 12])
        else
            ylim([-1 4])
        end
        xlim([-inf inf])
        line([0 0],get(gca,'ylim'),'LineStyle','--','color',[0 0 0])
        line([1.5 1.5],get(gca,'ylim'),'LineStyle','--','color',[0 0 0])

    end
end
legend([h(:).mainLine],AllKernelNames(vecvector))
saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels.fig']))
saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels.bmp']))

%%
UniqEVAllMice = reshape(UniqEVAllMice,length(timevec),length(AllKernels),length(AREAS),2,length(miceopt));
SumperTW = nan(length(TW),length(AllKernels),length(AREAS),2,length(miceopt));
for twid=1:length(TW)
    SumperTW(twid,:,:,:,:) = squeeze(nanmean(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:,:,:,:),1));
end

%%
kernelidx = vecvector;
figure('name',['barplots'])
% identify missing values; these hemispheres cannot be included for
% repeated measures
%     kernelidx = [2,3];
cols = colorcube(length(idx)+1);
cols(end,:)=[];
for twid=1:length(TW)
    subplot(1,length(TW),twid)
    tmp = reshape(SumperTW(twid,kernelidx,:,:,:).*100,length(kernelidx),length(AREAS),[]);
    h=bar(nanmean(tmp(:,:,idx),3));
    hold on
    for areaid=1:length(AREAS)
        h(areaid).EdgeColor = 'none';
        for kernelid=1:length(kernelidx)
            %         for id = 1:length(idx)
            %             plot(h(areaid).XEndPoints',squeeze(tmp(:,areaid,idx(id))),'.-','color',cols(idx(id),:))
            scatter(repmat(h(areaid).XEndPoints(kernelid),1,length(idx)),tmp(kernelid,areaid,idx),10,cols,'filled')
        end
    end
    set(gca,'XTick',1:length(kernelidx),'XTickLabel',AllKernelNames(kernelidx))
    ylabel(['\DeltaR^2 %'])
    FigureDefault

    % Statistics across hemispheres
    Ytmp = permute(tmp,[3,2,1]);
    tmpstr = arrayfun(@(X) arrayfun(@(Y) [AllKernelNames{X} '_AREA' num2str(Y)],1:length(AREAS),'Uni',0),kernelidx,'UniformOutput',0);
    tmpstr = [tmpstr{:}];
    %
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


    writetable(ttmp,fullfile(LocalFolder,['TABLE_' TimeWindowNames{twid} '.csv']))
end
legend(h,AREAGROUPNAMES)
linkaxes
saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BarPlot.fig']))
saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BarPlot.bmp']))

%% Test whether these are significantly more than 0 across mice
p = nan(length(TW),length(kernelidx),length(AREAS));
for twid=1:length(TW)
    for kernelid=1:length(kernelidx)
        for areaid=1:length(AREAS)
            %                for hemid=1:2
            %            SumperTW = nan(length(TW),length(AllKernels),length(AREAS),2,length(miceopt));
            tmp = reshape(squeeze(SumperTW(twid,kernelidx(kernelid),areaid,:,:)),1,[]).*100;
            h = kstest(tmp);
            if ~h
                [~,p(twid,kernelid,areaid)] = ttest(tmp,0,'tail','right');
            else
                p(twid,kernelid,areaid) = signrank(tmp,0,'tail','right');
            end
            %                end
        end
        p(twid,kernelid,:) = bonf_holm(p(twid,kernelid,:));

    end

end
%    p = bonf_holm(p);
% Information for figure 1B
for twid=1:length(TW)
    for areaid = 1:length(AREAS)
        for kernelid = 1:length(kernelidx)
            if any(p(twid,kernelid,areaid)<0.05)
                disp([TimeWindowNames{twid}  ' ' AREAGROUPNAMES{areaid} ' ' AllKernelNames{kernelidx(kernelid)} ' significant'])
            end
        end
    end
end


%% Summarize
MSEAllMice = nan(length(timevec),length(AllKernels),length(AREAS),2,length(miceopt));
UniqEVAllMice = MSEAllMice; % Difference between totalR and leaveoutR
for midx = 1:length(miceopt)
    % Load LM Data
    tmp = load(fullfile(LocalFolder,miceopt{midx},['LinearModels.mat']));

    MSEAllMice(:,:,:,:,midx) = tmp.MSE;
    UniqEVAllMice(:,:,:,:,midx) = tmp.UniqEV;% Difference between totalR and leaveoutR

end

cols = distinguishable_colors(length(AllKernels));
clear h
for areaid=1:length(AREAS)
    figure('name',['VarianceExplainedDifferentKernels ' AREAGROUPNAMES{areaid}])

    for hemid=1:2
        subplot(2,1,hemid)
        hold on
        for kernelid=1:length(AllKernels)
            h(kernelid) = shadedErrorBar(timevec./1000,nanmean(UniqEVAllMice(:,kernelid,areaid,hemid,:).*100,5),nanstd(UniqEVAllMice(:,kernelid,areaid,hemid,:).*100,[],5)./sqrt(length(miceopt)),{'-','color',cols(kernelid,:)},1);
        end
        FigureDefault
        xlabel('Time (s)')
        ylabel('\DeltaR^2 (%)')

        title(HemOpt{hemid})
    end
    legend([h(:).mainLine],AllKernelNames)
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels_' AREAGROUPNAMES{areaid} '.fig']))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_LinearModels_' AREAGROUPNAMES{areaid} '.bmp']))


end
%%
SumperTW = nan(length(TW),length(AllKernels),length(AREAS),2,length(miceopt));
for twid=1:length(TW)
    SumperTW(twid,:,:,:,:) = squeeze(nanmax(UniqEVAllMice(timevec>=TW{twid}(1)&timevec<=TW{twid}(2),:,:,:,:),[],1));
end


%%
kernelidx = [1,2,3,4,5];
hemid = 1;
for twid=2:length(TW)
    %
    figure('name',['Violins ' TimeWindowNames{twid}])
    h=Groupedviolinplot(permute(reshape(SumperTW(twid,kernelidx,:,hemid,:).*100,length(kernelidx),length(AREAS),[]),[3,1,2]),string(AllKernelNames(kernelidx)),string(AREAGROUPNAMES));

    ylabel(['\DeltaR^2 %'])
    FigureDefault

    % Statistics across hemispheres
    Ytmp = permute(reshape(squeeze(SumperTW(twid,kernelidx,:,:,:)),length(kernelidx),length(AREAS),[]),[3,1,2]);
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
        title([TimeWindowNames{twid} ', Intercept: p=' num2str(round(ranovatbl.pValue(1)*1000)./1000) ', Kernel: p=' num2str(round(ranovatbl.pValue(3)*1000)./1000) ', Area: p=' num2str(round(ranovatbl.pValue(5)*1000)./1000) ', Interaction: p=' num2str(round(ranovatbl.pValue(7)*1000)./1000)])
    catch
    end

    % see: help RepeatedMeasuresModel/multcompare
    multtbl = multcompare(Mdl,'Kernel','By','Area');
    [Sigid,id1,id2] = unique(multtbl.Area(multtbl.pValue<0.05));


    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BoxPlot_' TimeWindowNames{twid} '.fig']))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BoxPlot_' TimeWindowNames{twid} '.bmp']))
end

%%
%     kernelidx = [2,3];
for twid=2:length(TW)
    %
    figure('name',['barplots ' TimeWindowNames{twid}])
    tmp = reshape(SumperTW(twid,kernelidx,:,hemid,:).*100,length(kernelidx),length(AREAS),[]);
    h=barwitherr(nanstd(tmp,[],3)./sqrt(size(tmp,3)),nanmean(tmp,3));
    set(gca,'XTick',1:length(kernelidx),'XTickLabel',AllKernelNames(kernelidx))
    legend(h,AREAGROUPNAMES)
    ylabel(['\DeltaR^2 %'])
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
        title([TimeWindowNames{twid} ', Intercept: p=' num2str(round(ranovatbl.pValue(1)*1000)./1000) ', Kernel: p=' num2str(round(ranovatbl.pValue(3)*1000)./1000) ', Area: p=' num2str(round(ranovatbl.pValue(5)*1000)./1000) ', Interaction: p=' num2str(round(ranovatbl.pValue(7)*1000)./1000)])
    catch
    end



    writetable(ttmp,fullfile(LocalFolder,['TABLE_' TimeWindowNames{twid} '.csv']))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BarPlot_' TimeWindowNames{twid} '.fig']))
    saveas(gcf,fullfile(LocalFolder,['ExplainedKernelVariance_BarPlot_' TimeWindowNames{twid} '.bmp']))
end




%% Figure 3C/S2D
maxtrnr = 250;
for midx=1:length(miceopt)

    TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));

    % get model
    model = TMPDat.Model;

    % Get dFF
    disp('Loading in data, which will take a bit of time..')
    tmptrials = TMPDat.trialidentifier;
    tmptrialsall = [tmptrials{:}];
    tmpses = unique(tmptrialsall(2,:));

    % AllTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(ReactionOpt),length(SideOpt),length(HemOpt),length(OptoOpt));
    if midx==1
        timeline = TMPDat.timeline;
        AllMiceTC = cell(length(OptoOpt),2,2,length(AREAS),2,length(miceopt),maxtrnr);
        timevec=timeline;
        TimeLineHere = timeline;
        SR = nanmedian(diff(timevec));
    end
    figure('name',['LickDistribution ' miceopt{midx}])
    SideVec = nan(1,0);
    RewardVec = nan(1,0);
    for ridx=1%:2
        for sidx=1:size(AllDat,1)
            disp(['getting data for ridx = ' num2str(ridx) ', sidx=' num2str(sidx) '...'])
            OptoVec = nan(1,0);
            LeftActionVec = nan(length(timevec),0);
            RightActionVec = LeftActionVec;
            DifMotionVec = LeftActionVec;
            EyeDataVec = LeftActionVec;
            for optidx=1:2
                for sesid = 1:length(tmpses)
                    % Get correct log file
                    Log = TMPDat.LogsThisMouse;% matfile(fullfile(tmpdatFiles(tmpses(sesid)).folder,tmpdatFiles(tmpses(sesid)).name));
                    Log = Log{sesid};
                    %                     TimeLineHere = Log.timeline;
                    %                     Log = Log.Log;
                    takethesetrials = tmptrials{sidx,ridx,optidx}(1,tmptrials{sidx,ridx,optidx}(2,:)==tmpses(sesid));

                    Gavepassive = Log.Gavepassive(takethesetrials);


                    side = Log.Side(takethesetrials);
                    reaction = Log.Reaction(takethesetrials);
                    opto = Log.Opto(takethesetrials);

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
                    if length(Log.RTrightVec)~=length(Log.Opto)
                        keyboard
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
                    OptoVec = [OptoVec repmat(optidx,1,size(tmplicksLeft,2))];

                end
            end

            % Inspect lick distributions
            subplot(2,2,(sidx-1)*2+1)
            % Get equal lick distributions
            TrialIdx = StratifyTrials(timevec(1):100:timevec(end),LeftActionVec,OptoVec);
            title([AnalysisParameters.SideOpt{sidx} ' LeftLicks'])
            FigureDefault

            subplot(2,2,(sidx-1)*2+2)
            % Get equal lick distributions
            TrialIdx = StratifyTrials(timevec(1):100:timevec(end),RightActionVec,OptoVec);
            title([AnalysisParameters.SideOpt{sidx} ' RightLicks'])
            FigureDefault
        end
    end
    drawnow
end

%% Across mice - save out time course per area, keep trial information
if exist(fullfile(LocalFolder,'AllMiceTC.mat'))
    load(fullfile(LocalFolder,'AllMiceTC.mat'))
    % AllTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(ReactionOpt),length(SideOpt),length(HemOpt),length(OptoOpt));
    TMPDat = matfile(fullfile(LocalFolder,miceopt{1},'ProcessedData.mat'));
    timeline = TMPDat.timeline;
else
    for midx=1:length(miceopt)

        TMPDat = matfile(fullfile(LocalFolder,miceopt{midx},'ProcessedData.mat'));
        % get model
        model = TMPDat.Model;

        % Get dFF
        disp('Loading in data, which will take a bit of time..')
        tmptrials = TMPDat.trialidentifier;

        % AllTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(ReactionOpt),length(SideOpt),length(HemOpt),length(OptoOpt));
        if midx==1
            timeline = TMPDat.timeline;
            AllMiceTC = cell(length(OptoOpt),2,2,length(AREAS),2,length(miceopt));
            timevec=timeline;
            SR = nanmedian(diff(timevec));

        end

        % Average Across Grouped AreasTake
        for areaid=1:length(AREAS)
            for hemid=1:2
                % Find area pixels
                idx = find(ismember(model.Rnames,AREAS{areaid}));
                if isempty(idx)
                    idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),model.Rnames,'UniformOutput',0)));
                    areamask = false(400,400);
                    for i = 1:length(idx)
                        areamask(logical(round(model.Regions{idx(i)}))')=true;
                    end
                else
                    areamask = logical(any(cat(3,model.Regions{idx}),3))';
                end
                if hemid==2 %left hemisphere
                    areamask(1:round(model.Lambda(2)),:)=0;
                else
                    areamask(round(model.Lambda(2)):end,:)=0;
                end

                pixidx = find(areamask);
                %Extract data
                tmp = squeeze(AllDat(:,ismember(ReactionOpt,{'Hit','Miss'}),:));
                %reshape to index with pixidx
                tmp = cellfun(@(X) reshape(X,400*400,size(X,3),size(X,4)),tmp,'UniformOutput',0);
                %Average over the pixels
                tmp = cellfun(@(X) squeeze(nanmean(X(pixidx,:,:),1)),tmp,'UniformOutput',0);

                % Save out
                AllMiceTC(:,:,:,areaid,hemid,midx) = permute(tmp,[3,1,2]); %First opto, then side, the nreaction
            end
        end
        save(fullfile(LocalFolder,'AllMiceTC.mat'),'AllMiceTC')
    end
end
%% Check if response to stimulus increases over time
clear h
maxnrt = 700;
% optiidvec = {1, [2,3,4]};
figure('name',['StimulusResponseOverTrials '])
for midx=1:length(miceopt)
    subplot(ceil(sqrt(length(miceopt))),round(sqrt(length(miceopt))),midx)
    hold on
    for areaid=1:length(AREAGROUPNAMES)
        for hemid=1:2


            tmp = AllMiceTC{1,1,1,areaid,hemid,midx}; %only hits
            scatter(1:size(tmp,2),nanmean(tmp(timeline>0&timeline<500,:),1),10,repmat(areaid,1,size(tmp,2)),'filled')

        end
    end

    xlabel('Trial')
    ylabel('VisualResponse')
    title(miceopt{midx})
    makepretty
    offsetAxes
    drawnow



end
saveas(gcf,fullfile(LocalFolder,['StimulusResponseOverTrials.fig']))


%% Time courses baseline (no opto)
MouseGenOpt = unique(mousegen);
linestl = {'-','--'};
cols = [0 0 0; 0.25 0.25 0.25];
nrtri= 280;
FF = figure('name','AcrossMice_Control');

clear h
miny=[];
maxy=[];
for hemid = 1:length(HemOpt)


    for areaid = 1:length(AREAS)
        %     MouseGenOpt = unique(mousegen);
        legendnames = {};
        count=1;

        tmp2ana = nan(length(timeline),1,2,length(mousegen),nrtri);

        subplot(length(AREAS),2,(areaid-1)*2+hemid)
        for OPTOid = 1%:2%length(OptoOpt) %

            % Per Area
            for sidx = 1:2 %Contra/Ipsi (to visual stimulus)
                tmp = squeeze(AllMiceTC(OPTOid,sidx,1,areaid,hemid,:));
                for mid = 1:size(tmp,1)
                    % Z-score
                    tmpz = (tmp{mid}-nanmean(tmp{mid}(:)))./nanstd(tmp{mid}(:));
                    tmp{mid}(abs(tmpz)>3)=nan;
                    if size(tmp{mid},2)>nrtri
                        warning('more trials than anticipated')
                        keyboard
                    end
                    tmp2ana(:,OPTOid,sidx,mid,1:size(tmp{mid},2)) = tmp{mid};
                end
                tmp = squeeze(reshape(tmp2ana(:,OPTOid,sidx,:,:),length(timeline),[],size(tmp,1)*nrtri));
                tmp(:,sum(isnan(tmp),1)==length(timeline)) = [];
                tmpmean = nanmean(tmp,2);
                tmpsem = nanstd(tmp,[],2)./sqrt(size(tmp,2)-1);
                hold on
                h{count} = shadedErrorBar(timeline./1000,tmpmean,tmpsem,{'LineWidth',2,'color',cols(sidx,:),'LineStyle',linestl{1}},0);
                box off
                drawnow

                xlim([timeline(1)./1000+0.200 timeline(end)./1000])
                count=count+1;
            end
        end

        if any(~isnan(tmp(:)))
            g1 = repmat([1]',[1,2,size(tmp2ana,4),size(tmp2ana,5)]);
            g2 = repmat([1,2]',[1,1,size(tmp2ana,4),size(tmp2ana,5)]);
            g2 = permute(g2,[2,1,3,4]);
            g3 = repmat(1:size(tmp2ana,4),[1,1,2,size(tmp2ana,5)]);
            g3 = permute(g3,[2,3,1,4]);
            pvals = nan(1,length(timeline));
            for tid=1:length(timeline)
                tmp2 = tmp2ana(tid,:,:,:,:);
                idx = ~isnan(tmp2);
                tbl = table(double(tmp2(idx)),g1(idx),g2(idx),g3(idx),'VariableNames',{'dFF','Opto','Side','MOUSEID'});
                try
                    glme = fitglme(tbl,'dFF ~ 1 + Side + (1|MOUSEID)','Distribution','Normal');
                    tmp = anova(glme);
                    pvals(:,tid) = tmp.pValue(2:end);
                catch ME
                    disp(ME)
                end
            end

            ytmp = get(gca,'ylim');
            widthtmp = 0.05*(max(ytmp)-min(ytmp));

            for pid = find(pvals(1,:)<0.05)
                plot(timeline(pid)./1000,max(ytmp)-widthtmp*2,'.','color',[0 0 0],'MarkerSize',8)
            end
            if any(pvals(1,:)<0.05)
                disp([AREAGROUPNAMES{areaid} ' ' HemOpt{hemid} ' has Side effects!: ' num2str(timeline(find(pvals(1,:)<0.05)))])
            end


        end


        xlim([-0.600 3.5])

        line([-0.500 -0.500],ytmp,'color',[0 0 0],'LineStyle','--')
        line([1.500 1.500],ytmp,'color',[0 0 0],'LineStyle','--')
        line([0 0],ytmp,'color',[0 0 0],'LineStyle','--')
        ylabel(['\DeltaF/F'])

        %         if areaid==1
        %             title([HemOpt{hemid} ' ' AREAGROUPNAMES{areaid}])
        %         else
        %             title(AREAGROUPNAMES{areaid})
        %         end
        if areaid < length(AREAGROUPNAMES)
            set(gca,'XTickLabel',[])
        end

        makepretty
        %         offsetAxes

    end
end

h=cellfun(@(X) X.mainLine,h,'UniformOutput',0);

%
saveas(gcf,fullfile(LocalFolder,['AcrossMice_' [HemOpt{hemid}] 'Baseline.fig']))
saveas(gcf,fullfile(LocalFolder,['AcrossMice_'  [HemOpt{hemid}] 'Baseline.bmp']))


%% Time courses Fig 4/S2
MouseGenOpt = unique(mousegen);
linestl = {'-','--'};
nrtri= 280;
for hemid = 1:length(HemOpt)
    % Per Area
    for sidx = 1:2 %Contra/Ipsi (to visual stimulus)
        if sidx == 1
            cols = [0 0 0; 1 0 0; 0 0 1];

        else
            cols = [0.5 0.5 0.5; 1 0.2 0.2; 0.2 0.2 1];
        end
        FF(hemid) = figure('name',['AcrossMice_' [HemOpt{hemid}] ' ' SideOptNames{sidx}]);

        clear h
        miny=[];
        maxy=[];
        for areaid = 1:length(AREAS)
            %     MouseGenOpt = unique(mousegen);
            for gidx=1:length(MouseGenOpt)
                legendnames = {};
                count=1;

                tmp2ana = nan(length(timeline),2,sum(ismember(mousegen,MouseGenOpt{gidx})),nrtri);

                subplot(length(AREAS),2,(areaid-1)*2+gidx)
                for OPTOid = 1:2%length(OptoOpt) %
                    tmp = squeeze(AllMiceTC(OPTOid,sidx,1,areaid,hemid,ismember(mousegen,MouseGenOpt{gidx})));
                    for mid = 1:size(tmp,1)
                        % Z-score
                        tmpz = (tmp{mid}-nanmean(tmp{mid}(:)))./nanstd(tmp{mid}(:));
                        tmp{mid}(abs(tmpz)>3)=nan;
                        if size(tmp{mid},2)>nrtri
                            warning('more trials than anticipated')
                            keyboard
                        end
                        tmp2ana(:,OPTOid,mid,1:size(tmp{mid},2)) = tmp{mid};
                    end
                    tmp = squeeze(reshape(tmp2ana(:,OPTOid,:,:),length(timeline),[],size(tmp,1)*nrtri));
                    tmp(:,sum(isnan(tmp),1)==length(timeline)) = [];
                    %                                         tmp = arrayfun(@(k) cat(2,tmp{:,k}),1:size(tmp,2),'uniformoutput',0); %Concatenate over reactions

                    %                     tmp = cat(2,tmp{:});
                    if sum(ismember(mousegen,MouseGenOpt{gidx}))>1
                        tmpmean = nanmean(tmp,2);
                        tmpsem = nanstd(tmp,[],2)./sqrt(size(tmp,2)-1);
                        hold on
                        h{count} = shadedErrorBar(timeline./1000,tmpmean,tmpsem,{'LineWidth',2,'color',cols(1+(OPTOid-1).*gidx,:),'LineStyle',linestl{1}},0);
                    else
                        hold on
                        h{count} = plot(timeline./1000,tmp,'LineWidth',2,'color',cols(OPTOid,:),'LineStyle',linestl{1});
                    end
                    box off
                    drawnow

                    xlim([timeline(1)./1000+0.200 timeline(end)./1000])
                    %                     legendnames = {legendnames{:} [LegendInput{1} OptoOpt{OPTOid}]};
                    count=count+1;
                end

                if any(~isnan(tmp(:)))
                    g1=repmat([1,2]',[1,size(tmp2ana,3),size(tmp2ana,4)]);
                    g3 = repmat([1:size(tmp2ana,3)]',[1,2,size(tmp2ana,4)]);
                    g3 = permute(g3,[2,1,3]);
                    pvals = nan(1,length(timeline));
                    for tid=1:length(timeline)
                        tmp2 = tmp2ana(tid,:,:,:,:);
                        idx = ~isnan(tmp2);
                        tbl = table(double(tmp2(idx)),g1(idx),g3(idx),'VariableNames',{'dFF','Opto','MOUSEID'});
                        try
                            glme = fitglme(tbl,'dFF ~ 1 + Opto + (1|MOUSEID)','Distribution','Normal');
                            tmp = anova(glme);
                            pvals(:,tid) = tmp.pValue(2:end);
                        catch ME
                            disp(ME)
                        end
                    end
                    ytmp = get(gca,'ylim');
                    widthtmp = 0.05*(max(ytmp)-min(ytmp));
                    for pid = find(pvals(1,:)<0.05)
                        plot(timeline(pid)./1000,max(ytmp)-widthtmp,'.','color',[1 0 0],'MarkerSize',8)
                    end
                    if any(pvals(1,:)<0.05)
                        disp([AREAGROUPNAMES{areaid} ' ' HemOpt{hemid} ' ' SideOptNames{sidx} ' ' MouseGenOpt{gidx} ' has Opto effects!: ' num2str(timeline(find(pvals(1,:)<0.05)))])
                    end

                end



                if areaid <length(AREAGROUPNAMES)
                    set(gca,'XTickLabel',[])
                end
                line([-0.500 -0.500],ytmp,'color',[0 0 0],'LineStyle','--')
                line([1.500 1.500],ytmp,'color',[0 0 0],'LineStyle','--')
                line([0 0],ytmp,'color',[0 0 0],'LineStyle','--')
                %                 ylabel(['\DeltaF/F'])

                makepretty
            end

        end
        h=cellfun(@(X) X.mainLine,h,'UniformOutput',0);


        saveas(gcf,fullfile(LocalFolder,['AcrossMice_' [HemOpt{hemid}] '_' SideOptNames{sidx} '.fig']))
        saveas(gcf,fullfile(LocalFolder,['AcrossMice_'  [HemOpt{hemid}] '_' SideOptNames{sidx} '.bmp']))
    end

end

%% Time courses Fig 3/S1
MouseGenOpt = unique(mousegen);
linestl = {'-','--'};
cols = [0 0 0; 1 0 0; 0 0 1];
nrtri= 280;
for hemid = 1:length(HemOpt)

    FF(hemid) = figure('name',['AcrossMice_' [HemOpt{hemid}]]);

    clear h
    miny=[];
    maxy=[];
    for areaid = 1:length(AREAS)
        %     MouseGenOpt = unique(mousegen);
        for gidx=1:length(MouseGenOpt)
            legendnames = {};
            count=1;

            tmp2ana = nan(length(timeline),2,2,sum(ismember(mousegen,MouseGenOpt{gidx})),nrtri);

            subplot(length(AREAS),2,(areaid-1)*2+gidx)
            for OPTOid = 1:2%length(OptoOpt) %

                % Per Area
                for sidx = 1:2 %Contra/Ipsi (to visual stimulus)

                    tmp = squeeze(AllMiceTC(OPTOid,sidx,1,areaid,hemid,ismember(mousegen,MouseGenOpt{gidx})));
                    for mid = 1:size(tmp,1)
                        % Z-score
                        tmpz = (tmp{mid}-nanmean(tmp{mid}(:)))./nanstd(tmp{mid}(:));
                        tmp{mid}(abs(tmpz)>3)=nan;
                        if size(tmp{mid},2)>nrtri
                            warning('more trials than anticipated')
                            keyboard
                        end
                        tmp2ana(:,OPTOid,sidx,mid,1:size(tmp{mid},2)) = tmp{mid};
                    end
                    tmp = squeeze(reshape(tmp2ana(:,OPTOid,sidx,:,:),length(timeline),[],size(tmp,1)*nrtri));
                    tmp(:,sum(isnan(tmp),1)==length(timeline)) = [];
                    %                                         tmp = arrayfun(@(k) cat(2,tmp{:,k}),1:size(tmp,2),'uniformoutput',0); %Concatenate over reactions

                    %                     tmp = cat(2,tmp{:});
                    if sum(ismember(mousegen,MouseGenOpt{gidx}))>1
                        tmpmean = nanmean(tmp,2);
                        tmpsem = nanstd(tmp,[],2)./sqrt(size(tmp,2)-1);
                        hold on
                        h{count} = shadedErrorBar(timeline./1000,tmpmean,tmpsem,{'LineWidth',2,'color',cols(1+(OPTOid-1).*gidx,:),'LineStyle',linestl{sidx}},0);
                    else
                        hold on
                        h{count} = plot(timeline./1000,tmp,'LineWidth',2,'color',cols(OPTOid,:),'LineStyle',linestl{sidx});
                    end
                    box off
                    drawnow

                    xlim([timeline(1)./1000+0.200 timeline(end)./1000])
                    count=count+1;
                end
            end

            if any(~isnan(tmp(:)))
                g1=repmat([1,2]',[1,2,size(tmp2ana,4),size(tmp2ana,5)]);
                g2= repmat([1,2]',[1,2,size(tmp2ana,4),size(tmp2ana,5)]);
                g2 = permute(g2,[2,1,3,4]);
                g3 = repmat(1:size(tmp2ana,4),[1,2,2,size(tmp2ana,5)]);
                g3 = permute(g3,[2,3,1,4]);
                pvals = nan(3,length(timeline));
                for tid=1:length(timeline)
                    tmp2 = tmp2ana(tid,:,:,:,:);
                    idx = ~isnan(tmp2);
                    tbl = table(double(tmp2(idx)),g1(idx),g2(idx),g3(idx),'VariableNames',{'dFF','Opto','Side','MOUSEID'});
                    try
                        glme = fitglme(tbl,'dFF ~ 1 + Opto + Side + Opto:Side + (1|MOUSEID)','Distribution','Normal');
                        tmp = anova(glme);
                        pvals(:,tid) = tmp.pValue(2:end);
                    catch ME
                        disp(ME)
                    end
                end

                ytmp = get(gca,'ylim');
                widthtmp = 0.05*(max(ytmp)-min(ytmp));
                for pid = find(pvals(1,:)<0.05)
                    plot(timeline(pid)./1000,max(ytmp)-widthtmp,'.','color',[1 0 0],'MarkerSize',8)
                end
                if any(pvals(1,:)<0.05)
                    disp([AREAGROUPNAMES{areaid} ' ' HemOpt{hemid} ' ' MouseGenOpt{gidx} ' has Opto effects!: ' num2str(timeline(find(pvals(1,:)<0.05)))])
                end
                for pid = find(pvals(2,:)<0.05)
                    plot(timeline(pid)./1000,max(ytmp)-widthtmp*2,'.','color',[0 0 0],'MarkerSize',8)
                end
                if any(pvals(2,:)<0.05)
                    disp([AREAGROUPNAMES{areaid} ' ' HemOpt{hemid} ' ' MouseGenOpt{gidx} ' has Side effects!: ' num2str(timeline(find(pvals(2,:)<0.05)))])
                end
                for pid = find(pvals(3,:)<0.05)
                    plot(timeline(pid)./1000,max(ytmp)-widthtmp*3,'.','color',[0.5 0.5 0.5],'MarkerSize',8)
                end
                if any(pvals(3,:)<0.05)
                    disp([AREAGROUPNAMES{areaid} ' ' HemOpt{hemid} ' ' MouseGenOpt{gidx} ' has Interaction effects!: ' num2str(timeline(find(pvals(3,:)<0.05)))])
                end

            end



            line([-0.500 -0.500],ytmp,'color',[0 0 0],'LineStyle','--')
            line([1.500 1.500],ytmp,'color',[0 0 0],'LineStyle','--')
            line([0 0],ytmp,'color',[0 0 0],'LineStyle','--')
            ylabel(['\DeltaF/F'])

            if areaid==1
                title([MouseGenOpt{gidx} ', n='  num2str(sum(ismember(mousegen,MouseGenOpt{gidx}))) ', ' AREAGROUPNAMES{areaid}])
            else
                title(AREAGROUPNAMES{areaid})
            end
            FigureDefault

        end
        h=cellfun(@(X) X.mainLine,h,'UniformOutput',0);


    end
    saveas(gcf,fullfile(LocalFolder,['AcrossMice_' [HemOpt{hemid}] '.fig']))
    saveas(gcf,fullfile(LocalFolder,['AcrossMice_'  [HemOpt{hemid}] '.bmp']))

end

%% Time courses V1 Individual Mice
MouseGenOpt = unique(mousegen);
linestl = {'-','--'};
cols = [0 0 0; 1 0 0; 0 0 1];
nrtri= 280;
AREAS = AnalysisParameters.AREAS;
for areaid = 1%:length(AnalysisParameters.AREAS)
    FF(areaid) = figure('name',['EveryMouse_' [AnalysisParameters.AREAS{areaid}]]);

    clear h
    miny=[];
    maxy=[];
    %     MouseGenOpt = unique(mousegen);
    for hemid = 1:length(HemOpt)
        for gidx=2%1:length(MouseGenOpt)
            legendnames = {};

            for mid = 1:sum(ismember(mousegen,MouseGenOpt{gidx}))
                tmp2ana = nan(length(timeline),2,2,nrtri);
                count=1;

                subplot(sum(ismember(mousegen,MouseGenOpt{gidx})),2,(mid-1)*2+hemid)
                for OPTOid = 1:2%length(OptoOpt) %

                    % Per Area
                    for sidx = 1:2 %Contra/Ipsi (to visual stimulus)
                        tmp = squeeze(AllMiceTC(OPTOid,sidx,1,areaid,hemid,ismember(mousegen,MouseGenOpt{gidx})));
                        if size(tmp{mid},2)>nrtri
                            warning('more trials than anticipated')
                            keyboard
                        end
                        tmp2ana(:,OPTOid,sidx,1:size(tmp{mid},2)) = tmp{mid};
                        tmp = squeeze(reshape(tmp2ana(:,OPTOid,sidx,:),length(timeline),nrtri));
                        tmp(:,sum(isnan(tmp),1)==length(timeline)) = [];
                        %                     tmp = arrayfun(@(k) cat(2,tmp{:,k}),1:size(tmp,2),'uniformoutput',0); %Concatenate over reactions
                        if sum(ismember(mousegen,MouseGenOpt{gidx}))>1
                            tmpmean = nanmean(tmp,2);
                            tmpsem = nanstd(tmp,[],2)./sqrt(size(tmp,2)-1);
                            hold on
                            h{count} = shadedErrorBar(timeline,tmpmean,tmpsem,{'LineWidth',2,'color',cols(1+(OPTOid-1).*gidx,:),'LineStyle',linestl{sidx}},1);
                        else
                            hold on
                            h{count} = plot(timeline,tmp,'LineWidth',2,'color',cols(OPTOid,:),'LineStyle',linestl{sidx});
                        end
                        box off
                        drawnow

                        xlim([timeline(1)+200 timeline(end)])
                        count=count+1;
                    end
                end


                if any(~isnan(tmp(:)))
                    g1=repmat([1,2]',[1,2,size(tmp2ana,4)]);
                    g2= repmat([1,2]',[1,2,size(tmp2ana,4)]);
                    g2 = permute(g2,[2,1,3,4]);

                    pvals = nan(3,length(timeline));
                    for tid=1:length(timeline)
                        tmp2 = tmp2ana(tid,:,:,:);
                        idx = ~isnan(tmp2);
                        try
                            pvals(:,tid) = anovan(double(tmp2(idx)),{g1(idx),g2(idx)},'model',3,'VarNames',{'Opto','Side'},'Display','off');

                        catch ME
                            disp(ME)
                        end
                    end
                    ylim([-0.025 0.025])
                    ytmp = get(gca,'ylim');
                    widthtmp = 0.05*(max(ytmp)-min(ytmp));
                    for pid = find(pvals(1,:)<0.05)
                        plot(timeline(pid),max(ytmp)-widthtmp,'.','color',[1 0 0],'MarkerSize',8)
                    end
                    if any(pvals(1,:)<0.05)
                        disp([AREAS{areaid} ' ' HemOpt{hemid} ' ' MouseGenOpt{gidx} ' has Opto effects!: ' num2str(timeline(find(pvals(1,:)<0.05)))])
                    end
                    for pid = find(pvals(2,:)<0.05)
                        plot(timeline(pid),max(ytmp)-widthtmp*2,'.','color',[0 0 0],'MarkerSize',8)
                    end
                    if any(pvals(2,:)<0.05)
                        disp([AREAS{areaid} ' ' HemOpt{hemid} ' ' MouseGenOpt{gidx} ' has Side effects!: ' num2str(timeline(find(pvals(2,:)<0.05)))])
                    end
                    for pid = find(pvals(3,:)<0.05)
                        plot(timeline(pid),max(ytmp)-widthtmp*3,'.','color',[0.5 0.5 0.5],'MarkerSize',8)
                    end
                    if any(pvals(3,:)<0.05)
                        disp([AREAS{areaid} ' ' HemOpt{hemid} ' ' MouseGenOpt{gidx} ' has Interaction effects!: ' num2str(timeline(find(pvals(3,:)<0.05)))])
                    end

                end


                xlim([-600 2000])

                line([-500 -500],ytmp,'color',[0 0 0],'LineStyle','--')
                line([1500 1500],ytmp,'color',[0 0 0],'LineStyle','--')
                line([0 0],ytmp,'color',[0 0 0],'LineStyle','--')
                ylabel(['dFF'])


                title([MouseGenOpt{gidx} ', Mouse='  num2str(mid) ', Hem ' HemOpt{hemid}])
            end
        end
        h=cellfun(@(X) X.mainLine,h,'UniformOutput',0);


    end


end