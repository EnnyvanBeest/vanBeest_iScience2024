%% User Input
Datadir = '\Path\To\Data'
Scriptsdir = '\Path\To\Scripts\' %Path to all the (sub)scripts of this analysis
storepath = '\Path\To\Store\' %Main storage --> Will contain all processed data & model of the brain
miceopt ={'Faraday','Galilei','Haydn','Jemison','Karlsen','Lully','Mahler','Sebald'}
mousegen = {'Drd2','Drd2','Drd1','Drd2','Drd1','Drd1','Drd1','Drd1'}% 

AnalysisParameters.miceopt = miceopt;
AnalysisParameters.mousegen=mousegen;
AnalysisParameters.UsepRFMask = 0; % Use only those pixels that show up in the pRF mask
AnalysisParameters.UseEvokedMask = 0; %Use only pixels taht show up in the Evoked Response (to a figure on grey BG)
AnalysisParameters.WindowSize = 15; %Nr. trials around window
AnalysisParameters.baselinemethod = 1; %1: trial by trial baseline, 2: filtered baseline (F0), 3: weighted average baseline after detrending (over all trials)
AnalysisParameters.takeequalsample = 0;  %Eq nr of trials per condition
AnalysisParameters.dontthrowtrialsgenerally = 0; %Throw out trials based on motion
AnalysisParameters.smoothfact = 2; %Smoothing with neighbours (Gaussian); fill in standard deviation
AnalysisParameters.ScaleFact = 0.5; %between 0-1 is scale down
AnalysisParameters.normalize=0; %normalize to Figure on Grey condition
AnalysisParameters.AREAS = {'V1','Vpor','Vl','Val','Vrl','Va','Vam','Vpm','RSP','M1','M2'};%note: l includes li and l;%% %
AnalysisParameters.ReactionOpt = {'Passive'} %Which response types to include
AnalysisParameters.AREASGrouped = {{'V1'},{'Vl','Vli','Val','Vpor'},{'Vrl','Va'},{'Vam','Vpm'},{'RSP'},{'M1'},{'M2'}};%note: l includes li and l{'V1','HVLat_Vl_Val','HVAnterior_Vrl_Va','HVMedial_Vam_Vpm','RSP','M1','M2'};%% %
AREAGROUPNAMES = {'V1','Vlat','PPC',...
    'Vmed','RS','M1','M2'};

AnalysisParameters.OptoOpt = {'Off','On'};
AnalysisParameters.TimeWindowNames = {'Baseline','OptoEarly','OptoLate','After'};
AnalysisParameters.TW = {[-500 0],[0 500],[500 1000],[1000 1500]};
AnalysisParameters.templatemouse='Jemison'
AnalysisParameters.UseCluster = 0; %Make use of found clusters?

TW = AnalysisParameters.TW;
HemOpt = {'Stimulated','non-Stimulated'};
OPTOOpt = [0,1];
OptoOptNames = {'Off','On'};
MouseGenOpt = {'Drd1','Drd2'};
NewDim = 40;
[optimizer,metric]=imregconfig('Multimodal');

addpath(genpath(Scriptsdir))
global UserQuestions %Defines whether gui's show up to ask what to do with existing datasets. If 0, existing datasets are used and not overwritten
UserQuestions = 0; %put to 1 to do stuff manually

AREAS = AnalysisParameters.AREASGrouped;
%% Figure 2A
maxtr = 200;
if exist(fullfile(Datadir,'AllMiceTC.mat'))
    load(fullfile(Datadir,'AllMiceTC.mat'))
else
    for midx=1:length(miceopt)
        
        TMPDat = matfile(fullfile(Datadir,miceopt{midx},'ProcessedData.mat'));
        AllDat = TMPDat.AllDat;
        
        % get model
        model = load(fullfile(Datadir,miceopt{midx},[miceopt{midx} 'BrainModel.mat']));
                
        if midx==1
            timeline = TMPDat.timeline;
            AllMiceTC = nan(length(AnalysisParameters.AREASGrouped),length(timeline),length(HemOpt),length(OPTOOpt),maxtr,length(miceopt));
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
                
                tmp = AllDat;
                %reshape to index with pixidx
                tmp = cellfun(@(X) reshape(X,400*400,size(X,3),size(X,4)),tmp,'UniformOutput',0);
                %Average over the pixels
                tmp = cellfun(@(X) squeeze(nanmean(X(pixidx,:,:),1)),tmp,'UniformOutput',0);
                %                                 AllMiceTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(HemOpt),length(OPTOOpt),0,length(miceopt));
                for optid=1:2
                    % Save out
                    if size(tmp{optid},2)>maxtr
                        warning('Initialize with more trials!' )
                        keyboard
                    end
                    AllMiceTC(areaid,:,hemid,optid,1:size(tmp{optid},2),midx) = tmp{optid};
                end
            end
        end
        save(fullfile(Datadir,'AllMiceTC.mat'),'AllMiceTC')
    end
end


%% FIgure 2A -- Whole Brain
AvgWholeBrTC = nan(length(timeline),length(OPTOOpt),2,length(miceopt));
for midx=1:length(miceopt)
    for hemid=1:2
        % Load Brain outline
        load(fullfile(Datadir, miceopt{midx}, [miceopt{midx} 'BrainModel.mat']))
        Model = BrainModel;
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
        
        if hemid==1
            % Take only right hemisphere to be fair
            removepix(halfidx+1:end,:)=1;
        else
            removepix(1:halfidx,:)=1; %And left
        end
        % load data
        TMP = load(fullfile(Datadir,miceopt{midx},'ProcessedData.mat'));
        for optid = 1:2
            tmp = TMP.AllDat{optid};
            tmp(repmat(removepix,[1,1,length(timeline),size(tmp,4)]))=nan;
            AvgWholeBrTC(:,optid,hemid,midx)=squeeze(nanmean(nanmean(reshape(tmp,xpix*ypix,length(timeline),[]),1),3));
        end
    end
end
linestl = {'-','--'};
FF(areaid) = figure('name',['AcrossMice_WholeBrain']);
clear h
miny=[];
maxy=[];
showindividualmice=0;
cols = [0 0 0; 1 0 0; 0.5 0.5 0.5; 0 0 1];

%     MouseGenOpt = unique(mousegen);
for hemid=1:2;
    for gidx=1:length(MouseGenOpt)
        legendnames = {};
        count=1;
        subplot(2,length(MouseGenOpt)+1,[(hemid-1)*(length(MouseGenOpt)+1)+gidx])
        for OPTOid = 1:length(OPTOOpt) %
            %         AllTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(HemOpt),length(OPTOOpt),nrtrialspercond);
            
            
            %
            tmp = squeeze(AvgWholeBrTC(:,OPTOid,hemid,ismember(mousegen,MouseGenOpt{gidx})));
            
            if showindividualmice
                for i = 1:size(tmp,2)
                    plot(timeline,tmp(:,i),'LineWidth',1,'color',cols(OPTOid,:)*0.7.^(i-1),'LineStyle',linestl{hemid})
                    hold on
                end
                
                h{count} = plot(timeline,nanmean(tmp,2),'LineWidth',2,'color',cols(OPTOid,:),'LineStyle',linestl{hemid});
            else
                h{count} = shadedErrorBar(timeline,nanmean(tmp,2),nanstd(tmp,[],2)./sqrt(size(tmp,2)-1),{'LineWidth',2,'color',cols((gidx-1)*2+OPTOid,:),'LineStyle',linestl{1}},1);
                hold on
            end
            
            box off
            drawnow
            
            xlim([timeline(1)+100 timeline(end)])
            legendnames = {legendnames{:} ['OPTO=' num2str(OPTOOpt(OPTOid))]};
            count=count+1;
        end
        xlim([-200 2000])
        
        ytmp = get(gca,'ylim');
        widthtmp = 0.05*(max(ytmp)-min(ytmp));
        %     patch([0 Log.OptoDuration*1000 Log.OptoDuration*1000 0],[min(ytmp) min(ytmp) min(ytmp)+widthtmp min(ytmp)+widthtmp],[0 0 0])
        
        ylabel(['dFF'])
        
        
        title([MouseGenOpt{gidx} ', n='  num2str(sum(ismember(mousegen,MouseGenOpt{gidx})))])
        
    end
end
subplot(2,length(MouseGenOpt)+1,length(MouseGenOpt)+1)

hh=imagesc(~(smooth2a(removepix,5)>0.5));
colormap('jet')
hold on
plot(Model.AllX,Model.AllY,'k.','MarkerSize',3)
axis square
axis off
set(hh,'AlphaData',~(smooth2a(removepix,5)>0.5))
subh= subplot(2,length(MouseGenOpt)+1,2*(length(MouseGenOpt)+1));
Pos = get(gca,'Position');
axis square
axis off


if showindividualmice
    legend([h{:}],{legendnames{:}},'Position',[Pos(1) Pos(2) Pos(3) Pos(4)*0.2])
else
    h = cellfun(@(X) X.mainLine,h,'UniformOutput',0);
    legend([h{:}],{legendnames{:}},'Position',[Pos(1) Pos(2) Pos(3) Pos(4)*0.2])
end
saveas(gcf,fullfile(Datadir,['AcrossMice_WholeBrainTC.fig']))
saveas(gcf,fullfile(Datadir,['AcrossMice_WholeBrainTC.bmp']))

%% Figure 2A - per area
load(fullfile(Datadir, AnalysisParameters.templatemouse, [AnalysisParameters.templatemouse 'BrainModel.mat']))
Model = BrainModel;
mask = zeros(800*AnalysisParameters.ScaleFact,800*AnalysisParameters.ScaleFact);
xpix = 800*AnalysisParameters.ScaleFact;
ypix = 800*AnalysisParameters.ScaleFact;
halfidx = round(median(Model.AllY));
cols = cat(3,[0 0 0;1 0 0],[0 0 0;0 0 1]);
linestl = {'-','--'};
FF(areaid) = figure('name',['TCAcrossMice']);

for areaid = 1:length(AREAS)
    clear h
    miny=[];
    maxy=[];
    %     MouseGenOpt = unique(mousegen);
    
    for gidx=1:length(MouseGenOpt)
        legendnames = {};
        count=1;
        subplot(length(AREAS),3,(areaid-1)*3+gidx+1)
        for OPTOid = 1:length(OPTOOpt) %
            %         AllTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(HemOpt),length(OPTOOpt),nrtrialspercond);
            
            
            % Per Area
            for hemid = 1:length(HemOpt)
                
                if strcmp(HemOpt{hemid},'RightHem')
                    LegendInput={'IpsitoOPTO'};
                else
                    LegendInput={'ContratoOPTO'};
                end
                %            AllMiceTC = nan(length(AnalysisParameters.AREAS),length(timeline),length(HemOpt),length(OPTOOpt),maxtr,length(miceopt));
                
                tmp=squeeze(nanmean(AllMiceTC(areaid,:,hemid,OPTOid,:,ismember(mousegen,MouseGenOpt{gidx})),5));
                
                if showindividualmice
                    for i = 1:size(tmp,2)
                        plot(timeline,tmp(:,i),'LineWidth',1,'color',cols(OPTOid,:,gidx)*0.7.^(i-1),'LineStyle',linestl{hemid})
                        hold on
                    end
                    h{count} = plot(timeline,nanmean(tmp,2),'LineWidth',2,'color',cols(OPTOid,:,gidx),'LineStyle',linestl{hemid});
                else
                    h{count} = shadedErrorBar((timeline./1000),nanmean(tmp,2),nanstd(tmp,[],2)./sqrt(size(tmp,2)-1),{'LineWidth',2,'color',cols(OPTOid,:,gidx),'LineStyle',linestl{hemid}},0);
                    hold on
                end
                box off
                drawnow
                
                legendnames = {legendnames{:} [LegendInput{1} ' OPTO=' num2str(OPTOOpt(OPTOid))]};
                count=count+1;
            end
        end
        xlim([-0.200 2])
        ylim([-0.003 0.0045])
        
        ytmp = get(gca,'ylim');
        widthtmp = 0.05*(max(ytmp)-min(ytmp));
        line([0 0],ytmp,'color',[0 0 0],'LineStyle','--')
        line([Log.OptoDuration Log.OptoDuration],ytmp,'color',[0 0 0],'LineStyle','--')
        
        ylabel(['\DeltaF/F'])
        
        if areaid==1
            title([MouseGenOpt{gidx} ', n='  num2str(sum(ismember(mousegen,MouseGenOpt{gidx})))])
        end
        FigureDefault
        
    end
    
    subplot(length(AREAS),3,(areaid-1)*3+1)
    
    
    areamask = false(xpix,ypix);
    idx = find(ismember(Model.Rnames,AREAS{areaid}));
    if isempty(idx) %|| strcmp(AREAS{areaid},'Vl')
        idx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),Model.Rnames,'UniformOutput',0)));
        for i = 1:length(idx)
            areamask(logical(round(Model.Regions{idx(i)}))')=true;
        end
    else
        areamask = logical(any(cat(3,Model.Regions{idx}),3))';
    end
    
    areamask(halfidx+1:end,:) = 0;
    LegendInput={'IpsitoOPTO'};
    ha=imagesc(areamask);
    set(ha,'AlphaData',areamask==1)
    colormap(jet)
    hold on
    plot(Model.AllX,Model.AllY,'k.','MarkerSize',3)
    title(AREAGROUPNAMES{areaid})
    axis square
    axis off
    FigureDefault
    
  
    if ~exist(fullfile(Datadir,'Figures'))
        mkdir(fullfile(Datadir,'Figures'))
    end
end
print(gcf, '-dpdf', fullfile(Datadir,'Figures',['AcrossMiceTC.pdf']))
saveas(gcf,fullfile(Datadir,'Figures',['AcrossMiceTC.fig']))
saveas(gcf,fullfile(Datadir,'Figures',['AcrossMiceTC.bmp']))

% save(fullfile(localdir,'AllMiceTC.mat'),'AllMiceTC','AnalysisParameters','timeline','Model','Log','OPTOOpt','miceopt','mousegen','-v7.3')


%% Align mice
templatemouse = find(strcmp(miceopt,AnalysisParameters.templatemouse));

% Load template mouse first
load(fullfile(Datadir, miceopt{templatemouse}, [miceopt{templatemouse} 'BrainModel.mat']))
Model = BrainModel;

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
TemplateModel = Model;


% collapse onto one hemisphere
nrpix = max([halfidx xpix-halfidx]);
% halfbrain = flipud(brain(nrpix-1:end,:));
halfmask = removepix(1:nrpix,:);
templatepix=removepix;
templatepixnd = Nan_imresize(templatepix,NewDim./size(templatepix,1));

orimiceopt = miceopt;
if exist(fullfile(Datadir,'AlignToTemplateWithAllenFrame.mat'))
    load(fullfile(Datadir,'AlignToTemplateWithAllenFrame.mat'))
else
    AlignCell = cell(1,length(miceopt));
end
miceopt=orimiceopt;
for midx = 1:length(miceopt)
    load(fullfile(Datadir, miceopt{midx}, [miceopt{midx} 'BrainModel.mat']))
    Model = BrainModel;
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
    
    if length(AlignCell)<midx
        AlignCell{midx}=[];
    end
        
    if  ~isempty(AlignCell{midx})
        TM = AlignCell{midx};
        okay = 1;
    else
        %         if ~exist('rect_brain','var')
        %             H=figure
        %             imagesc(templateframe)
        %             [~,rect_brain] = imcrop();
        %             delete(H)
        %         end
        %
        rect_brain = [80 30 440 330];
        searchwin=50;
        sub_brain = imcrop(templatepix,rect_brain);
        sub_tmp = imcrop(removepix,[rect_brain(1)-searchwin./2 rect_brain(2)-searchwin./2 rect_brain(3)+searchwin rect_brain(4)+searchwin]);
        sub_brain(isnan(sub_brain))=0;
        
     
        c=normxcorr2(sub_brain,sub_tmp);
        
        % offset found by correlation
        [max_c, imax]=max(abs(c(:)));
        [ypeak,xpeak]=ind2sub(size(c),imax(1));
        corr_offset = [(xpeak-size(sub_brain,2)) (ypeak-size(sub_brain,1))];
        
        % relative offset of position of subimages
        rect_offset = [-searchwin./2 -searchwin./2];
        
        % total offset
        offset = corr_offset + rect_offset;
        
        TM = imregtform(uint8(removepix),uint8(templatepix),'similarity',optimizer,metric);
        TM.T = eye(3);
        TM.T(3,1)=-offset(1);
        TM.T(3,2)=-offset(2);
        
        okay = 0;
    end
    
    HH = figure('name',[miceopt{midx} '_Alignments']);
    subplot(1,2,1)
    h1 = imagesc(imfuse((templatepix),(removepix)));%,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
    if midx ==1
        title('Before')
    end
    axis square
    axis off
    % Now Use the reference to change pixels of
    subplot(1,2,2)
    h2 = imagesc(imfuse(templatepix,imwarp(removepix,TM,'OutputView',imref2d(size(templatepix)))));%,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
    if midx ==1
        title('After')
    end
    axis square
    axis off
    drawnow
    
    % Manually shift!
    key = '';
    while ~okay
        if exist('h1')
            delete(h1),delete(h2)
        end
        
        subplot(1,2,1)
        h1 = imagesc(imfuse((templatepix),(removepix)));%,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
        if midx ==1
            title('Before')
        end
        axis square
        axis off
        % Now Use the reference to change pixels of
        subplot(1,2,2)
        h2 = imagesc(imfuse(templatepix,imwarp(removepix,TM,'OutputView',imref2d(size(templatepix)))));%,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]));
        if midx ==1
            title('After')
        end
        axis square
        axis off
        drawnow
        
        if strcmp(key,'d')
            TM.T(3,1) = TM.T(3,1)+1;
            key = '0';
        elseif strcmp(key,'a')
            TM.T(3,1) = TM.T(3,1)-1;
            key = '0';
        elseif strcmp(key,'s')
            TM.T(3,2) = TM.T(3,2)+1;
            key = '0';
        elseif strcmp(key,'w')
            TM.T(3,2) = TM.T(3,2)-1;
            key = '0';
        elseif strcmp(key,'f')
            TM.T(1,1) = TM.T(1,1)*0.99;
            key = '0';
        elseif strcmp(key,'g')
            TM.T(1,1) = TM.T(1,1)*1.01;
            key = '0';
        elseif strcmp(key,'v')
            TM.T(2,2) = TM.T(2,2)*0.99;
            key = '0' ;
        elseif strcmp(key,'b')
            TM.T(2,2) = TM.T(2,2)*1.01;
            key = '0';
        elseif strcmp(key,'q')
            for i = 1:3
                for j = 1:3
                    try
                        TM.T(i,j) = 0;
                    catch
                        TM.T(i,j) = 1;
                    end
                end
            end
            key = '0';
        elseif strcmp(key,'k')
            okay = 1;
        else
            %         while ~waitforbuttonpress
            %             pause(0.01)
            %         end
            waitforbuttonpress
            key = get(gcf,'CurrentCharacter');
        end
    end
    AlignCell{midx} = TM;
    save(fullfile(Datadir,'AlignToTemplateWithAllenFrame'),'AlignCell','miceopt','TemplateModel','removepix')
    
    saveas(gcf,fullfile(Datadir,[miceopt{midx} '_AlignmentsFGM.fig']))
    saveas(gcf,fullfile(Datadir,[miceopt{midx} '_AlignmentsFGM.bmp']))
end
%% Figure 2B - Whole brain dprime
dpr = nan(400,400,length(miceopt));
alldFF = nan(400,400,2,200,length(miceopt));
for midx=1:length(miceopt)
    % Load Brain outline
    load(fullfile(Datadir, miceopt{midx}, [miceopt{midx} 'BrainModel.mat']))
    Model = BrainModel;
    
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
    removepix=smooth2a(removepix,3);
    removepix(removepix<1)=0;
    removepix = logical(removepix);
    % load data
    TMP = load(fullfile(Datadir,miceopt{midx},'ProcessedData.mat'));
    tmp2 = nan(400,400,2,200);
    
    TM = AlignCell{midx};
    
    for optid = 1:2
        tmp = TMP.AllDat{optid};
        tmp(repmat(removepix,[1,1,length(timeline),size(tmp,4)]))=nan;
        tmp2(:,:,optid,1:size(tmp,4))=nanmean(tmp(:,:,timeline>=0&timeline<=1000,:),3);
        alldFF(:,:,optid,1:size(tmp,4),midx) = imwarp(nanmean(tmp(:,:,timeline>=0&timeline<=1000,:),3),TM,'OutputView',imref2d([size(tmp,1),size(tmp,2),size(tmp,4)]));
    end
    tmp2 = imwarp(tmp2,TM,'OutputView',imref2d(size(tmp2)));
    dpr(:,:,midx) = (nanmean(tmp2(:,:,2,:),4)-nanmean(tmp2(:,:,1,:),4))./(0.5.*(sqrt(nanvar(tmp2(:,:,1,:),[],4)+nanvar(tmp2(:,:,2,:),[],4))));
    
end

dpr = reshape(dpr,400*400,[]);
p = nan(1,400*400);
parfor pixid = 1:400*400
    if sum(~isnan(dpr(pixid,:)))<length(miceopt)
        continue
    end
    p(pixid)= anovan(dpr(pixid,:),{mousegen(:)},'display','off');
end


figure('name','Wholebrain dprime')
for gidx=1:length(MouseGenOpt)
    tmp = nanmean(dpr(:,ismember(mousegen,MouseGenOpt{gidx})),2);
    subplot(1,3,gidx)
    h=imagesc(reshape(tmp,400,400));
    axis square
    axis off
    alphamap = ones(400,400);
    alphamap(isnan(reshape(tmp,400,400)))=0;
%     alphamap(reshape(p,400,400)>0.05)=0.5;
    set(h,'AlphaData',alphamap)
    colormap(redblue)
    set(gca,'clim',[-0.7 0.7])
    title(MouseGenOpt{gidx})
    
    hold on
    set(gca,'ydir','reverse')
%     ylim([0 200])
      plot(TemplateModel.AllX,TemplateModel.AllY,'.','MarkerSize',2,'color',[0.5 0.5 0.5])
end

saveas(gcf,fullfile(Datadir,['DprimeBrain.fig']))
%% Figure 2C
AreaCols = [15, 114, 186; 216 84 38; 236 177 32; 126 46 141; 117 172 66; 79 190 237; 128 130 133];
cols = AreaCols./255;
    
    tmp = squeeze(nanmean(AllMiceTC(:,timeline>=0&timeline<=1000,:,2,:,:),2)); %During Opto
    g1 = repmat(AREAGROUPNAMES',[1,length(HemOpt),size(tmp,3),length(miceopt)]); %Area
    g2 = repmat(HemOpt',[1,length(AREAGROUPNAMES),size(tmp,3),length(miceopt)]); %Hemisphere
    g2 = permute(g2,[2,1,3,4]);
    g4 = repmat(mousegen',[1,length(AREAGROUPNAMES),length(HemOpt),size(tmp,3)]); %Mouse gen
    g4 = permute(g4,[2,3,4,1]);
    g5 = repmat(miceopt',[1,length(AREAGROUPNAMES),length(HemOpt),size(tmp,3)]); %Mouse ID
    g5 = permute(g5,[2,3,4,1]);
    
    tmp = tmp(:);
    g1 = g1(:);
    g2 = g2(:);
%     g3 = g3(:);
    g4 = g4(:);
    idx = find(~isnan(tmp));

    
    % Mixed Linear Model Statistics
    tbl = table(double(tmp(idx)),g1(idx),g2(idx),g4(idx),g5(idx),'VariableNames',{'dFF','AREA','HEM','GEN','MOUSEID'});
    
    %% Figure 2 & Stats
    AREASHERE=AREAS;
    figure('name','dFF bar D1')
            pall = nan(length(AREASHERE),length(AREASHERE),2,2);

    for genid=1:2
        
        Directidx = find(ismember(g4(idx),MouseGenOpt{genid}));
        glme = fitglme(tbl(Directidx,:),'dFF ~ 1 + HEM + AREA + AREA:HEM + (1|MOUSEID)','Distribution','Normal');
        anova(glme)
        
        
        subplot(1,2,genid)
        tmp = squeeze(nanmean(AllMiceTC(:,timeline>=0&timeline<=1000,:,2,:,:),2));%-squeeze(nanmean(AllMiceTC(:,timeline>=TW{twid}(1)&timeline<=TW{twid}(2),1,1,:,:),2));
        
        clear dffavg
        clear dffsem
        dFFOptoNorm = reshape(tmp(:,:,:,ismember(mousegen,MouseGenOpt{genid})),length(AREAS),2,[]);
        dFFOptoNorm(:,:,sum(isnan(reshape(dFFOptoNorm,length(AREAS)*2,[])),1)==length(AREAS)*2)=[];
        dffavg = nanmean(dFFOptoNorm,3);
        dffsem = nanstd(dFFOptoNorm,[],3)./sqrt(size(dFFOptoNorm,3)-1);
        
        
        h=barwitherr(dffsem',dffavg');
        for areaid=1:length(AREAS)
            h(areaid).EdgeColor='none';
            h(areaid).FaceColor = cols(areaid,:);
            h(areaid).BarWidth = 1;
        end
        hold on
%         ylim([-0.003 0.003])
        box off
        axis square
        title(MouseGenOpt{genid})
        set(gca,'XTickLabel',HemOpt)
        ylabel(['\DeltaF/F'])
        FigureDefault
        emm = emmeans(glme, {'HEM','AREA'}, 'effects');
        
        AREASHERE = AREAGROUPNAMES;
        for hemid=1:2
            for areaid=1:length(AREASHERE)
                for areaid2=2:length(AREASHERE)
                    if areaid2<=areaid
                        continue
                    end
                    L_HL = (strcmp(emm.table.HEM,HemOpt{hemid})& strcmp(emm.table.AREA,AREASHERE{areaid}))'...
                        -(strcmp(emm.table.HEM,HemOpt{hemid})& strcmp(emm.table.AREA,AREASHERE{areaid2}))';
                    H0_Stim = contrasts_wald(glme,emm,L_HL);
                    p=H0_Stim.pVal;
                    if p>0.5 %Make a two sided test
                        p=1-p;
                    end
                    pall(areaid,areaid2,hemid,genid)=p;
                    
                    if p<0.05
                        sigstar({[h(areaid).XEndPoints(hemid) h(areaid2).XEndPoints(hemid)]},p)
                    end
                end
            end
        end
        
    end
linkaxes(get(gcf,'children'))

    %Random slope
    glme = fitglme(tbl,'dFF ~ 1 + HEM + GEN + AREA + GEN:AREA + HEM:GEN + GEN:HEM + GEN:AREA:HEM + (1|MOUSEID)','Distribution','Normal');
    VarExplained = glme.Rsquared.Ordinary*100;
    anova(glme)
    
    emm = emmeans(glme, {'HEM','GEN','AREA'}, 'effects');
    
    for hemid=1:2
        
        L_HL = (strcmp(emm.table.HEM,HemOpt{hemid})& strcmp(emm.table.GEN,MouseGenOpt{1}))'...
            -(strcmp(emm.table.HEM,HemOpt{hemid})& strcmp(emm.table.GEN,MouseGenOpt{2}))';  
        H0_Stim = contrasts_wald(glme,emm,L_HL);
        p=H0_Stim.pVal;
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        p=p*2 %Multiple corrections
    end

    
    
    figure('name','dFF bar D1D2')
for hemid=1:2
    subplot(1,2,hemid)
    tmp = squeeze(nanmean(AllMiceTC(:,timeline>=0&timeline<=1000,hemid,2,:,:),2));%-squeeze(nanmean(AllMiceTC(:,timeline>=TW{twid}(1)&timeline<=TW{twid}(2),1,1,:,:),2));
    
    clear dffavg
    clear dffsem
    for gidx=find(ismember(MouseGenOpt,{'Drd1','Drd2'}))%1:length(MouseGenOpt)
        dFFOptoNorm = reshape(tmp(:,:,ismember(mousegen,MouseGenOpt{gidx})),length(AREAS),[]);
        dFFOptoNorm(:,sum(isnan(dFFOptoNorm),1)==length(AREAS))=[];
        dffavg(:,gidx) = nanmean(dFFOptoNorm,2);
        dffsem(:,gidx) = nanstd(dFFOptoNorm,[],2)./sqrt(size(dFFOptoNorm,2)-1);
    end
    dffavg = reshape(dffavg,length(AREAS),[]);
    dffsem = reshape(dffsem,length(AREAS),[]);
    
    h=barwitherr(dffsem',dffavg');
    for areaid=1:length(AREAS)
        h(areaid).EdgeColor='none'
        h(areaid).FaceColor = cols(areaid,:);
        h(areaid).BarWidth = 1;
    end
    hold on
    ylim([-0.003 0.003])
    box off
    axis square
    title(HemOpt{hemid})
    set(gca,'XTickLabel',MouseGenOpt)
    ylabel('DF/F')
end


saveas(gcf,fullfile(Datadir,['D1D2_Bars.fig']))

%% Figure S2A
%Concatenate all data
AllDatMice = cell(1,length(miceopt)); %nan(xpix,ypix,length(TW),8,length(OPTOOpt),length(miceopt),'single');
AllOptoMice = cell(1,length(miceopt));
AllNaiveMice = cell(1,length(miceopt));
AllStateMice = cell(1,length(miceopt));
xpix=400;
ypix=400;
for midx=7:8
    TMP = load(fullfile(Datadir,miceopt{midx},'ProcessedData.mat'));
    for twid=2%1:length(TW)
        tmp2 = nan(xpix,ypix,0,'single');
        for optidx=1:2%length(OptoOpt)
            %        TMP.AlldFF(sidx,:,optidx);
            if ~exist('timeline','var')
                timeline = TMP.timeline;
            end
            tmp = cat(4,TMP.AllDat{optidx});
            tmp2 = cat(3,tmp2,squeeze(nanmean(tmp(:,:,timeline>=TW{twid}(1)&timeline<=TW{twid}(2),:),3)));
            AllOptoMice{midx} = cat(2,AllOptoMice{midx},repmat(optidx,1,size(tmp,4)));
            tmpnaive = TMP.AllNaive{optidx};
            
            if midx == 7
                tmpnaive = cat(2,repmat({'Naive'},1,100),repmat({'Expert'},1,50));
            elseif midx==8
                tmpnaive =  cat(2,repmat({'Naive'},1,50),repmat({'Expert'},1,50));
            end
            
            AllNaiveMice{midx} = cat(2,AllNaiveMice{midx},tmpnaive);
            AllStateMice{midx} = cat(2,AllStateMice{midx},TMP.AllState{optidx});
        end
        
        
%         TM = AlignCell{midx};
%         tmp2 = imwarp(tmp2,TM,'OutputView',imref2d(size(tmp2)));
        
        AllDatMice{midx}(:,:,twid,:) = tmp2;
        
    end
    
end


%% mixed Linear Model
g1 = [];%cat(2,AllOptoMice{:});
g2 = []; %Mouse
g3 = []; %NaiveExp
g4 = []; %Genotype
g5 = []; %State
for midx=7:8%1:length(miceopt)
 
    g1 = [g1 AllOptoMice{midx}];
    g2 = [g2 repmat(midx,1,length(AllOptoMice{midx}))];
    g3 = [g3 AllNaiveMice{midx}];
    g4 = [g4 repmat(find(strcmp(MouseGenOpt,mousegen{midx})),1,length(AllOptoMice{midx}))];
    g5 = [g5 AllStateMice{midx}];

end
g1=g1-1; %convert to 0 for off and 1 for on
g3(ismember(g3,'Naive')) = {0}; %convert to 0 for naive stim and 1 for expert
g3(cellfun(@ischar,g3)) = {1}; %convert to 0 for naive stim and 1 for expert
g3 = cell2mat(g3);
g6 = g4; %Dummy coding for Drd2
g6(g6==1| g6 ==2) =0;
g6(g6==3)=1; %Drd2
g4(g4==1 | g4 ==3)=0;
g4(g4==2)=1; %Drd1
g5(ismember(g5,'Anethetized')) = {0}; %convert to 0 for anesthetized and 1 for awake
g5(cellfun(@ischar,g5)) = {1}; %convert to 0 for naive stim and 1 for expert
g5 = cell2mat(g5);
inclmice = 1:length(miceopt);
twid=2;
TMP = cat(4,AllDatMice{:});
TMP = squeeze(TMP(:,:,twid,:));

%% dFF (Extra analysis naive versus expert)
% AREAS=AnalysisParameters.AREAS
% cols = cat(2,[108; 99; 91], [140; 20; 0], [110;81;76], [140;20;160], [189;105;41], [229;176;46], [154;192;209], [73;121;159],[111;60;151], [46;141;67], [209;211;212], [37;34;25])';
% % cols = cols./255;
% cols = cat(2,[108; 99; 91], [110;81;76], [189;105;41], [229;176;46], [154;192;209], [73;121;159],[111;60;151], [46;141;67], [209;211;212], [37;34;25])';
% cols = cols./255;
for midx=7:8
    figure('name',['dFFNaiveVsExpertMap ' miceopt{midx}])

      % Load Brain outline
    load(fullfile(Datadir, miceopt{midx}, [miceopt{midx} 'BrainModel.mat']))
    Model = BrainModel;
    
    idx = find(g2==midx);
    tmp = TMP(:,:,idx);
    tmp2 = nanmean(tmp,3);
    lims = quantile(abs(tmp2(:)),0.99);
    dffperarea = nan(2,length(AREAS),size(tmp,3));
    for j = 1:2 %Level
        for i = 2 %Opto
            subplot(1,2,(j))
            
            if isempty(tmp(g1(idx)==i-1&g3(idx)==j-1))
                continue
            end
            h=imagesc(nanmean(tmp(:,:,g1(idx)==i-1&g3(idx)==j-1),3),[-lims lims]);
            set(h,'AlphaData',~isnan(nanmean(tmp(:,:,g1(idx)==i-1&g3(idx)==j-1),3)))
            hold on
            plot(Model.AllX,Model.AllY,'k.','MarkerSize',4)
            if j==1
                title('Naive')
            elseif j==2
                title('Expert')
            end
            box off
            axis off
            hold on
            colormap redblue
            colorbar
            
            tmp2 = reshape(tmp(:,:,g1(idx)==i-1&g3(idx)==j-1),400*400,[]);

            for areaid=1:length(AREAS)
                areamask = false(xpix,ypix);
                aidx = find(ismember(Model.Rnames,AREAS{areaid}));
                if isempty(aidx) || strcmp(AREAS{areaid},'Vl')
                    aidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS{areaid}),Model.Rnames,'UniformOutput',0)));
                    for i = 1:length(aidx)
                        areamask(logical(round(Model.Regions{aidx(i)}))')=true;
                    end
                else
                    areamask(logical(round(Model.Regions{aidx}))')=true;
                end
                aridx= find(areamask);
                
                dffperarea(j,areaid,1:size(tmp2,2)) = nanmean(tmp2(aridx,:),1);
                
            end
           
            
            
        end
    end

    saveas(gcf,fullfile(Datadir,['dFFNaiveVsExpertMap_' miceopt{midx} '.fig']))
    
    figure('name',['BarExpertLevel ' miceopt{midx}])
    h=barwitherr(nanstd(dffperarea(:,~ismember(AREAS,'Vpor'),:),[],3)./sum(~isnan(dffperarea(:,~ismember(AREAS,'Vpor'),:)),3),nanmean(dffperarea(:,~ismember(AREAS,'Vpor'),:),3));
    
    for areaid=1:length(AREAS)-1
        h(areaid).EdgeColor='none'
        h(areaid).FaceColor = cols(areaid,:);
        h(areaid).BarWidth = 1;
    end
    set(gca,'XTickLabel',{'Naive','Trained'})
    
    y = dffperarea(:,~ismember(AREAS,'Vpor'),:);
    a1 = repmat({'Naive','Trained'}',[1,10,size(y,3)]);
    a2 = repmat(AREAS(~ismember(AREAS,'Vpor'))',[1,2,size(y,3)]);
    a2 = permute(a2,[2,1,3]);
    
    y=y(:);
    a1 = a1(:);
    a2= a2(:);
    
    iidx = find(~isnan(y));
    p=anovan(y(iidx),{a1(iidx),a2(iidx)},'model','interaction');
    
    box off
    title([miceopt{midx} ',  expertlevel p= ' num2str(p(1)) ', area p=' num2str(p(2)) ', interaction p=' num2str(p(3))])
        saveas(gcf,fullfile(Datadir,['dFFNaiveVsExpertBar_' miceopt{midx} '.fig']))

    
end


%% Figure S2B/C
% Now load in simple lick task data for mouse 7
midx = 7
TMP2 = load(fullfile(Datadir,miceopt{midx},'ProcessedData.mat'));
twid=2%1:length(TW)
optidx=1%OFF
%        TMP.AlldFF(sidx,:,optidx);
    timeline2 = TMP2.timeline;
tmp = cat(4,TMP2.AlldFF{:,:,1}); %Only take opto off
tmp = squeeze(nanmean(tmp(:,:,timeline2>=TW{twid}(1)&timeline2<=TW{twid}(2),:),3));
TaskDat = tmp;
TaskDat = reshape(TaskDat,400*400,[]);
ntr(3) = size(TaskDat,2);


TMP = cat(4,AllDatMice{:});
TMP = squeeze(TMP(:,:,twid,:));
idx = find(g2==7 & g1 ==1 & g3 == 0);
OptoNaiveDat = TMP(:,:,idx);
OptoNaiveDat = reshape(OptoNaiveDat,400*400,[]);
ntr(1) = size(OptoNaiveDat,2);

idx = find(g2==7 & g1 ==1& g3 == 1);
OptoTrainedDat = TMP(:,:,idx);
OptoTrainedDat = reshape(OptoTrainedDat,400*400,[]);
ntr(2) = size(OptoTrainedDat,2);


% Load Brain outline
load(fullfile(Datadir, miceopt{midx}, [miceopt{midx} 'BrainModel.mat']))
Model = BrainModel;
 
AREAS2 = Model.Rnames;
dFFPerArea = nan(length(AREAS2),sum(ntr)); %for opto naive, then opto trained, then the ntask
PixCol = nan(1,400*400);
for areaid=1:length(AREAS2)
    areamask = false(xpix,ypix);
    aidx = find(ismember(Model.Rnames,AREAS2{areaid}));
    if isempty(aidx) || strcmp(AREAS2{areaid},'Vl')
        aidx = find(~cellfun(@isempty,cellfun(@(X) strfind(X,AREAS2{areaid}),Model.Rnames,'UniformOutput',0)));
        for i = 1:length(aidx)
            areamask(logical(round(Model.Regions{aidx(i)}))')=true;
        end
    else
        areamask(logical(round(Model.Regions{aidx}))')=true;
    end
    aridx= find(areamask);
    PixCol(aridx)=areaid;
    
   % Opto OFF versus Task
   
   
   dFFPerArea(areaid,:) = cat(2,nanmean(OptoNaiveDat(aridx,:),1),nanmean(OptoTrainedDat(aridx,:),1),nanmean(TaskDat(aridx,:),1));

end
gn = cat(2,repmat(0,[1,ntr(1)]),repmat(1,[1,ntr(2)]),repmat(2,[1,ntr(3)]));
gn = repmat(gn,[length(AREAS2),1]);
ga = repmat(AREAS2',[1,sum(ntr)]);    
anovan(dFFPerArea(:),{gn(:),ga(:)},'model','interaction','varnames',{'TaskType','Area'});
      
avg = cat(2,nanmean(dFFPerArea(:,gn(1,:)==0),2),nanmean(dFFPerArea(:,gn(1,:)==1),2),nanmean(dFFPerArea(:,gn(1,:)==2),2))
stdv = cat(2,nanstd(dFFPerArea(:,gn(1,:)==0),[],2)./sqrt(sum(gn(1,:)==0)-1),nanstd(dFFPerArea(:,gn(1,:)==1),[],2)./sqrt(sum(gn(1,:)==1)-1),nanstd(dFFPerArea(:,gn(1,:)==2),[],2)./sqrt(sum(gn(1,:)==2)-1))
figure; 
h=barwitherr(stdv(~ismember(AREAS2,'Vpor'),:),avg(~ismember(AREAS2,'Vpor'),:));
set(gca,'XTick',1:length(AREAS2),'XTickLabel',AREAS2(~ismember(AREAS2,'Vpor')),'XTickLabelRotation',25)
for i = 1:length(h)
    h(i).BarWidth=1;
    h(i).EdgeColor = 'none' 
end
    
figure;
scatter3(nanmean(OptoNaiveDat,2),nanmean(OptoTrainedDat,2),nanmean(TaskDat,2),8,PixCol,'filled')
xlabel('Naive'); ylabel('Expert'); zlabel('In Task')
lims = [quantile(TaskDat(:),0.05) quantile(TaskDat(:),0.95)];
xlim(lims);ylim(lims);zlim(lims)

xvec = lims(1):0.001:lims(2);
figure;
subplot(1,3,1)
scatter(nanmean(OptoNaiveDat,2),nanmean(OptoTrainedDat,2),8,PixCol,'filled')
tmp1 = nanmean(OptoNaiveDat,2);
tmp2 = nanmean(OptoTrainedDat,2);
nidx = isnan(tmp1)|isnan(tmp2);
[r1,p] = corr(tmp1(~nidx),tmp2(~nidx))
title(['r = '  num2str(r1) ', p=' num2str(p)])
xlabel('Naive'); ylabel('Trained')
xlim(lims);ylim(lims);
axis square
%Regression
idx = ~isnan(tmp1)&~isnan(tmp2);
p = polyfit(tmp1(idx),tmp2(idx),1);
ypred = xvec*p(1)+p(2);
hold on;
plot(xvec,ypred,'r-','LineWidth',2);


subplot(1,3,2)
scatter(nanmean(OptoNaiveDat,2),nanmean(TaskDat,2),8,PixCol,'filled')
tmp1 = nanmean(OptoNaiveDat,2);
tmp2 = nanmean(TaskDat,2);
nidx = isnan(tmp1)|isnan(tmp2);
[r2,p] = corr(tmp1(~nidx),tmp2(~nidx))
title(['r = '  num2str(r2) ', p=' num2str(p)])
xlabel('Naive'); ylabel('In Task')
xlim(lims);ylim(lims);
axis square
%Regression
idx = ~isnan(tmp1)&~isnan(tmp2);
p = polyfit(tmp1(idx),tmp2(idx),1);
ypred = xvec*p(1)+p(2);
hold on;
plot(xvec,ypred,'r-','LineWidth',2);

subplot(1,3,3)
scatter(nanmean(OptoTrainedDat,2),nanmean(TaskDat,2),8,PixCol,'filled')
tmp1 = nanmean(OptoTrainedDat,2);
tmp2 = nanmean(TaskDat,2);
nidx = isnan(tmp1)|isnan(tmp2);
[r3,p] = corr(tmp1(~nidx),tmp2(~nidx))
title(['r = '  num2str(r3) ', p=' num2str(p)])
xlabel('Trained'); ylabel('In Task')
xlim(lims);ylim(lims);
axis square
%Regression
idx = ~isnan(tmp1)&~isnan(tmp2);
p = polyfit(tmp1(idx),tmp2(idx),1);
ypred = xvec*p(1)+p(2);
hold on;
plot(xvec,ypred,'r-','LineWidth',2);
saveas(gcf,fullfile(Datadir,['CorrelationOptoVsTask_' miceopt{midx} '.fig']))


[hyp,p,z] = mengz(r2, r1, r3, sum(~nidx))
[hyp,p,z] = mengz(r3, r2, r1, sum(~nidx))
[hyp,p,z] = mengz(r3, r1, r2, sum(~nidx))

figure;
imagesc(unique(PixCol(~isnan(PixCol))))
set(gca,'XTick',1:length(unique(PixCol(~isnan(PixCol)))),'XTickLabel',AREAS2,'XTickLabelRotation',90)


%% D1 mice naive versus expert
% g1 = [];%cat(2,AllOptoMice{:});
% g2 = []; %Mouse
% g3 = []; %NaiveExp
% g4 = []; %Genotype
% g5 = []; %State
twid=2;
TMP = cat(4,AllDatMice{:});
TMP = squeeze(TMP(:,:,twid,:));
TMP(repmat(templatepix,[1,1,size(TMP,3)]))=nan;
y = nanmean(reshape(TMP,400*400,[]),1);
idx = find(g2==9);

tmp = y(idx);
for j = 1:2 %Level
    [h,ppost(j)] = ttest2(tmp(g1(idx)==0&g3(idx)==j-1),tmp(g1(idx)==1&g3(idx)==j-1));
    for i = 1:2 %Opto
        meantmp(i,j) = nanmean(tmp(g1(idx)==i-1&g3(idx)==j-1));
        semtmp(i,j)=nanstd(tmp(g1(idx)==i-1&g3(idx)==j-1))./sqrt(sum(g1(idx)==i-1&g3(idx)==j-1)-1);
    end
end
[p, t] = anovan(y(idx),{g1(idx),g3(idx)},'model',2,'varnames',{'Opto','Level'});
figure; 
h=barwitherr(semtmp',meantmp');
set(gca,'XTickLabel',{'Naive','Expert'})
legend('Opto Off','Opto On')
title('D1')
box off
saveas(gcf,fullfile(storepath,'Figures',['ExpertVsNaive.fig']))
saveas(gcf,fullfile(storepath,'Figures',['ExpertVsNaive.bmp']))
