datapath = '\Path\To\Data\';
AllMice2Use = {'Alaska','Becquerel','Berlin','Cherry','China','Darwin','Druif','Faraday','Galilei','Grape','Haydn','Holland','HongKong','Jemison','Jenkins','Joule','Karlsen','Kirchoff','Leeuwenhoek','Lully','Ising','Mahler'};
Gentype = {'Control','D2','D2','D1','D1','Control','Control','D2','D2','Control','D1','D1','D2','D2','D1','Control','D1','D2','D2','D1','Control','D1'};
ActuallyUsed = {};
AllCondNames = {'Control','D1','D2'};
OptoOpt = {'Opto Off','Opto On','Opto 1st Lick'};
SideOptNames = {'left','right','Total'};
NormalizeLicks = 0; %normalize to median no opto cond
nback = 15;
TakeOnlyOptoBefore = 0; %Only take opto trails were laseronset <-500
takeoptoonly = 1; %Taking only opto trials (checkbox = 1)

%% Summarize Groups

AllMicePerGroup = cell(1,length(AllCondNames));
AllDatesPerGroup = cell(1,length(AllCondNames));
AllData = cell(1,length(AllCondNames));
maxnrmice = 0;


for dsidx =1:length(AllMicePerGroup)
    %% Find mice/days
    % addpath(genpath(datapath)) Directory file Mouse data
    
    MouseOpt2Use = AllMice2Use(ismember(Gentype,AllCondNames{dsidx}));
    
    
    % Mouse Name
    AllOpt = dir(fullfile(datapath,'**','*.mat'));
    AllOptSplit = cellfun(@(X) strsplit(X,'_'),{AllOpt(:).name},'UniformOutput',0);
    AllOptSplit(cell2mat(cellfun(@length,AllOptSplit,'UniformOutput',0))<3) = [];
    
    MouseOpt = unique(cellfun(@(X) X{1},AllOptSplit,'UniformOutput',0));
    midx = cell2mat(cellfun(@(X) find(ismember(MouseOpt,X)),MouseOpt2Use,'UniformOutput',false));
    clear DateOpt2Use
    if ~exist('DateOpt2Use','var')
        for m=1:length(midx)
            
            AllOpt = dir(fullfile(datapath,MouseOpt{midx(m)},[MouseOpt{midx(m)} '*'],'**',[MouseOpt{midx(m)} '*.mat']));
            AllOptSplit = cellfun(@(X) strsplit(X,'_'),{AllOpt(:).name},'UniformOutput',0);
            AllOptSplit(cell2mat(cellfun(@length,AllOptSplit,'UniformOutput',0))<3) = [];
            
            DateOpt = unique(cellfun(@(X) X{2},AllOptSplit,'UniformOutput',0));
            id = cell2mat(cellfun(@(X) strcmp(X{1},MouseOpt{midx(m)}),AllOptSplit,'UniformOutput',0));
            DateOpt2Use{m} = unique(cellfun(@(X) X{2},AllOptSplit(id),'UniformOutput',0));
        end
    end
    AllMicePerGroup{dsidx} = MouseOpt2Use(ismember(MouseOpt2Use,MouseOpt));
    AllDatesPerGroup{dsidx} = DateOpt2Use;
    
    if length(MouseOpt2Use)>maxnrmice
        maxnrmice = length(MouseOpt2Use);
    end
    
    % Mouse Name
    AllOpt = dir(fullfile(datapath,'**','*.mat'));
    AllOptSplit = cellfun(@(X) strsplit(X,'_'),{AllOpt(:).name},'UniformOutput',0);
    AllOptSplit(cell2mat(cellfun(@length,AllOptSplit,'UniformOutput',0))<3) = [];
    SideOpt = {'left','right'};
    DateOpt = unique(cellfun(@(X) X{2},AllOptSplit,'UniformOutput',0));
    
    didx = cellfun(@(X) find(ismember(DateOpt,X)),DateOpt2Use,'UniformOutput',0);
    %% Select Data
    
    maxdays = max(cell2mat(cellfun(@length,didx,'UniformOutput',0)));
    AllOptNames= {AllOpt(:).name};
    AllOptNames = unique(AllOptNames);
    % Preallocate variables
    AllData{dsidx}.AllReaction = cell(length(midx),maxdays);
    AllData{dsidx}.AllSides = cell(length(midx),maxdays);
    AllData{dsidx}.AllRT = cell(length(midx),maxdays);
    AllData{dsidx}.PassivePerc = cell(length(midx),maxdays);
    AllData{dsidx}.AllMoving = cell(length(midx),maxdays);
    AllData{dsidx}.BiasCell = cell(length(midx),maxdays);
    AllData{dsidx}.AllLeftLicks = cell(length(midx),maxdays);
    AllData{dsidx}.AllRightLicks = cell(length(midx),maxdays);
    AllData{dsidx}.AllOpto = cell(length(midx),maxdays);
    AllData{dsidx}.AllOptoOn = cell(length(midx),maxdays);
    AllData{dsidx}.AllOptoOff = cell(length(midx),maxdays);
    AllData{dsidx}.TakeIntoAccount = cell(length(midx),maxdays);
    AllData{dsidx}.LaserSetting = cell(length(midx),maxdays);
    AllData{dsidx}.TrialNumber = cell(length(midx),maxdays);

    %Importeer gewenste data
    for m = 1:length(midx)
        for d = 1:length(didx{m})
            fileidx = find(~cell2mat(cellfun(@isempty,(strfind(AllOptNames,[MouseOpt{midx(m)}])),'UniformOutput',0)) & ~cell2mat(cellfun(@isempty,(strfind(AllOptNames,[DateOpt{didx{m}(d)}])),'UniformOutput',0)));
            tmpreaction = [];
            tmpsides = [];
            tmprt = [];
            tmppassive=[];
            tmpbias = [];
            tmpmoving = [];
            tmpleftlicks = [];
            tmprightlicks = [];
            tmpopto = [];
            tmpoptoon = [];
            tmpoptooff = [];
            tmptakeintoaccount = [];
            tmplasersetting = [];
            tmptrialnr = [];
            for blockidx = 1:length(fileidx)
                blockname = strsplit(AllOptNames{fileidx(blockidx)},'_');
                blockname = strsplit(blockname{end},'.mat');
                if ~exist(fullfile(datapath,MouseOpt{midx(m)},[MouseOpt{midx(m)} DateOpt{didx{m}(d)}],blockname{1},AllOptNames{fileidx(blockidx)}))
                    continue
                end
                %                 CheckReactions(fullfile(datapath,MouseOpt{midx(m)},[MouseOpt{midx(m)} DateOpt{didx{m}(d)}],blockname{1},AllOptNames{fileidx(blockidx)}));
                try
                    tmp = load(fullfile(datapath,MouseOpt{midx(m)},[MouseOpt{midx(m)} DateOpt{didx{m}(d)}],blockname{1},AllOptNames{fileidx(blockidx)}));
                catch
                    try
                        tmp = load(fullfile(datapath,MouseOpt{midx(m)},[MouseOpt{midx(m)} DateOpt{didx{m}(d)}],blockname{1},AllOptNames{fileidx(blockidx)}),'-ascii');
                    catch
                        
                        disp([fullfile(datapath,MouseOpt{midx(m)},[MouseOpt{midx(m)} DateOpt{didx{m}(d)}],blockname{1},AllOptNames{fileidx(blockidx)}) ' doesn''t load'])
                        continue
                    end
                end
                ActuallyUsed = {ActuallyUsed{:} MouseOpt{midx(m)}};
                tmpmoving = [tmpmoving tmp.LOG.MovingGrating];
                tmpreaction = [tmpreaction tmp.LOG.Reaction];
                tmpsides = [tmpsides tmp.LOG.Side];
                tmprt = [tmprt tmp.LOG.RT];
                tmppassive = [tmppassive tmp.LOG.PassivPerc];
                [BiasMat,LickedBasedBias] = extractbiasidxBehavior(tmp.LOG,nback);
                tmpbias = [tmpbias BiasMat];
                tmpleftlicks = [tmpleftlicks tmp.LOG.RTleftVec];
                tmprightlicks = [tmprightlicks tmp.LOG.RTrightVec];
                tmptrialnr = [tmptrialnr tmp.LOG.Trial];
                %             if ismember(tmp.LOG.Date,DatesNoOpto)
                %                 tmp.LOG.Opto = repmat(0,1,length(tmp.LOG.Opto));
                %             end
                tmpopto = [tmpopto tmp.LOG.Opto];
                if ~isfield(tmp.LOG,'LaserSetting')
                    tmp.LOG.LaserSetting = repmat(0,1,length(tmp.LOG.Opto));
                end
                tmplasersetting = [tmplasersetting tmp.LOG.LaserSetting];
                if takeoptoonly && isfield(tmp.LOG,'OptoRegulator')
                    tmptakeintoaccount = [tmptakeintoaccount tmp.LOG.OptoRegulator];
                else
                    %                     keyboard
                    tmptakeintoaccount = [tmptakeintoaccount ones(1,length(tmp.LOG.Opto))];
                end
                tmptakeintoaccount = logical(tmptakeintoaccount);
                try
                    if length(tmp.LOG.OptoOnTime)<length(tmp.LOG.Opto)
                        tmp.LOG.OptoOnTime = [ tmp.LOG.OptoOnTime nan(1,length(tmp.LOG.Opto)-length(tmp.LOG.OptoOnTime))];
                    end
                    if length(tmp.LOG.OptoOffTime)<length(tmp.LOG.Opto)
                        tmp.LOG.OptoOffTime = [ tmp.LOG.OptoOffTime nan(1,length(tmp.LOG.Opto)-length(tmp.LOG.OptoOffTime))];
                    end
                    tmpoptoon = [tmpoptoon tmp.LOG.OptoOnTime];
                    tmpoptooff = [tmpoptooff tmp.LOG.OptoOffTime];
                catch
                    tmpoptoon = [tmpoptoon nan(1,length(tmp.LOG.Opto))];
                    tmpoptooff = [tmpoptooff nan(1,length(tmp.LOG.Opto))];
                    
                end
            end
            if length(tmplasersetting) == 1
                tmplasersetting = repmat(tmplasersetting,1,length(tmptakeintoaccount));
            end
            
            AllData{dsidx}.LaserSetting{m,d} = tmplasersetting(tmptakeintoaccount);
            AllData{dsidx}.AllReaction{m,d} = tmpreaction(tmptakeintoaccount);
            AllData{dsidx}.AllSides{m,d} = tmpsides(find(tmptakeintoaccount));
            AllData{dsidx}.AllRT{m,d} = tmprt(find(tmptakeintoaccount));
            AllData{dsidx}.PassivePerc{m,d} =tmppassive(find(tmptakeintoaccount));
            AllData{dsidx}.BiasCell{m,d} = tmpbias(find(tmptakeintoaccount));
            AllData{dsidx}.AllMoving{m,d} = tmpmoving(find(tmptakeintoaccount));
            AllData{dsidx}.AllLeftLicks{m,d} = tmpleftlicks(find(tmptakeintoaccount));
            AllData{dsidx}.AllRightLicks{m,d} = tmprightlicks(find(tmptakeintoaccount));
            AllData{dsidx}.AllOpto{m,d}=tmpopto(find(tmptakeintoaccount));
            AllData{dsidx}.AllOptoOn{m,d} = tmpoptoon(find(tmptakeintoaccount));
            AllData{dsidx}.AllOptoOff{m,d} = tmpoptooff(find(tmptakeintoaccount));
            AllData{dsidx}.TrialNumber{m,d} = tmptrialnr(find(tmptakeintoaccount));

            AllData{dsidx}.TakeIntoAccount{m,d} = tmptakeintoaccount(find(tmptakeintoaccount));
         
        end
    end
end


unique(ActuallyUsed)
optoopt = [0,1]; %optoopt = [0,1];

%% Laser setting
LaserOption2Use = cell(length(AllMicePerGroup),maxnrmice);
for dsidx = 1:length(AllMicePerGroup)

    MouseOpt = AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};

    for m = 1:length(MouseOpt)
        tmp = [AllData{dsidx}.LaserSetting{m,:}];
        tmp(isnan(tmp)| tmp == 0) = [];
        [UVals,id1,id2] = unique(tmp);
        [maxval,maxid] = max(arrayfun(@(X) sum(id2==X),1:length(UVals)));
        LaserOption2Use{dsidx,m} = [0,UVals(maxid)]; % Use only trials with this setting
    end
end


%% Licks
figure;
for dsidx =1:length(AllMicePerGroup)
    MouseOpt =AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    
    for m = 1:length(MouseOpt)
        subplot(maxnrmice,length(AllMicePerGroup),(m-1)*length(AllMicePerGroup)+dsidx)
        bigcount = 0;
        for d = 1:length(DateOpt{m})
            smallcount = 0;
            %Sort trials based on opto opt
            for p = 1:length(optoopt)
                
                trialvec = find(AllData{dsidx}.AllOpto{m,d}==optoopt(p));
                if isempty(trialvec)
                    continue
                end
                for tridx = 1:length(trialvec)
                    %OptoOnTime
                    h = plot( AllData{dsidx}.AllLeftLicks{m,d}{trialvec(tridx)},repmat(tridx+smallcount+bigcount,[1,length( AllData{dsidx}.AllLeftLicks{m,d}{trialvec(tridx)})]),'r.','MarkerSize',20);
                    hold on
                    i = plot( AllData{dsidx}.AllRightLicks{m,d}{trialvec(tridx)},repmat(tridx+smallcount+bigcount,[1,length( AllData{dsidx}.AllRightLicks{m,d}{trialvec(tridx)})]),'k.','MarkerSize',20);
                end
                
                if p == 2
                    tmpa = [ AllData{dsidx}.AllOptoOn{m,:}];
                    tmpb = [ AllData{dsidx}.AllOpto{m,:}];
                    tmpc = [ AllData{dsidx}.AllOptoOff{m,:}];
                    timinga = nanmean(tmpa(tmpb==p-1))*1000;
                    timingb = nanmean(tmpc(tmpb==p-1))*1000;
                    h=patch([timinga timingb timingb timinga]',[bigcount+smallcount bigcount+smallcount bigcount+smallcount+length(trialvec) bigcount+smallcount+length(trialvec)]','b');
                    alpha(0.3)
                    set(h,'EdgeAlpha',0)
                elseif p == 3
                    tmpa = [ AllData{dsidx}.AllOptoOn{m,:}];
                    tmpb = [ AllData{dsidx}.AllOpto{m,:}];
                    tmpc = [ AllData{dsidx}.AllOptoOff{m,:}];
                    timinga = nanmin(tmpa(tmpb==p-1))*1000;
                    tmpd = [ AllData{dsidx}.AllRT{m,:}];
                    if timinga>1.5*nanmean(tmpd(tmpb==p-1))
                        timinga = nanmin(tmpd(tmpb==p-1));
                    end
                    timingb = nanmean(tmpc(tmpb==p-1))*1000;
                    h=patch([timinga timingb timingb timinga]',[bigcount+smallcount bigcount+smallcount bigcount+smallcount+length(trialvec) bigcount+smallcount+length(trialvec)]','b');
                    alpha(0.3)
                    set(h,'EdgeAlpha',0)
                end
                
                smallcount = smallcount+length(trialvec);
                
            end
            bigcount = bigcount+smallcount;
            
            line([-1000 5000],[bigcount bigcount],'color',[0 0 0])
            
            xlim([-1000 5000])
            
            try
                ylim([1 bigcount])
            catch ME
                continue
            end
        end
        ylabel('Nr Trials')
        xlabel('Time (ms)')
        title(MouseOpt{m})
    end
end

%% Nr. Licks per mouse

boxsymb = ''; %if empty no outliers are shown
maxmouse = max(cellfun(@length,AllMicePerGroup));
LicksDuringBase = nan(length(AllMicePerGroup),maxmouse,3,3);
AllLicksDuringBase = cell(length(AllMicePerGroup),maxmouse,3,3);
maxdays = max(cell2mat(cellfun(@length,DateOpt,'UniformOutput',0)));
TrialNr =cell(length(AllMicePerGroup),maxmouse,3);
for dsidx =1:length(AllMicePerGroup)
    MouseOpt = AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    for m = 1:length(MouseOpt)
        AllLaserSetting = [AllData{dsidx}.LaserSetting{m,:}];
        trials2take = find(ismember(AllLaserSetting,LaserOption2Use{dsidx,m}));
        AllOptoOn = [AllData{dsidx}.AllOptoOn{m,:}];
        AllOptoOff = [AllData{dsidx}.AllOptoOff{m,:}];
        AllReaction =  [AllData{dsidx}.AllReaction{m,:}];
        AllOpto = [AllData{dsidx}.AllOpto{m,:}];
        TrialIdx = [AllData{dsidx}.TrialNumber{m,:}];% 
        OptoOnTime =  AllOptoOn;%(strcmp(AllReaction,'Hit'));
        OptoOffTime =  AllOptoOff;%(strcmp(AllReaction,'Hit'));
        OptoOn = AllOpto;%(strcmp(AllReaction,'Hit'));
        
        AllLeftLicks = [AllData{dsidx}.AllLeftLicks{m,:}];
        AllRightLicks = [AllData{dsidx}.AllRightLicks{m,:}];
        
        LeftLicks =AllLeftLicks;%(strcmp(AllReaction,'Hit'));
        RightLicks = AllRightLicks;%(strcmp(AllReaction,'Hit'));
        
        LeftLicksWithOpto1 = [];
        LeftLicksWithOpto2 = [];
        LeftLicksNoOpto = [];
        RightLicksWithOpto1 = [];
        RightLicksWithOpto2 = [];
        RightLicksNoOpto = [];
        TotalLicksNoOpto = [];
        TotalLicksWithOpto1 = [];
        TotalLicksWithOpto2 = [];
        TrialCount = zeros(3,1);
        
        for i =1:length(trials2take)
            % Save Left licks
            AllLicksDuringBase{dsidx,m,1,OptoOn(trials2take(i))+1}= [AllLicksDuringBase{dsidx,m,1,OptoOn(trials2take(i))+1} length(LeftLicks{trials2take(i)}((LeftLicks{trials2take(i)}>=-500 & LeftLicks{trials2take(i)} <= 0)))];
            % Save Right licks
            AllLicksDuringBase{dsidx,m,2,OptoOn(trials2take(i))+1}= [AllLicksDuringBase{dsidx,m,2,OptoOn(trials2take(i))+1} length(RightLicks{trials2take(i)}((RightLicks{trials2take(i)}>=-500 & RightLicks{trials2take(i)}<=0)))];
            % Save Total licks
            AllLicksDuringBase{dsidx,m,3,OptoOn(trials2take(i))+1}= [AllLicksDuringBase{dsidx,m,3,OptoOn(trials2take(i))+1} length(LeftLicks{trials2take(i)}((LeftLicks{trials2take(i)}>=-500 & LeftLicks{trials2take(i)} <=0)))+length(RightLicks{trials2take(i)}((RightLicks{trials2take(i)}>=-500 & RightLicks{trials2take(i)}<=0)))];
            % Save trial number        
            TrialNr{dsidx,m,OptoOn(trials2take(i))+1} = [TrialNr{dsidx,m,OptoOn(trials2take(i))+1} TrialIdx(trials2take(i))];

        end
    end
end
%% Example mouse
    MouseOpt = AllMicePerGroup{2};

midx = find(ismember(MouseOpt,'Lully'));
LeftLicks = [AllData{2}.AllLeftLicks{midx,:}];
RightLicks =  [AllData{2}.AllRightLicks{midx,:}];
AllOpto = [AllData{2}.AllOpto{midx,:}];

LeftLicks(AllOpto>0) = [];
RightLicks(AllOpto>0) = [];

figure('name','ExampleLickPlot')
hold on
arrayfun(@(X) scatter(LeftLicks{X}./1000,repmat(X,1,length(LeftLicks{X})),10,[1 0 0],'filled'),1:length(LeftLicks))
arrayfun(@(X) scatter(RightLicks{X}./1000,repmat(X,1,length(RightLicks{X})),10,[0 0 1],'filled'),1:length(RightLicks))
xlim([-0.5 3.20])
xlabel('time (s)')
ylabel('Trial')
makepretty
offsetAxes

%%
figure('name',['LicksDuringBaseline'])
cols = [0 0 0; 1 0 0; 0 0 1];
for dsidx =1:length(AllMicePerGroup)
    MouseOpt = AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    
    clear h
    for m = 1:length(MouseOpt)
        subplot(maxmouse,length(AllMicePerGroup),(m-1)*length(AllMicePerGroup)+dsidx)
        tmp = cat(2,AllLicksDuringBase{dsidx,m,:,:});
        normtmp = quantile(tmp,0.95);
        for sidx=1:3
            for optid=1:3
                % Licks
                hold on
                tmp =  AllLicksDuringBase{dsidx,m,sidx,optid};
                if NormalizeLicks
                    tmp = tmp./normtmp;
                end
                
                boxplot(tmp,'boxstyle','filled','positions',[(sidx-1).*1+(optid-1).*0.2+0.8],'colors',cols(optid,:),'symbol',boxsymb);
                h(optid) = line([(sidx-1).*1+(optid-1).*0.2+0.7 (sidx-1).*1+(optid-1).*0.2+0.9],[nanmean(tmp) nanmean(tmp)],'color',cols(optid,:),'LineWidth',2,'LineStyle','--');
                line([(sidx-1).*1+(optid-1).*0.2+0.7 (sidx-1).*1+(optid-1).*0.2+0.9],[nanmedian(tmp) nanmedian(tmp)],'color',cols(optid,:),'LineWidth',2,'LineStyle','--');
                
            end
        end
        set(gca,'XTick',[1:3],'XTickLabel',{'Left','Right','Total'},'xlim',[0.5 3.5])
        title(MouseOpt{m})
        if NormalizeLicks
            ylabel(['Normalized Lick Count'])
        else
            ylabel(['Lick Count'])
        end
        box off
    end
end
legend([h(:)],'No Opto','Opto Early','Opto 1st Lick')
%% During Visual
figure('name','AUCcurves')
hold on
colOpt = [0 0 0; 1 0 0; 0 0 1];
TimeCut = 1; %only take 1st second
LicksDuringVisual = nan(length(AllMicePerGroup),maxmouse,3,3);
AllLicksDuringVis = cell(length(AllMicePerGroup),maxmouse,3,3);
AllLicksDuringVisPerStim = cell(length(AllMicePerGroup),maxmouse,3,3,2);
FirstLicksDuringVisPerStim =  cell(length(AllMicePerGroup),maxmouse,3,2);
TrialNr =cell(length(AllMicePerGroup),maxmouse,3);
FractionMisses = nan(length(AllMicePerGroup),maxmouse,2); % Without/With opto 
AUCs = nan(length(AllMicePerGroup),maxmouse);
Switches = cell(length(AllMicePerGroup),maxmouse,3,2); %-1 right to left, 1 left to right
for dsidx =1:length(AllMicePerGroup)
    MouseOpt = AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    
    for m = 1:length(MouseOpt)
        AllLaserSetting = [AllData{dsidx}.LaserSetting{m,:}];
        AllReaction =  [AllData{dsidx}.AllReaction{m,:}];
        trials2take = find(ismember(AllLaserSetting,LaserOption2Use{dsidx,m}) & ~strcmp(AllReaction,'Miss'));

        AllOptoOn = [AllData{dsidx}.AllOptoOn{m,:}];
        AllOptoOff = [AllData{dsidx}.AllOptoOff{m,:}];
        AllOpto = [AllData{dsidx}.AllOpto{m,:}];
        TrialIdx = [AllData{dsidx}.TrialNumber{m,:}];%
        for optid=1:2
            FractionMisses(dsidx,m,optid) = sum(ismember(AllReaction,'Miss') & (AllOpto==optid-1))./sum(AllOpto==optid-1);
% 
%             [~,sortidx] = sort(TrialIdx,'ascend'); 
%             FractionMisses(dsidx,m,optid,1) = sum(ismember(AllReaction,'Miss') & (AllOpto==optid-1) & ismember(TrialIdx, TrialIdx(sortidx(1:round(0.1*length(sortidx))))))./sum(AllOpto==optid-1 & ismember(TrialIdx, TrialIdx(sortidx(1:round(0.1*length(sortidx))))));
% 
%              [~,sortidx] = sort(TrialIdx,'descend'); 
%             FractionMisses(dsidx,m,optid,2) = sum(ismember(AllReaction,'Miss') & (AllOpto==optid-1) & ismember(TrialIdx, TrialIdx(sortidx(1:round(0.1*length(sortidx))))))./sum(AllOpto==optid-1 & ismember(TrialIdx, TrialIdx(sortidx(1:round(0.1*length(sortidx))))));
%   
        end
        idx = logical(AllOpto<2);
        labels = ~double(ismember(AllReaction(idx),'Miss'));
        scores = AllOpto(idx);
        [x,y,~,AUCs(dsidx,m)] = perfcurve(labels,scores,1);
        plot(x,y,'color',colOpt(dsidx,:))


        AllRT = [AllData{dsidx}.AllRT{m,:}];
        OptoOnTime =  AllOptoOn(trials2take);%(AllRT>=100);%(strcmp(AllReaction,'Hit'));
        OptoOffTime =  AllOptoOff(trials2take);%(AllRT>=100);%(strcmp(AllReaction,'Hit'));
        OptoOn = AllOpto(trials2take);%(AllRT>=100);%(strcmp(AllReaction,'Hit'));
        AllStimSide = [AllData{dsidx}.AllSides{m,:}];
        AllStimSide = AllStimSide(trials2take);
        TrialIdx = TrialIdx(trials2take);

        
        AllLeftLicks = [AllData{dsidx}.AllLeftLicks{m,:}];
        AllRightLicks = [AllData{dsidx}.AllRightLicks{m,:}];
        RT=AllRT(trials2take);%(AllRT>=100);
        LeftLicks =AllLeftLicks(trials2take);%(AllRT>=100);%(strcmp(AllReaction,'Hit'));
        RightLicks = AllRightLicks(trials2take);%(AllRT>=100);%(strcmp(AllReaction,'Hit'));
        %         LeftLicksWithOpto1 = [];
        %         LeftLicksWithOpto2 = [];
        %         LeftLicksNoOpto = [];
        %         RightLicksWithOpto1 = [];
        %         RightLicksWithOpto2 = [];
        %         RightLicksNoOpto = [];
        %         TotalLicksNoOpto = [];
        %         TotalLicksWithOpto1 = [];
        %         TotalLicksWithOpto2 = [];
        
        TrialCount = zeros(3,1);
        
        for i = 1:length(OptoOn)
            % Save left licks
            AllLicksDuringVis{dsidx,m,1,OptoOn(i)+1}= [AllLicksDuringVis{dsidx,m,1,OptoOn(i)+1} length(LeftLicks{i}((LeftLicks{i}>0 & LeftLicks{i} < TimeCut*1000)))];
            
            % Save right licks
            AllLicksDuringVis{dsidx,m,2,OptoOn(i)+1}= [AllLicksDuringVis{dsidx,m,2,OptoOn(i)+1} length(RightLicks{i}((RightLicks{i}>0  & RightLicks{i}<TimeCut*1000)))];
            
            % Save all licks
            AllLicksDuringVis{dsidx,m,3,OptoOn(i)+1}= [AllLicksDuringVis{dsidx,m,3,OptoOn(i)+1} length(LeftLicks{i}((LeftLicks{i}>0 & LeftLicks{i} < TimeCut*1000)))+length(RightLicks{i}((RightLicks{i}>0 & RightLicks{i}<TimeCut*1000)))];
            
            % Save trial number
            TrialNr{dsidx,m,OptoOn(i)+1} = [TrialNr{dsidx,m,OptoOn(i)+1} TrialIdx(i)];

     
            %Save all licks for left and right stimulus seperately
            if strcmp(AllStimSide{i},'left')
                stimsid = 1;
            elseif strcmp(AllStimSide{i},'right')
                stimsid=2;
            end
            
            % Save left licks
            AllLicksDuringVisPerStim{dsidx,m,1,OptoOn(i)+1,stimsid}= [AllLicksDuringVisPerStim{dsidx,m,1,OptoOn(i)+1,stimsid} length(LeftLicks{i}((LeftLicks{i}>0 & LeftLicks{i} < TimeCut*1000)))];
            
            % Save right licks
            AllLicksDuringVisPerStim{dsidx,m,2,OptoOn(i)+1,stimsid}= [AllLicksDuringVisPerStim{dsidx,m,2,OptoOn(i)+1,stimsid} length(RightLicks{i}((RightLicks{i}>0  & RightLicks{i}<TimeCut*1000)))];
            
            % Save all licks
            AllLicksDuringVisPerStim{dsidx,m,3,OptoOn(i)+1,stimsid}= [AllLicksDuringVisPerStim{dsidx,m,3,OptoOn(i)+1,stimsid} length(LeftLicks{i}((LeftLicks{i}>0 & LeftLicks{i} < TimeCut*1000)))+length(RightLicks{i}((RightLicks{i}>0 & RightLicks{i}<TimeCut*1000)))];
            
          
            % indicate if the first lick (after 0) was a left (1) or a right(2) lick
            licktimeline = nan(1,1000);
            licktimeline(LeftLicks{i}(LeftLicks{i}>0&LeftLicks{i}<=1000))=1;
            licktimeline(RightLicks{i}(RightLicks{i}>0&RightLicks{i}<=1000))=2;
            licktimeline(isnan(licktimeline))=[];
            if ~isempty(licktimeline)
                Switches{dsidx,m,OptoOn(i)+1,stimsid} = [Switches{dsidx,m,OptoOn(i)+1,stimsid} sum(diff(licktimeline))];
                
                  % was first lick after stim onset left (1) or right (2)?
                  if LeftLicks{i}(find(LeftLicks{i}>0,1,'first')) < RightLicks{i}(find(RightLicks{i}>0,1,'first')) 
                      FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} = [FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} 1];
                  elseif LeftLicks{i}(find(LeftLicks{i}>0,1,'first')) > RightLicks{i}(find(RightLicks{i}>0,1,'first'))
                      FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} = [FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} 2];
                  elseif isempty(find(LeftLicks{i}>0)) && ~isempty(find(RightLicks{i}>0)) % Right first (or only)
                      FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} = [FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} 2];
                  elseif isempty(find(RightLicks{i}>0)) && ~isempty(find(LeftLicks{i}>0)) % Left first (or only)
                      FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} = [FirstLicksDuringVisPerStim{dsidx,m,OptoOn(i)+1,stimsid} 1];
                  else
                      keyboard
                  end
            end
            
        end
        
    end
end
makepretty
offsetAxes
title('p(Response|Opto)')
xlabel('False assumption a response was made')
ylabel('True assumption a repsonse was made')
%%
% FractionMisses = nan(length(AllMicePerGroup),maxmouse,2,2); % Without/With opto and early/late
figure('name','Fraction Omissions')
colOpt=[0 0 0; 1 0 0; 0 0 1];

nPerGroup = cellfun(@numel,AllMicePerGroup);
% hold on
% tmp = squeeze(FractionMisses(dsid,:,:));
% scatter(tmp(:,1),tmp(:,2),20,colOpt(dsid,:),'filled')
%     subplot(1,3,dsid)
h=barwitherr(squeeze(nanstd(FractionMisses,[],2))'./sqrt(nPerGroup-1),squeeze(nanmean(FractionMisses,2))');
hold on
set(gca,'XTick',1:2,'XTickLabel',{'Without','With'})
for hid = 1:length(h)
    h(hid).FaceColor = colOpt(hid,:);
    for optid = 1:2
        scatter(repmat(h(hid).XEndPoints(optid),1,size(FractionMisses,2)),squeeze(FractionMisses(hid,:,optid)),25,colOpt(hid,:),'filled')
    end
end
legend(AllCondNames)
makepretty
offsetAxes
ylabel('Fraction omissions')
%     offsetAxes
g1 = repmat(AllCondNames',[1,size(FractionMisses,2),2]);
g2 = repmat(OptoOpt(1:2)',[1,size(FractionMisses,2),3]);
g2 = permute(g2,[3,2,1]);
g3 = repmat(1:3*7,[1,2]);
g3 = reshape(g3,size(g1));
idx = ~isnan(FractionMisses(:));

t = table(FractionMisses(idx),g1(idx),g2(idx),g3(idx),'VariableNames',{'Omissions','Genotype','Opto','Mouse'});

result = fitglme(t,'Omissions ~ 1+ Genotype*Opto + (1|Mouse)','FitMethod','Laplace')
result.Rsquared.Adjusted
anova(result)
%
emm = emmeans(result, {'Opto','Genotype'}, 'effects');
for genid = 1:3

    % hypothesis testing of the interaction
    L_HL = zeros(1,height(emm.table));
    L_HL(ismember(emm.table.Opto,OptoOpt{2}) & ismember(emm.table.Genotype,AllCondNames{genid})) = 1;
    L_HL(ismember(emm.table.Opto,OptoOpt{1}) & ismember(emm.table.Genotype,AllCondNames{genid})) = -1;

    H0_Stim = contrasts_wald(result,emm,L_HL);
    p=H0_Stim.pVal;
    if p>0.5 %Make a two sided test
        p=1-p;
    end
    p=p*3; % correction multiple comparisons

    if p<0.05
        sigstar({[h(genid).XEndPoints(1) h(genid).XEndPoints(2)]},p)
    end
end
for optid = 1:2
    for genid = 1:3
        for genid2 = 2:3
            if genid2<=genid
                continue
            end
            % hypothesis testing of the interaction
            L_HL = zeros(1,height(emm.table));
            L_HL(ismember(emm.table.Opto,OptoOpt{optid}) & ismember(emm.table.Genotype,AllCondNames{genid})) = 1;
            L_HL(ismember(emm.table.Opto,OptoOpt{optid}) & ismember(emm.table.Genotype,AllCondNames{genid2})) = -1;

            H0_Stim = contrasts_wald(result,emm,L_HL);
            p=H0_Stim.pVal;
            if p>0.5 %Make a two sided test
                p=1-p;
            end
            p=p*2*3; % correction multiple comparisons

            if p<0.05
                sigstar({[h(genid).XEndPoints(1) h(genid).XEndPoints(2)]},p)
            end
        end
    end
end


%%
g1 = [];
figure('name','AUC values')
ps = nan(1,3);
for dsid = 1:3
    hold on
    scatter(dsid+rand(1,size(AUCs,2))*0.1-0.05,AUCs(dsid,:),10,colOpt(dsid,:),'filled')
    %     [~,ps(dsid)] = ttest(AUCs(dsid,:),0.5);
    [ps(dsid)] = signrank(AUCs(dsid,:),1);

    g1 = [g1 repmat(dsid,1,length(AUCs(dsid,:)))];
    line([dsid-0.1 dsid+0.1],[nanmedian(AUCs(dsid,:)) nanmedian(AUCs(dsid,:))],'color',colOpt(dsid,:))
end
AUCan = AUCs(:);
g1(isnan(AUCan)) = [];
AUCan(isnan(AUCan))= [];

p = anovan(AUCan,{g1})

set(gca,'XTick',1:3','XTickLabel',AllCondNames)
ylim([0 1])
ylabel('AUC')
title('AUC')
makepretty
offsetAxes


%% Decrease in licks with trial number?
ColOpt = [0 0 0; 1 0 0];

AvgLicks = figure('name','AvgLicksW/O Opto')
for dsidx = 1:length(AllMicePerGroup)
    MouseOpt = AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    LicksPerCond = nan(2,2,length(MouseOpt)); % OptoOff/ON, first/last 20 trials, mouse
    figure('name',['TrialNrVSLicks ' AllCondNames{dsidx}])

    for m = 1:length(MouseOpt)
        subplot(ceil(sqrt(length(MouseOpt))),round(sqrt(length(MouseOpt))),m)
        hold on

        for optid=1:2
            scatter(TrialNr{dsidx,m,optid},AllLicksDuringVis{dsidx,m,3,optid},10,ColOpt(optid,:),'filled')
            [sortedtr,sortidx] = sort(TrialNr{dsidx,m,optid},'ascend');
            LicksPerCond(optid,1,m) = nanmean(AllLicksDuringVis{dsidx,m,3,optid}(sortidx(1:30)));
            [sortedtr,sortidx] = sort(TrialNr{dsidx,m,optid},'descend');
            LicksPerCond(optid,2,m) = nanmean(AllLicksDuringVis{dsidx,m,3,optid}((sortidx(1:30))));
        end
    end

    figure(AvgLicks)
    subplot(1,3,dsidx)
    h = barwitherr(nanstd(LicksPerCond,[],3)./sqrt(length(MouseOpt)-1)',nanmean(LicksPerCond,3)');
    for hid = 1:length(h)
        h(hid).FaceColor = ColOpt(hid,:);
    end

    title(AllCondNames{dsidx})
    set(gca,'XTick',1:2,'XTickLabels',{'First 30','Last 30'})
    xlabel('Trial Number')
    ylabel('Licks/sec')
    makepretty
    offsetAxes

    legend('Opto Off','Opto On')



end

%% During Visual Plot

figure('name',['LicksDuringVisual'])
for dsidx =1:length(AllMicePerGroup)
    MouseOpt = AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    clear h
    for m = 1:length(MouseOpt)
        subplot(maxmouse,length(AllMicePerGroup),(m-1)*length(AllMicePerGroup)+dsidx)
        tmp = cat(2,AllLicksDuringVis{dsidx,m,:,:});
        normtmp = quantile(tmp,0.95);
        for sidx=1:3
            for optid=1:3
                % Licks
                hold on
                tmp =  AllLicksDuringVis{dsidx,m,sidx,optid};
                if NormalizeLicks
                    tmp = tmp./normtmp;
                end
                
                boxplot(tmp,'boxstyle','filled','positions',[(sidx-1).*1+(optid-1).*0.2+0.8],'colors',cols(optid,:),'symbol',boxsymb);
                h(optid) = line([(sidx-1).*1+(optid-1).*0.2+0.7 (sidx-1).*1+(optid-1).*0.2+0.9],[nanmean(tmp) nanmean(tmp)],'color',cols(optid,:),'LineWidth',2,'LineStyle','--');
                line([(sidx-1).*1+(optid-1).*0.2+0.7 (sidx-1).*1+(optid-1).*0.2+0.9],[nanmedian(tmp) nanmedian(tmp)],'color',cols(optid,:),'LineWidth',2,'LineStyle','--');
                
            end
        end
        set(gca,'XTick',[1:3],'XTickLabel',{'Left','Right','Total'},'xlim',[0.5 3.5])
        title(MouseOpt{m})
        if NormalizeLicks
            ylabel(['Normalized Licks (Hz)'])
        else
            ylabel(['Licks Count'])
        end
        box off
    end
end
legend([h(:)],'No Opto','Opto Early')


%% Lick errorbars
Cols = cols;
% distinguishable_colors(3);
OffAvgBaseline = nan(2,length(AllMice2Use));

clear h
scatterfigure = figure;
for dsidx =1:length(AllMicePerGroup)
    MouseOpt =AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    for sidx=1:2
        avgoff = [];
        avgon=[];
        for m= 1:length(MouseOpt)
            subplot(2,2,sidx)
            hold on
            off = nanmean(AllLicksDuringBase{dsidx,m,sidx,1});
            avgoff= [avgoff off];
            offneg = nanstd(AllLicksDuringBase{dsidx,m,sidx,1})./sqrt(length(AllLicksDuringBase{dsidx,m,sidx,1})-1);
            %             offneg = sqrt(nanmean(AllLicksDuringBase{dsidx,m,sidx,1}));
            offpos = offneg;
            OffAvgBaseline(sidx,ismember(AllMice2Use,MouseOpt{m})) = off;

            on =  nanmean(AllLicksDuringBase{dsidx,m,sidx,2});
            avgon = [avgon on];
            onneg = nanstd(AllLicksDuringBase{dsidx,m,sidx,2})./sqrt(length(AllLicksDuringBase{dsidx,m,sidx,2})-1);
            %             onneg = sqrt(nanmean(AllLicksDuringBase{dsidx,m,sidx,2}));
            onpos = onneg;
            
            
            h(dsidx) = errorbar(off,on,onneg,onpos,offneg,offpos,'MarkerSize',20,'color',Cols(dsidx,:),'Marker','.','LineWidth',0.5)
        end
        %         h(dsidx) = errorbar(nanmean(avgoff),nanmean(avgon),...
        %             nanstd(avgoff)./sqrt(length(avgoff)-1),nanstd(avgoff)./sqrt(length(avgoff)-1),...
        %             nanstd(onpos)./sqrt(length(onpos)-1),nanstd(onpos)./sqrt(length(onpos)-1),'MarkerSize',30,'color',Cols(dsidx,:),'Marker','.','LineWidth',0.5)
        
    end
end
for sidx=1:2
    
    subplot(2,2,sidx)
    hold on
    axis square
    line([-0.5 5],[-0.5 5],'LineWidth',1,'LineStyle','--')
    xlim([0 4])
    ylim([0 6])
    
    box off
    title(['Baseline Licks ' SideOptNames{sidx}])
    xlabel('Light Off')
    ylabel('Light On')
end

OptoOffFig = figure('name','OptoOff')
subplot(2,2,1)
hv = violinplot(OffAvgBaseline');
for hid = 1:length(hv)
hv(hid).ViolinAlpha = {[1]};
end
title('Baseline Licks')
set(gca,'XTick',1:2,'XTickLabels',SideOpt)
ylabel('Licks/sec')
makepretty
offsetAxes


%% Visual
figure(scatterfigure)
OffAvgVis = nan(2,length(AllMice2Use));
for dsidx =1:length(AllMicePerGroup)
    MouseOpt =AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    for sidx=1:2
        avgoff = [];
        avgon = [];
        for m= 1:length(MouseOpt)
            subplot(2,2,sidx+2)
            hold on
            off = nanmean(AllLicksDuringVis{dsidx,m,sidx,1});
            offneg = nanstd(AllLicksDuringVis{dsidx,m,sidx,1})./sqrt(length(AllLicksDuringVis{dsidx,m,sidx,1})-1);
            
            %             offneg = sqrt(nanmean(AllLicksDuringVis{dsidx,m,sidx,1}));
            offpos = offneg;
            avgoff = [avgoff off];
            OffAvgVis(sidx,ismember(AllMice2Use,MouseOpt{m})) = off;
            on =  nanmean(AllLicksDuringVis{dsidx,m,sidx,2});
            onneg = nanstd(AllLicksDuringVis{dsidx,m,sidx,2})./sqrt(length(AllLicksDuringVis{dsidx,m,sidx,2})-1);
            avgon = [avgon on];
            %             onneg = sqrt(nanmean(AllLicksDuringVis{dsidx,m,sidx,2}));
            onpos = onneg;
            
            h(dsidx) = errorbar(off,on,onneg,onpos,offneg,offpos,'MarkerSize',20,'color',Cols(dsidx,:),'Marker','.','LineWidth',0.5)
        end
    end
end
for sidx=1:2
    
    subplot(2,2,sidx+2)
    hold on
    axis square
    line([-0.5 8],[-0.5 8],'LineWidth',1,'LineStyle','--')
    xlim([0 4])
    ylim([0 6])
    
    box off
    title(['Visual Licks ' SideOptNames{sidx}])
    xlabel('Light Off')
    ylabel('Light On')
end

legend([h(:)],AllCondNames)


figure(OptoOffFig) 
subplot(2,2,2)
hv = violinplot(OffAvgVis');
for hid = 1:length(hv)
hv(hid).ViolinAlpha = {[1]};
end
title('Visual Licks')
set(gca,'XTick',1:2,'XTickLabels',SideOpt)
ylabel('Licks/sec')
makepretty
offsetAxes
%% Mixed linear model
%             off = nanmean(AllLicksDuringBase{dsidx,m,sidx,1});
licks = [];
momentvec = {};
genvec = {};
optvec = {};
sidevec = {};
mousevec = [];
mouseid=0;
for mid = 1:size(AllLicksDuringBase,2)
    for genid=1:3
        mouseid=1+mouseid;

        for optid=1:2
            for sidx=1:2
                % Base
                tmp = AllLicksDuringBase{genid,mid,sidx,optid};
                licks = [licks tmp];
                
                momentvec = cat(2,momentvec{:},repmat({'Base'},1,length(tmp)));
                genvec = cat(2,genvec{:},repmat(AllCondNames(genid),1,length(tmp)));
                optvec = cat(2,optvec{:},repmat(OptoOpt(optid),1,length(tmp)));
                sidevec = cat(2,sidevec{:},repmat(SideOpt(sidx),1,length(tmp)));
                mousevec = cat(2,mousevec,repmat(mouseid,1,length(tmp)));
                
                %Visual
                tmp = AllLicksDuringVis{genid,mid,sidx,optid};
                licks = [licks tmp];
                
                momentvec = cat(2,momentvec{:},repmat({'Visual'},1,length(tmp)));
                genvec = cat(2,genvec{:},repmat(AllCondNames(genid),1,length(tmp)));
                optvec = cat(2,optvec{:},repmat(OptoOpt(optid),1,length(tmp)));
                sidevec = cat(2,sidevec{:},repmat(SideOpt(sidx),1,length(tmp)));
                mousevec = cat(2,mousevec,repmat(mouseid,1,length(tmp)));
                
            end
        end
    end
end


%% Make bar-plots of the same data

% %for specific licks:
pidnames = {'Base','Visual'}
AvgMat = nan(3,2,2);
StdMat = nan(3,2,2);
for genid=1:3
    for sidx=1:2
        for pid=1:2
            AvgMat(genid,sidx,pid)=nanmean(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genid})&ismember(optvec,OptoOpt{2})))-...
                nanmean(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genid})&ismember(optvec,OptoOpt{1})));
         
            StdMat(genid,sidx,pid)=sqrt(nanvar(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genid})&ismember(optvec,OptoOpt{1})))./...
                sum(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})& ismember(genvec,AllCondNames{genid})&ismember(optvec,OptoOpt{1}))+...
                nanvar(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genid})&ismember(optvec,OptoOpt{2})))./...
                sum(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})& ismember(genvec,AllCondNames{genid})&ismember(optvec,OptoOpt{2})));
        end
    end
end

%% Draw
figure('name',[OptoOpt{2} 'BarPlots'])
h=barwitherr(reshape(StdMat,3,[])',...
    reshape(AvgMat,3,[])');
for hid = 1:length(h)
    set(h(hid),'EdgeColor','none','FaceColor',cols(hid,:))
end
set(gca,'XTickLabel',{'ContraLicks Baseline','IpsiLicks Baseline','ContraLicks Visual','IpsiLicks Visual'})
box off
ylabel([OptoOpt{2} ' - ' OptoOpt{1}])
for pid = 1:2
    % Epoch baseline
    idx = ismember(momentvec,pidnames{pid}) & ismember(optvec,OptoOpt(1:2));
    t = table(licks(idx)',genvec(idx)',optvec(idx)',sidevec(idx)',mousevec(idx)','VariableNames',{'Licks','Genotype','Opto','LickSide','Mouse'});

    result = fitglme(t,'Licks ~ 1+ Genotype*Opto*LickSide + (1|Mouse)','FitMethod','Laplace','Distribution','Poisson')
    result.Rsquared.Adjusted
    disp([pidnames{pid}])
    anova(result)

    emm = emmeans(result, {'Opto','Genotype'}, 'effects');

    for genid = 1:3
        for genid2 = 2:3
            if genid2<=genid
                continue
            end
            % hypothesis testing of the interaction
            L_HL = ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid}))'...
                -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid}))') - ...;
                ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid2}) )'...
                -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid2}))');
            H0_Stim = contrasts_wald(result,emm,L_HL);
            p=H0_Stim.pVal;
            if p>0.5 %Make a two sided test
                p=1-p;
            end
            p=p*2*2*2*3

            if p<0.05
                sigstar({[h(genid).XEndPoints((pid-1)*2+sidx) h(genid2).XEndPoints((pid-1)*2+sidx)]},p)
            end

        end
    end


    emm = emmeans(result, {'Opto','Genotype','LickSide'},'effects')
    for genid = 1:3
        % hypothesis testing of the interaction
        L_HL = ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid})&strcmp(emm.table.LickSide,SideOpt{1}))'...
            -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid})&strcmp(emm.table.LickSide,SideOpt{1}))') - ...;
            ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid})&strcmp(emm.table.LickSide,SideOpt{2}))'...
            -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid})&strcmp(emm.table.LickSide,SideOpt{2}))');
        H0_Stim = contrasts_wald(result,emm,L_HL);
        p=H0_Stim.pVal;
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        if p<0.05
            sigstar({[h(genid).XEndPoints((pid-1)*2+1) h(genid).XEndPoints((pid-1)*2+2)]},p)
        end
    end

end


%% Draw 
figure('name',[OptoOpt{optid} 'BarPlots'])
for genid=1:length(AllCondNames)
    hsub(genid) = subplot(1,length(AllCondNames),genid);
    h=barwitherr(squeeze(StdMat(genid,:,1))',...
        squeeze(AvgMat(genid,:,1))');
    for hid = 1:length(h)
        set(h(hid),'EdgeColor','none','FaceColor',cols(genid,:))
    end
    set(gca,'XTickLabel',{'ContraLicks Baseline','IpsiLicks Baseline'})
    box off
    ylabel([OptoOpt{2} ' - ' OptoOpt{1}])
    for pid = 1 %only baseline
        
        % Epoch baseline
        idx = ismember(momentvec,pidnames{pid}) & ismember(genvec,AllCondNames{genid});
        t = table(licks(idx)',optvec(idx)',sidevec(idx)',mousevec(idx)','VariableNames',{'Licks','Opto','Side','Mouse'});
        
        result = fitglme(t,'Licks ~ 1+ Opto*Side + (1|Mouse)','FitMethod','Laplace','Distribution','Poisson')
        result.Rsquared.Adjusted
        disp([pidnames{pid} ' ' AllCondNames{genid}])
        p=anova(result)
        
        if p.pValue(3)*3<0.05
            sigstar({[1 2]}, p.pValue(3)*3)
        end
        
    end
    title(AllCondNames{genid})
end
linkaxes(hsub(:))
     
%% Baseline vs Visual window effect?
opt2take = [1,2]

AvgMat = nan(3,2,2);
StdMat = nan(3,2,2);
figure('name',['BasevsVisual BarPlots'])
clear emm
for genidx=1:length(AllCondNames)
    
    idx = ismember(genvec,AllCondNames{genidx});
    t = table(licks(idx)',momentvec(idx)',sidevec(idx)',optvec(idx)',mousevec(idx)','VariableNames',{'Licks','Timewindow','LickSide','Opto','Mouse'});
    
    AllCondNames{genidx}
    result = fitglme(t,'Licks ~ 1+ Timewindow*LickSide*Opto + (1|Mouse)','FitMethod','Laplace','Distribution','Poisson')
    result.Rsquared.Adjusted
    anova(result)
    
    %% Graphics
    emm = emmeans(result, {'Opto','LickSide','Timewindow'}, 'effects');
    for sidx=1:2
        
        for pid=1:2
            AvgMat(genidx,sidx,pid)=nanmean(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{opt2take(2)})))-...
                nanmean(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{1})));
       
            StdMat(genidx,sidx,pid)=sqrt(nanvar(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{1})))./...
                sum(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})& ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{1}))+...
                nanvar(licks(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{opt2take(2)})))./...
                sum(ismember(momentvec,pidnames{pid})&ismember(sidevec,SideOpt{sidx})& ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{opt2take(2)})));
            
            
        end
    end
    hsub(genidx) = subplot(1,length(AllCondNames),genidx);
    h=barwitherr(squeeze(StdMat(genidx,:,:))',...
        squeeze(AvgMat(genidx,:,:))');
    for hid = 1:length(h)
        tmpcol = cols(genidx,:)+[0.5 0.5 0.5].*(hid-1);
        tmpcol(tmpcol>1)=1;
        set(h(hid),'EdgeColor','none','FaceColor',tmpcol)
    end
    set(gca,'XTickLabel',pidnames)
    box off
    ylabel([OptoOpt{opt2take(2)} ' - ' OptoOpt{1}])
    for pid = 1:2
        for genid = 1:2
            for genid2 = 2
                if genid2<=genid
                    continue
                end
                % hypothesis testing (use emm that has not been releveled)
                % joint test that USA is different
                L_HL = ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.Timewindow,pidnames{pid})&strcmp(emm.table.LickSide,SideOpt{genid}) )'...
                    -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.Timewindow,pidnames{pid})&strcmp(emm.table.LickSide,SideOpt{genid}) )') - ...;
                    ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.Timewindow,pidnames{pid})&strcmp(emm.table.LickSide,SideOpt{genid2}))'...
                    -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.Timewindow,pidnames{pid})&strcmp(emm.table.LickSide,SideOpt{genid2}))');
                H0_Stim = contrasts_wald(result,emm,L_HL);
                p=H0_Stim.pVal;
                if p>0.5 %Make a two sided test
                    p=1-p;
                end
                p=p*2*2*2*3 %Multiple comparisons
                
                if p<0.05
                    sigstar({[h(genid).XEndPoints(pid) h(genid2).XEndPoints(pid)]},p)
                end

            end
        end
    end
    
    for sidx = 1:2
        
        % hypothesis testing (use emm that has not been releveled)
        % joint test that USA is different
        L_HL = ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.Timewindow,pidnames{1})&strcmp(emm.table.LickSide,SideOpt{sidx}))'...
            -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.Timewindow,pidnames{1})&strcmp(emm.table.LickSide,SideOpt{sidx}))') - ...;
            ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.Timewindow,pidnames{2})&strcmp(emm.table.LickSide,SideOpt{sidx}))'...
            -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.Timewindow,pidnames{2})&strcmp(emm.table.LickSide,SideOpt{sidx}))');
        H0_Stim = contrasts_wald(result,emm,L_HL);
        p=H0_Stim.pVal;
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        p=p*2*2*3 %multiple comparisons
        
        if p<0.05
            sigstar({[h(sidx).XEndPoints(1) h(sidx).XEndPoints(2)]},p)
        end
        
    end
    title([AllCondNames{genidx} ' Licks'])
    legend(h,cellfun(@(X) ['Licks ' X],SideOpt,'UniformOutput',0))
end
linkaxes([hsub(:)])



%% Mixed linear model - Visual Epoch
opt2take = [1,2]
licks = [];
visualvec = {};
genvec = {};
optvec = {};
sidevec = {};
mousevec = [];
mouseid=1;
for mid = 1:size(AllLicksDuringVisPerStim,2)
    for genid=1:3
        
        for optid=1:2
            for sidx=1:2
                for vidx = 1:2
                    
                    %Visual
                    tmp = AllLicksDuringVisPerStim{genid,mid,sidx,opt2take(optid),vidx};
                    licks = [licks tmp];
                    
                    genvec = cat(2,genvec{:},repmat(AllCondNames(genid),1,length(tmp)));
                    optvec = cat(2,optvec{:},repmat(OptoOpt(opt2take(optid)),1,length(tmp)));
                    sidevec = cat(2,sidevec{:},repmat(SideOpt(sidx),1,length(tmp)));
                    mousevec = cat(2,mousevec,repmat(mouseid,1,length(tmp)));
                    visualvec = cat(2,visualvec{:},repmat(SideOpt(vidx),1,length(tmp)));
                end
            end
        end
        mouseid=1+mouseid;
    end
    mouseid-1
end

% %for specific licks:
%
AvgMat = nan(3,2,2);
StdMat = nan(3,2,2);
figure('name',['Visual BarPlots'])

for genidx=1:length(AllCondNames)
    
    idx = ismember(genvec,AllCondNames{genidx});
    t = table(licks(idx)',visualvec(idx)',sidevec(idx)',optvec(idx)',mousevec(idx)','VariableNames',{'Licks','VisualSide','LickSide','Opto','Mouse'});
    
    AllCondNames{genidx}
    result = fitglme(t,'Licks ~ 1+ VisualSide*LickSide*Opto + (1|Mouse)','FitMethod','Laplace','Distribution','Poisson')
    result.Rsquared.Adjusted
    anova(result)
    
    %% Graphics
    emm = emmeans(result, {'Opto','LickSide','VisualSide'}, 'effects');
    for sidx=1:2
        
        for pid=1:2
            AvgMat(genidx,sidx,pid)=nanmean(licks(ismember(visualvec,SideOpt{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{opt2take(2)})))-...
                nanmean(licks(ismember(visualvec,SideOpt{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{1})));
        
            StdMat(genidx,sidx,pid)=sqrt(nanvar(licks(ismember(visualvec,SideOpt{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{1})))./...
                sum(ismember(visualvec,SideOpt{pid})&ismember(sidevec,SideOpt{sidx})& ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{1}))+...
                nanvar(licks(ismember(visualvec,SideOpt{pid})&ismember(sidevec,SideOpt{sidx})&...
                ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{opt2take(2)})))./...
                sum(ismember(visualvec,SideOpt{pid})&ismember(sidevec,SideOpt{sidx})& ismember(genvec,AllCondNames{genidx})&ismember(optvec,OptoOpt{opt2take(2)})));
            
            
        end
    end
    hsub(genidx) = subplot(1,length(AllCondNames),genidx);
    h=barwitherr(squeeze(StdMat(genidx,:,:))',...
        squeeze(AvgMat(genidx,:,:))');
    for hid = 1:length(h)
        tmpcol = cols(genidx,:)+[0.5 0.5 0.5].*(hid-1);
        tmpcol(tmpcol>1)=1;
        set(h(hid),'EdgeColor','none','FaceColor',tmpcol)
    end
    set(gca,'XTickLabel',{'Contralateral Stimulus','Ipsilateral Stimulus'})
    box off
    ylabel([OptoOpt{opt2take(2)} ' - ' OptoOpt{1}])
    for pid = 1:2
        for genid = 1:2
            for genid2 = 2
                if genid2<=genid
                    continue
                end
                % hypothesis testing (use emm that has not been releveled)
                % joint test that USA is different
                L_HL = ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.VisualSide,SideOpt{pid})&strcmp(emm.table.LickSide,SideOpt{genid}) )'...
                    -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.VisualSide,SideOpt{pid})&strcmp(emm.table.LickSide,SideOpt{genid}) )') - ...;
                    ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.VisualSide,SideOpt{pid})&strcmp(emm.table.LickSide,SideOpt{genid2}))'...
                    -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.VisualSide,SideOpt{pid})&strcmp(emm.table.LickSide,SideOpt{genid2}))');
                H0_Stim = contrasts_wald(result,emm,L_HL);
                p=H0_Stim.pVal;
                if p>0.5 %Make a two sided test
                    p=1-p;
                end
                p=p*2*2*2*3
                
                if p<0.05
                    sigstar({[h(genid).XEndPoints(pid) h(genid2).XEndPoints(pid)]},p)
                end
                
                
                
                
            end
        end
    end
    
    for sidx = 1:2
        
        % hypothesis testing (use emm that has not been releveled)
        % joint test that USA is different
        L_HL = ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.VisualSide,SideOpt{1})&strcmp(emm.table.LickSide,SideOpt{sidx}))'...
            -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.VisualSide,SideOpt{1})&strcmp(emm.table.LickSide,SideOpt{sidx}))') - ...;
            ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.VisualSide,SideOpt{2})&strcmp(emm.table.LickSide,SideOpt{sidx}))'...
            -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.VisualSide,SideOpt{2})&strcmp(emm.table.LickSide,SideOpt{sidx}))');
        H0_Stim = contrasts_wald(result,emm,L_HL);
        p=H0_Stim.pVal;
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        p=p*2*2*3 %multiple comparisons
        
        if p<0.05
            sigstar({[h(sidx).XEndPoints(1) h(sidx).XEndPoints(2)]},p)
        end
        
    end
    title([AllCondNames{genidx} ' Licks'])
    legend(h,cellfun(@(X) ['Licks ' X],SideOpt,'UniformOutput',0))
end
linkaxes([hsub(:)])


%% 
figure('name',['Visual BarPlots'])
for genidx = 1:length(AllCondNames)
     hsub(genidx) = subplot(1,length(AllCondNames),genidx);
    h=barwitherr(squeeze(StdMat(genidx,:,:))',...
        squeeze(AvgMat(genidx,:,:))');
    for hid = 1:length(h)
        tmpcol = cols(genidx,:)+[0.5 0.5 0.5].*(hid-1);
        tmpcol(tmpcol>1)=1;
        set(h(hid),'EdgeColor','none','FaceColor',tmpcol)
    end
    set(gca,'XTickLabel',{'Contralateral Stimulus','Ipsilateral Stimulus'})
    box off
    ylabel([OptoOpt{opt2take(2)} ' - ' OptoOpt{1}])
end
linkaxes
for sidx=1:length(SideOpt)
    
    idx = ismember(sidevec,SideOpt{sidx});
    t = table(licks(idx)',visualvec(idx)',genvec(idx)',optvec(idx)',mousevec(idx)','VariableNames',{'Licks','VisualSide','Genotype','Opto','Mouse'});
    
    SideOpt{sidx}
    result = fitglme(t,'Licks ~ 1+ VisualSide*Genotype*Opto + (1|Mouse)','FitMethod','Laplace','Distribution','Poisson')
    result.Rsquared.Adjusted
    anova(result)
    
    % posthoc
    emm = emmeans(result, {'Opto','Genotype','VisualSide'}, 'effects');

    for visidx = 1:2
        for genid = 1:3
            L_HL = zeros(1,height(emm.table));
            L_HL(ismember(emm.table.VisualSide,SideOpt{visidx}) & ismember(emm.table.Genotype,AllCondNames{genid}) & ismember(emm.table.Opto,OptoOpt{1})) = -1;
            L_HL(ismember(emm.table.VisualSide,SideOpt{visidx}) & ismember(emm.table.Genotype,AllCondNames{genid}) & ismember(emm.table.Opto,OptoOpt{2})) = 1;
            H0_Stim = contrasts_wald(result,emm,L_HL);
            p=H0_Stim.pVal;
            if p>0.5 %Make a two sided test
                p=1-p;
            end
            p=p*2*2*3;

            if p<0.05
                subplot(1,3,genid)
                sigstar({[h(sidx).XEndPoints(visidx) h(sidx).XEndPoints(visidx)]},p)
            end

        end
    end  
end



%% Separates for licks
figure('name',['Visual BarPlots'])
for genidx = 1:length(AllCondNames)
     hsub(genidx) = subplot(1,length(AllCondNames),genidx);
    h=barwitherr(squeeze(StdMat(genidx,:,:)),...
        squeeze(AvgMat(genidx,:,:)));
    for hid = 1:length(h)
        tmpcol = cols(genidx,:)+[0.5 0.5 0.5].*(hid-1);
        tmpcol(tmpcol>1)=1;
        set(h(hid),'EdgeColor','none','FaceColor',tmpcol)
    end
    set(gca,'XTickLabel',{'Contralateral Licks','Ipsilateral Licks'})
    box off
    ylabel([OptoOpt{opt2take(2)} ' - ' OptoOpt{1}])
end
legend('Contra Visual','Ipsi Visual')
linkaxes

for sidx = 1:2 % Lick side
    idx = ismember(sidevec,SideOpt{sidx});
    t = table(licks(idx)',visualvec(idx)',genvec(idx)',optvec(idx)',mousevec(idx)','VariableNames',{'Licks','VisualSide','Genotype','Opto','Mouse'});

    SideOpt{sidx}
    result = fitglme(t,'Licks ~ 1+ VisualSide*Genotype*Opto + (1|Mouse)','FitMethod','Laplace','Distribution','Poisson')
    result.Rsquared.Adjusted
    anova(result)

    % posthoc
    emm = emmeans(result, {'Opto','Genotype','VisualSide'}, 'effects');
    for genid = 1:3
        subplot(1,3,genid)
        % hypothesis testing (use emm that has not been releveled)
        % joint test that USA is different
        L_HL = ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.VisualSide,SideOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid}))'...
            -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.VisualSide,SideOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid}))') - ...;
            ((strcmp(emm.table.Opto,OptoOpt{1})& strcmp(emm.table.VisualSide,SideOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid}))'...
            -(strcmp(emm.table.Opto,OptoOpt{opt2take(2)})& strcmp(emm.table.VisualSide,SideOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid}))');
        H0_Stim = contrasts_wald(result,emm,L_HL);
        p=H0_Stim.pVal;
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        p=p*3; %Multiple comparisons

        if p<0.05
            sigstar({[h(1).XEndPoints(sidx) h(2).XEndPoints(sidx)]},p)
        end

    end
end
title([AllCondNames{genidx} ' Licks'])
legend(h,cellfun(@(X) ['Licks ' X],SideOpt,'UniformOutput',0))

%% Split up right licks increase for D1 mice in Right and Left Stimulus
figure('name','LicksPerStimSide')
OffAvgVisPerSide = nan(2,2,length(AllMice2Use));

for dsidx = 1:3
    MouseOpt =AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    for m= 1:length(MouseOpt)
        for stimsideidx=1:2
            for licksidx=1:2
                subplot(2,2,(stimsideidx-1)*2+licksidx)
                hold on
                off = nanmean(AllLicksDuringVisPerStim{dsidx,m,licksidx,1,stimsideidx});
                offneg = nanstd(AllLicksDuringVisPerStim{dsidx,m,licksidx,1,stimsideidx})./sqrt(length(AllLicksDuringVisPerStim{dsidx,m,licksidx,1,stimsideidx})-1);
                
                offpos = offneg;
                OffAvgVisPerSide(stimsideidx,licksidx,ismember(AllMice2Use,MouseOpt{m})) = off;

                on =  nanmean(AllLicksDuringVisPerStim{dsidx,m,licksidx,2,stimsideidx});
                onneg = nanstd(AllLicksDuringVisPerStim{dsidx,m,licksidx,2,stimsideidx})./sqrt(length(AllLicksDuringVisPerStim{dsidx,m,licksidx,2,stimsideidx})-1);
                
                onpos = onneg;
                
                h(dsidx) = errorbar(off,on,onneg,onpos,offneg,offpos,'MarkerSize',20,'color',Cols(dsidx,:),'Marker','.','LineWidth',0.5)
            end
        end
    end
end
for stimsideidx=1:2
    for licksidx=1:2
        subplot(2,2,(stimsideidx-1)*2+licksidx)
        hold on
        axis square
        line([-0.5 8],[-0.5 8],'LineWidth',1,'LineStyle','--')
        xlim([0 4])
        ylim([0 7])
        
        box off
        title(['Visual Stimulus ' SideOptNames{stimsideidx} ', Licks ' SideOptNames{licksidx}])
        xlabel('Light Off')
        ylabel('Light On')
    end
end


figure(OptoOffFig) 
subplot(2,2,3)
hv = violinplot(squeeze(OffAvgVisPerSide(1,:,:))');
for hid = 1:length(hv)
hv(hid).ViolinAlpha = {[1]};
end
title('Left Visual Licks')
set(gca,'XTick',1:2,'XTickLabels',SideOpt)
ylabel('Licks/sec')

makepretty
offsetAxes

subplot(2,2,4)
hv = violinplot(squeeze(OffAvgVisPerSide(2,:,:))');
for hid = 1:length(hv)
hv(hid).ViolinAlpha = {[1]};
end
title('Right Visual Licks')
set(gca,'XTick',1:2,'XTickLabels',SideOpt)
ylabel('Licks/sec')

makepretty
offsetAxes
linkaxes
%% Make bar-plots of the same data
OnOff = squeeze(cell2mat(cellfun(@(X,Y) nanmean(X)-nanmean(Y),AllLicksDuringVisPerStim(:,:,:,2,:),AllLicksDuringVisPerStim(:,:,:,1,:),'UniformOutput',0)));
figure('name','BarPlots Per VisStim')
OnOff = reshape(OnOff(:,:,1:2,:),size(OnOff,1),size(OnOff,2),[]);
h=barwitherr(squeeze(nanstd(OnOff(:,:,:),[],2)./sqrt(sum(~isnan(OnOff(:,:,:)),2)-1))',...
    squeeze(nanmean(OnOff(:,:,:),2))');
for hid = 1:length(h)
    set(h(hid),'EdgeColor','none','FaceColor',cols(hid,:))
end
set(gca,'XTickLabel',{'ContraLicks ContraGrating','IpsiLicks ContraGrating','ContraLicks IpsiGrating','IpsiLicks IpsiGrating'})
box off
ylabel('Difference nr. Licks Opto On - Opto Off')


%ANOVA
OnOff = squeeze(cell2mat(cellfun(@(X,Y) nanmean(X)-nanmean(Y),AllLicksDuringVisPerStim(:,:,:,2,:),AllLicksDuringVisPerStim(:,:,:,1,:),'UniformOutput',0)));
OnOff = OnOff(:,:,1:2,:);
g1 = repmat(AllCondNames',[1,size(OnOff,2),length(SideOpt),2]);
g2 = repmat(SideOpt',[1,length(AllCondNames),size(OnOff,2),2]);
g2 = permute(g2,[2,3,1,4]);
g3 = repmat(SideOpt',[1,length(AllCondNames),size(OnOff,2),2]);
g3 = permute(g3,[2,3,4,1]);
g4 = repmat([1:size(OnOff,2)*length(AllCondNames)]',[1,length(SideOpt),2]);

OnOff = OnOff(:);
g1 = g1(:);
g2 = g2(:);
g3= g3(:);
g4=g4(:);

t = table(OnOff,g1,g2,g3,g4,'VariableNames',{'LicksDifference','Genotype','LeftvsRightLick','LeftvsRightGrating','Mouse'});
result = fitglme(t,'LicksDifference ~ 1+ Genotype*LeftvsRightLick*LeftvsRightGrating + (1|Mouse)','FitMethod','Laplace')
result.Rsquared.Adjusted
anova(result)
% result = fitglm(t,'LicksDifference~Genotype*LeftvsRightLick*LeftvsRightGrating');

%post-hoc test
%D1 versus D2
[h,p] = ttest2(OnOff(ismember(g1,'D1')),OnOff(ismember(g1,'D2')))
[h,p] = ttest2(OnOff(ismember(g1,'D1')),OnOff(ismember(g1,'Control')))
[h,p] = ttest2(OnOff(ismember(g1,'D2')),OnOff(ismember(g1,'Control')))

%for specific licks:
for pid = 1:2
    for sidx = 1:2
        [~,p(pid,sidx,1)] = ttest2(OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})),OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})))
        
        [~,p(pid,sidx,2)] = ttest2(OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})),OnOff(ismember(g1,'Control')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})))
        [~,p(pid,sidx,3)] = ttest2(OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})),OnOff(ismember(g1,'Control')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})))
    end
end

% D1 mice Contra versus IPsi:
[~,p] = ttest(OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{1})),OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{2})))
[~,p] = ttest(OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{1})),OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{2})))

[~,p] = ttest(OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{1})&ismember(g3,SideOpt{1})),OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{1})&ismember(g3,SideOpt{2})))
%% Look at fist lick after vis on
figure('name','FirstLickLeft percentage after VisOn')
Percvec = [];
genvec = {};
optvec = {};
sidevec = {};
mousevec = {};
for stimsideidx=1:2
    for dsidx = 1:3
        subplot(2,3,(stimsideidx-1)*3+dsidx)

        MouseOpt =AllMicePerGroup{dsidx};
        DateOpt = AllDatesPerGroup{dsidx};
        nFirstLicks = nan(2,2,length(MouseOpt)); % Left/Right X Off/On X mice
        for m= 1:length(MouseOpt)
            hold on
            for optid=1:2
                for lickid=1:2
                    nFirstLicks(lickid,optid,m) = (sum(FirstLicksDuringVisPerStim{dsidx,m,optid,stimsideidx}==lickid));
                end
                optvec = cat(1,optvec,OptoOpt{optid});
                mousevec = cat(1,mousevec, MouseOpt{m});
                Percvec = [Percvec  nFirstLicks(1,optid,m)./nansum(nFirstLicks(:,optid,m),1)];
                sidevec = cat(1,sidevec, SideOpt{stimsideidx});
                genvec = cat(1,genvec, AllCondNames{dsidx});

            end
        end
        PercFirstLicksLeft = squeeze(nFirstLicks(1,:,:)./(nansum(nFirstLicks,1)));
        h=violinplot(PercFirstLicksLeft',{'Off','On'},'ViolinColor',[0.5 0.5 0.5; Cols(dsidx,:)]);
        %         scatter(PercFirstLicksLeft(1,:),PercFirstLicksLeft(2,:),20,Cols(dsidx,:),'filled');
        for hid=1:length(h)
            h(hid).ViolinAlpha={[1]}; % For copying to inkscape get rid of transparency, fix there
            h(hid).EdgeColor='none'
        end
        
        %         h(dsidx) = plot(off,on,'MarkerSize',20,'color',Cols(dsidx,:),'Marker','.','LineWidth',0.5)
        axis square
        line(get(gca,'xlim'),[0.5 0.5],'LineWidth',1,'LineStyle','--','color',[0 0 0])
        
        
        %     xlim([0 1])
            ylim([0 1])
        %
        box off
        if dsidx==1
        ylabel(['Visual Stimulus ' SideOptNames{stimsideidx}])
          else
            ylabel('Proportion of first licks being left')
        end
        if stimsideidx==1
            title(AllCondNames{dsidx})
      
        end      
        FigureDefault
    end
     
    
end

% Epoch baseline
t = table(Percvec',genvec,optvec,sidevec,mousevec,'VariableNames',{'PercFirstLeft','Genotype','Opto','VisualSide','Mouse'});

result = fitglme(t,'PercFirstLeft ~ 1+ Genotype*Opto*VisualSide + (1|Mouse)','FitMethod','Laplace')
result.Rsquared.Adjusted
anova(result)

emm = emmeans(result, {'Opto','Genotype'}, 'effects');
%%
for genid = 1:3
    for genid2 = 2:3
        if genid2<=genid
            continue
        end
        % hypothesis testing of the interaction
        L_HL = ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid}))'...
            -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid}))') - ...;
            ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid2}) )'...
            -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid2}))');
        H0_Stim = contrasts_wald(result,emm,L_HL);
        p=H0_Stim.pVal;
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        p=p*3; % 3 tests

        if p<0.05
            disp(['Difference between opto on and opto off significantly different between ' AllCondNames{genid} ' and ' AllCondNames{genid2} ', p = ' num2str(round(p*100)/100)])
        end



    end
end

%% Differene opto on/off
emm = emmeans(result, {'Opto','Genotype'},'effects')
for genid = 1:3
    % hypothesis testing of the interaction
    L_HL = ((strcmp(emm.table.Opto,OptoOpt{2})&strcmp(emm.table.Genotype,AllCondNames{genid}))'...
        -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid}))');
    H0_Stim = contrasts_wald(result,emm,L_HL);
    p=H0_Stim.pVal;
    if p>0.5 %Make a two sided test
        p=1-p;
    end
    p=p*3; % 3 tests
    if p<0.05
        disp(['Difference between opto on and opto off significantly different for' AllCondNames{genid} ', p = ' num2str(round(p*100)/100)])
    end
   
end

%% Differene opto off between genotypes
for genid = 1:3
    for genid2 = 2:3
        if genid2<=genid
            continue
        end

        % hypothesis testing of the interaction
        L_HL = ((strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid}))'...
            -(strcmp(emm.table.Opto,OptoOpt{1})&strcmp(emm.table.Genotype,AllCondNames{genid2}))');
        H0_Stim = contrasts_wald(result,emm,L_HL);
        p=H0_Stim.pVal
        if p>0.5 %Make a two sided test
            p=1-p;
        end
        p=p*3; % 3 tests
        if p<0.05
            disp(['Difference between ' AllCondNames{genid} ' and ' AllCondNames{genid2} 'for opto off, p = ' num2str(round(p*100)/100)])
        end
    end
end

%% Look at SwitchRates
figure('name','Switches')
clear h
for dsidx = 1:3
    MouseOpt =AllMicePerGroup{dsidx};
    DateOpt = AllDatesPerGroup{dsidx};
    for m= 1:length(MouseOpt)
        for stimsideidx=1:2
            subplot(1,2,stimsideidx)
            hold on
            off = (sum(Switches{dsidx,m,1,stimsideidx}==1)-sum(Switches{dsidx,m,1,stimsideidx}==-1))./...
                (sum(Switches{dsidx,m,1,stimsideidx}==1)+sum(Switches{dsidx,m,1,stimsideidx}==-1));
            
            on = (sum(Switches{dsidx,m,2,stimsideidx}==1)-sum(Switches{dsidx,m,2,stimsideidx}==-1))./...
                (sum(Switches{dsidx,m,2,stimsideidx}==1)+sum(Switches{dsidx,m,2,stimsideidx}==-1));
            
            
            
            h(dsidx) = plot(off,on,'MarkerSize',20,'color',Cols(dsidx,:),'Marker','.','LineWidth',0.5)
        end
    end
end
for stimsideidx=1:2
    subplot(1,2,stimsideidx)
    hold on
    axis square
    line([-2 2],[-2 2],'LineWidth',1,'LineStyle','--')
    line([-2 2],[0 0],'LineWidth',1,'LineStyle','--')
    line([0 0],[-2 2],'LineWidth',1,'LineStyle','--')
    
    xlim([-2 2])
    ylim([-2 2])
    
    box off
    title(['Visual Stimulus ' SideOptNames{stimsideidx}])
    xlabel('Light Off Switch Index')
    ylabel('Light On Switch Index')
end

%% Same in barplots
figure('name','Switches BarPlots')
SwitchIndex = cell2mat(cellfun(@(X) (sum(X==1)-sum(X==-1))./(sum(X==1)+sum(X==-1)),Switches(:,:,1:2,:),'UniformOutput',0));
SwitchIndex = permute(SwitchIndex,[3,2,1,4]);
SwitchIndex = reshape(SwitchIndex,2,size(SwitchIndex,2),[]);
h=barwitherr(squeeze(nanstd(SwitchIndex(:,:,:),[],2)./sqrt(sum(~isnan(SwitchIndex(:,:,:)),2)-1))',...
    squeeze(nanmean(SwitchIndex(:,:,:),2))');

set(gca,'XTickLabel',{'Control ContraGrating','D1 ContraGrating','D2 ContraGrating','Control IpsiGrating','D1 IpsiGrating','D2 IpsiGrating',})
box off
ylabel('Switch Index')

%ANOVA
SwitchIndex = cell2mat(cellfun(@(X) (sum(X==1)-sum(X==-1))./(sum(X==1)+sum(X==-1)),Switches(:,:,1:2,:),'UniformOutput',0));

g1 = repmat(AllCondNames',[1,size(SwitchIndex,2),length(SideOpt),2]);
g2 = repmat(OptoOpt(1:2)',[1,length(AllCondNames),size(SwitchIndex,2),2]);
g2 = permute(g2,[2,3,1,4]);
g3 = repmat(SideOpt',[1,length(AllCondNames),size(SwitchIndex,2),2]);
g3 = permute(g3,[2,3,4,1]);
g4 = repmat([1:size(SwitchIndex,2)*length(AllCondNames)]',[1,length(SideOpt),2]);

SwitchIndex = SwitchIndex(:);
g1 = g1(:);
g2 = g2(:);
g3= g3(:);
g4=g4(:);
idx = ~isnan(SwitchIndex);

% just an ANOVA is suitable here
[p,st,tbl] = anovan(SwitchIndex(idx),{g1(idx),g2(idx),g3(idx)},'Model',4,'Varnames',{'Genotype','Opto','GratingSide'})
t = table(SwitchIndex,g1,g2,g3,g4,'VariableNames',{'SwitchIndex','Genotype','Opto','LeftvsRightGrating','Mouse'});

result = fitglme(t,'SwitchIndex ~ 1+ Genotype*Opto*LeftvsRightGrating + (1|Mouse)','FitMethod','Laplace')
result.Rsquared.Adjusted
anova(result)
% result = fitglm(t,'SwitchIndex~Genotype*Opto*LeftvsRightGrating');

%post-hoc test
%D1 versus D2
[h,p] = ttest2(OnOff(ismember(g1,'D1')),OnOff(ismember(g1,'D2')))
[h,p] = ttest2(OnOff(ismember(g1,'D1')),OnOff(ismember(g1,'Control')))
[h,p] = ttest2(OnOff(ismember(g1,'D2')),OnOff(ismember(g1,'Control')))

%for specific licks:
for pid = 1:2
    for sidx = 1:2
        [~,p(pid,sidx,1)] = ttest2(OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})),OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})))
        
        [~,p(pid,sidx,2)] = ttest2(OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})),OnOff(ismember(g1,'Control')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})))
        [~,p(pid,sidx,3)] = ttest2(OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})),OnOff(ismember(g1,'Control')&ismember(g2,SideOpt{sidx})&ismember(g3,SideOpt{pid})))
    end
end

% D1 mice Contra versus IPsi:
[~,p] = ttest(OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{1})),OnOff(ismember(g1,'D1')&ismember(g2,SideOpt{2})))
[~,p] = ttest(OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{1})),OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{2})))

[~,p] = ttest(OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{1})&ismember(g3,SideOpt{1})),OnOff(ismember(g1,'D2')&ismember(g2,SideOpt{1})&ismember(g3,SideOpt{2})))



%% Poisson GLM Baseline

TMP = []; %length(AllMicePerGroup),length(midx),Sides,OptoCond;
g1=[]; %Genotype
g2=[]; %Mouse
g3=[]; %Left/Right/Total Licks
g4=[]; %No Opto, Opto EARLY, (Opto Lick excluded)
for dsid=1:length(AllMicePerGroup)
    for midx=1:length(AllMicePerGroup{dsid})
        for sidx=1:2
            for optid=1:3
                TMP = [TMP AllLicksDuringBase{dsid,midx,sidx,optid}];
                g1 = [g1  repmat(AllCondNames(dsid),[1,length(AllLicksDuringBase{dsid,midx,sidx,optid})])];
                g2 = [g2 repmat(AllMicePerGroup{dsid}(midx),[1,length(AllLicksDuringBase{dsid,midx,sidx,optid})])];
                g3 = [g3 repmat(SideOptNames(sidx),[1,length(AllLicksDuringBase{dsid,midx,sidx,optid})])];
                g4 = [g4 repmat(OptoOpt(optid),[1,length(AllLicksDuringBase{dsid,midx,sidx,optid})])];
            end
        end
    end
end
%form to table
TMP(ismember(g4,'Opto 1st Lick'))=[];
g1(ismember(g4,'Opto 1st Lick'))=[];
g2(ismember(g4,'Opto 1st Lick'))=[];
g3(ismember(g4,'Opto 1st Lick'))=[];
g4(ismember(g4,'Opto 1st Lick'))=[];
tbl = table(TMP',g1',g2',g3',g4','VariableNames',{'LickRate','Genotype','Mouse','Side','Opto'});
%% Fit model without fixed effects:

glmerand = fitglme(tbl,'LickRate ~ 1+(1|Mouse)','Distribution','Poisson','Link','log','FitMethod','Laplace')
anova(glmerand)
%Plot residuals
figure; glmerand.plotResiduals

%eBLUPS
eBLUPS = exp(glmerand.randomEffects);
disp(eBLUPS)

%Graph with grouped data and eBLUPS and observed means
pencol = jet(length(AllMice2Use));
figure; hold on
for j = 1:length(AllMice2Use)
    %get licks from each mouse)
    vals = TMP(ismember(g2,AllMice2Use{j}));
    %Scatter the data with mouse on the x-axis
    scatter(j+linspace(-0.1,0.1,length(vals)),vals,[],pencol(j,:),'filled');
    %add the eBLUPS (black line) and observed means (red line)
    plot([j-0.35 j+0.35],[glmerand.Coefficients.Estimate+eBLUPS(j) glmerand.Coefficients.Estimate+eBLUPS(j)],'Color',[0 0 0],'LineWidth',2)
    plot([j-0.35 j+0.35],[nanmean(vals) nanmean(vals)],'Color',[1 0 0],'LineWidth',2)
end
xlabel('Mouse'),ylabel('Lick Counts')
%Plot the fixed intercept term, which is an estimate of teh grand mean, and
%teh observed grand mean
h1 = plot([0 length(AllMice2Use)],[exp(glmerand.fixedEffects) exp(glmerand.fixedEffects)],'k:');
h2 = plot([0 length(AllMice2Use)],[nanmean(TMP) nanmean(TMP)],'r:');
legend([h1 h2],'Estimated','Observed'),legend('boxoff')
%% Calculate the ICC

%Between mice variance is given by the random intercept variance term
%Within mice is given by the residual variance (accessed as lmemod.MSE).
betwmic = glmerand.covarianceParameters{1};
withinmice = (glmerand.SSE)./(glmerand.NumObservations-2);
ICC = betwmic./(betwmic+withinmice);
disp(ICC)

% Does Opto effect nr. Licks?
clear h
figure('name','LicksPerMouse&Condition')
satval = [0 0.25 0.5];
for sidx=1:2
    subplot(2,1,sidx)
    hold on
    
    mcount = 1;
    for dsid=1:length(AllMicePerGroup)
        for midx=1:length(AllMicePerGroup{dsid})
            clear outvals
            for optid=1:3
                %Get the vals for this condition
                vals = TMP(ismember(g2,AllMicePerGroup{dsid}{midx})&ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid}));
                outvals{optid}=vals; %for connecting the points
                %Scatter
                satcol = (1-pencol).*satval(optid);
                h(optid)=scatter(mcount+ones(numel(vals),1).*(optid-2).*0.25,vals,[],pencol(mcount,:)+satcol(mcount,:),'filled');
                h(4)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmedian(vals) nanmedian(vals)],'Color',[0 0 0],'LineWidth',2);
                h(5)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmean(vals) nanmean(vals)],'Color',[1 0 0],'LineWidth',2);
                
            end
            mcount = mcount+1;
        end
        line([mcount-0.5 mcount-0.5],[0 5],'color',[0 0 0],'LineStyle','--')
    end
    xlabel('Mouse'),ylabel('Lick count')
    set(gca,'XTick',[length(AllMicePerGroup{1})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})+length(AllMicePerGroup{3})./2],'XTickLabel',AllCondNames)
    title(SideOptNames{sidx})
    xlim([0 length(AllMice2Use)+0.5])
end
legend([h(:)],{OptoOpt{:},'Median','Mean'},'Location','best')
%% Include fixed effects

glme = fitglme(tbl,'LickRate ~ 1+ Genotype + Side + Opto + Genotype:Opto:Side + (1|Mouse)','Distribution','Poisson','Link','log','FitMethod','Laplace')
anova(glme)
%Why are there only two levels for Opto and one for Side?
%Let's take a look at the design matrix:
figure,imagesc(glme.designMatrix)
set(gca,'Xtick',[1:length(glme.CoefficientNames)]); set(gca,'XTickLabel',glme.CoefficientNames,'XTickLabelRotation',25); colormap gray
title('Design Matrix')
%We leave out the regressor for Opto Off because otherwise we are
%over-specifying the model. This means that our 'intercept' term becomes
%the response to Opto Off Right Licks

figure('name','Mixed Effects Visual')
coefftable = glme.Coefficients;
hold on
for i = 1:size(coefftable,1)
    line([i i],[coefftable(i,7) coefftable(i,8)],'color',[0 0 0])
    h(i) = plot(i,coefftable(i,2),'ko');
    pval = coefftable(i,6);
    if double(pval)<0.001
        text(i-0.18,double(coefftable(i,8))+0.1,'***','color','k')
    elseif double(pval)<0.01
        text(i-0.12,double(coefftable(i,8))+0.1,'**','color','k')
    else
        text(i-0.06,double(coefftable(i,8))+0.1,'*','color','k')
    end
    
end
set(gca,'xlim',[0.5 size(coefftable,1)+0.5], 'XTick',1:size(coefftable,1),'XTickLabel',glme.CoefficientNames,'XTickLabelRotation',25)
ylabel('Beta-value')
title(['Licks Modelled, R2=' num2str(round(glme.Rsquared.Ordinary*100)) '%'])
line([0.5 size(coefftable,1)+0.5],[0 0],'color', [0 0 0],'LineStyle','--')

difmicecols = distinguishable_colors(maxmouse);
OptoCols = [0 0 0; 1 0 0; 0.7 0.3 0];
figure('name','baseline');
clear h
for sidx=1:2
    subplot(1,2,sidx)
    %Plot
    for dsid=1:length(AllMicePerGroup)
        for optid=1:3
            
            tmp = TMP(ismember(g1,AllCondNames{dsid})&ismember(g4,OptoOpt(optid))&ismember(g3,SideOptNames(sidx)));
            hold on
            boxplot(tmp,'boxstyle','filled','positions',[(dsid-1).*1+(optid-1).*0.2+0.8],'colors',OptoCols(optid,:),'symbol','');
            h(optid) = line([(dsid-1).*1+(optid-1).*0.2+0.7 (dsid-1).*1+(optid-1).*0.2+0.9],[nanmedian(tmp) nanmedian(tmp)],'color',OptoCols(optid,:),'LineWidth',2,'LineStyle','-');
            
        end
        
    end
    
    %         set(gca,'ylim',[-1 3])
    set(gca,'xlim',[0.5 3.5],'XTick',[1:3],'XTickLabel',AllCondNames)
    ylim([-1 5])
    
    %divide sections
    line([1.5 1.5],get(gca,'ylim'),'LineWidth',1,'LineStyle','--','color',[0 0 0])
    line([2.5 2.5],get(gca,'ylim'),'LineWidth',1,'LineStyle','--','color',[0 0 0])
    box off
    ylabel('Licks (hz)')
    title([SideOptNames{sidx} ' licks'])
end
legend([h(:)],OptoOpt)

%If we just want to know the main-effect of Opto we can use the ANOVA method
glme.anova('DFMethod','none')
%You may want use the more conservative Satterthwaite method for estimating the degrees of
%freedom, gives you better control over Type I errors, but has lower power.
glme.anova('DFMethod','residual')

%What is going wrong here?
figure,
r = glme.residuals;
for sidx=1:2
    for optoid = 1:length(OptoOpt)
        subplot(2,3,(sidx-1).*length(OptoOpt)+optoid)
        tmp = cell2mat(cellfun(@(X) find(ismember(AllMice2Use,X)),g2(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})),'UniformOutput',0));
        scatter(tmp,r(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})))
        xlabel('MouseID'),ylabel('Residual')
        title([OptoOpt{optoid} ' ' SideOptNames{sidx}])
    end
end

%The random intercept model can only shift data along the y-axis
%Each mouse gets its own y-offset (or intercept).
satval = [0 0.25 0.5];
p = glme.predict;
clear h
figure('name','PredictedLicksPerMouse&Condition')
satval = [0 0.25 0.5];
for sidx=1:2
    subplot(2,1,sidx)
    hold on
    mcount = 1;
    for dsid=1:length(AllMicePerGroup)
        for midx=1:length(AllMicePerGroup{dsid})
            for optid=1:3
                %Get the predicted output from the model for each pen
                predvals = p(ismember(g2,AllMicePerGroup{dsid}{midx})&ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid}));
                %Scatter the data with pen on teh x-axis
                satcol = (1-pencol).*satval(optid);
                h(optid) = scatter(mcount+ones(numel(predvals),1).*(optid-2).*0.25,predvals,[],pencol(mcount,:)+satcol(mcount,:),'filled');
                
                
                h(4)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmedian(predvals) nanmedian(predvals)],'Color',[0 0 0],'LineWidth',2);
                h(5)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmean(predvals) nanmean(predvals)],'Color',[1 0 0],'LineWidth',2);
                
            end
            mcount = mcount+1;
        end
        line([mcount-0.5 mcount-0.5],[0 2],'color',[0 0 0],'LineStyle','--')
    end
    xlabel('Mouse'),ylabel('Predicted Lick count')
    set(gca,'XTick',[length(AllMicePerGroup{1})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})+length(AllMicePerGroup{3})./2],'XTickLabel',AllCondNames)
    title(SideOptNames{sidx})
    xlim([0 length(AllMice2Use)+0.5])
end
legend([h(:)],{OptoOpt{:},'Median','Mean'},'Location','best')
%% Random slope model

%We will now allow each mouse to have it's own trajectory (i.e. the
%response to the three different opto will be allowed to vary for each
%mouse).
%Note that we should also now specify the covariance pattern for the random
%effect of mouse SEE PPT
lme_cue_slope_diag = fitglme(tbl,'LickRate ~ 1+ Genotype + Side + Opto + Genotype:Opto:Side + (1+Opto|Mouse)','Distribution','Poisson','Link','log','FitMethod','Laplace')
anova(lme_cue_slope_diag)
%Look at our fixed effects as main effects
lme_cue_slope_diag.anova
lme_cue_slope_diag.anova('DFMethod','none')

%What is going wrong here?
figure('name','residuals random slope'),
r = lme_cue_slope_diag.residuals;
for sidx=1:2
    for optoid = 1:length(OptoOpt)
        subplot(2,3,(sidx-1).*length(OptoOpt)+optoid)
        tmp = cell2mat(cellfun(@(X) find(ismember(AllMice2Use,X)),g2(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})),'UniformOutput',0));
        scatter(tmp,r(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})))
        xlabel('MouseID'),ylabel('Residual')
        title([OptoOpt{optoid} ' ' SideOptNames{sidx}])
    end
end

%The random intercept model can only shift data along the y-axis
%Each mouse gets its own y-offset (or intercept).
satval = [0 0.25 0.5];
p = lme_cue_slope_diag.predict;
clear h
figure('name','PredictedLicksPerMouse&Condition random slope')
satval = [0 0.25 0.5];
for sidx=1:2
    subplot(2,1,sidx)
    hold on
    mcount = 1;
    for dsid=1:length(AllMicePerGroup)
        for midx=1:length(AllMicePerGroup{dsid})
            for optid=1:3
                %Get the predicted output from the model for each pen
                predvals = p(ismember(g2,AllMicePerGroup{dsid}{midx})&ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid}));
                %Scatter the data with pen on teh x-axis
                satcol = (1-pencol).*satval(optid);
                h(optid) = scatter(mcount+ones(numel(predvals),1).*(optid-2).*0.25,predvals,[],pencol(mcount,:)+satcol(mcount,:),'filled');
                
                
                h(4)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmedian(predvals) nanmedian(predvals)],'Color',[0 0 0],'LineWidth',2);
                h(5)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmean(predvals) nanmean(predvals)],'Color',[1 0 0],'LineWidth',2);
                
            end
            mcount = mcount+1;
        end
        line([mcount-0.5 mcount-0.5],[0 2],'color',[0 0 0],'LineStyle','--')
    end
    xlabel('Mouse'),ylabel('Predicted Lick count')
    set(gca,'XTick',[length(AllMicePerGroup{1})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})+length(AllMicePerGroup{3})./2],'XTickLabel',AllCondNames)
    title(SideOptNames{sidx})
    xlim([0 length(AllMice2Use)+0.5])
end
legend([h(:)],{OptoOpt{:},'Median','Mean'},'Location','best')

figure('name','Mixed Effects Baseline')
coefftable = lme_cue_slope_diag.Coefficients;
hold on
for i = 1:size(coefftable,1)
    line([i i],[coefftable(i,7) coefftable(i,8)],'color',[0 0 0])
    h(i) = plot(i,coefftable(i,2),'ko');
    pval = coefftable(i,6);
    if double(pval)<0.001
        text(i-0.18,double(coefftable(i,8))+0.1,'***','color','k')
    elseif double(pval)<0.01
        text(i-0.12,double(coefftable(i,8))+0.1,'**','color','k')
    else
        text(i-0.06,double(coefftable(i,8))+0.1,'*','color','k')
    end
    
end
set(gca,'xlim',[0.5 size(coefftable,1)+0.5], 'XTick',1:size(coefftable,1),'XTickLabel',glme.CoefficientNames,'XTickLabelRotation',25)
ylabel('Beta-value')
title(['Licks Modelled, R2=' num2str(round(glme.Rsquared.Ordinary*100)) '%'])
line([0.5 size(coefftable,1)+0.5],[0 0],'color', [0 0 0],'LineStyle','--')
%% Poisson GLM Visual

TMP = []; %length(AllMicePerGroup),length(midx),Sides,OptoCond;
g1=[]; %Genotype
g2=[]; %Mouse
g3=[]; %Left/Right/Total Licks
g4=[]; %No Opto, Opto EARLY, Opto Lick
for dsid=1:length(AllMicePerGroup)
    for midx=1:length(AllMicePerGroup{dsid})
        for sidx=1:2
            for optid=1:3
                TMP = [TMP AllLicksDuringVis{dsid,midx,sidx,optid}];
                g1 = [g1  repmat(AllCondNames(dsid),[1,length(AllLicksDuringVis{dsid,midx,sidx,optid})])];
                g2 = [g2 repmat(AllMicePerGroup{dsid}(midx),[1,length(AllLicksDuringVis{dsid,midx,sidx,optid})])];
                g3 = [g3 repmat(SideOptNames(sidx),[1,length(AllLicksDuringVis{dsid,midx,sidx,optid})])];
                g4 = [g4 repmat(OptoOpt(optid),[1,length(AllLicksDuringVis{dsid,midx,sidx,optid})])];
            end
        end
    end
end
TMP(ismember(g4,'Opto 1st Lick'))=[];
g1(ismember(g4,'Opto 1st Lick'))=[];
g2(ismember(g4,'Opto 1st Lick'))=[];
g3(ismember(g4,'Opto 1st Lick'))=[];
g4(ismember(g4,'Opto 1st Lick'))=[]; % Don't use opto1stlick
%form to table
tbl = table(TMP',g1',g2',g3',g4','VariableNames',{'LickRate','Genotype','Mouse','Side','Opto'});
%% Fit model without fixed effects:

glmerand = fitglme(tbl,'LickRate ~ 1+(1|Mouse)','Distribution','Poisson','Link','log','FitMethod','Laplace')
anova(glmerand)
%Plot residuals
figure; glmerand.plotResiduals

%eBLUPS
eBLUPS = exp(glmerand.randomEffects);
disp(eBLUPS)

%Graph with grouped data and eBLUPS and observed means
pencol = jet(length(AllMice2Use));
figure; hold on
for j = 1:length(AllMice2Use)
    %get licks from each mouse)
    vals = TMP(ismember(g2,AllMice2Use{j}));
    %Scatter the data with mouse on the x-axis
    scatter(j+linspace(-0.1,0.1,length(vals)),vals,[],pencol(j,:),'filled');
    %add the eBLUPS (black line) and observed means (red line)
    plot([j-0.35 j+0.35],[glmerand.Coefficients.Estimate+eBLUPS(j) glmerand.Coefficients.Estimate+eBLUPS(j)],'Color',[0 0 0],'LineWidth',2)
    plot([j-0.35 j+0.35],[nanmean(vals) nanmean(vals)],'Color',[1 0 0],'LineWidth',2)
end
xlabel('Mouse'),ylabel('Lick Counts')
%Plot the fixed intercept term, which is an estimate of teh grand mean, and
%teh observed grand mean
h1 = plot([0 length(AllMice2Use)],[exp(glmerand.fixedEffects) exp(glmerand.fixedEffects)],'k:');
h2 = plot([0 length(AllMice2Use)],[nanmean(TMP) nanmean(TMP)],'r:');
legend([h1 h2],'Estimated','Observed'),legend('boxoff')
%% Calculate the ICC

%Between mice variance is given by the random intercept variance term
%Within mice is given by the residual variance (accessed as lmemod.MSE).
betwmic = glmerand.covarianceParameters{1};
withinmice = (glmerand.SSE)./(glmerand.NumObservations-2);
ICC = betwmic./(betwmic+withinmice);
disp(ICC)

% Does Opto effect nr. Licks?
clear h
figure('name','LicksPerMouse&Condition')
satval = [0 0.25 0.5];
for sidx=1:2
    subplot(2,1,sidx)
    hold on
    
    mcount = 1;
    for dsid=1:length(AllMicePerGroup)
        for midx=1:length(AllMicePerGroup{dsid})
            clear outvals
            for optid=1:3
                %Get the vals for this condition
                vals = TMP(ismember(g2,AllMicePerGroup{dsid}{midx})&ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid}));
                outvals{optid}=vals; %for connecting the points
                %Scatter
                satcol = (1-pencol).*satval(optid);
                h(optid)=scatter(mcount+ones(numel(vals),1).*(optid-2).*0.25,vals,[],pencol(mcount,:)+satcol(mcount,:),'filled');
                h(4)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmedian(vals) nanmedian(vals)],'Color',[0 0 0],'LineWidth',2);
                h(5)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmean(vals) nanmean(vals)],'Color',[1 0 0],'LineWidth',2);
                
            end
            mcount = mcount+1;
        end
        line([mcount-0.5 mcount-0.5],[0 5],'color',[0 0 0],'LineStyle','--')
    end
    xlabel('Mouse'),ylabel('Lick count')
    set(gca,'XTick',[length(AllMicePerGroup{1})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})+length(AllMicePerGroup{3})./2],'XTickLabel',AllCondNames)
    title(SideOptNames{sidx})
    xlim([0 length(AllMice2Use)+0.5])
end
legend([h(:)],{OptoOpt{:},'Median','Mean'},'Location','best')
%% Include fixed effects

glme = fitglme(tbl,'LickRate ~ 1+ Genotype + Side + Opto + Genotype:Opto:Side + (1|Mouse)','Distribution','Poisson','Link','log','FitMethod','Laplace')
anova(glme)
%Why are there only two levels for Opto and one for Side?
%Let's take a look at the design matrix:
figure,imagesc(glme.designMatrix)
set(gca,'Xtick',[1:length(glme.CoefficientNames)]); set(gca,'XTickLabel',glme.CoefficientNames,'XTickLabelRotation',25); colormap gray
title('Design Matrix')
%We leave out the regressor for Opto Off because otherwise we are
%over-specifying the model. This means that our 'intercept' term becomes
%the response to Opto Off Right Licks

figure('name','Mixed Effects Visual')
coefftable = glme.Coefficients;
hold on
for i = 1:size(coefftable,1)
    line([i i],[coefftable(i,7) coefftable(i,8)],'color',[0 0 0])
    h(i) = plot(i,coefftable(i,2),'ko');
    pval = coefftable(i,6);
    if double(pval)<0.001
        text(i-0.18,double(coefftable(i,8))+0.1,'***','color','k')
    elseif double(pval)<0.01
        text(i-0.12,double(coefftable(i,8))+0.1,'**','color','k')
    else
        text(i-0.06,double(coefftable(i,8))+0.1,'*','color','k')
    end
    
end
set(gca,'xlim',[0.5 size(coefftable,1)+0.5], 'XTick',1:size(coefftable,1),'XTickLabel',glme.CoefficientNames,'XTickLabelRotation',25)
ylabel('Beta-value')
title(['Licks Modelled, R2=' num2str(round(glme.Rsquared.Ordinary*100)) '%'])
line([0.5 size(coefftable,1)+0.5],[0 0],'color', [0 0 0],'LineStyle','--')

difmicecols = distinguishable_colors(maxmouse);
OptoCols = [0 0 0; 1 0 0; 0.7 0.3 0];
figure('name','Visual');
clear h
for sidx=1:2
    subplot(1,2,sidx)
    %Plot
    for dsid=1:length(AllMicePerGroup)
        for optid=1:3
            
            tmp = TMP(ismember(g1,AllCondNames{dsid})&ismember(g4,OptoOpt(optid))&ismember(g3,SideOptNames(sidx)));
            hold on
            boxplot(tmp,'boxstyle','filled','positions',[(dsid-1).*1+(optid-1).*0.2+0.8],'colors',OptoCols(optid,:),'symbol','');
            %         for mid=1:length(AllMicePerGroup{dsid})
            %             tmp = TMP(ismember(g1,AllCondNames{dsid})&ismember(g2,AllMicePerGroup{dsid}(mid))&ismember(g4,OptoOpt(optid)));
            %             % PLOT
            %             hold on
            %             scatter(repmat((dsid-1).*1+(optid-1).*0.2+0.8+(0.01.*mid)-0.005.*length(AllMicePerGroup{dsid}),1,length(tmp)),tmp,10,difmicecols(mid,:))
            %         end
            h(optid) = line([(dsid-1).*1+(optid-1).*0.2+0.7 (dsid-1).*1+(optid-1).*0.2+0.9],[nanmedian(tmp) nanmedian(tmp)],'color',OptoCols(optid,:),'LineWidth',2,'LineStyle','-');
            
        end
        
    end
    
    %         set(gca,'ylim',[-1 3])
    set(gca,'xlim',[0.5 3.5],'XTick',[1:3],'XTickLabel',AllCondNames)
    ylim([-1 5])
    
    %divide sections
    line([1.5 1.5],get(gca,'ylim'),'LineWidth',1,'LineStyle','--','color',[0 0 0])
    line([2.5 2.5],get(gca,'ylim'),'LineWidth',1,'LineStyle','--','color',[0 0 0])
    box off
    ylabel('Licks (hz)')
    title([SideOptNames{sidx} ' licks'])
end
legend([h(:)],OptoOpt)

%If we just want to know the main-effect of Opto we can use the ANOVA method
glme.anova('DFMethod','none')
%You may want use the more conservative Satterthwaite method for estimating the degrees of
%freedom, gives you better control over Type I errors, but has lower power.
glme.anova('DFMethod','residual')

%What is going wrong here?
figure,
r = glme.residuals;
for sidx=1:2
    for optoid = 1:length(OptoOpt)
        subplot(2,3,(sidx-1).*length(OptoOpt)+optoid)
        tmp = cell2mat(cellfun(@(X) find(ismember(AllMice2Use,X)),g2(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})),'UniformOutput',0));
        scatter(tmp,r(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})))
        xlabel('MouseID'),ylabel('Residual')
        title([OptoOpt{optoid} ' ' SideOptNames{sidx}])
    end
end

%The random intercept model can only shift data along the y-axis
%Each mouse gets its own y-offset (or intercept).
satval = [0 0.25 0.5];
p = glme.predict;
clear h
figure('name','PredictedLicksPerMouse&Condition')
satval = [0 0.25 0.5];
for sidx=1:2
    subplot(2,1,sidx)
    hold on
    mcount = 1;
    for dsid=1:length(AllMicePerGroup)
        for midx=1:length(AllMicePerGroup{dsid})
            for optid=1:3
                %Get the predicted output from the model for each pen
                predvals = p(ismember(g2,AllMicePerGroup{dsid}{midx})&ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid}));
                %Scatter the data with pen on teh x-axis
                satcol = (1-pencol).*satval(optid);
                h(optid) = scatter(mcount+ones(numel(predvals),1).*(optid-2).*0.25,predvals,[],pencol(mcount,:)+satcol(mcount,:),'filled');
                
                
                h(4)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmedian(predvals) nanmedian(predvals)],'Color',[0 0 0],'LineWidth',2);
                h(5)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmean(predvals) nanmean(predvals)],'Color',[1 0 0],'LineWidth',2);
                
            end
            mcount = mcount+1;
        end
        line([mcount-0.5 mcount-0.5],[0 2],'color',[0 0 0],'LineStyle','--')
    end
    xlabel('Mouse'),ylabel('Predicted Lick count')
    set(gca,'XTick',[length(AllMicePerGroup{1})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})+length(AllMicePerGroup{3})./2],'XTickLabel',AllCondNames)
    title(SideOptNames{sidx})
    xlim([0 length(AllMice2Use)+0.5])
end
legend([h(:)],{OptoOpt{:},'Median','Mean'},'Location','best')
%% Random slope model

%We will now allow each mouse to have it's own trajectory (i.e. the
%response to the three different opto will be allowed to vary for each
%mouse).
%Note that we should also now specify the covariance pattern for the random
%effect of mouse SEE PPT
lme_cue_slope_diag = fitglme(tbl,'LickRate ~ 1+ Genotype + Side + Opto + Genotype:Opto:Side + (1+Opto|Mouse)','Distribution','Poisson','Link','log','FitMethod','Laplace')
anova(lme_cue_slope_diag)

%Look at our fixed effects as main effects
lme_cue_slope_diag.anova
lme_cue_slope_diag.anova('DFMethod','none')

%What is going wrong here?
figure('name','residuals random slope'),
r = lme_cue_slope_diag.residuals;
for sidx=1:2
    for optoid = 1:length(OptoOpt)
        subplot(2,3,(sidx-1).*length(OptoOpt)+optoid)
        tmp = cell2mat(cellfun(@(X) find(ismember(AllMice2Use,X)),g2(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})),'UniformOutput',0));
        scatter(tmp,r(ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid})))
        xlabel('MouseID'),ylabel('Residual')
        title([OptoOpt{optoid} ' ' SideOptNames{sidx}])
    end
end

%The random intercept model can only shift data along the y-axis
%Each mouse gets its own y-offset (or intercept).
satval = [0 0.25 0.5];
p = lme_cue_slope_diag.predict;
clear h
figure('name','PredictedLicksPerMouse&Condition random slope')
satval = [0 0.25 0.5];
for sidx=1:2
    subplot(2,1,sidx)
    hold on
    mcount = 1;
    for dsid=1:length(AllMicePerGroup)
        for midx=1:length(AllMicePerGroup{dsid})
            for optid=1:3
                %Get the predicted output from the model for each pen
                predvals = p(ismember(g2,AllMicePerGroup{dsid}{midx})&ismember(g3,SideOptNames{sidx})&ismember(g4,OptoOpt{optid}));
                %Scatter the data with pen on teh x-axis
                satcol = (1-pencol).*satval(optid);
                h(optid) = scatter(mcount+ones(numel(predvals),1).*(optid-2).*0.25,predvals,[],pencol(mcount,:)+satcol(mcount,:),'filled');
                
                
                h(4)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmedian(predvals) nanmedian(predvals)],'Color',[0 0 0],'LineWidth',2);
                h(5)=plot([mcount+1.*(optid-2).*0.25-0.1 mcount+1.*(optid-2).*0.25+0.1],[nanmean(predvals) nanmean(predvals)],'Color',[1 0 0],'LineWidth',2);
                
            end
            mcount = mcount+1;
        end
        line([mcount-0.5 mcount-0.5],[0 2],'color',[0 0 0],'LineStyle','--')
    end
    xlabel('Mouse'),ylabel('Predicted Lick count')
    set(gca,'XTick',[length(AllMicePerGroup{1})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})./2 length(AllMicePerGroup{1})+length(AllMicePerGroup{2})+length(AllMicePerGroup{3})./2],'XTickLabel',AllCondNames)
    title(SideOptNames{sidx})
    xlim([0 length(AllMice2Use)+0.5])
end
legend([h(:)],{OptoOpt{:},'Median','Mean'},'Location','best')

figure('name','Mixed Effects Visual')
coefftable = lme_cue_slope_diag.Coefficients;
hold on
for i = 1:size(coefftable,1)
    line([i i],[coefftable(i,7) coefftable(i,8)],'color',[0 0 0])
    h(i) = plot(i,coefftable(i,2),'ko');
    pval = coefftable(i,6);
    if double(pval)<0.001
        text(i-0.18,double(coefftable(i,8))+0.1,'***','color','k')
    elseif double(pval)<0.01
        text(i-0.12,double(coefftable(i,8))+0.1,'**','color','k')
    else
        text(i-0.06,double(coefftable(i,8))+0.1,'*','color','k')
    end
    
end
set(gca,'xlim',[0.5 size(coefftable,1)+0.5], 'XTick',1:size(coefftable,1),'XTickLabel',glme.CoefficientNames,'XTickLabelRotation',25)
ylabel('Beta-value')
title(['Licks Modelled, R2=' num2str(round(glme.Rsquared.Ordinary*100)) '%'])
line([0.5 size(coefftable,1)+0.5],[0 0],'color', [0 0 0],'LineStyle','--')