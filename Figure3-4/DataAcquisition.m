%Training Script Enny van Beest for Mice: EasyLickTask. Start with
%water-restriction. Habituation
%to head-restriction before start training.
% Moving Stim is rewarded (80%) and still stim is not rewarded (20%)
% (eventually, but first  trai n on 100% moving)
%% Manual User Settings
ProcMoving = 100; %Moving stimulus = rewarded (in percentage)
arduino = 1; %Arduino connected?
MetaData = 0; %Metafile (Jason) saving?
debug = 0;%Debug mode (i.e. Destkop computer is cogent screen)
savefig = 0; %Save figure at the end?
accuracyThrs = 65; %Accuracy threshold (in performace figure)
nback = 15; %check how many last perfor mances for live scores? (per side x ori combination!)
MovingGratings = 1; %if 1; stimuli move. If 0; stimuli don't move
lraltern = 2; %how many of the same stimuli behind eachother? (only for phase 1 & 2)
punishEarlyLicking =0; % if 1, errorbeep when licks are early. If 0, no punishment
rng('shuffle'); % avoid pseudorandomization
propR = 1; %proportion of trials that is left (rest is right)
cleanbaseline = 1; % Clean baseline before trial starts; helps getting rid of random licks
GraceP = 0;%inms; time after stim onset that licks do not count for reward
if cleanbaseline ==1
    cleanbasetime = 2; %Time there should not have been a lick
else
    cleanbasetime =0;
end
OnlineFeedback = 0;
%% Paths + Current Directory
user_name =getenv('USERNAME');
tmpname = cd;
if ~isempty(strfind(tmpname,[user_name '.HERSEN']))
    user_name = regexp(tmpname,'\','split');
    user_name = user_name{3};
end
LocalFolder = mfilename('fullpath');
WMpath = strsplit(LocalFolder,'\');

DomainFolder = fullfile(WMpath{1:end-3});
LocalFolder = WMpath{1};
WMpath = fullfile(WMpath{1:end-1});
vcTrainingLogs = '\Path\To\LogFiles\';

try
    addpath(genpath(WMpath))
    cd(WMpath)
catch ME
    if debug
        WMpath = 'I:\GitHub\Mouse\WM';
        addpath(genpath(WMpath))
        cd(WMpath)
        addpath(genpath('\Path\To\MATLAB\CogGph'))
    else
        disp(ME)
        WMpath = ['C:\Users\' user_name '\Documents\Github\Mouse\WM\'];
        addpath(genpath(WMpath))
        cd(WMpath)
    end
end



%% Upload arduino
if arduino
    disp('Sure you uploaded correct arduino script? Press key to continue or abort...')
    pause
end

%% Training setup
setup = {'1', '2', '3', '4', '5', '6', 'Testing', 'Recording','Imaging','Opto','Imaging+Opto'};
setup = menu('Choose the setup', setup);

%% Start Serial connection with Arduino
if arduino
    if ~exist('sport', 'var')
        sport = serial('com3');
    else
        clear LOG
    end
    set(sport,'InputBufferSize', 10240)
    try
        if strcmp(get(sport, 'status'), 'closed'); fopen(sport); end;
    catch
        if setup==10
            sport = serial('com20'); %changed this to port 20 was 21 - Moh
        elseif setup==11
            sport = serial('com9');
        end        
        if strcmp(get(sport, 'status'), 'closed'); fopen(sport); end;
    end
    set(sport, 'baudrate', 250000);
    set(sport, 'timeout', 0.5);
end


if setup == 9 || setup==11
    addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github'],'OpenAccessStorage')))
    MetaData=1;
    addpath('Z:\OpenAccessStorageParameters')

end
%     if setup==10
%     addpath(genpath(fullfile(DomainFolder,'OpenAccessStorage')))
%     addpath(genpath('E:\CogGph\CogGphTB'))
%     addpath(genpath('E:\das64bit'))
% elseif debug
%     addpath(genpath('I:\GitHub\OpenAccessStorage'))
% else
%     addpath(genpath(fullfile(['C:\Users\' user_name '\Dropbox\'],'OpenAccessStorage')))
% end



%% Box/Imaging setup
if setup == 9 || setup==11
    prerecordingtime = 0.5;
    TotalRecordingtime = 5.5; %??
    eval(['run ' fullfile(['C:\Users\' user_name '\Documents\Github'],'Mouse','imagingStim','runSettings.m')])
    D2ScreenWidth = 51./2; %??
    ITI = 9;
    LickPos = 70;    %motor
    BreakPos = 64;   %motor
    displayscript = 'correctDisplayWF'
    gammaconversion = 'gammaconWF'
elseif setup == 3
    screenx = 1280;
    screeny = 720;
    refresh = 60;
    ScreenDistance = 13;
    D2ScreenWidth = 51./2;
    ITI = 4;
    LickPos = 60;    %motor
    BreakPos = 70;   %motor
    prefsf = 0.075;
    prefspeed = 24; %deg/s
    displayscript = 'correctDisplay2'
    gammaconversion = 'gammacon'
elseif setup==10
    screenx = 1920;
    screeny = 1200;
    refresh = 60;
    gammaconversion = 'gammaconOpto'
    displayscript = 'correctDisplay2Opto'
    ITI = 4;
    LickPos = 60;    %motor
    BreakPos = 70;   %motor
    prefsf = 0.075;
    prefspeed = 24; %deg/s
    ScreenDistance = 11;
    D2ScreenWidth = 51/2;     % width of screen, in cm
    TotalRecordingtime = 3.5; %??

else
    screenx = 1280;
    screeny = 720;
    refresh = 60;
    ScreenDistance = 13;
    D2ScreenWidth = 51./2;
    ITI = 4;
    LickPos = 60;    %motor
    BreakPos = 70;   %motor
    prefsf = 0.075;
    prefspeed = 24; %deg/s
    displayscript = 'correctDisplay2'
    gammaconversion = 'gammacon'
end

%% Session name (nts)
warning on backtrace
logon = 0;
prompt = {'Mouse Name', 'Exp Nr', 'Date','Laser Setting'};
dlg_title = 'Please enter parameters';
num_lines = 1;
def = {'Name', '1', datestr(datenum(date), 'yyyymmdd'), '220', '1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~strcmp(answer{1}, 'Name')
    
    mouse = answer{1};
    expnum = answer{2};
    date = answer{3};
    expname = [mouse expnum];
    CurLaserSetting = answer{4};
    LOG.Stimulus = 'SimpleLick';
    logon = 1;
    
    %     % add dasbit stuff
    %     addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github\Mouse\DIO24'])))
else
    mouse = 'Unknown';
    leftside = round(rand(1,1))+1; % if 1, left is 0 degrees, if 2, left is 90 degrees
    CurPhase = 1;
    ezvalue = 1;
    spatialoffset = screenx/4;
    EarlyLickP = round((Par.maxDelay + Par.StimDurMin)); %EarlyLickP is the period in which the mouse is allowed to lick without it being punished
    
end
nts = [mouse '_' date '_B' expnum];
expname = [mouse num2str(expnum)];

while exist(fullfile(vcTrainingLogs,mouse,[mouse date],['B' expnum],[nts '.mat']), 'file')
    disp(['Log-file expnum ' num2str(expnum) ' already used, expnum number changes...']);
    expnum = expnum+1;
    nts = [mouse '_' date '_B' num2str(expnum)]
    expname = [mouse num2str(expnum)];    
end
expname = [mouse expnum];

%% JSON FILE
% the metadata in this schema is stored in the database,
% but much more can be added to the json files
if MetaData
    if setup == 9 || setup == 10 || setup == 11
        fields.project = 'WorkingMemory_EB';
        fields.dataset = 'SimpleLickTask';
        fields.investigator = 'EnnyvB';
        fields.date = date;
        fields.subject = mouse;
        if setup==9
            fields.setup = 'WF';
        elseif setup==10
            fields.setup = 'Opto';
        elseif setup==11
            fields.setup = 'OptoAndWF';
        else
            fields.setup = 'Behavior';
        end
        fields.stimulus = 'SpatialWM';
        fields.condition = 'awake';
        
        fields = getdbfields('VC',fields); %retrieves info from mysql tables using a GUI, fill in already.
        
        %not necessary if you are extremely consistent in filling the schema below
        %schema for structuring fields in json file
        %the php script that parses the json files requires that the data in the json file
        %is structured accordingly, to fill the database with records
        %
        json = fields;
        json.logfile = [expname '.mat']; %name for your logfile
        json.version = '1.0'
        
    else
        MetaData = 0;
    end
end

%% initialize PAR-variable: contains general Settings.
global Par
Par.ItiRandMin = 0;
Par.ItiRandMax = 2;
Par.StimDurMin = 500; %min. ms stim on
Par.Orientations = [45 135]; %With Old mice : [0 90]),new mice [ 45 135];
Par.BeepFreq = 2000; %Hz
Par.BeepTime = 0.5;  %Enter time constant beep
Par.StartleFreq = 5000; %Hz :Switched to WHITENOISE! %Tone for preliminary licking
Par.StartleTime = 1.5;

Par.SpatFreq = 0.08;%0.05; 12.5 degree sper cycle m ?
Par.RewardTime = 3000; %max reward time (Opto is on till halfway this time!)
Par.FigSize = 35; % in visual degrees

%% Initialize default values for other variables
%controlbox (overwritten later on by 'load CurrentPhase.m' if existing)
global RTrightVec
global RTleftVec
ezvalue = 1;
shuffleval = 1;
CurPassivPerc = 0;
freeze = 0;
drumm = 0;
semidrumm = 0;
passL = 1; %passives ON (1)/OFF (0)
passR = 1; %passives ON (1)/OFF (0)
spatialoffset = 40; %in visual degrees

%Trial-specific
passfirst = 0;
ESC = 0;
Miss = 0;
trialcount = 0;
gavepassive = 0;
wentthrough = 0;
details = [];
prefspeed = 24; %deg/s
changepassives = 0;
GaveExtraAlready = 0;
RewardCount = 0; %Counts rewards given for callibration purposes
Reaction = [];

%visual stimuli
centercontrast = 100;
whichway = 'rgb2lum';
maxlum = eval([gammaconversion '(0.4,''rgb2lum'')']);
minlum = eval([gammaconversion '(0,''rgb2lum'')']);
greylum = (maxlum+minlum)./2;
Par.greylum = greylum;
grey = eval([gammaconversion '(greylum,''lum2rgb'')']);
flashstim = grey;
Par.grey = grey;
fullscrfig = 1;

%% Initialize other pars Computer/cabinet dependent, so after loading of PAR for mouse
%Screen
Par.Screenx = screenx;
Par.Screeny = screeny;
Par.Refresh = refresh;
Par.DistanceToScreen = ScreenDistance;%
Par.Centerx = screenx/2;
Par.ScreenWidthD2 = D2ScreenWidth;
Par.ITI = ITI;
Par.HW = Par.Screenx/2;
Par.HH = Par.Screeny/2;
Par.PixPerDeg = Par.HW/atand(Par.ScreenWidthD2/Par.DistanceToScreen);
Par.DegPerPix = 1/Par.PixPerDeg;
Par.RadPerPix = Par.DegPerPix*pi/180;
Ang2Pix = Par.PixPerDeg;
display(['Pixel per Degree: ' num2str(Ang2Pix)]);


%% Initialize  Cogent

cgloadlib

cgshut %shut cg if one's open
if debug
    warning('Debug Mode enabled')
    cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 0)
else
    if setup == 9 || setup == 11
        try
            cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 1) %Second window?
        catch
            cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 0)
        end
    elseif setup==10
        cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 1)
    else
        cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 2)
    end
end

cgpenwid(1)
cgpencol(0,0,0)
cgfont('Arial',100)
cogstd('sPriority','high')
cgflip(.5, .5, .5)
cgflip(.5, .5, .5)


%% Check whether screen details match
gsd = cggetdata('gsd');
if gsd.ScreenWidth ~= Par.Screenx || gsd.ScreenHeight ~= Par.Screeny
    error('screen details do not match')
end

AllPerfVec = nan(1,500);
AccVec = nan(1,500);

%% Ccreate visual stimuli
SafeStimPath = fullfile(LocalFolder,'WM_OriGratings');
if ~exist(SafeStimPath,'dir')
    mkdir(SafeStimPath)
end

% pixpercycle = Par.PixPerDeg * 1/Par.SpatFreq;
secpercycle = 1/Par.SpatFreq/prefspeed;
framespercycle = round(Par.Refresh * secpercycle);
phasesteps = 0:2*pi/framespercycle:2*pi-(2*pi/framespercycle);

for j = 1:length(Par.Orientations)
    for i = 1:length(phasesteps)        
        CURR.Orient = Par.Orientations(j);
        CURR.Phase = phasesteps(i);
        CURR.Contrast = 100; %Full contrast for center
        CURR.pedestal = greylum;
        if (setup==9) || (setup==11)
            a1 = correctDisplayWF(makegrating2GammaconV(160,CURR.Orient,Par.PixPerDeg * 1/prefsf, CURR.Phase, minlum, maxlum,gammaconversion));
        elseif setup==10
            a1 = correctDisplay2Opto(makegrating2GammaconV(160,CURR.Orient,Par.PixPerDeg * 1/prefsf, CURR.Phase, minlum, maxlum,gammaconversion));
        else
            a1 = correctDisplay2(makegrating2GammaconV(160,CURR.Orient,Par.PixPerDeg * 1/prefsf, CURR.Phase, minlum, maxlum,gammaconversion));
        end
        imwrite(a1,[SafeStimPath '\' num2str(Par.Orientations(j)) '_' num2str(phasesteps(i)) '.bmp'])
    end
end

%% Dasbit
if setup == 9 || setup == 11
    dasinit(22)
    dasbit(0, 0) %initializes camera (set to 0)
    dasbit(3, 0) %stimbit
    dasbit(1, 1) %Shutter open
    disp('Shutter Opened')
    dasbit(2, 0)  %starts recording
    if setup == 11
        optoport = 4;
        dasbit(optoport,0) %opto
        OPTOREG = 0;
    end
end

if setup == 10
%     dasinit(23);
%     optoport = 0;
%     dasbit(optoport,0);
    OPTOREG = 0;
end

optotrialvec = [0]; %0 = off, 1 = before trial onset, 2 = with first lick (after stimonset)
optoOn = 0;

movementtrialvec = [repmat(1,1,ProcMoving/10) repmat(0,1,(100-ProcMoving)/10)];

figure(1) %show control-window
switch propR
    case 1
        oriopt = [Par.Orientations(2),Par.Orientations(2),Par.Orientations(1),Par.Orientations(1)];
    case 2
        oriopt = [Par.Orientations(2) Par.Orientations(1) Par.Orientations(1) Par.Orientations(1)];
    case 3
        oriopt = [Par.Orientations(2) Par.Orientations(2) Par.Orientations(2) Par.Orientations(1)];
    otherwise
        error('unknown PropR value')
end


%% Open Little window
OpenContrWindowSimpleLickTask

%% Pause until start
disp('Pausing. Press any key to start session')
pause

%% Final settings
ititime = tic;
start = tic;
left = 'left';
right = 'right';

%% Make orientation matrix: 'details'
clear details
z = 1;
for p = optotrialvec
    for m = movementtrialvec %moving stim is rewarded
        for o =oriopt  %orientations
            details(z,:) = [o,p,m];
            z = z+1;
        end
    end
end
%% Shuffling: create (random) 'order' variable
[ntrials, nvars] = size(details);

%Shuffle trials?
if get(shflbox,'value') == 0
    trials = details;
    shuffleval = 0;
elseif get(shflbox,'value') == 1
    % Order such that it's trying to have a maximum of 3 trials of the same
    % side after each other
    trials = [];
    flag = 0;
    tmpdetails = details;
    while (size(trials,1))<ntrials && ~flag
        order = randsample(oriopt,4,0);
        for i = 1:length(order)
            vec = find(tmpdetails(:,1) == order(i));
            if length(vec)>0
                pick = randsample(vec,1);
            else
                %More from 1 condition --> Probably because of propl.
                %Put these extra oin random locations
                flag=1;
                break
            end
            trials = [trials; tmpdetails(pick,:)];
            tmpdetails(pick,:)=[];
        end
    end
    %More from 1 condition --> Probably because of propl.
    %Put these extra oin random locations
    while ~isempty(tmpdetails)
        for sid = unique(tmpdetails(:,1))'
            RandIdxVec = find(trials(:,1)==sid);
            RandIdxVec = RandIdxVec(diff(RandIdxVec') == lraltern);
            if isempty(RandIdxVec)
                RandIdx = size(trials,1)+1;
                trials = cat(1,trials,tmpdetails(find(tmpdetails(:,1)==sid,1),:));
                tmpdetails(find(tmpdetails(:,1)==sid,1),:)=[];
                break
            end
            count = 0;
            for RandIdx = RandIdxVec'
                if isempty(find(tmpdetails(:,1)==sid,1))
                    break
                else
                    trials = cat(1,trials(1:RandIdx+count,:),tmpdetails(find(tmpdetails(:,1)==sid,1),:),trials(RandIdx+count+1:end,:));
                    tmpdetails(find(tmpdetails(:,1)==sid,1),:) = [];
                end
                count = count+1;
            end
        end
    end
    
    shuffleval = 1;
end


maxcount =countrepeats(trials(:,1));
if maxcount>4
    warning('More than 4 trials of the same side after each other in trial struct')
end

%% Set passive percentage and define passive trial vector
% set(PercPassBox,'String',num2str(CurPassivPerc))
GivePercTrial = false(1,size(trials,1));
if ~shuffleval
    ii = 1;
    curPerc = 0;
    totalPerc = CurPassivPerc/100;
    GivePercTrial(ii:lraltern:length(GivePercTrial)) = 1;
    while sum(GivePercTrial)./length(GivePercTrial) < totalPerc %increase number of passives
        ii = ii+1;
        GivePercTrial(ii:lraltern:length(GivePercTrial)) = 1;
    end
    if sum(GivePercTrial)./length(GivePercTrial) > totalPerc %decrease number of passives
        onevec = find(GivePercTrial == 1);
        onevec = onevec(randperm(length(onevec)));
        GivePercTrial(onevec(1:round((sum(GivePercTrial)./length(GivePercTrial) - totalPerc)*length(GivePercTrial))))=0;
    end
else %shuffle
    GivePercTrial(1:(round(size(trials,1)*CurPassivPerc/100))) = 1;
    GivePercTrial = GivePercTrial(randperm(size(trials,1)));
end

while ~ESC
    %% Read in GUI
    CurPassivPerc = str2double(get(PercPassBox,'string'));    
    passrew = str2double(get(rewbox,'string'));
    passL = get(plbox,'value');
    passR = get(prbox,'value');
    pdelay = str2double(get(pdbox,'string'));
    cleanbasetime = str2double(get(cleanbaselinebox,'string'));
    shuffleval = get(shflbox,'value');
    ezvalue = get(ezbox,'value');
    drumm = get(drummbox,'value');
    semidrumm = get(semidrummbox,'value');
    
    threshold = str2double(get(thbox, 'string'));
    rdurL = str2double(get(rdurleft, 'string'));
    rdurR = str2double(get(rdurright, 'string'));
    spatialoffset = str2double(get(xdispbox,'String'));
    
    if OPTOREG == 1
        CurLaserSetting = str2double(get(optoValBox,'string'));
    end    
    OPTOREG = get(optobox,'value');
    
    MovingGratings = get(movingBox,'value');
    fullscrfig = get(fcheckbox,'Value');
    propR = get(PropBox,'value');
    
    optoOnTime = nan;
    optoOFFTime = nan;
    
    %% Update Arduino settings
    if arduino
        sendardwm(sport, ['IL ' num2str(rdurR)]) %This seems strange, but fits with the arduino script which has been build the other way around
        sendardwm(sport, ['IR ' num2str(rdurL)])
        sendardwm(sport, ['IF ' get(thbox, 'string')])
        sendardwm(sport, ['IM ' num2str(get(ezbox, 'Value'))])
    end
    
    %% Set current trial parameters
    ori = trials(1,1);                                 %Which orientation?
    opto = trials(1,2);
    MovingGratings = trials(1,3);    
    trialcount = trialcount + 1;                        %Trial number
    
    %% Reset:
    gavepassive = 0;
    ValveOpenTime = 0;
    didflip = 0;
    Reaction = [];
    RTrightVec = [];
    RTleftVec = [];
    
    %%  Passive from percentage?
    if GivePercTrial(1) == 1
        passrew = passrew+1;
        set(rewbox,'string',num2str(passrew));
    end
    if setup == 9 || setup==11
        caminitinITI = 0;
    end
    
    %%  check for licks
    if arduino
        set(rightlck,'Background',[1 1 1])
        set(leftlck,'Background',[1 1 1])
        checkLicks
        if llick
            set(leftlck,'Background',[0 0 0])
        end
        if rlick
            set(rightlck,'Background',[0 0 0])
        end
        drawnow
    end
    
    %%  Define which side should be licked for this orientation
    if ori == Par.Orientations(1)
        side = left;
    else
        side = right;
    end
    
    %%  Give passives?
    if passrew && strcmp(side,'left') && passL
        passrew = passrew-1;
        set(rewbox,'string',num2str(passrew));
        passivereward = 1;
    elseif passrew && strcmp(side,'right') && passR
        passrew = passrew-1;
        set(rewbox,'string',num2str(passrew));
        passivereward = 1;
    else
        passivereward = 0;
    end
    
    %%  check for licks
    if arduino
        set(rightlck,'Background',[1 1 1])
        set(leftlck,'Background',[1 1 1])
        checkLicks
        if llick
            set(leftlck,'Background',[0 0 0])
        end
        if rlick
            set(rightlck,'Background',[0 0 0])
        end
        drawnow
    end
    
    %%  Create stimuli
    if strcmp(side,'right')
        x = Par.PixPerDeg*spatialoffset;
    elseif strcmp(side,'left')
        x =  -spatialoffset*Par.PixPerDeg;
    end
    y = Par.PixPerDeg*str2double(get(ydispbox,'String'));
    Par.FigSize = str2double(get(fsizebox, 'string'));
    
    if get(fcheckbox,'value')
        matb1 = repmat(grey, Par.Screeny, Par.Screenx);
        matb2 = repmat(grey, Par.Screeny, Par.Screenx);
        if (setup==9 || setup == 11)
            circle = logical(correctDisplayWF(selectcircle(160, Par.FigSize * Par.PixPerDeg,x, y), 1));
        elseif setup==10
            circle = logical(correctDisplay2Opto(selectcircle(160, Par.FigSize * Par.PixPerDeg,x, y), 1));            
        else
            circle = logical(correctDisplay2(selectcircle(160, Par.FigSize * Par.PixPerDeg,x, y), 1));
        end
        matb1(circle) = 1;
        matb2(circle) = 0;
        mat(:,:,1) = matb1;
        mat(:,:,2) = matb2;
        mat(:,:,3) = matb2;
        written = 0;
        while ~written
            try
                imwrite(mat, [SafeStimPath,'\tmpA_' num2str(setup) '.bmp'])
                written = 1;
            catch
                disp('Trouble writing tmp img')
            end
        end
        cgloadbmp(50, [SafeStimPath '\tmpA_' num2str(setup) '.bmp'])
        cgtrncol(50,'r')
    end
    
    if MovingGratings
        for phaseidx = 1:length(phasesteps)
            cgloadbmp(phaseidx, [SafeStimPath '\' num2str(ori) '_' num2str(phasesteps(phaseidx)) '.bmp'])
            if get(fcheckbox,'value')
                cgsetsprite(phaseidx)
                cgdrawsprite(50,0,0)
                cgsetsprite(0)
            end
        end
    else%still image
        phasestep1 = randsample(1:length(phasesteps),1);
        cgloadbmp(1, [SafeStimPath '\' num2str(ori) '_' num2str(phasesteps(phasestep1)) '.bmp'])
        if get(fcheckbox,'value')
            cgsetsprite(1)
            cgdrawsprite(50, 0, 0)
            cgsetsprite(0)
        end
    end    
    CENT.Contrast = 100; %???
    CENT.Orient = ori;
    if MovingGratings
        CENT.Phase(trialcount) = nan;
    else
        CENT.Phase(trialcount) = phasesteps(phasestep1);
    end
    CENT.Size = Par.FigSize;
    
    %%  ITI (check licks)
    I = '';
    cleanbasetimenow = cleanbasetime + (Par.ItiRandMax -Par.ItiRandMin) * rand;      
    if setup == 9 || setup==11 %imaging
        % Most of the ITI already done
        while toc(ititime) < Par.ITI - caminit - prerecordingtime - cleanbasetimenow - 0.500 %Last 500ms is for opto to start          
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            checkLicks
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
        end
        
        % Clean baseline (no licking period) if enabled
        if cleanbaseline
            cbasetimer = tic;
            while toc(cbasetimer) < cleanbasetimenow
                set(rightlck,'Background',[1 1 1])
                set(leftlck,'Background',[1 1 1])
                checkLicks
                if llick
                    set(leftlck,'Background',[0 0 0])
                    cbasetimer = tic;
                end
                if rlick
                    set(rightlck,'Background',[0 0 0])
                    cbasetimer = tic;
                end
                drawnow
            end
        end
        
        % Turn on microscope, wait for camera to initialize               
        dasbit(0, 1)
        camtime = tic;
        while toc(camtime) < caminit
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            checkLicks
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
        end
        
        % Start recording
        dasbit(2, 1)
        recordtime = tic;
        while toc(recordtime) < prerecordingtime
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            checkLicks
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
        end
    else %standard training:
        while toc(ititime) <  Par.ITI - 0.5 - cleanbasetimenow
            if arduino         %check licks during ITI
                set(rightlck,'Background',[1 1 1])
                set(leftlck,'Background',[1 1 1])
                checkLicks
                if llick
                    set(leftlck,'Background',[0 0 0])
                end
                if rlick
                    set(rightlck,'Background',[0 0 0])
                end
                drawnow
            end
        end
        
        %Clean Baseline
        if cleanbaseline
            cbasetimer = tic;
            while toc(cbasetimer) < cleanbasetimenow
                set(rightlck,'Background',[1 1 1])
                set(leftlck,'Background',[1 1 1])
                checkLicks
                if llick
                    set(leftlck,'Background',[0 0 0])
                    cbasetimer = tic;
                end
                if rlick
                    set(rightlck,'Background',[0 0 0])
                    cbasetimer = tic;
                end
                drawnow
            end
        end        
    end   
    
    %Turn on Opto (if enabled)
    OptoOnOnce = 0;    
    %Optotrial where opto start before trial?
    if opto==1 && ~optoOn && ~OptoOnOnce
        disp('Sending TTL for Opto ON')
        dasbit(optoport,1);
        optoOnTime = tic;
        optoOn = 1;
        OptoOnOnce = 1;
    end
    
    OptoITItime = tic;
    while toc(OptoITItime) < 0.500
        if arduino         %check licks during ITI
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            checkLicks
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
        end
    end
    itiduration = toc(ititime);    
    
    %%  Stimulus Presentation
    if setup==9 || setup==11
        dasbit(3, 1) %Send stimbit of stimulus
    end
    
    %%  Read relative trialtime from arduino
    if arduino %
        fprintf(sport, 'IS');
        I = '';
        to = tic;
        while ~strcmp(I, 'R') && toc(to) < 0.2
            I = fscanf(sport, '%s'); % break test
            if strcmp(I, 'R'); break; end;
        end
        trialtime = str2num(fscanf(sport, '%s'));
        comtime = toc(to);
    end
    
    if MovingGratings        
        if strcmp(side,'right') %
            sendardwm(sport, 'IE 1');
            disp('Enabling 1')
        elseif strcmp(side,'left')
            sendardwm(sport, 'IE 2');
            disp('Enabling 2')
        end
        enable = 1;
    end
    
    stimonset = tic;
    f = 1;                                  %sprite number
    cgdrawsprite(f, 0, 0)
    cgflip(grey,grey,grey)
    if debug
        cgscrdmp
    end
    I = '';                                                         %initialize serial port message
    stimtime = toc(start);                                          %Relative to recordtime    % ????? kan weg?
    justonce =0;    
    while (toc(stimonset) < Par.RewardTime/1000)
        if arduino
            if sport.BytesAvailable
                while ~strcmp(I, 'O') && sport.BytesAvailable
                    I = fscanf(sport, '%s'); % break test
                    if strcmp(I, 'O'); break; end;
                end
                pause(0.001)
                llick = 0;
                rlick = 0;
                Inext = fscanf(sport, '%s');
                switch Inext
                    case 'X' % regular respons
                        %Optotrial where opto start before trial?
                        
                        if opto==2 && ~optoOn && ~OptoOnOnce
                            disp('Sending TTL for Opto ON')
                            dasbit(optoport,1);
                            optoOnTime = tic;
                            optoOn=1;
                            
                            OptoOnOnce = 1;
                        end
                        if sport.BytesAvailable
                            Reaction = fscanf(sport, '%s');
                        end
                        if sport.BytesAvailable
                            RT = fscanf(sport, '%s');
                        end
                        if sport.BytesAvailable
                            passfirst = str2double(fscanf(sport, '%s'));
                        end
                        if sport.BytesAvailable
                            wentthrough = str2double(fscanf(sport, '%s'));
                            set(wt, 'string', num2str(wentthrough));
                        end
                        if sport.BytesAvailable
                            thres1 = fscanf(sport, '%s');
                            disp(['R ' thres1])
                        end
                        if sport.BytesAvailable
                            thres2 = fscanf(sport, '%s');
                            disp(['L ' thres2])
                        end
                        if sport.BytesAvailable
                            lickSide = fscanf(sport, '%s');
                        end
                        if strcmp(Reaction, '0')
                            sendardwm(sport, 'ID');
                            ititime = tic;
                            %                             disp('Disabling')
                            enable = 0;
                            justonce = 1;
                        end
                    case 'Y' % all licks right
                        %Optotrial where opto start before trial?
                        
                        if opto==2 && ~optoOn && ~OptoOnOnce
                            disp('Sending TTL for Opto ON')
                            dasbit(optoport,1);
                            optoOnTime = tic;
                            optoOn=1;
                            
                            OptoOnOnce = 1;
                        end
                        if sport.BytesAvailable
                            RTright = str2num(fscanf(sport, '%s'));
                        end
                        try
                            RTrightVec = [RTrightVec, RTright];
                        catch ME
                            disp(ME)
                        end
                        rlick = 1;
                        
                    case 'Z' % all licks left
                        %Optotrial where opto start before trial?
                        
                        if opto==2 && ~optoOn && ~OptoOnOnce
                            disp('Sending TTL for Opto ON')
                            dasbit(optoport,1);
                            optoOnTime = tic;
                            optoOn=1;
                            
                            OptoOnOnce = 1;
                        end
                        if sport.BytesAvailable
                            RTleft = str2num(fscanf(sport, '%s'));
                        end
                        RTleftVec = [RTleftVec, RTleft];
                        llick = 1;
                    case 'Q' %Timing for opening rewardports
                        if sport.BytesAvailable
                            ValveOpenTime = ValveOpenTime+str2num(fscanf(sport, '%s'));
                        end
                end
                
                pause(0.001)
                
                %% show licks
                set(rightlck,'Background',[1 1 1])
                set(leftlck,'Background',[1 1 1])
                if llick
                    set(leftlck,'Background',[0 0 0])
                end
                if rlick
                    set(rightlck,'Background',[0 0 0])
                end
                drawnow
            end
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
            
            if (setup == 9 || setup == 11) && OnlineFeedback
                if dasreadbit(8) == 1
                    soundsc(ErrorBeep)
                    disp(['Errorbeep... ABORT TRIAL'])
                    pause(Par.StartleTime)
                    aborttrial = 1;
                    dasbit(2, 0);
                    toc(recordtime)
                    camoff = 1;
                    pause(0.005)
                    dasbit(0, 0);
                    if optoOn
                        disp('Sending TTL for Opto OFF')
                        dasbit(optoport,0);
                        optoOFFTime = tic;
                        optoOn=0;
                    end
                    break                    
                end                
            end
        end
        
        %Move gratings
        if MovingGratings && ~didflip
            cgdrawsprite(f,0,0)
            cgflip(grey,grey,grey)
            f = f+1;
            if f > length(phasesteps)
                f = 1;
            end
        end
        
        if passivereward && ~gavepassive && toc(stimonset) > (pdelay/1000 + str2double(get(gpbox, 'string'))/1000) && enable == 1
            sendardwm(sport, 'IP');
            disp('Giving passive')
            gavepassive = 1;
            RewardCount = RewardCount + 0.5;
        end
        
        
        %% END OPTO
        if optoOn && ((Par.RewardTime/1000 - toc(stimonset))<((Par.RewardTime/1000)/2)) % Half of the time opto on, other half not
            disp('Sending TTL for Opto OFF')
            dasbit(optoport,0);
            optoOFFTime = tic;
            optoOn=0;
        end   
    end    
     
    %% Turn of stimulus
    if toc(stimonset) >= Par.RewardTime/1000 && ~justonce %Responsetime over
        ititime = tic;
        sendardwm(sport, 'ID');
        %                     disp('Disabling')
        drawnow
        enable = 0;
        justonce = 1;
    end
    
    if ~didflip
        cgflip(grey,grey,grey)
        didflip = 1;
        if setup==9||setup==11
            dasbit(3, 0) %Send stimbit of stimulus off
        end
    end
    
    %% Process response
    if strcmp(Reaction,'1') %right
        Reaction = 'Hit';
        RewardCount = RewardCount +1;
    elseif strcmp(Reaction, '2') %left
        Reaction = 'Hit';
        RewardCount = RewardCount +1;
    elseif strcmp(Reaction,'0')
        Reaction = 'Error';
        if get(ezbox, 'Value')
            RewardCount = RewardCount +0.5;
        end
    elseif strcmp(Reaction,'4')
        Reaction = 'TooFast';
    else
        Reaction = 'Miss';
        Miss = Miss + 1;
        RT = NaN;
    end
    
    if exist('lickSide','var') && strcmp(lickSide,'1')
        lickSide = 'right';
    elseif exist('lickSide','var') && strcmp(lickSide,'2')
        lickSide = 'left';
    else
        lickSide = 'unknown';
    end
    
    %% Update Hits/Errors in figure window 1
    if strcmp(side,'left') && strcmp(lickSide,'left')
        set(LeftText(2),'string',num2str(str2double(get(LeftText(2),'string'))+1))
    elseif strcmp(side,'right') && strcmp(lickSide,'left')
        set(RightText(2),'string',num2str(str2double(get(RightText(2),'string'))+1))
    elseif strcmp(side,'left') && strcmp(lickSide,'right')
        set(LeftText(1),'string',num2str(str2double(get(LeftText(1),'string'))+1))
    elseif strcmp(side,'right') && strcmp(lickSide,'right')
        set(RightText(1),'string',num2str(str2double(get(RightText(1),'string'))+1))
    end
    
    
    %% Checklicks
    if arduino
        set(rightlck,'Background',[1 1 1])
        set(leftlck,'Background',[1 1 1])
        checkLicks
        if llick
            set(leftlck,'Background',[0 0 0])
        end
        if rlick
            set(rightlck,'Background',[0 0 0])
        end
        drawnow
    end
    
    %% Finish Trial (read keyboard input)
    [kd,kp] = cgkeymap;
    if length(find(kp)) == 1
        if find(kp) == 1;
            ESC = 1;
        end
    end
    
    %% shared LOG data
    if ~debug
        LOG.Par = Par;
        LOG.Setup = setup;
        LOG.Type = 'SimpleLickTask';
        LOG.Mouse = mouse;
        LOG.Date = date;
        LOG.ITI = Par.ITI;
        LOG.shuffleval(trialcount) = shuffleval;
        LOG.Trial(trialcount) = trialcount;
        LOG.Orientation(trialcount) = CENT.Orient;
        LOG.Opto(trialcount) = opto;
        if ~isnan(optoOnTime)
            LOG.OptoOnTime(trialcount) = double(toc(stimonset) - toc(optoOnTime));
        end
        if ~isnan(optoOFFTime)
            LOG.OptoOffTime(trialcount) = double(toc(stimonset)-toc(optoOFFTime));
        end
        LOG.Side(trialcount) = {side};
        if ischar(RT)
            LOG.RT(trialcount) = str2double(RT);
        else
            LOG.RT(trialcount) = RT;
        end
        LOG.FigContrast(trialcount) = CENT.Contrast;
        LOG.Reaction(trialcount) = {Reaction};
        LOG.lickSide(trialcount) = {lickSide};
        LOG.GracePeriod(trialcount) = str2double(get(gpbox, 'string'));
        LOG.Itiduration(trialcount) = itiduration;
        LOG.Gavepassive(trialcount) = gavepassive;
        LOG.Passivefirst(trialcount) = passfirst;
        LOG.PassivPerc(trialcount) = CurPassivPerc;
        LOG.OptoRegulator(trialcount) = OPTOREG;
        LOG.LaserSetting(trialcount) = str2double(get(optoValBox,'string'));
        LOG.Ezbox(trialcount) = ezvalue;
        LOG.Passivedelay(trialcount) = pdelay;
        LOG.Drum(trialcount) = drumm;
        LOG.SemiDrumm(trialcount) = semidrumm;
        LOG.StimPosX(trialcount) = spatialoffset; %degrees
        LOG.StimPosY(trialcount) = str2double(get(ydispbox,'String'));%degrees
        LOG.FullScreen(trialcount) = fullscrfig;
        LOG.MovingGrating(trialcount) = MovingGratings;
        LOG.propl(trialcount) = propR;
        LOG.cleanbasetime(trialcount) = cleanbasetime;
        % Reward times
        LOG.RewardTime(trialcount) = Par.RewardTime;
        
        LOG.spatialoffset = spatialoffset;
        LOG.RewardCount(trialcount) = RewardCount;
        LOG.rdurright(trialcount) = rdurR;
        LOG.rdurleft(trialcount) = rdurL;
        LOG.ValveOpenTime(trialcount) = ValveOpenTime;
        disp(['Valve open for ' num2str(ValveOpenTime)])
        
        %%  Print trial results
        if strcmp(Reaction, 'Error')
            cprintf( 'Red', [ 'Trial: ' num2str(trialcount) ,' Ori: ' num2str(CENT.Orient), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
        elseif strcmp(Reaction, 'Hit')
            cprintf( [0, 0.5, 0], [ 'Trial: ' num2str(trialcount) ,' Ori: ' num2str(CENT.Orient), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
        elseif  strcmp(Reaction, 'Too Early') && (strcmp(lickSide, 'left') || strcmp(lickSide, 'right'))
            cprintf( [0, 0, 1],[ 'Trial: ' num2str(trialcount) ,' Ori: ' num2str(CENT.Orient), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
        elseif strcmp(Reaction,'Too Early') && (strcmp(lickSide,'Unknown'))
            cprintf('SystemCommands',['Trial ' num2str(trialcount), ' Ori: ' num2str(CENT.Orient), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
        else
            disp([ 'Trial: ' num2str(trialcount) , ' Side: ', side ,', Orientation: ' , num2str(CENT.Orient), ' RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string')]);
        end
        
        %%  Update figure and performance
        if strcmp(Reaction,'Hit')
            AllPerfVec(trialcount) = 1;
        elseif strcmp(Reaction,'Error')
            AllPerfVec(trialcount) = 0;
        end
        
        if (sum(~isnan(AllPerfVec)) >= nback)
            LastPerf = [];
            nbackidx = 1:nback;
            stopthisloop = 1;            
            while stopthisloop
                
                LastPerf = AllPerfVec(trialcount-nbackidx+1);
                
                if sum(~isnan(LastPerf)) >= nback || any(trialcount - nbackidx +1 < 1)
                    stopthisloop = 0;
                else
                    nbackidx = [nbackidx nbackidx(end)+1];
                end
            end
            AccVec(trialcount) = sum(LastPerf)./nback;        
        end
        
        %% Checklicks
        if arduino
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            checkLicks
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
        end
        %% Update trial structure
        if strcmp(Reaction, 'Hit') || (~drumm && ~semidrumm)
            trials(1,:) = [];
            GivePercTrial(1) = [];
        elseif semidrumm
            if size(trials,1)>=4
                randplace = randperm(4,1); %Find a random place in the stimulus info matrix
            else
                randplace = randperm(size(trials,1),1);
            end
            trials(randplace+1:end+1,:) = trials(randplace:end,:); %Put the trialinformation of the randomplace at the end of the vector (to not lose this trial)
            trials(randplace,:) = trials(1,:); %Put the information of current trial on the random place
            trials(1,:) = [];
            GivePercTrial(randplace+1:end+1) = GivePercTrial(randplace:end);
            GivePercTrial(randplace) = GivePercTrial(1);
            GivePercTrial(1) = [];
        end
        if  trialcount>1 && LOG.OptoRegulator(trialcount) ~= LOG.OptoRegulator(trialcount-1)
            trials = [];
            if  LOG.OptoRegulator(trialcount) ==1
                optotrialvec = [0,0,0,1,2];
                set(optoValBox,'string',CurLaserSetting)
            else
                optotrialvec = [0];
                set(optoValBox,'string',NaN)
            end
        end
        if isempty(trials) || (trialcount > 1 &&  LOG.shuffleval(trialcount) ~= LOG.shuffleval(trialcount-1)  || (trialcount>1 && LOG.propl(trialcount) ~= LOG.propl(trialcount-1)))
            switch propR
                case 1
                    oriopt = [Par.Orientations(2),Par.Orientations(2),Par.Orientations(1),Par.Orientations(1)];
                case 2
                    oriopt = [Par.Orientations(2) Par.Orientations(1) Par.Orientations(1) Par.Orientations(1)];
                case 3
                    oriopt = [Par.Orientations(2) Par.Orientations(2) Par.Orientations(2) Par.Orientations(1)];
                otherwise
                    error('unknown PropR value')
            end
            
            %% Make orientation matrix: 'details'
            clear details
            z = 1;
            for p = optotrialvec
                for m = movementtrialvec %moving stim is rewarded
                    for o =oriopt  %orientations
                        details(z,:) = [o,p,m];
                        z = z+1;
                    end
                end
            end
                        
            %% Shuffling: create (random) 'order' variable
            [ntrials, nvars] = size(details);
            
            %Shuffle trials?
            if get(shflbox,'value') == 0
                trials = details;
                shuffleval = 0;
            elseif get(shflbox,'value') == 1
                % Order such that it's trying to have a maximum of 3 trials of the same
                % side after each other
                trials = [];
                flag = 0;
                tmpdetails = details;
                while (size(trials,1))<ntrials && ~flag
                    order = randsample(oriopt,4,0);
                    for i = 1:length(order)
                        vec = find(tmpdetails(:,1) == order(i));
                        if length(vec)>0
                            pick = randsample(vec,1);
                        else
                            %More from 1 condition --> Probably because of propl.
                            %Put these extra oin random locations
                            flag=1;
                            break
                        end
                        trials = [trials; tmpdetails(pick,:)];
                        tmpdetails(pick,:)=[];
                    end
                end
                %More from 1 condition --> Probably because of propl.
                %Put these extra oin random locations
                while ~isempty(tmpdetails)
                    for sid = unique(tmpdetails(:,1))'
                        RandIdxVec = find(trials(:,1)==sid);
                        RandIdxVec = RandIdxVec(diff(RandIdxVec') == lraltern);
                        if isempty(RandIdxVec)
                            RandIdx = size(trials,1)+1;
                            trials = cat(1,trials,tmpdetails(find(tmpdetails(:,1)==sid,1),:));
                            tmpdetails(find(tmpdetails(:,1)==sid,1),:)=[];
                            break
                        end
                        count = 0;
                        for RandIdx = RandIdxVec'
                            if isempty(find(tmpdetails(:,1)==sid,1))
                                break
                            else
                                trials = cat(1,trials(1:RandIdx+count,:),tmpdetails(find(tmpdetails(:,1)==sid,1),:),trials(RandIdx+count+1:end,:));
                                tmpdetails(find(tmpdetails(:,1)==sid,1),:) = [];
                            end
                            count = count+1;
                        end
                    end
                end
                
                shuffleval = 1;
            end           
            
            maxcount =countrepeats(trials(:,1));
            if maxcount>4
                warning('More than 4 trials of the same side after each other in trial struct')
            end
            
        end
        
        %Need to read in all phases:
        if trialcount > 1 && LOG.PassivPerc(trialcount) ~= LOG.PassivPerc(trialcount-1)
            changepassives = 1;
        end
        
        if isempty(GivePercTrial) || changepassives
            changepassives = 0;
            GaveExtraAlready = 0;
            GivePercTrial = false(1,size(trials,1));
            if ~get(shflbox,'value') %Ordered trials
                ii = 1;
                curPerc = 0;
                totalPerc = CurPassivPerc/100;
                GivePercTrial(ii:lraltern:length(GivePercTrial)) = 1;
                while sum(GivePercTrial)./length(GivePercTrial) < totalPerc
                    ii = ii+1;
                    GivePercTrial(ii:lraltern:length(GivePercTrial)) = 1;
                end
                if sum(GivePercTrial)./length(GivePercTrial) > totalPerc
                    onevec = find(GivePercTrial == 1);
                    onevec = onevec(randperm(length(onevec)));
                    GivePercTrial(onevec(1:round((sum(GivePercTrial)./length(GivePercTrial) - totalPerc)*length(GivePercTrial))))=0;
                end
            else
                GivePercTrial(1:(round(size(trials,1)*CurPassivPerc/100))) = 1;
                GivePercTrial = GivePercTrial(randperm(size(trials,1)));
            end
        end
        
        %% Checklicks
        if arduino
            set(rightlck,'Background',[1 1 1])
            set(leftlck,'Background',[1 1 1])
            checkLicks
            if llick
                set(leftlck,'Background',[0 0 0])
            end
            if rlick
                set(rightlck,'Background',[0 0 0])
            end
            drawnow
        end
        %add isempty check
        
        %% store lick data
        try
            LOG.RTrightVec{trialcount} = RTrightVec - trialtime;
            LOG.RTleftVec{trialcount} = RTleftVec - trialtime;
        catch
            disp('no trialtime recorded' )
        end
        
        %% Save logs
        if logon
            try
                if ~isdir(fullfile(vcTrainingLogs,mouse,[mouse date],['B' expnum]))
                    mkdir(fullfile(vcTrainingLogs,mouse,[mouse date],['B' expnum]))
                end
                save(fullfile(vcTrainingLogs,mouse,[mouse date],['B' expnum], nts),'LOG');
                if MetaData
                    savejson('', json, fullfile(vcTrainingLogs,mouse,[mouse date],['B' expnum],[expname '_' fields.date '_session.json']));
                end
                
            catch
                save(fullfile(LocalFolder, nts),'LOG');
                if MetaData
                    
                    savejson('', json, fullfile(LocalFolder,[expname '_' fields.date '_session.json']));
                end
            end
        end
        
        %% Check ESC in order to Pause        
        if ESC == 1
            
            PauseButton = questdlg('Want to Pause or to Exit?','Pause Menu','Pause','Exit','Pause');
            
            switch PauseButton
                case 'Pause'
                    
                    %Put lickspout in lickposition
                    %                     if  servoavailable
                    %                         sendardwm(sport, ['IO ' num2str(LickPos)])
                    %                     end
                    
                    ContinueFigure = figure;
                    ContinueButton = uicontrol('Position',[20 20 200 40],'String','Continue',...
                        'Callback','uiresume(gcbf)');
                    uiwait(gcf);
                    close(ContinueFigure);
                    
                    %                     %Put lickspout in breakpos
                    %                     if  servoavailable
                    %                         sendardwm(sport, ['IO ' num2str(BreakPos)])
                    %                     end
                    %
                    ESC = 0;
                    
                    if arduino
                        sendardwm(sport,'IC')
                        if sport.BytesAvailable
                            thres1 = fscanf(sport, '%s');
                            disp(['R ' thres1])
                        end
                        if sport.BytesAvailable
                            thres2 = fscanf(sport, '%s');
                            disp(['L ' thres2])
                        end
                        
                        reconnectard = input('Want to reconnect arduino to get better values? (y/n)\n','s');
                        if strcmp(reconnectard,'y') || strcmp(reconnectard,'Y')
                            if strcmp(get(sport, 'status'), 'open'); fclose(sport); end;
                            if ~exist('sport', 'var')
                                sport = serial('com3');
                            end
                            set(sport,'InputBufferSize', 10240)
                            if strcmp(get(sport, 'status'), 'closed'); fopen(sport); end;
                            
                            set(sport, 'baudrate', 250000);
                            set(sport, 'timeout', 0.5);
                        end
                        
                        
                    end
                    
                case 'Exit'
                    ESC = 1;
                    
            end
        end
        
        if setup==9 || setup==11
            while toc(recordtime)<TotalRecordingtime
                if arduino         %check licks during ITI
                    set(rightlck,'Background',[1 1 1])
                    set(leftlck,'Background',[1 1 1])
                    checkLicks
                    if llick
                        set(leftlck,'Background',[0 0 0])
                    end
                    if rlick
                        set(rightlck,'Background',[0 0 0])
                    end
                    drawnow
                end
            end
            
            % Microscope off
            cgflip(grey,grey,grey)
            dasbit(2, 0);
            toc(recordtime)
            camoff = 1;
            pause(0.005)
            dasbit(0, 0);
            recofftime = toc(start);
        end
    end
end

%% Close it all
try
    if setup==9 || setup==11
        dasbit(1, 0) %Closes shutter
    end
    pause(.1)
    if strcmp(get(sport, 'status'), 'open'); fclose(sport); end;
    pause(.1)
    clear all
    pause(.1)
    cogstd('sPriority','normal')
    pause(.1)
    cgshut
    pause(.1)
catch
    pause(.1)
    cgshut
    clear all
    
end

