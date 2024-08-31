%% Working Memory Trainin g + Opto + Imaging Script
%Training Script Enny van Beest for Mice: Working memory. Start
%wa ter-restriction. Habituation
%to head -restriction before start training.
%Increase level of paradigm pvo (motor) l ickspout. At the end of an
%imaging session, passive licker mouse. This script has 5 different
%phases, and an option for sers will be given and responses recorded.
opengl('save', 'software') 
%% Manual User Settin gs
arduino=1; 
debug = 0;%Debug mode (i.e. Destkop comput er is cogent screen)
savefig = 0; %Save figure at the end?
accuracyThrs = 65; %Accuracy threshold (in performace figure)
nback = 15; %check how last performances for live scores? (per side x ori combination!)
MovingGratings = 1; %if 1; stimuli move. If 0; stimuli don't move
lraltern = 2; %how many of the same stimuli behind eachother? (only for phase 1 & 2)
punishEarlyLicking =0; % if 1, errorbeep when licks are early. If 0, no punishment
% passivesatend = 30; %Number of passive trials at the end (only when imaging)
rng('shuffle'); % avoid pseudorandomization
orientationclearness = 100; %in percentage, how clear is it which orientation it is? (if 100, orientation discrimination possible, otherwise not)
propR = 1; %proportion of trials that is left (rest is right)
cleanbaseline = 1; % Clean baseline before trial starts; helps getting rid of random licks
optoduration = 0.5; % in seconds
if cleanbaseline ==1
    cleanbasetime = 3; %Time there should not have been a lick before starting new trial
else
    cleanbasetime = 0;
end
OnlineFeedback = 1; %Online motion feedback with imaging
FigOpt = [1]; %0 = no figure present (catch), 1 = normal 1 figure present, 2 = both figures present (catch)
MetaData = 1;

%motor
dlgTitle    = 'User Question';
dlgQuestion = 'Are you using the motor?';
choiceM = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'No'); %No = default
servoavailable = 0;
if strcmp(choiceM,'Yes')
    servoavailable = 1
end

if punishEarlyLicking
    disp('DO NOT MAKE SOUND >10%!!')
    pause
end

% LickPos = 70 %70 for imaging set-up, 120 for trainingboxes
% BreakPos = 64 %64 for imaging set-up, 126 for trainingboxes
%% Paths + Current Directory
user_name =getenv('USERNAME');
tmpname = cd;
if ~isempty(strfind(tmpname,[user_name '.HERSEN']))
    user_name = regexp(tmpname,'\','split');
    user_name = user_name{3};
end
tmppath=mfilename('fullpath');
tmpfile = mfilename;

WMpath = strsplit(tmppath,tmpfile);
WMpath = WMpath{1};

LocalFolder = ['C:\Users\' user_name '\Documents\MATLAB\WorkingMemoryTask\'];
try
    addpath(genpath(WMpath))
    cd(WMpath)
catch ME
    disp(ME)
    WMpath = ['C:\Users\' user_name '\Dropbox\Mouse\WM\'];
    addpath(genpath(WMpath))
    cd(WMpath)
end
%Make sure that subscripts is the most current folder (so that scripts will
%be taken from this subfolder)
addpath(genpath([WMpath 'subscripts']))
%% Upload arduino
if arduino
    disp('Sure you uploaded correct arduino script? Press key to continue or abort...')
    pause
end

%% Better safe than sorry!
if servoavailable
    disp('Paused until made sure that lickspout is far away from the mouse ...')
    disp('Press key to continue...')
    pause
end

%% Training setup
setup = {'1', '2', '3', '4', '5', '6', 'Testing', 'Recording','Imaging','Opto','Imaging+Opto'};
setup = menu('Choose the setup', setup);

if setup == 9 || setup==11
    addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github'],'OpenAccessStorage')))
else
    addpath(genpath(fullfile(['C:\Users\' user_name '\Dropbox\'],'OpenAccessStorage')))
end

%% Start Serial connection with Arduino
if arduino
    if ~exist('sport', 'var')
        if setup==9 || setup==11
            sport = serial('com9');
        elseif setup==10
            sport=serial('com20');            
        else
            try
                sport = serial('com3');
            catch
                sport = serial('com4');
            end
        end
    else
        clear LOG
    end
    set(sport,'InputBufferSize', 10240)% was 10240

    try
        if strcmp(get(sport, 'status'), 'closed'); fopen(sport); end;
    catch ME
        disp(ME)
       keyboard
    end
    set(sport, 'baudrate', 250000); %was 250000
    set(sport, 'timeout', 0.5);
end

%% Box/Imaging setup
if setup == 9 || setup==11
    eval(['run ' fullfile(['C:\Users\' user_name '\Documents\Github'],'Mouse','imagingStim','runSettings.m')])
    D2ScreenWidth = 51./2; %??
    ITI = 12;
    LickPos = 70;    %motor
    BreakPos = 64;   %motor
    OnlineFeedback = 1
    displayscript = 'correctDisplayWF'
    gammaconversion = 'gammaconWF'
    addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github'],'OpenAccessStorage')))%Json files
    addpath('Z:\OpenAccessStorageParameters')
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
    addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github'],'OpenAccessStorage')))%Json files
    addpath('Z:\OpenAccessStorageParameters')
elseif setup==10
    screenx = 1920;
    screeny = 1200;
    refresh = 60;
    gammaconversion = 'gammaconOpto'
    displayscript = 'correctDisplay2Opto'
    ITI = 7;
    LickPos = 70;    %motor
    BreakPos = 60;   %motor
    prefsf = 0.075;
    prefspeed = 24; %deg/s
    ScreenDistance = 11;
    D2ScreenWidth = 51/2;     % width of screen, in cm
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
    addpath(genpath(fullfile(['C:\Users\' user_name '\Dropbox\OpenAccessStorage'],'OpenAccessStorage')))%Json files
    addpath('Z:\OpenAccessStorageParameters')
end
caminit = 1.5;
prerecordingtime = 0.5;

vcTrainingLogs = ['\Path\To\LogFiles'];

%% Session name (nts)
warning on backtrace

logon = 0;
if setup == 9 || setup==11   % Imaging training: Don't change; imaging has same way of giving names etc.
    waitfor(msgbox({'Final check:' '' '- Zoom = 10.0x' '' '- Arduino cables' '' '- Motion recorder ON' '' '- Microscope position (lowered?)' '' '- Nr of Stimuli (reset?)'},'Final Check:'))
    prompt = {'Mouse Name', 'Exp Nr', 'Date', 'Exposure'};
    dlg_title = 'Please enter parameters';
    num_lines = 1;
    def = {'Name', '1', datestr(datenum(date), 'yyyymmdd'), '50', '1'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~strcmp(answer{1}, 'Name')
        
        mouse = answer{1};
        expnum = answer{2};
        date = answer{3};
        exposure = str2double(answer{4});
        
        expname = [mouse expnum];
        folder = ['\\Nin518\e\Imaging\' [mouse date]];
        logfile = [folder '\' expname '\' expname '.mat'];
        
        LOG.Start = -prerecordingtime;
        LOG.Stimulus = 'SpatialWM';
        LOG.Exposure = exposure;
        logon = 1;
        
        if exist(logfile, 'file')
            error('Log-file already used, please pick a different name.');
        end
        
        % add dasbit stuff
        addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github\Mouse\DIO24'])))
    else
        mouse = 'Unknown';
        leftside = round(rand(1,1))+1; % if 1, left is 0 degrees, if 2, left is 90 degrees
        CurPhase = 1;
        ezvalue = 1;
        spatialoffset = screenx/4;
    end
    nts = [mouse '_' date '_B' expnum];
    expname = [mouse num2str(expnum)];
    
else    % standard training
    
    logon = 1;   %Save LOG
    choice = menu('Which Mouse?','Illy','Jemison','Galilei','Karlsen','Jenkins','Lully','Other','Test (no log)');
    switch choice
        case 1
            mouse = 'Illy';
        case 2
            mouse = 'Jemison' ;
        case 3
            mouse = 'Galilei';
        case 4
            mouse = 'Karlsen';
        case 5
            mouse = 'Jenkins';
        case 6
            mouse = 'Lully';            
        case 7
            mouse = input('Enter mouse name:' ,'s');            
        case 8
            mouse = 'Test';
            logon = 0; %don't save LOG
    end
    
    nts = [mouse '_' datestr(date,'yyyymmdd') '_B'];
    expnum = '1' ;
    
    nts = sprintf('%s',nts,num2str(expnum));
    expname = [mouse num2str(expnum)];
    
    fprintf('%s\n',nts);
end

% choice = menu('Which Memory version?','All trials delay','shuffled no delay and delay');
% switch choice
%     case 1
memoryshuffle = 0;
% %     case 2
%         memoryshuffle = 1;
% end

while exist(fullfile(vcTrainingLogs,[nts '.mat']), 'file')
    disp(['Log-file expnum ' num2str(expnum) ' already used, expnum number changes...']);
    if isstr(expnum)
        expnum = str2num(expnum);
    end
    expnum = expnum+1;
    nts = [mouse '_' datestr(date,'yyyymmdd') '_B' num2str(expnum)]
    expname = [mouse num2str(expnum)];
    
end
expname = [mouse expnum];

%% JSON FILE
if MetaData
    if setup == 11
        fields.project = 'WorkingMemory_EB';
        fields.dataset = 'WM_WithOpto';
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
Par.StimDurMin = 500; %min. ms stim on, rest depenent on delay
Par.Orientations = [45 135]; %With Old mice : [0 90]),new mice [ 45 135];
Par.BeepFreq = 2000; %Hz
Par.BeepTime = 0.5;  %Enter time constant beep
Par.StartleFreq = 5000; %Hz :Switched to WHITENOISE! %Tone for preliminary licking
Par.StartleTime = 1.5;

Par.SpatFreq = 0.08;%0.05; 12.5 degree sper cycle m ?
Par.minDelay = -500; %ms
Par.maxDelay = 1500; %ms --> Was 2000, now 1500
Par.DelayStep = 200; %Only used when mice do complete task
Par.RewardTime =2000; %max reward time
Par.DelayOptions = Par.minDelay:Par.DelayStep: Par.maxDelay;
Par.FigSize = 35;

Par.Yoffset = 15;


TotalRecordingtime = 2+caminit + prerecordingtime + optoduration + Par.RewardTime/1000 + Par.maxDelay/1000 + Par.StimDurMin/1000;
Par.TotalRecordingtime = TotalRecordingtime;
LOG.Duration = TotalRecordingtime;

%% Initialize default values for other variables
%controlbox (overwritten later on by 'load CurrentPhase.m' if existing)
global RTrightVec
global RTleftVec
NewTrialMat = 0;
CurPhase = 1;
ezvalue = 1;
withflash = 0;
flashlum = 100;
shuffleval = 0;
curdelay = Par.DelayOptions(1);
CurPassivPerc = 50;
freeze = 0;
drumm = 0;
semidrumm = 1;
passL = 1; %passives ON (1)/OFF (0)
passR = 1; %passives ON (1)/OFF (0)
RepError = 0; %Counts Repetition of same error
ErrorTimeOut = 0; %Time out of 5 seconds when error
%Trial-specific
passfirst = 0;
ESC = 0;
Miss = 0;
trialcount = 0;
gavepassive = 0;
wentthrough = 0;
details = [];
ThisContrTrialcount = 1;
TimePunish = 1.5; %Time punishment if too early response
prefspeed = 24; %deg/s
changepassives = 0;
GaveExtraAlready = 0;
RewardCount = 0; %Counts rewards given for callibration purposes
Reaction = [];
OptoAllowed = 0;
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

%% Load CurrentPhase.m or initialize new variables
if exist(fullfile(vcTrainingLogs,[mouse '_CurrentPhase.mat']),'file')
    load(fullfile(vcTrainingLogs,[mouse '_CurrentPhase.mat'])) % In this file it says per mouse whether 0 is coupled to right or left and what the current learning phase is
elseif  exist(fullfile(LocalFolder,[mouse '_CurrentPhase.mat']),'file')
        load(fullfile(LocalFolder,[mouse '_CurrentPhase.mat'])) % In this file it says per mouse whether 0 is coupled to right or left and what the current learning phase is
else %Create new variables if necessary
    leftside = 1;
    spatialoffset=40;
    fullscrfig = 1;
end

%Re-initalize
Par.minDelay = -500; %ms
Par.maxDelay = 1500; %ms --> Was 2000, now 1500
Par.DelayStep = 200; %Only used when mice do complete task

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

if spatialoffset > 100
    spatialoffset = round(Par.DegPerPix*spatialoffset);
end

Par.Yoffset = 15;

%% Open control box GUI
OpenContrDiscWindow_EnnyV2
set(ydispbox,'string',num2str(Par.Yoffset));

%% Initialize  Cogent
cgloadlib

cgshut %shut cg if one's open
if debug
    warning('Debug Mode enabled')
    cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 0)
else
    if setup == 9 || setup==11
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

%% Performance Figure Online
AllPerfVec = nan(1,500);
TooEarlyVec = nan(1,500);
Perf0 = nan(1,nback);
Perf0NoOpto = Perf0;
Perf90 = nan(1,nback);
Perf90NoOpto = Perf90;
idx0 = 1;
idx90 = 1;
Par.DelayOptions = Par.minDelay:Par.DelayStep: Par.maxDelay;
DelayOptions4Plot = Par.minDelay:Par.DelayStep/2:Par.maxDelay;
Delayvec = [];
Delayidx = 1;
PerfMat = nan(length(Par.Orientations),length(DelayOptions4Plot));
if  exist('HP')==0
    HP = figure('Name',['Online Performance ' nts]);
    HAX = axes('parent',HP);
    hold on
    line([-520 2020],[accuracyThrs accuracyThrs],'Color',[1 0 0])
    
    xlabel('ResponseDelay')
    ylabel('Performance (%)')
    axis([-520 2020 0 100])
end

%% Create sound stimuli
% create reward sound
rfreq=2*pi*Par.BeepFreq;                            %Convert to radian frequency
Ts=[linspace(0,Par.BeepTime,8192*Par.BeepTime)];    %Specify sample times over 4 seconds
envelope=exp(-Ts/Par.BeepTime);                     %Calculate the amplitude envelope
signal=cos(rfreq*Ts);                               %calculate the cosine
ResponseBeep=envelope.*signal;                      %apply the envelope to the signal

%Create startle-stimulus: sound + display
%sound
rfreq = 2*pi*Par.StartleFreq;
Ts=[linspace(0,Par.StartleTime,8192*Par.StartleTime)];                    %Specify sample times over 4 seconds
envelope=exp(-Ts/Par.StartleTime);                            %Calculate the amplitude envelope
signal=cos(rfreq*Ts);                             %calculate the cosine
ErrorBeep=envelope.*signal;                       %apply the envelope to the signal
%display
cgflip(grey,grey,grey)

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
        if (setup==9 || setup==11)
            a1 = correctDisplayWF(makegrating2GammaconV(160,CURR.Orient,Par.PixPerDeg * 1/prefsf, CURR.Phase, minlum, maxlum,gammaconversion));
        elseif setup==10
            a1 = correctDisplay2Opto(makegrating2GammaconV(160,CURR.Orient,Par.PixPerDeg * 1/prefsf, CURR.Phase, minlum, maxlum,gammaconversion));
        else
            a1 = correctDisplay2(makegrating2GammaconV(160,CURR.Orient,Par.PixPerDeg * 1/prefsf, CURR.Phase, minlum, maxlum,gammaconversion));
        end
        imwrite(a1,[SafeStimPath '\' num2str(Par.Orientations(j)) '_' num2str(phasesteps(i)) '.bmp'])
    end
end

if servoavailable
    sendardwm(sport, ['IO ' num2str(LickPos)])
    disp('Paused... Put lickspout in correct position now!')
    pause
end
if servoavailable
    %% Test LickPos+BreakPos
    while 1
        dlgTitle    = 'User Question';
        dlgQuestion = 'Test motor?';
        choiceR = questdlg(dlgQuestion,dlgTitle,'Yes','No', 'Yes'); %No = default
        if strcmp(choiceR,'No')
            break
        end
        prompt = {'BreakPos', 'LickPos','# of test runs'};
        dlg_title = 'Motor settings';
        num_lines = 1;
        def = {num2str(BreakPos), num2str(LickPos), '1'};
        m_answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if ~isempty(m_answer)%
            
            BreakPos = str2num(m_answer{1});
            LickPos = str2num(m_answer{2});
            n=1;
            fprintf('\n')
            while n <= str2num(m_answer{3})
                sendardwm(sport, ['IO ' num2str(BreakPos)])
                fprintf('B')
                pause(1)
                sendardwm(sport, ['IO ' num2str(LickPos)]);
                fprintf('L')
                pause(1)
                n=n+1;
            end
        end
        %
    end
    fprintf('\n')
end

%% Final settings
ititime = tic;
start = tic;

if leftside == 1
    left = 'left';
    right = 'right';
elseif leftside == 2
    left = 'right';
    right = 'left';
end


%% Dasbit
if setup == 9 || setup==11
    dasinit(22)
    dasbit(0, 0) %initializes camera (set to 0)
    dasbit(3, 0) %stimbit
    dasbit(1, 1) %Shutter open
    disp('Shutter Opened')
    dasbit(2, 0)  %starts recording
end

%Opto choices possible, but not now:
% 0 = off
% 1 = baseline (-500 to 0ms)
% 2 = visual (-100 to 600ms)
% 3 = delay (600 to 1800ms)
% 4 = response (1800-2400ms)
% 5 = post-trail (>4000)
% if setup == 10
%     dasinit(23);
%     optoport = 0;
%     dasbit(optoport,0);
%     optotrialvec = [0,0,0,0,0,1,2,3,3,4,5]; %
%     set(OptoBox,'value',0);
% else
if setup==11
    optoport = 4;
    dasbit(optoport,0);
    optotrialvec = [0,0,0,0,0,1,2,3,3,4,5]; %
    set(OptoBox,'value',0);
else
    optotrialvec = [0]; %0 = off
    set(OptoBox,'value',0);
end

% if any(optotrialvec)
%     punishEarlyLicking=0;
% end

figure(1) %show control-window
%motor
if servoavailable
    sendardwm(sport, ['IO ' num2str(BreakPos)])
end

%% Outer Trialloop
while ~ESC
    %% Update settings
    CurPhase = str2double(get(PhaseBox,'string'));
    curdelay = str2double(get(DelayBox,'string'));
    pdelay = str2double(get(pdbox,'string'));
    CurPassivPerc = str2double(get(PercPassBox,'string'));
    memoryshuffle = get(memoryBox,'value');
    OptoAllowed = get(OptoBox,'value');
    
    passrew = str2double(get(rewbox,'string'));
    passL = get(plbox,'value');
    passR = get(prbox,'value');
    
    freeze = get(frzbox,'value');
    shuffleval = get(shflbox,'value');
    ezvalue = get(ezbox,'value');
    drumm = get(drummbox,'value');
    semidrumm = get(semidrummbox,'value');
    threshold = str2double(get(thbox,'string'));
    rdurL = str2double(get(rdurleft,'string'));
    rdurR = str2double(get(rdurright,'string'));
    spatialoffset = str2num(get(xdispbox,'String'));
    %     orientationclearness = str2double(get(OriClearbox,'String'));
    
    MovingGratings = get(movingBox,'value');
    fullscrfig = get(fcheckbox,'Value');
    propR = get(PropBox,'value');
    optoOnTime = nan;
    optoOFFTime = nan;
    
    %% Phase specific settings
    switch CurPhase
        case 1
            set(pdbox, 'string','0')
            if exist('NewPhase','var')
                CurPassivPerc = 50; %fill in maximum percentage
            end
            accuracyThrs = 60;
            disp('Threshold for this phase is 60% accurate')
            set(gpbox, 'string',num2str(200))
        case 2
            set(pdbox,'string','0') %passive delay
            CurPassivPerc = 0;
            set(PercPassBox,'string', num2str(CurPassivPerc))
            accuracyThrs = 65;
            if curdelay>0
                set(gpbox, 'string',num2str(0))
            else
                set(gpbox,'string',num2str(200))
            end
        case 3
            set(pdbox,'string','0') %passive delay
            CurPassivPerc = 0;
            set(PercPassBox,'string', num2str(CurPassivPerc))
            accuracyThrs = 65;
            if curdelay>0
                set(gpbox, 'string',num2str(100))
            else
                set(gpbox,'string',num2str(200))
            end
            set(shflbox,'value',1)
    end
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
    if memoryshuffle
        for m = [0,1,1,1] %25% of trials lower memory
            for p = optotrialvec %Optogenetics
                for o = oriopt  %orientations
                    details(z:z+lraltern-1,:) = repmat([o,p,m],[lraltern,1]);
                    z = z+lraltern;
                end
            end
        end
    else
        for p = optotrialvec
            for o =oriopt  %orientations
                details(z:z+lraltern-1,:) = repmat([o,p],[lraltern,1]);
                z = z+lraltern;
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
        OptoflagOff = 0;
        tmpdetails = details;
        while (size(trials,1))<ntrials && ~OptoflagOff
            order = randsample(oriopt,4,0);
            for i = 1:length(order)
                vec = find(tmpdetails(:,1) == order(i));
                if length(vec)>0
                    pick = randsample(vec,1);
                else
                    %More from 1 condition --> Probably because of propl.
                    %Put these extra oin random locations
                    OptoflagOff=1;
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
    
    %% INNER TRIAL LOOP
    NewPhase = 0;
    while ~ESC && ~NewPhase
        %% Update Controlbox settings
        CurPhase = str2double(get(PhaseBox,'string'));
        curdelay = str2double(get(DelayBox,'string'));
        pdelay = str2double(get(pdbox,'string'));
        CurPassivPerc = str2double(get(PercPassBox,'string'));
        memoryshuffle = get(memoryBox,'value');
        OptoAllowed = get(OptoBox,'value');
        if OptoAllowed
            set(semidrummbox,'value',0); 
        end
        OptoflagOff = 0;
        OptoflagOn = 0;
        
        passrew = str2double(get(rewbox,'string'));
        passL = get(plbox,'value');
        passR = get(prbox,'value');
        
        freeze = get(frzbox,'value');
        shuffleval = get(shflbox,'value');
        ezvalue = get(ezbox,'value');
        drumm = get(drummbox,'value');
        semidrumm = get(semidrummbox,'value');
        
        threshold = str2double(get(thbox, 'string'));
        rdurL = str2double(get(rdurleft, 'string'));
        rdurR = str2double(get(rdurright, 'string'));
        spatialoffset = str2double(get(xdispbox,'String'));
        
        MovingGratings = get(movingBox,'value');
        fullscrfig = get(fcheckbox,'Value');
        propR = get(PropBox,'value');
        
        %% Update Arduino settings
        if arduino
            sendardwm(sport, ['IL ' num2str(rdurR)]) %This seems strange, but fits with the arduino script which has been build the other way around
            sendardwm(sport, ['IR ' num2str(rdurL)])
            sendardwm(sport, ['IF ' get(thbox, 'string')])
            sendardwm(sport, ['IM ' num2str(get(ezbox, 'Value'))])
        end
        
        %%  Initialize Trial
        aborttrial = 0;
        ori = trials(1,1);                                 %Which orientation?
        opto = trials(1,2); %which opto condition?
        if ~OptoAllowed
            opto = 0;
            OptoflagOn =1;
            OptoflagOff = 1;
        end
        
        if ori==0
            keyboard
        end
        try
            if memoryshuffle
                memory = trials(1,3);
            else
                memory = 1;
            end
        catch ME
            disp(ME)
        end
        savecurdelay = curdelay;
        
        trialcount = trialcount + 1;                        %Trial number
        %         if ~exist('LOG','var') || ~isfield(LOG,'Reaction')
        %             curdelay = Par.DelayOptions(1);
        %             opto = 0;
        %             OptoflagOn =1;
        %             OptoflagOff = 1;
        %             OptoAllowed = 0;
        %         elseif sum(strcmp(LOG.Reaction,'Hit'))<20 && curdelay>Par.minDelay%First trials make sure curdelay is shorter to get used to trials again (4 hits per time delay)
        %             curdelay = floor(sum(strcmp(LOG.Reaction,'Hit'))/4)*(((str2double(get(DelayBox,'string'))/5)-Par.minDelay)+Par.minDelay);
        %             opto = 0;
        %             OptoflagOn =1;
        %             OptoflagOff = 1;
        %             OptoAllowed = 0;
        %         elseif memory == 0 %In case of memory = 0, change curdelay
        %             savecurdelay = curdelay;
        %             curdelay = 0;
        %             OptoAllowed = 0;
        %             opto = 0;
        %             OptoflagOn =1;
        %             OptoflagOff = 1;
        %         end
        
        %         if trialcount==1 || (trialcount>1 && sum(strcmp(LOG.Reaction,'Hit'))<10)%First trials make sure curdelay is shorter to get used to trials again (4 hits per time delay)
        %             shuffleval=0;
        %             drumm=1;
        %             if trialcount<6
        %                 passrew = 1;
        %                 pdelay =0;
        %             end
        %         end
        
        %Give passive if more than 3 errors on same side
        if RepError>=3
            if length(unique(LOG.Side(trialcount-3:trialcount-1)))==1
                passrew=1;
                pdelay=1;
                RepError=0;
            end
        end
        
        stimdur = Par.StimDurMin;     %stimulus duration
        % Reset:
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
        
        
        %%  Spatial offset (x/y)
        if strcmp(side,'right')
            x = Par.PixPerDeg*spatialoffset;
        elseif strcmp(side,'left')
            x =  -spatialoffset*Par.PixPerDeg;
        end
        y = Par.PixPerDeg*str2double(get(ydispbox,'String'));
        
        %%  Create (circle) figure (if enabled)
        Par.FigSize = str2double(get(fsizebox, 'string'));
        if get(fcheckbox,'value')
            matb1 = repmat(grey, Par.Screeny, Par.Screenx);
            matb2 = repmat(grey, Par.Screeny, Par.Screenx);
            if (setup==9||setup==11)
                circle = logical(correctDisplayWF(selectcircle(160, Par.FigSize * Par.PixPerDeg,x, y), 1));
            elseif setup == 10
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
        
        %%  load moving gratings OR still image
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
        
        %%  save settings for log (CENT)
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
        randomITI = Par.ItiRandMax * rand + ErrorTimeOut;
        
        % Clean baseline
        TotalCleanbasetime = tic;
        if cleanbaseline
            cbasetimer = tic;
            while toc(cbasetimer) < cleanbasetime
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
                if llick || rlick
                    cbasetimer = tic;
%                     soundsc(ErrorBeep)
                end
            end
        end
        TotalCleanbasetime = toc(TotalCleanbasetime);
        
        if setup == 9 || setup==11 %imaging
            while toc(ititime) < Par.ITI - caminit - prerecordingtime - TotalCleanbasetime-optoduration%No Par --> Loads training variables
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
            recordtime = tic;
            
            while toc(ititime) < Par.ITI - cleanbasetime - optoduration
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
        end
        set(rightlck,'Background',[1 1 1])
        set(leftlck,'Background',[1 1 1])
        checkLicks
        if llick
            set(leftlck,'Background',[0 0 0])
        end
        if rlick
            set(rightlck,'Background',[0 0 0])
        end
        
        %Optotrial where opto start before trial?
        if opto==1 && ~OptoflagOn
            disp('Sending TTL for Opto ON')
            dasbit(optoport,1);
            optoOnTime = tic;
            OptoflagOn = 1;
        end
        
        OptoITItime = tic;
        while toc(OptoITItime) < optoduration
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
        
        %Optotrial where opto during baseline only?
        if opto==1 && OptoflagOn
            disp('Sending TTL for Opto OFF')
            dasbit(optoport,0);
            optoOFFTime = tic;
            OptoflagOff=1;
        end
        
        itiduration = toc(ititime);
        
        %%  Stimulus Presentation
        if setup==9 || setup==11
            dasbit(3, 1) %Send message of stimulus
        end
        
        stimonset = tic;
        f = 1;                                  %sprite number
        cgdrawsprite(f, 0, 0)
        cgflip(grey,grey,grey)
        if debug
            cgscrdmp
        end
        stimtime = toc(start);                                          %Relative to recordtime    % ????? kan weg?
        
        earlylicks = 0;
        TooEarly = 0;                                                   %When response is premature
        
        didflip = 0;
        I = '';
        enable = 0;
        justonce = 0;
        motorjustonce =0;
        rewardonset = uint16(0);
        while toc(stimonset) < (stimdur/1000  + curdelay/1000 + Par.RewardTime/1000)  % while stimdur+curdelay+time punishmentm,
            if arduino
                if sport.BytesAvailable
                    while ~strcmp(I, 'O') && sport.BytesAvailable
                        I = fscanf(sport, '%s'); % break test
                        if strcmp(I, 'O'); break; end;
                    end
                    %             disp('Escaped')
                    pause(0.001)
                    llick = 0;
                    rlick = 0;
                    Inext = fscanf(sport, '%s');
                    switch Inext
                        case 'X' % regular respons
                            %Optotrial where opto start before trial?
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
                            if strcmp(Reaction, '0')
                                sendardwm(sport, 'ID');
                                ititime = tic;
                                %                             disp('Disabling')
                                enable = 0;
                            end
                        case 'Y' % all licks right
                            if sport.BytesAvailable
                                RTright = str2num(fscanf(sport, '%s'));
                            end
                            RTrightVec = [RTrightVec, RTright];
                            rlick = 1;
                        case 'Z' % all licks left
                            if sport.BytesAvailable
                                RTleft = str2num(fscanf(sport, '%s'));
                            end
                            try
                                RTleftVec = [RTleftVec, RTleft];
                            catch ME
                                disp(ME)
                            end
                            llick = 1;
                        case 'Q' %Timing for opening rewardports
                            if sport.BytesAvailable
                                ValveOpenTime = ValveOpenTime+str2num(fscanf(sport, '%s'));
                            end
                    end
                    
                    pause(0.001)
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
                
                if (setup == 9 || setup == 11) && OnlineFeedback && toc(stimonset) < (stimdur/1000  + curdelay/1000)
                    if dasreadbit(8) == 1
                        soundsc(ErrorBeep)
                        disp(['Errorbeep... ABORT TRIAL'])
                        if OptoflagOn && ~OptoflagOff
                            disp('Sending TTL for Opto OFF')
                            dasbit(optoport,0);
                            optoOFFTime = tic;
                            OptoflagOff=1;
                        end
                        pause(Par.StartleTime)
                        aborttrial = 1;
                        dasbit(2, 0);
                        recordtime = tic;
                        camoff = 1;
                        pause(0.005)
                        dasbit(0, 0);
                        break
                    end
                end
                if punishEarlyLicking && toc(stimonset) < (stimdur/1000  + curdelay/1000)
                    if ~isempty(RTrightVec) || ~isempty(RTleftVec)
                        soundsc(ErrorBeep)
                        disp(['Errorbeep... ABORT TRIAL'])
                        if OptoflagOn  && ~OptoflagOff
                            disp('Sending TTL for Opto OFF')
                            dasbit(optoport,0);
                            optoOFFTime = tic;
                            OptoflagOff=1;
                        end
                        TooEarly = 1;
                        pause(Par.StartleTime)
                        aborttrial = 1;
                        recordtime = tic;
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
            
            %Turn off stimulus when stimulus is on the screen longer than
            %stimdur
            if ~didflip && toc(stimonset) >= (stimdur)/1000
                cgflip(grey,grey,grey)
                didflip = 1;
                if setup==9 || setup==11
                    dasbit(3,0)
                end
            end
            
            % Start rewardtime if stimdur + delay is done
            
             %         disp(['Communication with Arduino took ' num2str(toc(to))])
             if (toc(stimonset) >= (stimdur/1000  + curdelay/1000)) && ~motorjustonce
                 if servoavailable && ~aborttrial
                     sendardwm(sport, ['IO ' num2str(LickPos)])
                     motorjustonce=1;
                 end
             end
            
            if (toc(stimonset) >= (stimdur/1000  + curdelay/1000 + str2double(get(gpbox, 'string'))/1000)) && ~justonce
                justonce = 1;
                startRewardPeriod = toc(stimonset);
%                 soundsc(ResponseBeep);                                    % play the reward cue sound
                
                if strcmp(side,'left')
                    set(cueboxL,'Background',[0 0 0]) %black
                elseif strcmp(side,'right')
                    set(cueboxR,'Background',[0 0 0]) %black
                end
                drawnow
                if setup == 9 || setup==11
                    dasbit(3, 1) %StimSync with soundcue
                end
                %  Read relative trialtime from arduino
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
               
                
                if arduino
                    if enable == 0
                        if strcmp(side,'right') %
                            sendardwm(sport, 'IE 1');
                            %                         disp('Enabling 1')
                        elseif strcmp(side,'left')
                            sendardwm(sport, 'IE 2');
                            %                         disp('Enabling 2')
                        end
                        enable = 1;
                        rewardonset = tic;
                    end
                end
            end % of reward condition
            
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
            if enable && toc(rewardonset) >= Par.RewardTime/1000  %Responsetime over
                ititime = tic;
                sendardwm(sport, 'ID');
                %                     disp('Disabling')
                drawnow
                enable = 0;
            end
            
            if enable == 1 && passivereward && ~gavepassive && toc(rewardonset) > (pdelay + str2double(get(gpbox, 'string')))/1000
                sendardwm(sport, 'IP');
                disp('Giving passive')
                gavepassive = 1;
                RewardCount = RewardCount + 0.5;
            end
            
            %% Start or end OPTO
            %Opto choices possible, but not now:
            % 0 = off
            % 1 = baseline (-500 to -0ms)
            % 2 = visual (0 to 500ms)
            % 3 = delay (1000 to 1500ms)
            % 4 = response (2000 to 2500ms)
            % 5 = post-trial (3500 to 4000)
            if OptoflagOn && ~OptoflagOff
                if toc(optoOnTime)>=optoduration
                    OptoflagOff =1;
                    disp('Sending TTL for Opto OFF')
                    dasbit(optoport,0);
                    optoOFFTime = tic;
                end
            elseif ~OptoflagOff && ~OptoflagOn
                switch opto
                    case 2
                        if toc(stimonset)>=0
                            OptoflagOn = 1;
                        end
                    case 3
                        if toc(stimonset)>=1
                            OptoflagOn = 1;
                        end
                    case 4
                        if toc(stimonset)>=2
                            OptoflagOn = 1;
                        end
                end
                if OptoflagOn
                    disp('Sending TTL for Opto ON')
                    dasbit(optoport,1);
                    optoOnTime = tic;
                end
            end
        end
        
        if  enable
            ititime = tic;
            sendardwm(sport, 'ID');
            %                     disp('Disabling')
            drawnow
            enable = 0;
            justonce = 1;
        end
        
        if ~didflip && (toc(stimonset) >= (stimdur)/1000)
            cgflip(grey,grey,grey)
            didflip = 1;
        end
        
        set(cueboxL,'Background',[1 1 1])
        set(cueboxR,'Background',[1 1 1])
        drawnow
        
        cgflip(grey,grey,grey)
        
        if setup==9|| setup==11
            dasbit(3, 0) %Send message of stimulus           
        end
        
        % In case of early break out of the previous loop, turn off opto
        if TooEarly == 1 || aborttrial==1          
            if OptoflagOn && ~OptoflagOff
                disp('Sending TTL for Opto OFF')
                dasbit(optoport,0);
                optoOFFTime = tic;
                OptoflagOff=1;
            end
        end
        
        %motor
        if servoavailable
            sendardwm(sport, ['IO ' num2str(BreakPos)])
        end
        
        %Make recordingtime full
        while toc(recordtime) < TotalRecordingtime
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
            
            %% Opto post-trial
             % 0 = off
            % 1 = baseline (-500 to -0ms)
            % 2 = visual (0 to 500ms)
            % 3 = delay (1000 to 1500ms)
            % 4 = response (2000 to 2500ms)
            % 5 = post-trial (3500 to 4000)
            if ~TooEarly && ~aborttrial
                if OptoflagOn && ~OptoflagOff
                    if toc(optoOnTime)>=optoduration
                        OptoflagOff =1;
                        disp('Sending TTL for Opto OFF')
                        dasbit(optoport,0);
                        optoOFFTime = tic;
                    end
                elseif ~OptoflagOff && ~OptoflagOn
                    if opto==5
                        if toc(stimonset)>=3.5
                            OptoflagOn = 1;
                        end
                    end
                    if OptoflagOn
                        disp('Sending TTL for Opto ON')
                        dasbit(optoport,1);
                        optoOnTime = tic;
                    end
                end
            end
        end
        
        if setup == 9 || setup == 11
            % Microscope off
            cgflip(grey,grey,grey)
            dasbit(2, 0);
            recofftime = toc(start);
            
            pause(0.05);
            dasbit(0, 0);
        end
        
        %% Process response
        lickSide = 'unknown';
        if strcmp(Reaction,'1') %right
            lickSide = 'right';
            Reaction = 'Hit';
            RewardCount = RewardCount +1;
            ErrorTimeOut = 0;
        elseif strcmp(Reaction, '2') %left
            Reaction = 'Hit';
            lickSide = 'left';
            RewardCount = RewardCount +1;
            ErrorTimeOut = 0;
        elseif strcmp(Reaction,'0')
            Reaction = 'Error';
            if get(ezbox, 'Value')
                RewardCount = RewardCount +0.5;
            end
            ErrorTimeOut = 2.5;
        elseif strcmp(Reaction,'4')
            Reaction = 'TooFast';
            ErrorTimeOut = 0;
        else
            Reaction = 'Miss';
            Miss = Miss + 1;
            RT = NaN;
            ErrorTimeOut = 0;
        end
        
        %% Update Hits/Errors in figure window 1
        if strcmp(side,'left') && strcmp(Reaction,'Hit')
            set(HitText(2),'string',num2str(str2double(get(HitText(2),'string'))+1))
        elseif strcmp(side,'right') && strcmp(Reaction,'Hit')
            set(HitText(1),'string',num2str(str2double(get(HitText(1),'string'))+1))
        elseif strcmp(side,'left') && strcmp(Reaction,'Error') 
            set(ErrText(2),'string',num2str(str2double(get(ErrText(2),'string'))+1))
        elseif strcmp(side,'right') && strcmp(Reaction,'Error') 
            set(ErrText(1),'string',num2str(str2double(get(ErrText(1),'string'))+1))
        end
        
        %% Update figure window 2
        switch ori
            case Par.Orientations(1)
                if strcmp(Reaction,'Hit') && opto ==0
                    Perf0(idx0) = 1;
                elseif strcmp(Reaction,'Error') || (strcmp(lickSide,'unknown') && TooEarly) && opto ==0
                    Perf0(idx0) = 0;
                end
                PerfMat(1,Delayidx) = nansum(Perf0)./sum(~isnan(Perf0)).*100;
                idx0 = mod(idx0,nback)+1;
                
            case Par.Orientations(2)
                if strcmp(Reaction,'Hit') && opto ==0
                    Perf90(idx90) = 1;
                elseif strcmp(Reaction,'Error') || (strcmp(lickSide,'unknown') && TooEarly) && opto ==0
                    Perf90(idx90) = 0;
                end
                PerfMat(2,Delayidx) =  nansum(Perf90)./sum(~isnan(Perf90)).*100;
                idx90 = mod(idx90,nback) + 1;
        end
        
        %% if 'Too Early'
        if TooEarly
            rtemp = nan;
            if strcmp(Reaction,'Hit')
                if punishEarlyLicking
                    RewardCount = RewardCount - 0.5; %Only half of the reward received
                end
                rtemp = 1;
            elseif strcmp(Reaction,'Error')
                rtemp = 0;
            end
            if strcmp(lickSide,'unknown')
                rtemp = 0;
            end
            Reaction = 'Too Early';
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
            LOG.Setup = setup;
            LOG.Type = 'OriDiscrWM';
            LOG.Mouse = mouse;
            LOG.Date = date;
            LOG.ITI = Par.ITI;
            if  ~exist('comtime','var')
                comtime=nan;
            end
            LOG.comtime(trialcount) = comtime;
            LOG.aborttrial(trialcount) = aborttrial;
            LOG.earlylicks(trialcount) = earlylicks;
            LOG.RealTask(trialcount) = 1; %Real task or no visual sitmulus?
            LOG.shuffleval(trialcount) = shuffleval;
            LOG.Trial(trialcount) = trialcount;
            LOG.Orientation(trialcount) = CENT.Orient;
            LOG.Opto(trialcount) = opto;
            if ~isnan(optoOnTime)
                LOG.OptoOnTime(trialcount) = toc(optoOnTime)-toc(stimonset);
            else
                LOG.OptoOnTime(trialcount) = nan;
            end
            if ~isnan(optoOFFTime)
                LOG.OptoOffTime(trialcount) = toc(optoOFFTime)-toc(stimonset);
            else
                LOG.OptoOffTime(trialcount)=nan;
            end
            LOG.Side(trialcount) = {side};
            LOG.RepError(trialcount) = RepError;
            
            if ischar(RT)
                LOG.RT(trialcount) = str2double(RT);
            else
                LOG.RT(trialcount) = RT;
            end
            LOG.FigContrast(trialcount) = CENT.Contrast;
            LOG.Reaction(trialcount) = {Reaction};
            LOG.Flash(trialcount) = flashstim;
            LOG.EarlyLickP(trialcount) = nan;
            LOG.GracePeriod(trialcount) = str2double(get(gpbox, 'string'));
            LOG.OptoAllowed(trialcount) = OptoAllowed;
            LOG.OptoDurationPlanned(trialcount) = optoduration;
            LOG.Itiduration(trialcount) = itiduration;
            LOG.Gavepassive(trialcount) = gavepassive;
            LOG.Passivefirst(trialcount) = passfirst;
            LOG.PassivPerc(trialcount) = CurPassivPerc;
            LOG.Ezbox(trialcount) = ezvalue;
            LOG.Passivedelay(trialcount) = pdelay;
            LOG.Drum(trialcount) = drumm;
            LOG.SemiDrumm(trialcount) = semidrumm;
            LOG.StimPosX(trialcount) = spatialoffset; %degrees
            LOG.StimPosY(trialcount) = str2double(get(ydispbox,'String'));%degrees
            LOG.FullScreen(trialcount) = fullscrfig;
            LOG.MovingGrating(trialcount) = MovingGratings;
            LOG.orientationclearness(trialcount) = orientationclearness;
            LOG.propl(trialcount) = propR;
            LOG.RadioChannel = nan;
            % Reward times
            LOG.RewardTime(trialcount) = Par.RewardTime;
            LOG.Stimdur(trialcount) = stimdur;
            LOG.currentdelay(trialcount) = curdelay;
            if ~exist('startRewardPeriod','var')
                startRewardPeriod=nan;
            end
            LOG.ActualRewPeriod(trialcount) = startRewardPeriod;
            LOG.CurrentPhase(trialcount) = CurPhase;
            LOG.leftside = leftside;
            LOG.spatialoffset = spatialoffset;
            LOG.RewardCount(trialcount) = RewardCount;
            LOG.memory(trialcount) = memory;
            LOG.rdurright(trialcount) = rdurR;
            LOG.rdurleft(trialcount) = rdurL;
            LOG.ValveOpenTime(trialcount) = ValveOpenTime;
            LOG.memoryshuffle(trialcount) = memoryshuffle;
            LOG.TotalCleanbasetime(trialcount) = TotalCleanbasetime;
            disp(['Valve open for ' num2str(ValveOpenTime)])
            
            
            %%  Print trial results
            if strcmp(Reaction, 'Error')
                cprintf( 'Red', [ 'Trial: ' num2str(trialcount) ,' Ori: ' num2str(CENT.Orient), ' Opto: ' num2str(opto), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
            elseif strcmp(Reaction, 'Hit')
                cprintf( [0, 0.5, 0], [ 'Trial: ' num2str(trialcount) ,' Ori: ' num2str(CENT.Orient), ' Opto: ' num2str(opto), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
            elseif  strcmp(Reaction, 'Too Early') && (strcmp(lickSide, 'left') || strcmp(lickSide, 'right'))
                cprintf( [0, 0, 1],[ 'Trial: ' num2str(trialcount) ,' Ori: ' num2str(CENT.Orient), ' Opto: ' num2str(opto), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
            elseif strcmp(Reaction,'Too Early') && (strcmp(lickSide,'Unknown'))
                cprintf('SystemCommands',['Trial ' num2str(trialcount), ' Ori: ' num2str(CENT.Orient), ' Side: ', side ,', RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string') '\n']);
            else
                disp([ 'Trial: ' num2str(trialcount) , ' Side: ', side ,', Orientation: ' ,   num2str(CENT.Orient), ' Opto: ' num2str(opto), ' RT: ' num2str( LOG.RT(trialcount)) ', Reaction: ' Reaction , ', on ' lickSide, ' ,Passive: ', num2str(gavepassive), ' ,Passfirst: ', num2str(passfirst), ' pdelay: ', get(pdbox, 'string')]);
            end
            
            %%  Update figure and performance
            if strcmp(Reaction,'Hit')
                AllPerfVec(trialcount) = 1;
                TooEarlyVec(trialcount) = 1;
            elseif strcmp(Reaction,'Error')
                AllPerfVec(trialcount) = 0;
                TooEarlyVec(trialcount) = 1;
            elseif strcmp(Reaction,'Too Early')
                TooEarlyVec(trialcount) = 0;
                AllPerfVec(trialcount) = rtemp;
                %         elseif strcmp(Reaction,'Too Early')
                %             AllPerfVec(trialcount) = -1;
            end
            Delayvec = [Delayvec curdelay];
            Delayidx = round((curdelay-Par.minDelay)./(Par.DelayStep/2))+1; %ms
            if trialcount > 1 && Delayvec(trialcount) ~= Delayvec(trialcount-1)
                Perf0 = nan(1,nback);
                Perf90 = nan(1,nback);
                idx0 = 1;
                idx90 = 1;
                ThisContrTrialcount = 1;
            else
                ThisContrTrialcount = ThisContrTrialcount+1;
            end
            
            %% Calculate total performance and biasindex
            if sum(~isnan(AllPerfVec)) >= nback
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
                
                %%  Phases
                %Phase 1: decrease % passive licks on full stimdur, phase 2:
                %increase delay
                if ~freeze && sum(strcmp(LOG.Reaction,'Hit'))>20
                    switch CurPhase
                        case 1 %Passive Percentage
                            if nansum(LastPerf)/nback > accuracyThrs/100 && strcmp(Reaction,'Hit') && PerfMat(1,Delayidx)>accuracyThrs && PerfMat(2,Delayidx) > accuracyThrs%&& nback/length(nbackidx) > 0.20
                                disp(['Performance over ' num2str(accuracyThrs) '% and hit plus enough responses and no bias! Decreasing number of passives and increasing passive delay'])
                                CurPassivPerc = CurPassivPerc - 5;
                                set(PercPassBox, 'string',num2str(CurPassivPerc))
                                changepassives = 1;
                                
                                if CurPassivPerc < 0
                                    CurPassivPerc = 0;
                                    set(PercPassBox, 'string',num2str(CurPassivPerc))
                                    CurPhase = CurPhase + 1;
                                    set(PhaseBox,'String',num2str(CurPhase));
                                    disp('Entering Phase 2: no more passives, decreasing extra stim duration')
                                    NewPhase = 1;
                                end
                            elseif (nansum(LastPerf)/nback < (accuracyThrs-10)/100 && strcmp(Reaction,'Error'))
                                disp(['Performance below ' num2str((accuracyThrs-10)) '% and error. Increasing number of passives and decreasing passive delay'])
                                CurPassivPerc = CurPassivPerc + 2;
                                if CurPassivPerc > 75
                                    CurPassivPerc = 75;
                                end
                                set(PercPassBox,'String',num2str(CurPassivPerc)) %Set passive percentage
                                changepassives = 1;
                            else
                                %                                                             disp('No change in number of passives/delay')
                            end
                            
                        case 2 %Start increase delay
                            %fprintf('Ratio early versus normal respond: %.2f\n',nansum(TooEarlyVec(trialcount-nbackidx+1))/nback) % Shows ratio between early licks & normal licks
                            if nansum(LastPerf)/nback > accuracyThrs/100 && strcmp(Reaction,'Hit') && PerfMat(1,Delayidx)>accuracyThrs && PerfMat(2,Delayidx) > accuracyThrs
                                disp(['Performance over ' num2str(accuracyThrs) '% and hit and no bias! Increasing delay...'])
                                curdelay = curdelay + 50;
                                if curdelay > Par.maxDelay
                                    curdelay = Par.maxDelay;
                                    CurPhase = CurPhase + 1;
                                    set(PhaseBox,'String',num2str(CurPhase));
                                    set(DelayBox,'string',curdelay)
                                    disp('Full WM Reached')
                                    NewPhase = 1;
                                    
                                end
                            elseif nansum(LastPerf)/nback < (accuracyThrs-20)/100 && strcmp(Reaction,'Error')
                                disp(['Performance below ' num2str((accuracyThrs-20)) '% and error. Decrease delay...'])
                                curdelay = curdelay - 20;
                                if curdelay < Par.minDelay
                                    curdelay = Par.minDelay;
                                end
                            else
                                %                             disp('No change in delay')
                            end
                            
                            set(DelayBox,'string',curdelay)
                        case 3
                            if setup>=10 && ~OptoAllowed && nansum(LastPerf)/nback > accuracyThrs/100 && strcmp(Reaction,'Hit') && PerfMat(1,Delayidx)>accuracyThrs && PerfMat(2,Delayidx) > accuracyThrs
                                disp(['Performance over ' num2str(accuracyThrs) '% and hit and no bias! Opto On...'])
                                set(OptoBox,'value',1);
                                NewTrialMat = 1;
                                set(semidrummbox,'value',0); %drumming off
                                set(drummbox,'value',0); % in case immediate drumming on; turn also off
                                %                                 memoryshuffle = 0; %turn on shuffle
                                %                                 set(memoryBox,'value',0);
                                
                            elseif setup >= 10 && OptoAllowed && nansum(LastPerf)/nback < (accuracyThrs-20)/100 && strcmp(Reaction,'Error')
                                disp(['Performance below ' num2str((accuracyThrs-10)) '% and error. Opto Off...'])
                                set(OptoBox,'value',0);
                                NewTrialMat = 1;
                                set(semidrummbox,'value',1); %drumming on
                                %                                 memoryshuffle = 1;
                                %                                 set(memoryBox,'value',1);
                            end
                    end
                end
            end
        end
        %% Update performance figure
        if trialcount >= nback
%             figure(HP)
            if exist('hperf')
                delete(hperf)
            end
            hperf = bar(HAX,DelayOptions4Plot,PerfMat',0.8);
            set(HAX, 'XTick', [-500:500:2000])
            set(HAX,'XTickLabel',[-500:500:2000])
            %             hold on
            %             hline = plot(HAX,DelayOptions4Plot,LastPerf);
            if leftside == 1
                legend(hperf,['Orientation ' num2str(Par.Orientations(1)) ' , left'],['Orientation '  num2str(Par.Orientations(2)) ' , right'])
            else
                legend(hperf,['Orientation ' num2str(Par.Orientations(1)) ' , right'],['Orientation ' num2str(Par.Orientations(2)) ' , left'])
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
        
        if  strcmp(Reaction, 'Error')
            RepError = RepError+1;
        end
        %% Update trial structure
        if strcmp(Reaction, 'Hit') || (~drumm && ~semidrumm)
            trials(1,:) = [];
            GivePercTrial(1) = [];
            RepError = 0;
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
        
        if NewTrialMat || isempty(trials) || (trialcount > 1 &&  LOG.shuffleval(trialcount) ~= LOG.shuffleval(trialcount-1)  || (trialcount>1 && LOG.propl(trialcount) ~= LOG.propl(trialcount-1)) || (trialcount>1 && LOG.memoryshuffle(trialcount) ~= LOG.memoryshuffle(trialcount-1)))
            NewTrialMat = 0;
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
            
            clear details
            z = 1;
            if memoryshuffle
                for m = [0,1,1,1] %25% of trials lower memory
                    for p = optotrialvec %Optogenetics
                        for o = oriopt  %orientations
                            
                            details(z:z+lraltern-1,:) = repmat([o,p,m],[lraltern,1]);
                            z = z+lraltern;
                        end
                    end
                end
            else
                for p = optotrialvec
                    for o =oriopt  %orientations
                        details(z:z+lraltern-1,:) = repmat([o,p],[lraltern,1]);
                        z = z+lraltern;
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
                OptoflagOff = 0;
                tmpdetails = details;
                while (size(trials,1))<ntrials && ~OptoflagOff
                    order = randsample(oriopt,4,0);
                    for i = 1:length(order)
                        vec = find(tmpdetails(:,1) == order(i));
                        if length(vec)>0
                            pick = randsample(vec,1);
                        else
                            %More from 1 condition --> Probably because of propl.
                            %Put these extra oin random locations
                            OptoflagOff=1;
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
            
            disp('new trial struct')
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
                save(fullfile(vcTrainingLogs, nts),'LOG');
                if MetaData
                    savejson('', json, fullfile(vcTrainingLogs,[expname '_' fields.date '_session.json']));
                end
            catch
                save(fullfile(LocalFolder, nts),'LOG');
                if MetaData
                    savejson('', json, fullfile(LocalFolder,[expname '_' fields.date '_session.json']));
                end
            end
            
            if  memory == 0
                curdelay = savecurdelay;
            end
            % Save CurrentPhase settings
            try
                save(fullfile(vcTrainingLogs, [mouse '_CurrentPhase.mat']),'CurPhase','leftside','Par','CurPassivPerc','ezvalue','spatialoffset','flashstim','shuffleval','curdelay','fullscrfig','MovingGratings','freeze','orientationclearness')
            catch
                save(fullfile(LocalFolder, [mouse '_CurrentPhase.mat']),'CurPhase','leftside','Par','CurPassivPerc','ezvalue','spatialoffset','flashstim','shuffleval','curdelay','fullscrfig','MovingGratings','freeze','orientationclearness')
                if  trialcount == 1
                    disp('Warning, could not save on network. Saved in local directory')
                end
            end
        end
        
        %Check for newphase
        if CurPhase ~= str2double(get(PhaseBox,'String'));
            NewPhase = 1;
        end
        
        %% Check ESC in order to Pause
        
        if ESC == 1
            
            PauseButton = questdlg('Want to Pause or to Exit?','Pause Menu','Pause','Exit','Pause');
            
            switch PauseButton
                case 'Pause'
                    
                    %Put lickspout in lickposition
                    if  servoavailable
                        sendardwm(sport, ['IO ' num2str(LickPos)])
                    end
                    
                    ContinueFigure = figure;
                    ContinueButton = uicontrol('Position',[20 20 200 40],'String','Continue',...
                        'Callback','uiresume(gcbf)');
                    uiwait(gcf);
                    close(ContinueFigure);
                    
                    %Put lickspout in breakpos
                    if  servoavailable
                        sendardwm(sport, ['IO ' num2str(BreakPos)])
                    end
                    
                    ESC = 0;
                    
                case 'Exit'
                    ESC = 1;
                    
            end
        end
        
    end
end

%% Save figure
if savefig
    try
        saveas(HP,fullfile(vcTrainingLogs,[nts '_Perf']), 'bmp')
        saveas(HP,fullfile(vcTrainingLogs,[nts '_Perf']),'fig')
    catch
        saveas(HP,fullfile(LocalFolder,[nts '_Perf']), 'bmp')
        saveas(HP,fullfile(LocalFolder,[nts '_Perf']),'fig')
    end
end

try
    save(fullfile(LocalFolder, [mouse '_CurrentPhase.mat']),'CurPhase','leftside','Par','CurPassivPerc','ezvalue','spatialoffset','flashstim','shuffleval','curdelay','fullscrfig','MovingGratings','freeze','orientationclearness')
    save(fullfile(vcTrainingLogs, [mouse '_CurrentPhase.mat']),'CurPhase','leftside','Par','CurPassivPerc','ezvalue','spatialoffset','flashstim','shuffleval','curdelay','fullscrfig','MovingGratings','freeze','orientationclearness')
catch
    disp('Warning, could not save on network. Saved in local directory') 
end
%% Show Rewards given:
Totalhits = str2double(get(HitText(1),'string'))+str2double(get(HitText(2),'string'));
Totalerrors = str2double(get(ErrText(1),'string'))+str2double(get(ErrText(2),'string'));
fprintf('Total hits: %g\n',Totalhits)
fprintf('Total errors: %g\n',Totalerrors)
fprintf('Reward count for this session is %.1f\n',RewardCount)

%% Close it all
try
    if setup==9 || setup == 11
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

