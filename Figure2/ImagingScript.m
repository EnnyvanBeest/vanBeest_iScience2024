%% Passive Imaging Opto Script
% Enny van Beest, January 2019. Script made for passive Imaging + Opto.

%% Manual settings
arduino = 1; %Arduino connected?
MetaData = 1; %Metafile (Jason) saving?
debug = 0;%Debug mode (i.e. Destkop computer is cogent screen)
MovingGratings = 1; %if 1; stimuli move. If 0; stimuli don't move
rng('shuffle'); % avoid pseudorandomization
MaxTrials = 120 %15 trial repetitions, 2 sides, 2 orientations, 2 opto (no and yes)
setup = 9;
%% Paths + Current Directory
user_name =getenv('USERNAME');
tmpname = cd;
if ~isempty(strfind(tmpname,[user_name '.HERSEN']))
    user_name = regexp(tmpname,'\','split');
    user_name = user_name{3};
end
LocalFolder = mfilename('fullpath');
WMpath = strsplit(LocalFolder,'PassiveImagingOpto');
DomainFolder = fullfile('C:\Users\',user_name);
WMpath = fullfile(WMpath{1});
vcLogs = '\\vs02\VandC\WorkingMemory_EB\PassiveImgOpto\';
DasPath = strsplit(WMpath,'Mouse');
DasPath = fullfile(DasPath{1},'Mouse','DIO24');
addpath(genpath(DasPath))
try
    addpath(genpath(WMpath))
    cd(WMpath)
catch ME
    if debug
        WMpath = 'I:\GitHub\Mouse\WM';
        addpath(genpath(WMpath))
        cd(WMpath)
        addpath(genpath('C:\Users\beest\Documents\MATLAB\CogGph'))
    else
        disp(ME)
        WMpath = ['C:\Users\' user_name '\Documents\Github\Mouse\WM\'];
        addpath(genpath(WMpath))
        cd(WMpath)
    end
end
if MetaData
    addpath(genpath(fullfile(['C:\Users\' user_name '\Documents\Github'],'OpenAccessStorage')))
    addpath('Z:\OpenAccessStorageParameters')
end


%% Upload arduino
if arduino
    disp('Sure you uploaded correct arduino script? Press key to continue or abort...')
    pause
end

% Start Serial connection with Arduino
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
        sport = serial('com4');
        if strcmp(get(sport, 'status'), 'closed'); fopen(sport); end;
    end
    set(sport, 'baudrate', 250000);
    set(sport, 'timeout', 0.5);
end

%% Settings for Set-up (always imaging + opto)
caminit = 1.2; %1.2 seconds minimally
camstoretime = 2; %Seconds it takes software to store images  %For exposure 50 ms, take minimally 2 seconds, for exposure 20ms, take 3.5!
run(fullfile(['C:\Users\' user_name '\Documents\Github'],'Mouse','imagingStim','runSettings.m'))
D2ScreenWidth = 51./2; %??
ITI = caminit+camstoretime;%1.5; %N.B. IF SHORTER THAN CAMINIT + CAMSTORETIME IT WILL BE MORE THAN THIS VALUE;
displayscript = 'correctDisplayWF'
gammaconversion = 'gammaconWF'

%% Session name (nts)
warning on backtrace
logon = 0;
prompt = {'Mouse Name', 'Exp Nr', 'Date'};
dlg_title = 'Please enter parameters';
num_lines = 1;
def = {'Name', '1', datestr(datenum(date), 'yyyymmdd'), '50', '1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if ~strcmp(answer{1}, 'Name')
    
    mouse = answer{1};
    expnum = answer{2};
    date = answer{3};
    expname = [mouse expnum];
    LOG.Stimulus = 'SimpleLick';
    logon = 1;
else
    mouse = 'Unknown';
    leftside = round(rand(1,1))+1; % if 1, left is 0 degrees, if 2, left is 90 degrees
    CurPhase = 1;
    ezvalue = 1;
    EarlyLickP = round((Par.maxDelay + Par.StimDurMin)); %EarlyLickP is the period in which the mouse is allowed to lick without it being punished
    
end
nts = [mouse '_' date '_B' expnum];
expname = [mouse num2str(expnum)];

while exist(fullfile(vcLogs,mouse,[mouse date],['B' expnum],[nts '.mat']), 'file')
    disp(['Log-file expnum ' num2str(expnum) ' already used, expnum number changes...']);
    expnum = expnum+1;
    nts = [mouse '_' datestr(date,'yyyymmdd') '_B' num2str(expnum)]
    expname = [mouse num2str(expnum)];
    
end
expname = [mouse expnum];

%% JSON FILE
% the metadata in this schema is stored in the database,
% but much more can be added to the json files
if MetaData
    try
        fields.project = 'WorkingMemory_EB';
        fields.dataset = 'PassiveImagingOpto';
        fields.investigator = 'EnnyvB';
        fields.subject = mouse;
        fields.setup = 'WF';
        fields.stimulus = 'SpatialWM';
        fields.condition = 'anesthetized';
        fields = getdbfields('VC',fields); %retrieves info from mysql tables using a GUI, fill in already.
        
    catch ME
        disp(ME)
        fields = getdbfields('VC'); %retrieves info from mysql tables using a GUI, fill in already.\
    end

    
    %not necessary if you are extremely consistent in filling the schema below
    %schema for structuring fields in json file
    %the php script that parses the json files requires that the data in the json file
    %is structured accordingly, to fill the database with records
    %
    json = fields;
    json.logfile = [expname '.mat']; %name for your logfile
    
end

%% initialize PAR- and other variable: contains general Settings.
global Par
Par.ItiRandMin = 0;
Par.ItiRandMax = 2;
Par.StimDur = 1000; %Stimulus on
Par.Orientations = [45 135]; %With Old mice : [0 90]),new mice [ 45 135];
Par.SideOpt = [1 2]; %1 = left, 2 = right
Par.SpatFreq = 0.08;%0.05; 12.5 degree sper cycle m ?
Par.FigSize = 35; % in visual degrees
Par.prefspeed = 24; %deg/s
spatialoffset = 40; %in visual degrees
yoffset = 15;
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

ESC = 0;
details = [];

%% Initialize  Cogent
cgloadlib

cgshut %shut cg if one's open
if debug
    warning('Debug Mode enabled')
    cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 0)
else
    
    try
        cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 1) %Second window?
    catch
        cgopen(Par.Screenx, Par.Screeny, 32,Par.Refresh , 0)
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

%% Ccreate visual stimuli
SafeStimPath = fullfile(DomainFolder,'WM_OriGratings');
if ~exist(SafeStimPath,'dir')
    mkdir(SafeStimPath)
end

% pixpercycle = Par.PixPerDeg * 1/Par.SpatFreq;
secpercycle = 1/Par.SpatFreq/prefspeed;
framespercycle = round(Par.Refresh * secpercycle);
phasesteps = 0:2*pi/framespercycle:2*pi-(2*pi/framespercycle);
disp('Preparing stimuli...')

%% Final settings
ititime = tic;
start = tic;
left = 'left';
right = 'right';

dasinit(22);
optoport = 4;
dasbit(optoport,0);
optotrialvec = [0,1]; %0 = off, 1 = before trial onset,
optoOn = 0;
OPTOREG = 1;
lraltern = 2;
figure(1) %show control-window
oriopt = [Par.Orientations(2),Par.Orientations(2),Par.Orientations(1),Par.Orientations(1)];
sideopt = [Par.SideOpt(1),Par.SideOpt(2),Par.SideOpt(1),Par.SideOpt(2)];


%% Make trial matrix: 'details'
clear details
z = 1;
for p = optotrialvec
    for o = 1:length(oriopt)  %orientations
        details(z,:) = [oriopt(o),sideopt(o),p];
        z = z+1;
    end
end

%% Shuffling: create (random) 'order' variable
[ntrials, nvars] = size(details);

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

maxcount =countrepeats(trials(:,1));
if maxcount>4
    warning('More than 4 trials of the same side after each other in trial struct')
end
trialcount = 0;

dasbit(1,1) %Shutter open

start = tic;
%% Trial Loops
while ~ESC && trialcount<MaxTrials
    
    
    optoOnTime = nan;
    optoOFFTime = nan;
    
    %% Update Arduino settings
    ori = trials(1,1);                                 %Which orientation?
    sidenr = trials(1,2); %Which side?
    opto = trials(1,3); % opto?
    trialcount = trialcount + 1;                        %Trial number
    
    % Reset:
    didflip = 0;
    caminitinITI = 0;
    
    %%  Define which side should the stimulus be on for this orientation
    if sidenr==1
        side = left;
    else
        side = right;
    end
    
    %%  Spatial offset (x/y)
    if strcmp(side,'right')
        x = Par.PixPerDeg*spatialoffset;
    elseif strcmp(side,'left')
        x =  -spatialoffset*Par.PixPerDeg;
    end
    y = Par.PixPerDeg*yoffset;
    
    
    %%  Create (circle) figure (if enabled)
    if fullscrfig
        matb1 = repmat(grey, Par.Screeny, Par.Screenx);
        matb2 = repmat(grey, Par.Screeny, Par.Screenx);
        circle = logical(correctDisplayWF(selectcircle(160, Par.FigSize * Par.PixPerDeg,x, y), 1));
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
            if fullscrfig
                cgsetsprite(phaseidx)
                cgdrawsprite(50,0,0)
                cgsetsprite(0)
            end
        end
    else%still image
        phasestep1 = randsample(1:length(phasesteps),1);
        cgloadbmp(1, [SafeStimPath '\' num2str(ori) '_' num2str(phasesteps(phasestep1)) '.bmp'])
        if fullscrfig
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
    randomITI = Par.ItiRandMax * rand;
    while toc(ititime) < Par.ITI + randomITI - 0.500 %500ms opto on before hand
        cgflip(grey,grey,grey)
    end
    
    if ~debug
        dasbit(0, 1) %Initialises camera
        pause(caminit)
        dasbit(2, 1) %Starts image acquisition
    end
    %Time of acquisition onset
    recordtime = toc(start);
    
    % Wait 500ms for baseline data
    pause(0.5)
    
    %Optotrial where opto start before trial?
    if opto==1 && ~optoOn
        disp('Sending TTL for Opto ON')
        dasbit(optoport,1);
        optoOnTime = tic;
        optoOn = 1;
    end
    
    OptoITItime = tic;
    while toc(OptoITItime) < 0.500
        cgflip(grey,grey,grey)
    end
    itiduration = toc(ititime);
    
    %%  Stimulus Presentation
    dasbit(3, 1) %Send message of stimulus
    
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
    
    while (toc(stimonset) < Par.StimDur/1000)
        %Move gratings
        if MovingGratings && ~didflip
            cgdrawsprite(f,0,0)
            cgflip(grey,grey,grey)
            f = f+1;
            if f > length(phasesteps)
                f = 1;
            end
        end
    end
    
    %% End stimulus, start ITI
    
    %Turn off stimulus
    if ~didflip
        cgflip(grey,grey,grey)
        didflip = 1;
    end
    
    if ~debug
        dasbit(3,0)%Reset stimulus bit
        dasbit(2,0)%Turns off image acquisition
        %Store images Here
        dasbit(0,0) %Deinitalizes and stores, this takes some time.
    end
    
    ititime = tic;
    % Opto 500ms longer than stimulus
    while ititime<0.500
        cgflip(grey,grey,grey)
    end
    
    %% END OPTO
    if optoOn
        disp('Sending TTL for Opto OFF')
        dasbit(optoport,0);
        optoOFFTime = tic;
        optoOn=0;
    end
    
    %Need more time to store?
    while ititime<(0.5-camstoretime)
        cgflip(grey,grey,grey)
    end
    
    %% Finish session?? (read keyboard input)
    [kd,kp] = cgkeymap;
    if length(find(kp)) == 1
        if find(kp) == 1;
            ESC = 1;
        end
    end
    
    %% shared LOG data
    LOG.Setup = 'WF';
    LOG.Type = 'PassiveImagingOpto';
    LOG.Mouse = mouse;
    LOG.Date = date;
    LOG.ITI = Par.ITI;
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
    LOG.FigContrast(trialcount) = CENT.Contrast;
    LOG.Itiduration(trialcount) = itiduration;
    LOG.OptoRegulator(trialcount) = OPTOREG;
    LOG.StimPosX(trialcount) = spatialoffset; %degrees
    LOG.StimPosY(trialcount) = str2double(yoffset);%degrees
    LOG.FullScreen(trialcount) = fullscrfig;
    LOG.MovingGrating(trialcount) = MovingGratings;
    
    disp([ 'Trial: ' num2str(trialcount) , ' Side: ', side ,', Orientation: ' , num2str(CENT.Orient), ' Opto: ' num2str(opto)]);
    
    %% Update Trial Struct
    trials(1,:) = [];
    
    if  trialcount>1 && LOG.OptoRegulator(trialcount) ~= LOG.OptoRegulator(trialcount-1)
        trials = [];
        if  LOG.OptoRegulator(trialcount) ==1
            optotrialvec = [0,1];
        else
            optotrialvec = [0];
        end
    end
    if isempty(trials)
        clear details
        z = 1;
        for p = optotrialvec
            for o = 1:length(oriopt)  %orientations
                details(z,:) = [oriopt(o),sideopt(o),p];
                z = z+1;
            end
        end
        % Shuffling: create (random) 'order' variable
        [ntrials, nvars] = size(details);
        
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
        
        maxcount =countrepeats(trials(:,1));
        if maxcount>4
            warning('More than 4 trials of the same side after each other in trial struct')
        end
    end
    %% Save logs
    if logon
        try
            if ~isdir(fullfile(vcLogs,mouse,[mouse date],['B' expnum]))
                mkdir(fullfile(vcLogs,mouse,[mouse date],['B' expnum]))
            end
            save(fullfile(vcLogs,mouse,[mouse date],['B' expnum], nts),'LOG');
            if MetaData
                savejson('', json, fullfile(vcLogs,mouse,[mouse date],['B' expnum],[expname '_' fields.date '_session.json']));
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
                ContinueFigure = figure;
                ContinueButton = uicontrol('Position',[20 20 200 40],'String','Continue',...
                    'Callback','uiresume(gcbf)');
                uiwait(gcf);
                close(ContinueFigure);
                
                ESC = 0;
                
            case 'Exit'
                ESC = 1;
                
        end
    end
end

%% Close it all
try
    dasbit(1, 0) %Closes shutter
    dasbit(optoport,0)
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
