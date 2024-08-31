ImagePath = '\Path\To\Zenodo\Data'
allen_atlas_path = '\Path\To\GitHub\allenCCF\';
Redo = 0;


AllMice2Use = {'Becquerel','Berlijn','Cherry','China','Chopin','Faraday','Galilei','Haydn',...
   'Holland','HongKong','Jemison','Jenkins','Karlsen','Kirchoff','Leeuwenhoek','Lully','Mahler','Sebald'};
Genotype =  {'D2','D2','D1','D1','D2','D2','D2','D1','D1','D2','D2','D1','D1','D2','D2','D1','D1','D1'}
ChannelPerMouse = [2,1,2,2,1,1,1,1,2,2,1,1,1,2,2,1,1,1];
%% Add toolbox directories to path - you need thse toolboxes below
addpath(genpath('\Path\To\GitHub\npy-matlab')) 
addpath(genpath('\Path\To\GitHub\AP_histology'))
% addpath(genpath('\Path\To\BrewerMap'))

%% Add data path
addpath(genpath(ImagePath))
addpath(genpath(allen_atlas_path))

%% Load CCF and set paths for slide and slice images
% Load CCF atlas
tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']);
av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']);
st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']);

%% Loop over individual mice
AllFolders = dir(ImagePath);
AllFolders(~[AllFolders(:).isdir]) = [];

%% Manual part
for fid = 1:length(AllFolders)
    if ~ismember(AllFolders(fid).name,AllMice2Use)
        continue
    end
    
    % Set paths for histology images and directory to save slice/alignment
    im_path = fullfile(AllFolders(fid).folder,AllFolders(fid).name);
    slice_path = [im_path filesep 'slices'];
    
    if exist(fullfile(slice_path,'histology_ccf.mat')) && ~Redo
        disp('Already aligned manually... skip')
        continue
    end
    
    % Preprocess slide images to produce slice images
    % Set white balance and resize slide images, extract slice images
    % (Note: this resizes the images purely for file size reasons - the CCF can
    % be aligned to histology no matter what the scaling. If pixel size is
    % available in metadata then automatically scales to CCF resolution,
    % otherwise user can specify th e resize factor as a second argument)
    
    % Set resize factor
    resize_factor = 0.1; % (slides ome.tiff: auto-resize ~CCF size 10um/px)
    %     resize_factor = 1; % (slides tiff: resize factor)
    
    % Set slide or slice images
    % slice_images = false; % (images are slides - extract individual slices)
    slice_images = true; % (images are already individual slices)
    
    % Preprocess images
    AP_process_histology(im_path,resize_factor,slice_images);
    
    % (optional) Rotate, center, pad, flip slice images
    AP_rotate_histology(slice_path);
    pause()
    
    % Align CCF to slices
    % Find CCF slices corresponding to each histology slice
    AP_grab_histology_ccf(tv,av,st,slice_path);
    pause()
    
end

%% Automatic part
for fid = 1:length(AllFolders)
    if ~ismember(AllFolders(fid).name,AllMice2Use)
        continue
    end
    
    % Set paths for histology images and directory to save slice/alignment
    im_path = fullfile(AllFolders(fid).folder,AllFolders(fid).name);
    slice_path = [im_path filesep 'slices'];
    %
    % Align CCF slices and histology slices
    % (first: automatically, by outline)
    AP_auto_align_histology_ccf(slice_path);
end

%% continue manual part
for fid = 1:length(AllFolders)
    if ~ismember(AllFolders(fid).name,AllMice2Use)
        continue
    end
    
    % Set paths for histology images and directory to save slice/alignment
    im_path = fullfile(AllFolders(fid).folder,AllFolders(fid).name);
    slice_path = [im_path filesep 'slices'];
    %
    % (second: curate manually)
    AP_manual_align_histology_ccf(tv,av,st,slice_path);
    pause()
end

%% Utilize aligned CCF
for fid = 1:length(AllFolders)
    if ~ismember(AllFolders(fid).name,AllMice2Use)
        continue
    end
    % Set paths for histology images and directory to save slice/alignment
    im_path = fullfile(AllFolders(fid).folder,AllFolders(fid).name);
    slice_path = [im_path filesep 'slices'];
    %
    
    
    %     % Display histology within 3D CCF
    %     AP_view_aligned_histology_volume(tv,av,st,slice_path,channel);
    %     pause()
    %
    % Get probe trajectory from histology, convert to CCF coordinates
    AP_get_probe_histology(tv,av,st,slice_path);
    pause()
end

%% Projection - from here for all mice
for fid = 1:length(AllFolders)
    if ~ismember(AllFolders(fid).name,AllMice2Use)
        continue
    end
    AllFolders(fid).name
    % Set paths for histology images and directory to save slice/alignment
    im_path = fullfile(AllFolders(fid).folder,AllFolders(fid).name);
    slice_path = [im_path filesep 'slices'];
    %
    % Display aligned CCF over histology slices
    AP_view_aligned_histology(st,slice_path);
    
    channel = input('Channel: 1=Red, 2=Green, 3=Blue');
    HistPlusFiber
    ChannelPerMouse(fid) = channel;
end

%% 2D projection across mice
GenOpt = unique(Genotype);
figure;

MiceIncl = {};
for genid=1:length(GenOpt)
    colors = jet(sum(ismember(Genotype,GenOpt{genid}))+1);
    clear allsliceallmice
    flag = 0;
    countmouse = 1;
    for fid = 1:length(AllFolders)
        if ~ismember(AllFolders(fid).name,AllMice2Use)
            continue
        end
        if ~strcmp(Genotype(find(ismember(AllMice2Use,AllFolders(fid).name))),GenOpt{genid})
            continue
        end
        clear allslice
        % Set paths for histology images and directory to save slice/alignment
        im_path = fullfile(AllFolders(fid).folder,AllFolders(fid).name);
        slice_path = [im_path filesep 'slices'];
        channel = ChannelPerMouse(find(ismember(AllMice2Use,AllFolders(fid).name)));
        
        gui_data = struct;
        gui_data.st = st;
        
        % Load in slice images
        gui_data.slice_im_path = slice_path;
        slice_im_dir = dir([slice_path filesep '*.tif']);
        slice_im_fn = natsortfiles(cellfun(@(path,fn) [path filesep fn], ...
            {slice_im_dir.folder},{slice_im_dir.name},'uni',false));
        gui_data.slice_im = cell(length(slice_im_fn),1);
        for curr_slice = 1:length(slice_im_fn)
            gui_data.slice_im{curr_slice} = imread(slice_im_fn{curr_slice});
        end
        
        % Load corresponding CCF slices
        ccf_slice_fn = [slice_path filesep 'histology_ccf.mat'];
        load(ccf_slice_fn);
        gui_data.histology_ccf = histology_ccf;
        
        % Load histology/CCF alignment
        ccf_alignment_fn = [slice_path filesep 'atlas2histology_tform.mat'];
        load(ccf_alignment_fn);
        gui_data.histology_ccf_alignment = atlas2histology_tform;
        
        % Warp area labels by histology alignment
        gui_data.histology_aligned_av_slices = cell(length(gui_data.slice_im),1);
        for curr_slice = 1:length(gui_data.slice_im)
            curr_av_slice = gui_data.histology_ccf(curr_slice).av_slices;
            curr_av_slice(isnan(curr_av_slice)) = 1;
            curr_slice_im = gui_data.slice_im{curr_slice};
            
            tform = affine2d;
            tform.T = gui_data.histology_ccf_alignment{curr_slice};
            tform_size = imref2d([size(curr_slice_im,1),size(curr_slice_im,2)]);
            gui_data.histology_aligned_av_slices{curr_slice} = ...
                imwarp(curr_av_slice,tform,'nearest','OutputView',tform_size);
        end
        
        
        % Warp histology to CCF
        gui_data.atlas_aligned_histology = cell(length(gui_data.slice_im),1);
        for curr_slice = 1:length(gui_data.slice_im)
            curr_av_slice = gui_data.histology_ccf(curr_slice).av_slices;
            curr_av_slice(isnan(curr_av_slice)) = 1;
            curr_slice_im = gui_data.slice_im{curr_slice};
            
            tform = affine2d;
            tform.T = gui_data.histology_ccf_alignment{curr_slice};
            % (transform is CCF -> histology, invert for other direction)
            tform = invert(tform);
            
            tform_size = imref2d([size(gui_data.histology_ccf(curr_slice).av_slices,1), ...
                size(gui_data.histology_ccf(curr_slice).av_slices,2)]);
            
            gui_data.atlas_aligned_histology{curr_slice} = ...
                imwarp(curr_slice_im,tform,'nearest','OutputView',tform_size);
            
        end
        
        
        tmp = cellfun(@(X) double(X(:,:,channel)),gui_data.atlas_aligned_histology,'UniformOutput',0);
        tmp = cat(3,tmp{:});
        value_thresh = double(quantile(tmp(:),0.99));%minimum value)
        
        for curr_slice = 1:length(gui_data.slice_im)
            
            % Get thresholded image
            curr_slice_im = double(gui_data.atlas_aligned_histology{curr_slice}(:,:,channel));%./...
            nanmean(nanmean(gui_data.atlas_aligned_histology{curr_slice}(:,:,channel),2),1);
            slice_alpha = zeros(size(curr_slice_im));
            if ~exist('allsliceallmice')
                allsliceallmice = zeros(size(curr_slice_im));
                dv = gui_data.histology_ccf(curr_slice).plane_dv;
                ml = gui_data.histology_ccf(curr_slice).plane_ml;
            end
            if ~exist('allslice')
                allslice = zeros(size(curr_slice_im));
            end
            % Draw if thresholded pixels (ignore if not)
            if any(curr_slice_im(:) > value_thresh)
                tmpslice = curr_slice_im;%gui_data.atlas_aligned_histology{curr_slice}(:,:,channel);
                tmpslice(tmpslice<value_thresh)=0;
                
                tmpslice(tmpslice>=value_thresh) = countmouse;
                [rid,cid] = find(tmpslice~=0);
                % project in z
                for idx = 1:length(rid)
                    [ridx,cidx] = find(dv==rid(idx)&ml==cid(idx));
                    allslice(rid(idx),cid(idx)) = allslice(rid(idx),cid(idx))+1;
                end
            end
        end
        
        
        
        %         freezeColors
        load([gui_data.slice_im_path filesep 'probe_ccf.mat'])
        clear probe_lines
        for curr_probe = 1:length(probe_ccf)
            probe_lines{curr_probe} = probe_ccf.trajectory_coords;
        end
        probe_lines = cat(1,probe_lines{:});
        probe_lines = probe_lines(:,2:3);
        AllProbeFits{countmouse} = probe_lines;
        allsliceallmice(allslice>=1) = allsliceallmice(allslice>=1)+1;
        countmouse = countmouse+1;
        MiceIncl = {MiceIncl{:} AllFolders(fid).name};
    end
    
    subplot(1,2,genid)
    h=imagesc(allsliceallmice,[0 countmouse])
    
    set(h,'alphadata',allsliceallmice~=0)
    colormap hot
    freezeColors
    hold on
    
    
    for midx = 1:length(AllProbeFits)
        line(AllProbeFits{midx}(:,2),AllProbeFits{midx}(:,1), ...
            'color',colors(midx,:),'linewidth',2)
    end
    box off
    axis off
    axis square
    title(GenOpt{genid})
    colorbar
    
end

saveas(gcf,fullfile(ImagePath,['2DProjectionHistoWithFiber.fig']))
saveas(gcf,fullfile(ImagePath,['2DProjectionHistoWithFiber.bmp']))

%Get histology
curr_av_slice_warp = gui_data.histology_aligned_av_slices{5};
av_warp_boundaries = round(conv2(curr_av_slice_warp,ones(3)./9,'same')) ~= curr_av_slice_warp;

figure; imagesc((av_warp_boundaries.*-1)+1)
colormap gray