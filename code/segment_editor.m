function segment_editor_new()

    %%%%% Set up the starting folder %%%%%
    starting_dir = '../data/';
    if ~exist(starting_dir, 'dir')
        % If there is no data folder just start back one level
        cur_dir = pwd;
        idcs   = strfind(pwd,'/');
        starting_dir = cur_dir(1:idcs(end)-1);
    end
    dir_data = what(starting_dir);
    starting_dir = dir_data.path;

    % Get the name of the file that the user wants to use.
    %defaultFileName = fullfile(startingFolder, '*.*');
    dname = uigetdir(starting_dir, 'Select a data directory');
    if dname == 0
      % User clicked the Cancel button.
      fprintf('Process terminated by user. No file was selected. \n')
      return;
    end

    % Select the data series and segmentation seed frame for the data
    [segs, seg_select] = initialize_dialogue(dname,starting_dir);

    % Check for the segmentation fields and either load them or begin segmentations from scratch
    working_seg = segs.(seg_select);

    % Field structure 
    % seg.'series'.bounds.lungs.left/right
    % seg.'series'.bounds.ext
    % seg.'series'.z_sep

    % If bounds do not exist segment the files 
    if ~isfield(working_seg,'bounds')
        % Segment the data from directory 
        [imgs,working_seg,dcm_info] = segment_dicom([dname, '/', seg_select], working_seg.initial_frame);

        % Downsample the segmentations 
        frames = 5; % TODO 2 above and 2 below the intercostal space
        for i=1:frames
            seg = working_seg{i};
         % Downsample to 20... 
            if isfield(seg,'exterior')
                o_seg{i}.exterior = seg.exterior(round(linspace(1,length(seg.exterior),20)),:);
                o_seg{i}.exterior(end,:) = [];
                %h(2) = plot(app.UIAxes,o_seg.exterior(:,2),o_seg.exterior(:,1),'-o','LineWidth',2,'Color',[217,95,2]/256);
            end 
            if isfield(seg, 'l_lung')
                o_seg{i}.l_lung = seg.l_lung(round(linspace(1,length(seg.l_lung),20)),:);
                o_seg{i}.l_lung(end,:) = [];
                %h(3) = plot(app.UIAxes,o_seg.l_lung(:,2),o_seg.l_lung(:,1),'-o','LineWidth',2,'Color',[27,158,119]/256);
            end
            if isfield(seg, 'r_lung')
                o_seg{i}.r_lung = seg.r_lung(round(linspace(1,length(seg.r_lung),20)),:);
                o_seg{i}.r_lung(end,:) = [];
                %h(4) = plot(app.UIAxes,o_seg.r_lung(:,2),o_seg.r_lung(:,1),'-o','LineWidth',2,'Color',[117,112,179]/256);
            end
        end
    else
        [imgs,~,dcm_info] = segment_dicom([dname, '/', seg_select], working_seg.initial_frame); % Still have to laod the images...
        o_seg = working_seg.bounds;
    end
    segs.(seg_select).bounds = o_seg;
    
    % Sett up the app and buttons
    app.segmentApp = uifigure('Name','Figure 1','Position',[50 50 1100 700]); 
    app.UIAxes = uiaxes(app.segmentApp,'Position',[5 5 900 600]);
    axes(app.UIAxes);

    % Initialize app variables
    pathparts=strsplit(dname, filesep);

    setappdata(app.segmentApp, 'saveName', pathparts{end})
    setappdata(app.segmentApp, 'segName', seg_select)
    setappdata(app.segmentApp, 'ogSeg', o_seg)
    setappdata(app.segmentApp, 'newSeg', o_seg)
    setappdata(app.segmentApp, 'mainSeg', segs)
    setappdata(app.segmentApp, 'imgs', imgs)
    setappdata(app.segmentApp, 'dcm_info', dcm_info)
    %dd = uidropdown(app.segmentApp,'Items',{'All Regions','Exterior','Left Lung','Right Lung'},'Position',[800 500 150 20]);
    dd.Value = 0;
    % Listen for keybindings and use them to triger specific callbacks
    set(app.segmentApp,'KeyPressFcn',{@keybind, app, dd});

    % Set the first frame so the initial screen does not appear blank
    btn_load_callback(0, 0, app, 1); % 0s are place holders 

    % Modify the new points and save in an array
    uibutton(app.segmentApp,'Text','Edit Segment', ...
                            'Position',[800 550 150 20], ...
                            'ButtonPushedFcn',{@edit_callback, app, dd.Value});

    uibutton(app.segmentApp,'Text','Undo Last Point', ...
                            'Position',[800 500 150 20], ...
                            'ButtonPushedFcn',{@undo_callback, app});
                               
    uibutton(app.segmentApp,'Text','Confirm Segment', ...
                            'Position',[800 450 150 20], ...
                            'ButtonPushedFcn',{@save_callback, app});

    uibutton(app.segmentApp,'Text', '+2cm', 'Position',[800 400 150 20],'ButtonPushedFcn',{@btn_load_callback, app, 1}); 
    uibutton(app.segmentApp,'Text', '+1cm', 'Position',[800 350 150 20],'ButtonPushedFcn',{@btn_load_callback, app, 2});  
    uibutton(app.segmentApp,'Text', '4th Intercostal','Position',[800 300 150 20],'ButtonPushedFcn',{@btn_load_callback, app, 3});  
    uibutton(app.segmentApp,'Text', '-1cm', 'Position',[800 250 150 20],'ButtonPushedFcn',{@btn_load_callback, app, 4});  
    uibutton(app.segmentApp,'Text', '-2cm', 'Position',[800 200 150 20],'ButtonPushedFcn',{@btn_load_callback, app, 5}); 
    uibutton(app.segmentApp,'Text', 'Generate Mesh', 'Position',[800 150 150 20],'ButtonPushedFcn',{@gen_mesh, app});
end

function [seg_data, seg_select] = initialize_dialogue(dname, starting_dir)
    folder = [];
    ref_frame = [];
    data = [];
    % Check for dicom files...
    dcm_check = dir(fullfile(dname, '**/*.DCM'));
    if isempty(dcm_check)
        fprintf('No dicom files were detected in the selected folder, please select a folder with dicom images. \n');
        return
    end
    all_folders = {dcm_check.folder};
    unique_folders = unique(all_folders);
    last_dir = [];
    for i=1:length(unique_folders)
        split_dir = strsplit(unique_folders{i},filesep);
        last_dir{i} = split_dir{end};
        sec_last_dir{i} = split_dir{end-1};
    end
    unique_last_dirs = unique(last_dir);
    if (length(unique_last_dirs) < length(unique_folders))
        msg = strcat('\n ___________________________________________________________________________________\n', ...
                       '|                                                                                   |\n', ...
                       '| ERROR: The folder you selected contains multiple Dicom series with the same name. |\n', ...
                       '| Please ensure only one patient is selected.                                       |\n', ...
                       '|___________________________________________________________________________________|\n');
        error('usr:err',msg)
    end
    patient_folder = unique(sec_last_dir);
    if (length(patient_folder) > 1)
        msg = strcat('\n _______________________________________________________________________________________\n', ...
                       '|                                                                                       |\n', ...
                       '| ERROR: It looks like there is Dicom data belonging to more than one patient selected. |\n', ...
                       '| Please ensure only one patient is selected.                                           |\n', ...
                       '|_______________________________________________________________________________________|\n');
        error('usr:err',msg)
        return
    end

    % Select the file from list
    if ismember('SRS00000', unique_last_dirs)
        [~, idx] = ismember('SRS00000', unique_last_dirs);
        unique_last_dirs(idx) = [];
    end
    if ismember('SRS00002', unique_last_dirs)
        [~, idx] = ismember('SRS00002', unique_last_dirs);
    end    
    [indx,~] = listdlg('ListString',unique_last_dirs,'SelectionMode','single','InitialValue',idx);
    var_select = unique_last_dirs{indx};

    % Check for any saved data for the patient folder
    try_load = dir([starting_dir,'/**/', patient_folder{1} ,'.mat']);
    if isempty(try_load)
        save_name = [starting_dir, '/' , patient_folder{1} ,'.mat'];
        default_frame = '';
        alt_seg_level = '';
    else
        save_name = [try_load.folder,'/',try_load.name];
        w = who('-file',[try_load.folder,'/',try_load.name]);
        if ismember('segs',w)
            load([try_load.folder,'/',try_load.name],'segs');

            if isfield(segs,var_select)
                default_frame = segs.(var_select).initial_frame;
                alt_seg_level = '';
            else
                msg = strcat('\n _______________________________________________________________________________\n', ...
                             '|                                                                                 |\n', ...
                             '| WARNING: The saved data for this patient does not contain the expected fields.  |\n', ...
                             '| The saved file will be overwritten.                                             |\n', ...
                             '|_________________________________________________________________________________|\n');
                warning('usr:err',msg)
                default_frame = '';
                alt_seg_level = '';
            end
        else
            msg = strcat('\n _______________________________________________________________________________\n', ...
                         '|                                                                                 |\n', ...
                         '| WARNING: The saved data for this patient does not contain the expected fields.  |\n', ...
                         '| The saved file will be overwritten.                                             |\n', ...
                         '|_________________________________________________________________________________|\n');
            warning('usr:err',msg)
            default_frame = '';
            alt_seg_level = '';
        end
    end

    prompt = {'Enter the frame with the 4^{th} intercostal space:','Change the segmentation seed value [0-1] (Optional):'};
    dlg_title = 'Input';
    num_lines = 1;
    defaultans = {default_frame,alt_seg_level};
    response = inputdlg(prompt,dlg_title,num_lines,defaultans);
    segs.(var_select).initial_frame = response{1}; 
    % Update variables in save file
    % If bounds exist ask before overwriting the files... 
    if isfield(segs.(var_select),'bounds')
        % Segment the data from directory (10 frames)
        answer = questdlg('A segmentation already exists for this series in the file selected, what would you like to do?', ...
                          'Confirm segmentation', ...
	                  'Load segmentation','Overwrite segmentation','Cancel','Load segmentation');
        % Handle response
        switch answer
            case 'Load segmentation'
                % do nothing
            case 'Overwrite segmentation'
                segs.(var_select) = rmfield(segs.(var_select),'bounds');
            case 'Cancel'
                msg = strcat('\n ______________________________________________\n', ...
                               '|                                              |\n', ...
                               '| ERROR: Process terminated by user.           |\n', ...
                               '|______________________________________________|\n');
                error('usr:err',msg)
        end
    end
    save(save_name,'segs');
    seg_data = segs;
    seg_select = var_select;
end

function keybind(src, event, app, dd)
   switch event.Key
       case 'e'
           % Edit segment
           edit_callback(src, event, app, dd);
       case 'z'
           % revert last modified point to original location
           undo_callback(src, event, app);
       case 's'
           % save the entire segmentation 
           save_callback(src, event, app);
       case '1' 
           % Frame 1
           btn_load_callback(src, event, app, 1);
       case '2' 
           % Frame 2
           btn_load_callback(src, event, app, 2);
       case '3' 
           % Frame 3
           btn_load_callback(src, event, app, 3);
       case '4' 
           % Frame 4
           btn_load_callback(src, event, app, 4);
       case '5' 
           % Frame 5
           btn_load_callback(src, event, app, 5);
   end
end

function [new_points] = selectDatapoints(ax, region, app)
    % Get coordinates of all data points in the axes
    %xyobj = findall(ax.Children, '-Property','xData','LineStyle','-','Type','line'); 

    % Check the app object to see what has been segmented  
    % Determine how many objects are in the frame
    frame = getappdata(app.segmentApp , 'frame');
    update_seg = getappdata(app.segmentApp , 'newSeg');
    cur_frame = update_seg{frame};

    combined_data = [];
    if isfield(cur_frame,'exterior')
        ext_index = length(cur_frame.exterior);
        new_points.exterior = flipud(cur_frame.exterior(:,1:2)');
        combined_data = [combined_data, new_points.exterior];
    end 
    if isfield(cur_frame, 'l_lung')
        l_lung_index = length(cur_frame.l_lung);
        new_points.l_lung = flipud(cur_frame.l_lung(:,1:2)');
        combined_data = [combined_data, new_points.l_lung];
    end
    if isfield(cur_frame, 'r_lung')
        r_lung_index = length(cur_frame.r_lung);
        new_points.r_lung = flipud(cur_frame.r_lung(:,1:2)');
        combined_data = [combined_data, new_points.r_lung];
    end
     
    % change title of axes to instructions, in red
    originalTitle = get(ax.Title, {'String', 'Color'}); 
    set(ax.Title, 'String', 'Select the point you would like to move...', 'Color', 'r')
    
    pan(ax, 'off') %turn off panning so the interaction doesn't drag the data.
    roi = drawpoint(ax); 
    % dx, dy to scale the distance
    dx = (ax.XLim(2) - ax.XLim(1));
    dy = (ax.YLim(2) - ax.YLim(1));
    xy = roi.Position;
    xdata = combined_data(1,:);
    ydata = combined_data(2,:);
    
    % what point was closest?
    [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy);
    
    h(1) = plot(ax,xdata(pointslist),ydata(pointslist),'go','LineWidth',2,'Color',[0.2 0.57 0.06]);
    delete(roi);
    
    set(ax.Title, 'String', 'Select the new location of the point...', 'Color', 'r')
    roi = drawpoint(ax);
    xy = roi.Position;
    % Update the points array
    % Determine which shape was edited
    if pointslist <= ext_index
        idx = pointslist;
        if isfield(cur_frame,'exterior')
            new_points.exterior(:,idx) = [xy(1);xy(2)];
            update_seg{frame}.exterior(:,1:2) = fliplr(new_points.exterior'); 
            shape = 'exterior';
        elseif isfield(cur_frame, 'l_lung')
            new_points.l_lung(:,idx) = [xy(1);xy(2)];
            update_seg{frame}.l_lung(:,1:2) = fliplr(new_points.l_lung');
            shape = 'l_lung';
        elseif isfield(cur_frame, 'r_lung')
            new_points.r_lung(:,idx) = [xy(1);xy(2)];
            update_seg{frame}.r_lung(:,1:2) = fliplr(new_points.r_lung');
            shape = 'r_lung';
        end
    elseif pointslist <= ext_index + l_lung_index 
        idx = pointslist-ext_index;
        if isfield(cur_frame, 'l_lung')
            new_points.l_lung(:,idx) = [xy(1);xy(2)];
            update_seg{frame}.l_lung(:,1:2) = fliplr(new_points.l_lung');
            shape = 'l_lung';
        elseif isfield(cur_frame, 'r_lung')
            new_points.r_lung(:,idx) = [xy(1);xy(2)];
            update_seg{frame}.r_lung(:,1:2) = fliplr(new_points.r_lung');
            shape = 'r_lung';
        end  
    elseif pointslist <= ext_index + l_lung_index+ r_lung_index
        idx = pointslist - ext_index - l_lung_index;
        if isfield(cur_frame, 'r_lung')
            new_points.r_lung(:,idx) = [xy(1);xy(2)];
            update_seg{frame}.r_lung(:,1:2) = fliplr(new_points.r_lung');
            shape = 'r_lung';
        end
    else
        disp('Error: the point does not correspond to a shape...')
        return
    end
    
    %Set point last edited
    setappdata(app.segmentApp, 'last_idx', idx)
    setappdata(app.segmentApp, 'last_shape', shape)
    setappdata(app.segmentApp, 'last_frame', frame)
    
    delete(roi)
    setappdata(app.segmentApp, 'newSeg', update_seg)
    
    % Plot the new segment 
    btn_load_callback(0, 0, app, frame)
    
    % Show the current change
    h(2) = quiver(ax, xselect, yselect, xy(1)-xselect, xy(2)-yselect,'LineWidth',1.5,'MaxHeadSize',10);

    % Return original title
    set(ax.Title, 'String', originalTitle{1}, 'Color', originalTitle{2})   
end

function [pointslist,xselect,yselect] = closestpoint(xy,xdata,ydata,dx,dy)
% find the single closest point to xy, in scaled units
if ~iscell(xdata)
  % just one set of points to consider
  D = sqrt(((xdata - xy(1))/dx).^2 + ((ydata - xy(2))/dy).^2);
  [junk,pointslist] = min(D(:)); %#ok
  xselect = xdata(pointslist);
  yselect = ydata(pointslist);
else
  % there is more than one set of points
  Dmin = inf;
  pointslist = cell(size(xdata));
  for i = 1:numel(xdata)
    D = sqrt(((xdata{i} - xy(1))/dx).^2 + ((ydata{i} - xy(2))/dy).^2);
    [mind,ind] = min(D(:)); 
    
    if mind < Dmin
      % searching for the closest point
      Dmin = mind;
      
      pointslist = cell(size(xdata));
      xselect = cell(size(xdata));
      yselect = cell(size(xdata));
      
      pointslist{i} = ind;
      xselect{i} = xdata{i}(ind);
      yselect{i} = ydata{i}(ind);
    end
  end
end
end % subfunction end

function edit_callback(~, ~, app, region)
    selectDatapoints(app.UIAxes, region, app);
end

function save_callback(~, ~, app)
% The variable names are a bit out of control
    seg = getappdata(app.segmentApp, 'newSeg'); 
    segs = getappdata(app.segmentApp, 'mainSeg');
    seg_select = getappdata(app.segmentApp, 'segName');
    segs.(seg_select).bounds = seg;
    savename = getappdata(app.segmentApp, 'saveName');
    save(['../data/',savename, '.mat'],'segs')
end

function undo_callback(src, event, app)
    %Load point last edited
    idx = getappdata(app.segmentApp, 'last_idx');
    shape = getappdata(app.segmentApp, 'last_shape');
    frame = getappdata(app.segmentApp, 'last_frame');
    
    % Revert to original
    working_seg = getappdata(app.segmentApp, 'newSeg');
    og_seg = getappdata(app.segmentApp, 'ogSeg');
    
    og_data      = og_seg{frame}.(shape);
    working_data = working_seg{frame}.(shape);
    working_data(idx,:) = og_data(idx,:);
    working_seg{frame}.(shape) = working_data;
    setappdata(app.segmentApp, 'newSeg', working_seg)
    % display new data
    btn_load_callback(src, event, app, frame)
end

function btn_load_callback(~, ~, app, frame)
    cla(app.UIAxes) % Remove the old objects
    % Load the selected frame and segmentation
    % Plot them on the GUI
    img = getappdata(app.segmentApp, 'imgs');
    I = img(:,:,frame);
    %seg = segs{frame};
    seg = getappdata(app.segmentApp , 'newSeg');
    seg = seg{frame};
    colormap(app.UIAxes, gray(256));
    setappdata(app.segmentApp, 'frame', frame)
    %colormap(app.UIAxes,hsv(256));
    %colormap(app.UIAxes,[hsv(64);hsv(64);hsv(64);hsv(64)])
    h(1) = imagesc(app.UIAxes, I);
    hold(app.UIAxes,'on')
    if isfield(seg,'exterior')
        h(2) = plot(app.UIAxes,seg.exterior(:,2),seg.exterior(:,1),'-o','LineWidth',2,'Color',[217,95,2]/256);
    end 
    if isfield(seg, 'l_lung')
        h(3) = plot(app.UIAxes,seg.l_lung(:,2),seg.l_lung(:,1),'-o','LineWidth',2,'Color',[27,158,119]/256);
    end
    if isfield(seg, 'r_lung')
        h(4) = plot(app.UIAxes,seg.r_lung(:,2),seg.r_lung(:,1),'-o','LineWidth',2,'Color',[117,112,179]/256);
    end
    axis(app.UIAxes,'equal')
    title(app.UIAxes, 'Segmentation Results')
    axis(app.UIAxes,'off')   
end

function [imgs,segs,dcm_info] = segment_dicom(dname,frame)
    frame = str2num(frame);
    % Pick the correct folder and segment the data
    fileinfo = dir(fullfile(dname, '**', '*.DCM'));
    if isempty(fileinfo)
        fileinfo = dir(fullfile(dname, '**', '*.dcm')); % After moving files this can change
    end
    filenames = fullfile({fileinfo.folder}, {fileinfo.name});
    frame = size(filenames,2)-frame;
    
    dcm_info = dicominfo(filenames{1});
    spacing = abs(dcm_info.SliceThickness); % Starts at the top
    frames_per_cm = ((20/spacing)/2);
    segment_size = round([frame-frames_per_cm*2, frame-frames_per_cm, frame, frame+frames_per_cm, frame+frames_per_cm*2]);
    c=0;
    
    for i=1:length(filenames)
        all_imgs(:,:,i) = mat2gray(dicomread(filenames{i}));
    end
    
    for i = segment_size
        c = c+1;
        imgs(:,:,c) = mat2gray(dicomread(filenames{i}));
        info = dicominfo(filenames{i});
        z_dim(c) = (info.SliceThickness/info.PixelSpacing(1))*(i-1); % This is a bit crude - gives z dims
    end

    NumberSlices = 5; 

    for i = 1:length(segment_size)
        [disp_slice(:,:,i), ext, ~, ~, total_lung] = seg_thorax(all_imgs, segment_size(i));

        l_bounds = bwboundaries(total_lung);
        ext_bound = bwboundaries(ext);
        [~,idx] = max(cellfun('length', ext_bound));
        ext_bound = ext_bound{idx};
        ext_bound(:,3) = z_dim(i);
        segs{i}.exterior = ext_bound;
        %plot(ext_bound(:,2), ext_bound(:,1), 'g', 'LineWidth', 2);
        for j = 1:length(l_bounds)
            if length(l_bounds{j}) < 50
                % do nothing with small boundaries
            elseif mean(l_bounds{j}(:,2)) < 230
                %plot(l_bounds{j}(:,2), l_bounds{j}(:,1), 'r', 'LineWidth', 2)
                segs{i}.r_lung = l_bounds{j}; 
                segs{i}.r_lung(:,3) = z_dim(i);
            elseif mean(l_bounds{j}(:,2)) > 280
                %plot(l_bounds{j}(:,2), l_bounds{j}(:,1), 'y', 'LineWidth', 2)
                segs{i}.l_lung = l_bounds{j};
                segs{i}.l_lung(:,3) = z_dim(i);
            end
        end
    end
end

function [adjusted_slice, external, chest_cavity, lung_healthy, total_lung]=seg_thorax(IMG,frame, NumberSlices,thresh_adjust)
% IMG is a series of loaded CT images (512,512, number of frames)
% Frame - a frame selected from the thorax
% Number of slices (should be an odd number, used to define the rib region surrounding the slice)
    lungs = 0;
    heart = 0;
    if nargin<3
        NumberSlices=9;
    end
    if nargin<3
        thresh_adjust = 0;
    end
    % Pixel Dimensions- Required for calculations
    Pixel=0.68164;  %%ENTER THE PIXEL DIMENSION
    ribs=zeros(512,512,NumberSlices);
    % CHECK IMAGE SIZE!!
    if (size(IMG,1) ~= 512) || (size(IMG,2) ~= 512)
        msg = strcat('\n __________________________________________________________________________\n', ...
                       '|                                                                          |\n', ...
                       '| ERROR: These images is not the correct size. Expecing 512 x 512.         |\n', ...
                       '| Ending segmentation...                                                   |\n', ...
                       '|__________________________________________________________________________|\n');
        error('usr:err',msg)
    end
    for i=1:NumberSlices
        % SEGMENTION OF THE CHEST CAVITY
        % getting the ribs to define the chest cavity... 
        Ic = IMG(:,:,frame+(NumberSlices-1)/2-i+1);%%ENTER THE NEW IMAGE USED AS REFERENCE FOR RIBS  
        Iw = wiener2(Ic);           % noise-removal
        J1 = imadjust(Iw);            % basic increase contrast
        % Improve contrast level based on the image
        dist = imhist(J1); 
        dist(1:10) = 0;
        [~,locs] = findpeaks(-dist,'MinPeakProminence',100,'MinPeakDistance',15);
        locs = locs/255;
        level = 0.5; % Make sure it is somewhat close to the original code...
        [~,loc]=min(abs(locs-level));
        level = locs(loc); % This variable naming is weird - just finding threshold for lung tissue
        J1 = imadjust(J1,[level,1]); % increase contrast but better
        %SEGMENT THE BODY BOUNDARY
        se = strel('disk', 20);
        Ie = imerode(J1, se);
        se = strel('disk', 10);
        Iobr = imreconstruct(Ie, J1);
        Iobrd = imdilate(Iobr, se);
        Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
        Iobrcbr = imcomplement(Iobrcbr); % This has a pretty defined heart region %%%%%% TODO!!
        bw = imbinarize(Iobrcbr, 0.5*max(max(Iobrcbr)));
        bw = imfill(bw,'holes');
        bw = imclose(bw,strel('disk',2));
        body_final = imopen(bw,strel('disk', 5));
        if ((NumberSlices-1)/2-i+1) == 0 
            % Save the external boundary
            external = bwareafilt(body_final,1); % Largest object
            adjusted_slice = J1;
            lung_heart = Iobrcbr;
        end
        body_small = imerode(body_final, strel('disk', 8)); % Mask to find ribs and discard skin layer

        %SEGMENT THE RIBS  
        [~,locs] = findpeaks(dist,'MinPeakProminence',1000,'MinPeakDistance',20);
        locs = locs/255;
        if isempty(locs)
            thresh = 0.5;
        else
            level = 0.5;     %TODO can this be better?  %%THRESHOLD-1 LEVEL BETWEEN  0-1
            [~,loc]=min(abs(locs-level));
            thresh = locs(loc); % This variable naming is weird
        end
        t= imbinarize(J1, thresh-0.05-thresh_adjust); % This should show the bones... RIGHT?!?!?!?
        %sum(sum(t))
        c = 0;
        while(sum(sum(t)) > 15000)
            t=imbinarize(J1, (thresh-thresh_adjust+c));
            c=c+0.01;
        end
        ta = (t & body_small);
        tb = imfill(ta,'holes');
        % Clear small stuff in the back 2 thirds
        t_temp = tb;
        t_temp(1:170,:) = 0; % Remove top of image 
        t_temp = bwareafilt(t_temp,20);
        t_temp2 = t;
        t_temp2(171:end,:) = 0;
        tb = t_temp | t_temp2; 

        % remove really small things from total image
        tc = bwareaopen(tb,5);
        td = imclose(tc,strel('disk',5));
        te = imopen(td,strel('disk',2));
        ribs(:,:,i) = te;  
    end
    % Combine to try and make a solid ribcage
    r      = sum(ribs,3);
    rCom   = imbinarize(r,1); % things that occur more than once
    if sum(sum(rCom)) < 12000
        rCom = imbinarize(r,0);
    end
    % Try to close the top
    topEnc = imclose(rCom,strel('rectangle',[15 80]));
    topEnc(256:end,:) = 0; % Remove bottom of image
    topEnc = imclose(topEnc,strel('disk',30));
    rCont   = topEnc | rCom;
    rEnc = imclose(rCont,strel('disk',15));
    rEnc = bwmorph(rEnc,'thicken'); % Just to be sure
    % Keep biggest object
    ribcage = imcomplement(bwareafilt(rEnc,1));
    cavity = imclearborder(ribcage);
    cc = bwconncomp(cavity,4);
    number  = cc.NumObjects;
    if number > 1 % Make sure there is only one object in the ribs
        c = 0;
        while number > 1 % Make sure it is very closed
            cavity = imclose(cavity,strel('rectangle',[10,60+c]));
            cavity = imclose(cavity,strel('disk',15+c));
            c=c+1;
            cc = bwconncomp(cavity,4);
            number  = cc.NumObjects;
        end
    elseif number < 1
        msg = strcat('\n __________________________________________________________________________\n', ...
                       '|                                                                          |\n', ...
                       '| Warning: The lung segmentation was unable to locate any ribs.            |\n', ...
                       '| Trying to segment without the ribs...                                    |\n', ...
                       '|__________________________________________________________________________|\n');
        warning('usr:err',msg)
    end
    c_c  = imdilate(cavity,strel('disk',5));
    c_c2 = c_c | rCom;
    c_c3 = imfill(c_c2,'holes'); % 
    c_c4 = imclose(c_c3,strel('disk',25));
    c_c5 = c_c4 - rCom;
    c_c6 = imopen(c_c5,strel('disk',5));
    c_c7 = bwareafilt(imbinarize(c_c6),1);
    
    windowSize = 41;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(c_c7), kernel, 'same');
    chest_cavity = blurryImage > 0.50; % Rethreshold
    [~,loc] = findpeaks(smoothdata(-sum(chest_cavity),'gaussian',5),'MinPeakProminence',20,'MinPeakDistance',20);
    %%%%%%%%%% NOTE if this does not work we can maybe do a watershed?
    %findpeaks(smoothdata(-sum(chest_cavity),'gaussian',5),'MinPeakProminence',20)
%     test = chest_cavity;
%  keyboard
    % Set up the heart oval...
    imageSizeX = 512;
    imageSizeY = 512;
    R = sum(chest_cavity(:,loc))/2+10;
    span = find(chest_cavity(:,loc) == 1);
    a = 0.8;
    [columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
    centerY = round(mean(span))+10;
    centerX = loc+10; % Shift it to the left...
    heart_region = (rowsInImage - centerY).^2 + ((columnsInImage - centerX).^2)/a^2 <= R.^2;
    test = (heart_region + chest_cavity);
    lung_region = chest_cavity-heart_region;
    %C1 = imfuse(lung_region, adjusted_slice,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    %imresize(C1,0.5);
    %subplot(311)
    %imshow(C1)
    lung_region = chest_cavity-heart_region;
   %keyboard
    % DO THE HEALTHY LUNG SEGMENTATION
    % Get healthy lungs
    la = imcomplement(adjusted_slice);
    lb = imclearborder(la);
    lb = lb .* external; 
    lc = imbinarize(imbinarize(lb)-ribs(:,:,5)); % Binarize and remove the bones
    ld = imopen(lc,strel('disk',3));
    le = bwareaopen(ld,250);%% Just keep the biggest 2 objects instead?
    le = bwareafilt(le,2);
    le = imclose(le,strel('disk',2));
    le = imfill(le,'holes');
    %C2 = imfuse(la,le,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    %imresize(C2,0.5);
    %subplot(312)
    %imshow(C2)
    %subplot(313)
    lung_estimate = imbinarize(lung_region - ribs(:,:,5));
    % Smooth the lung estimate a bit?
    temp = lung_estimate - le;
    temp2 = imopen(temp,strel('disk',10));
    temp3 = temp2 | le;
    temp4 = imclose(temp3,strel('disk',5));
    %imshow(temp4)
    %imshow(le + lung_estimate)
    %C3 = imfuse(adjusted_slice,temp4,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    %imresize(C3,0.5);
    %subplot(313)
    %imshow(C3)
    lb = imclearborder(la);
    lung_rough = imbinarize(lb);
    lung_healthy = le;
    
    % Smooth the lungs a little bit...
    windowSize = 5;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(temp4), kernel, 'same');
    total_lung = blurryImage > 0.50; % Rethreshold
    % Make sure there are two lungs
    cc = bwconncomp(total_lung,4);
    num_objects  = cc.NumObjects;
    c=0;
    while num_objects == 1 % If there is only one lung shape - FIX IT!!
        c = c+1;
        total_lung = imopen(total_lung,strel('disk',2+c));
        total_lung = imclose(total_lung,strel('disk',2));
        total_lung = imfill(total_lung,'holes');
        %lung_healthy = total_lung;
        cc = bwconncomp(total_lung,4);
        num_objects  = cc.NumObjects;
    end
return    
end

function gen_mesh(~, ~, app)
run C:\Users\symon\Documents\eidors-v3.10-ng\eidors-v3.10-ng\eidors\startup
% Save before meshing since it should be the final result that is meshed
    save_callback(0,0,app);
    % We want to give a 2d extruded option and also a 3D based on enclosure
    seg = getappdata(app.segmentApp, 'newSeg');
    dcm_info = getappdata(app.segmentApp, 'dcm_info');
    savename = getappdata(app.segmentApp, 'saveName');
    x_spacing = dcm_info.PixelSpacing(1); 
    y_spacing = dcm_info.PixelSpacing(2); 
    % Select the middle segment
    slice = 3;
    ext   = seg{slice}.exterior(:,1:2);
    l_lng = seg{slice}.l_lung(:,1:2);
    r_lng = seg{slice}.r_lung(:,1:2);
    % Center and scale the points
%     ext(:,1)   = round(((ext(:,1) -   round(double(dcm_info.Width/2)) )*x_spacing)/100,1);
%     ext(:,2)   = round(((ext(:,2) -   round(double(dcm_info.Height/2)))*y_spacing)/100,1);
%     l_lng(:,1) = round(((l_lng(:,1) - round(double(dcm_info.Width/2)) )*x_spacing)/100,1);
%     l_lng(:,2) = round(((l_lng(:,2) - round(double(dcm_info.Height/2)))*y_spacing)/100,1);
%     r_lng(:,1) = round(((r_lng(:,1) - round(double(dcm_info.Width/2)) )*x_spacing)/100,1);
%     r_lng(:,2) = round(((r_lng(:,2) - round(double(dcm_info.Height/2)))*y_spacing)/100,1);
    
    ext(:,1)   = round(((ext(:,1)  )*x_spacing)/100,1);
    ext(:,2)   = round(((ext(:,2)  )*y_spacing)/100,1);
    l_lng(:,1) = round(((l_lng(:,1))*x_spacing)/100,1);
    l_lng(:,2) = round(((l_lng(:,2))*y_spacing)/100,1);
    r_lng(:,1) = round(((r_lng(:,1))*x_spacing)/100,1);
    r_lng(:,2) = round(((r_lng(:,2))*y_spacing)/100,1);
    
    
    outline = seg{slice}.exterior;
    lung_1  = seg{slice}.l_lung(:,1:2);
    lung_2  = seg{slice}.r_lung(:,1:2);
    trunk  = outline(:,1:2)/256;
    l_a = lung_1(:,1:2)/256;
    l_b = lung_2(:,1:2)/256; 
    %  Height (cm) 
    h = 5/10;
    elec_radius = 0.3/10;
    fmdl = ng_mk_extruded_model({0.2,{trunk, l_a, l_b}, [4,40], 0.02},[16,0,0.1], [0.01]);
    %keyboard
    %fmdl = ng_mk_extruded_model({h,{ext, l_lng, r_lng}, [4,40], 0.08},[16,0,h/2], [elec_radius]);
    %fmdl = ng_mk_extruded_model({h,{ext, r_lng}, [4,41], 0.08},[16,0,h/2], [elec_radius]);
    %fmdl = ng_mk_extruded_model({h,{ext}, [4,40], 0.08},[16,0,h/2], [elec_radius]);
    img = mk_image(fmdl,1);
    img.elem_data( fmdl.mat_idx{2} ) = 0.9;
    img.elem_data( fmdl.mat_idx{3} ) = 0.9;
    fmdl= img;
    show_fem(fmdl)
    % Save the FMDL with the filename from before...
    save([savename, '_fmdl.mat'],'fmdl');
    return
    
    % MAYBE LATER WHEN TRYING TO DO A 3D MESH OF THE LUNGS
    if 0
        % generate the point clouds
        ext_cloud    = [];
        l_lung_cloud = [];
        r_lung_cloud = [];
        for i = 1:length(seg)
            ext_cloud    = [ext_cloud;    seg{i}.exterior];
            l_lung_cloud = [l_lung_cloud; seg{i}.l_lung];
            r_lung_cloud = [r_lung_cloud; seg{i}.r_lung];
        end
        % Offset everything to 0 in the z axis
        offset = min(ext_cloud(:,3));
        ext_cloud(:,3) = (ext_cloud(:,3) - offset);
        l_lung_cloud(:,3) = (l_lung_cloud(:,3) - offset);
        r_lung_cloud(:,3) = (r_lung_cloud(:,3) - offset);
        ext_bounds    = boundary(ext_cloud,0);
        l_lung_bounds = boundary(l_lung_cloud,0.6);
        r_lung_bounds = boundary(l_lung_cloud,0.6);
        figure(1); clf;
        trisurf(ext_bounds,ext_cloud(:,1),ext_cloud(:,2),ext_cloud(:,3),'Facecolor','blue','FaceAlpha',0.05)
        hold on
        trisurf(l_lung_bounds,l_lung_cloud(:,1),l_lung_cloud(:,2),l_lung_cloud(:,3),'Facecolor','red','FaceAlpha',0.1)
        trisurf(r_lung_bounds,r_lung_cloud(:,1),r_lung_cloud(:,2),r_lung_cloud(:,3),'Facecolor','green','FaceAlpha',0.1)
        axis equal
        hold off
    end
    keyboard
end