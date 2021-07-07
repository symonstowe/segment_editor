function reconstruct_breaths(fname,seg_fname,CT_fname,ref_frame)
	% Reconstruct all individual breaths in a file 
	% and the EA breath in a seperate figure
	[dd,auxdata,stim] = eidors_readdata(fname,"DRAEGER-EIT");
	pp.FR = stim(1).framerate;
	pp.LPF = 0.8;
	pp = set_parameters(pp,dd);
	dd = preproc_data(dd,pp);
	% TODO select complete breaths to reconstruct...
	breaths = detect_breaths(dd);
	% Just do the first 5 breaths
	num_breaths = 5;
	seg_sel = [breaths(1).trgh(1):breaths(num_breaths).trgh(2)];
	clf;
	tiledlayout(4,5, 'Padding', 'none', 'TileSpacing', 'compact');
	nexttile([1 5])
	plot(sum(dd(:,seg_sel)))
	hold on 
	for i=1:num_breaths
		xline(breaths(i).pk-seg_sel(1),'r')
		xline(breaths(i).trgh(1)-seg_sel(1),'b')
		xline(breaths(i).trgh(2)-seg_sel(1),'b')
	end
	% Get the custom data to generate the plots
	seg_data = load(seg_fname);
	bounds = seg_data.segs.SRS00002.bounds;
	lung_masks = get_lung_masks(CT_fname,ref_frame);
	% Reconstruct the images 
	for j = 1:3 % use the 3 different types of models
		[imdl,ROIs,pp] = mk_imdl(pp,j,bounds,lung_masks); 
		pp.ROIs = ROIs;
		for i=1:num_breaths
			img = inv_solve(imdl, dd(:,breaths(i).trgh(1)), dd(:,breaths(i).pk));
			nexttile 
			show_slices(img)
		end
		keyboard
	end
	%TODO - make a 2x2 figure to show boundary differences 

	%TODO - make a 2x3 figure to show the lung shape differences 

keyboard
	

function  pp = set_parameters(pp,dd);
	if ~isfield(pp,'Nel'); 
		switch size(dd,1);
		case  208; pp.Nel = 16;
		case 1024; pp.Nel = 32;
		otherwise; error 'huh?';
		end
	end
	if ~isfield(pp,'FR'); 
		pp.FR = 50;
	end
	if ~isfield(pp,'LPF');
		pp.LPF = 1.2;
	end

function dd = preproc_data(di,pp)
	di = real(freq_filt(di,@(f) f<pp.LPF,pp.FR,2));
	
	[~, msel] = mk_stim_patterns(pp.Nel,1,[0,1],[0,1],{'no_meas_current'},1);
	dd = zeros(pp.Nel^2,size(di,2));
	dd(msel, :) = di; 
 
function [imdl,ROIs,pp] = mk_imdl(pp,mdl_sel,bnds,msks)
	%   fmdl.electrode = fmdl.electrode(fliplr([9:16,1:8]));
	% bnds should have fields:
	% bnds.ext
	% bnds.l_lung
	% bnds.r_lung
	% bnds.mask.l_lung % these are to make point clouds of the lungs and find points within it... should be very large 
	% bnds.mask.r_lung % these are to make point clouds of the lungs and find points within it... should be very large 
	switch mdl_sel;
	case 0; % Don't need this 
		fmdl= mk_library_model('adult_male_32el_lungs');
		fmdl.electrode(2:2:32)=[];
		evol = get_elem_volume(fmdl);
		pp.lmag   = sum(evol(fmdl.mat_idx{3})) / ...
			sum(evol(fmdl.mat_idx{2}));
	case 1; % Use the library model
		fmdl = mk_library_model({'adult_male','boundary'},...
			[16 1.015 0.5],[0.05],0.08);
		pp.lmag = 1.2;
	case 2; % Use the extruded model with custom bounds
		fmdl = ng_mk_extruded_model({model_height,{trunk, l_l, l_r}, [4,num_points], 0.02},[16,elec_spacing,elec_height], [elec_size]);
		pp.lmag = 1.5;
	case 3; % Use the lung point clouds with alphashapes
		fmdl = ng_mk_extruded_model({model_height,{trunk}, [4,num_points], 0.02},[16,elec_spacing,elec_height], [elec_size]);
        
	end
	%pp = get_boundary(fmdl,pp);

	fmdl = mdl_normalize(fmdl,1);
	[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current_next1'},1);
	
	[img,vh,vi] = set_elem_background(fmdl);
	opt.imgsz = [100 100];
	opt.noise_figure = 0.8;
	opt.square_pixels = 1;
	imdl=mk_GREIT_model(img, 0.25, [], opt);
	rmdl = imdl.rec_model;
	xy = rmdl.coarse2fine'*interp_mesh(rmdl);
	ROI = getfield(inv_solve(imdl,vh,vh.meas*2),'elem_data')<-10;
	%ROIl=ROI.*(xy(:,1)<0);
	%ROIr=ROI.*(xy(:,1)>0);
	ROIl=(xy(:,1)<0);
	ROIr=(xy(:,1)>0);
	ROIs = [ROIl(:), ROIr(:), ROI(:)];
	ROI = getfield(inv_solve(imdl,vh,vi),'elem_data')>500;
	ROIs = [ROIs, ROI(:)];

%function pp = get_boundary(fmdl,pp);
	%fmdl2 = mdl2d_from3d(fmdl);
	%b = fmdl2.boundary; bs = size(b);
	%pp.boundaryx = reshape(fmdl2.nodes(b,1),bs)';
	%pp.boundaryy = reshape(fmdl2.nodes(b,2),bs)';

function [img,vh,vi] = set_elem_background(fmdl)
	img = mk_image(fmdl,1);
	vh = fwd_solve(img);
	%   img.elem_data(vertcat(fmdl.mat_idx{2:3}))= 2;
	m_frac = elem_select(fmdl, @(x,y,z) ...
		(x.^2)/0.6 + (y-0.25).^2/0.55 < 1);
	img.elem_data = img.elem_data + m_frac*2;
	vi = fwd_solve(img);

function [breath] = detect_breaths(dd);
	breath_sig = sum(dd);
	breath_sig = rescale(breath_sig,0,1);
	prominence = 0.25; 
    
        pks = islocalmax(breath_sig,'MinProminence',prominence,'FlatSelection','first');
        p = prominence;
        while isempty(pks)
            p = p-0.05;
            pks = islocalmax(breath_sig,'MinProminence',p,'FlatSelection','first');
        end
        pks = find(pks);
        % Find the breath troughs using findchangepts
        idx = findchangepts(breath_sig,'Statistic','linear','MinThreshold',0.2); % TODO needs to be adjusted base on Fs eventually
        dx = diff(idx);
        dy = diff(breath_sig(idx));
        slope = dy./dx;
        breathTroughs = idx(slope > 0.0001); % remove negative slopes
        % discard peaks that don't have at least one trough before or after
        pks = pks(pks<max(breathTroughs) & pks>min(breathTroughs));
	breath = [];
        if isempty(pks)
            keyboard
        end
        for i=1:numel(pks)
            breath(i).pk   = pks(i);%   = pks(pks>breathTroughs(i) & pks<breathTroughs(i+1));
            trghs_diff = breathTroughs-pks(i);
            low_diffs  = trghs_diff;
            high_diffs = trghs_diff;
            low_diffs(low_diffs>0)    = inf;
            high_diffs(high_diffs<0)  = inf;
            [~,T1] = min(abs(low_diffs));
            [~,T2] = min(high_diffs);
            breath(i).trgh(:) = [breathTroughs(T1), breathTroughs(T2)];
            length_inhale(i) = breath(i).pk - breath(i).trgh(1);
            length_exhale(i) = breath(i).trgh(2) - breath(i).pk; 
        end
        % Remove high outliers in breath lengths - assume something weird is going on...
        while max(length_exhale)>1.25*mean(length_exhale) 
            [~,reject] = max(length_exhale);
            length_exhale(reject) = [];
            length_inhale(reject) = [];
            breath(reject) = [];
            pks(reject) = [];
        end

function mask = get_lung_masks(CT_dir,ref_frame)
	% TODO could we have just used poly2mask?!?!?!?!?!?
	% Segment the lung frame for everywhere except the frames too close to the edge
	% create a meshgrid of the masks to return a point cloud for alphashape generation
	fileinfo = dir(fullfile(CT_dir, '**', '*.DCM'));
	filenames = fullfile({fileinfo.folder}, {fileinfo.name});
	for i=1:length(filenames)
		imgs(:,:,i) = mat2gray(dicomread(filenames{i}));
	end
	c=0;
	for i=ref_frame-10:ref_frame+10
		c=c+1;
		se = strel('disk',25); % just some extra smoothing to make everything go smoother...
		[~, ~, ~, lung_healthy(:,:,c), total_lung(:,:,c)]=seg_thorax(imgs,i);
		%lung_healthy(:,:,c) = imopen(lung_healthy(:,:,c));
		total_lung(:,:,c) = imopen(total_lung(:,:,c), se);
		s = regionprops(total_lung(:,:,c),'centroid');
		x1 = s(1).Centroid(1);
		x2 = s(2).Centroid(1);
		[~,I_min] = min([x1,x2]); 
		b = regionprops(total_lung(:,:,c),'BoundingBox');
		for k=1:2 % If there are more than 2 lungs we have bigger issues
			x = [floor(b(k).BoundingBox(1)),floor(b(k).BoundingBox(1)), ...
			     floor(b(k).BoundingBox(1))+floor(b(k).BoundingBox(3)),floor(b(k).BoundingBox(1))+floor(b(k).BoundingBox(3))];
			y = [floor(b(k).BoundingBox(2)),floor(b(k).BoundingBox(2))+floor(b(k).BoundingBox(4)), ...
			     floor(b(k).BoundingBox(2))+floor(b(k).BoundingBox(4)),floor(b(k).BoundingBox(2))];
			l_msk = poly2mask(x,y,512,512);
			if k == I_min % Right lung
				right_lung(:,:,c) = total_lung(:,:,c).*l_msk;
			else
				left_lung(:,:,c) = total_lung(:,:,c).*l_msk;
			end
		end
	end
	[r,c,v] = ind2sub(size(right_lung),find(right_lung));
	mask.A.r = r;
	mask.A.c = c;
	mask.A.v = v;
	[r,c,v] = ind2sub(size(left_lung),find(left_lung));
	mask.B.r = r;
	mask.B.c = c;
	mask.B.v = v;
	