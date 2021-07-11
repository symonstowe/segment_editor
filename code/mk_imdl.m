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
		model_height = 0.8;
		elec_height = model_height/2;
		elec_size = 0.04;
		elec_spacing = 1.25;
		num_points = 40;
		trunk = bnds{3}.exterior(:,1:2)/256;
		l_l   = bnds{3}.l_lung(:,1:2)/256;
		l_r   = bnds{3}.r_lung(:,1:2)/256;
		fmdl = ng_mk_extruded_model({model_height,{trunk, l_l, l_r}, [4,num_points], 0.02},[16,elec_spacing,elec_height], [elec_size]);
		pp.lmag = 1.5;
	case 3; % Use the lung point clouds with alphashapes
		model_height = 0.8;
		elec_height = model_height/2;
		elec_size = 0.04;
		elec_spacing = 1.25;
		num_points = 40;
		trunk = bnds{3}.exterior(:,1:2)/256;
		% Do lung A
		layers = unique(msks.A.v);
		new_layers = linspace(0,model_height,numel(layers));
		zpts = msks.A.v;
		for i=1:size(layers)
			zpts(zpts == layers(i)) = new_layers(i);
		end
		a_shp = alphaShape([msks.A.r/256,msks.A.c/256,flipud(zpts)],0.1);
		% Do lung B
		layers = unique(msks.B.v);
		new_layers = linspace(0,model_height,numel(layers));
		zpts = msks.B.v;
		for i=1:size(layers)
			zpts(zpts == layers(i)) = new_layers(i);
		end
		b_shp = alphaShape([msks.B.r/256,msks.B.c/256,flipud(zpts)],0.1);
		fmdl = ng_mk_extruded_model({model_height,{trunk}, [4,num_points], 0.02},[16,elec_spacing,elec_height], [elec_size]);
		mesh_pts = fmdl.nodes;
		idx = inShape(a_shp,mesh_pts(:,1),mesh_pts(:,2),mesh_pts(:,3));
		I = find(idx == 1);
        	% Find elements made up of only selected nodes 
        	elem_sel = ismember(fmdl.elems,I);
        	I2 = find(sum(elem_sel,2) == 4); 
        	fmdl.mat_idx{2} = I2;
        	idx = inShape(b_shp,mesh_pts); % Compare 
        	I = find(idx == 1);
        	% Find elements made up of only selected nodes 
        	elem_sel = ismember(fmdl.elems,I);
        	I2 = find(sum(elem_sel,2) == 4); 
        	fmdl.mat_idx{3} = I2;
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
end

function [img,vh,vi] = set_elem_background(fmdl)
	img = mk_image(fmdl,1);
	vh = fwd_solve(img);
	%   img.elem_data(vertcat(fmdl.mat_idx{2:3}))= 2;
	if numel(fmdl.mat_idx)>1
		img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 2; % lungs 
	end
	%m_frac = elem_select(fmdl, @(x,y,z) ...
	%	(x.^2)/0.6 + (y-0.25).^2/0.55 < 1);
	%img.elem_data = img.elem_data + m_frac*2;
	%keyboard
	vi = fwd_solve(img);
end