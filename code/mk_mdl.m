function [fmdl] = mk_mdl(pp,mdl_sel,bnds,msks,pt_num)
	%   fmdl.electrode = fmdl.electrode(fliplr([9:16,2:8]));
	% bnds should have fields:
	% bnds.ext
	% bnds.l_lung
	% bnds.r_lung
	% bnds.mask.l_lung % these are to make point clouds of the lungs and find points within it... should be very large 
	% bnds.mask.r_lung % these are to make point clouds of the lungs and find points within it... should be very large 
	if pt_num == 2
		num_points = 41; % needs to be played with A LOT for this method to work!!!!!!!
		%elec_spacing = 1.25; % This is perfect - don't ever touch it
		%elec_spacing = 1.265; % This is better makes the electrode nugget go away
		elec_spacing = 1.25; % This is better makes the electrode nugget go away
		% WHY IS AN EXTRA ELECTRODE NUGGET APPEARING?!?!?!
		%elec_order = [8,7,6,5,4,3,2,1,16,15,14,13,12,11,10,9];
		elec_order = [1:16];
		%elec_order = fliplr([1:16]);
	elseif pt_num == 3
		%num_points = 50; % needs to be played with A LOT for this method to work!!!!!!!
		%elec_spacing = 1.31; % SET don't change from 1.31!!
		num_points = 41; % needs to be played with A LOT for this method to work!!!!!!!
		elec_spacing = 1.30; % SET don't change from 1.31!!
		elec_order = [1:16];
	elseif pt_num == 4
		num_points = 41; % needs to be played with A LOT for this method to work!!!!!!!
		elec_spacing = 1.28;
		elec_order = [1:16];
	elseif pt_num == 5
		num_points = 41; % needs to be played with A LOT for this method to work!!!!!!!
		elec_spacing = 1.28;
		% 1.28 works but is not quite aligned
		% Electrode 6 is a bit stretched on a protrusion - will have to work
		elec_order = [16,1:15];
		%elec_order = [8,7,6,5,4,3,2,1,16,15,14,13,12,11,10,9];
	end
	switch mdl_sel;
	case 0; 

		elec_pos = [16,1.015,.5]; elec_shape=[0.05]; maxsz=0.05; nfft=27;
		fmdl = mk_library_model({'neonate','boundary','left_lung','right_lung'}, ...
		      elec_pos, elec_shape, maxsz,nfft);
		%evol = get_elem_volume(fmdl);
		%pp.lmag   = sum(evol(fmdl.mat_idx{3})) / ...
		%	sum(evol(fmdl.mat_idx{2}));
	case 1; % Use the library model
		%fmdl = mk_library_model({'adult_male_32el_lungs','boundary'},...
		%	[16 1.015 0.5],[0.05],0.02);
		fmdl = mk_library_model({'adult_male','boundary','left_lung','right_lung'},...
			[32 1.0 0.5],[0.05],0.05);
		fmdl.electrode(1:2:32)=[];
		pp.lmag = 1.2;
	case 2; % Use the extruded model with custom bounds
		model_height = 0.8;
		elec_height = model_height/2;
		elec_size = 0.04;
		% for num_points usually something around 40 works 
		trunk = bnds{3}.exterior(:,1:2)/256;
		l_l   = bnds{3}.l_lung(:,1:2)/256;
		l_r   = bnds{3}.r_lung(:,1:2)/256;
		fmdl = ng_mk_extruded_model({model_height,{trunk, l_l, l_r}, [4,num_points], 0.05},[16,elec_spacing,elec_height], [elec_size]);
		pp.lmag = 1.5;
	case 3; % Use the lung point clouds with alphashapes
		model_height = 0.8;
		elec_height = model_height/2;
		elec_size = 0.04;
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
		b_shp = alphaShape([msks.B.r/256,msks.B.c/256,flipud(zpts)],0.2);
		fmdl = ng_mk_extruded_model({model_height,{trunk}, [4,num_points], 0.05},[16,elec_spacing,elec_height], [elec_size]);
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
	fmdl.electrode = fmdl.electrode(elec_order);