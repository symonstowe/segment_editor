function [imdl,ROIs,pp] = mk_imdl(pp,mdl_sel,bnds,msks,pt_num)
	%   fmdl.electrode = fmdl.electrode(fliplr([9:16,2:8]));
	% bnds should have fields:
	% bnds.ext
	% bnds.l_lung
	% bnds.r_lung
	% bnds.mask.l_lung % these are to make point clouds of the lungs and find points within it... should be very large 
	% bnds.mask.r_lung % these are to make point clouds of the lungs and find points within it... should be very large 
	fmdl = mk_mdl(pp,mdl_sel,bnds,msks,pt_num);
	%pp = get_boundary(fmdl,pp);

	%fmdl = mdl_normalize(fmdl,1);
	[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current_next1'},1);
	
	[img,vh,vi] = set_elem_background(fmdl);
	%opt.imgsz = [100 100];
	%opt.noise_figure = 0.8;
	%opt.square_pixels = 1;
	opt.imgsz = [50 50];
	opt.distr = 3; % non-random, uniform
	opt.Nsim = 500; % 500 hundred targets to train on, seems enough
	opt.target_size = 0.01; %small targets % original was 0.3
	opt.target_offset = 0;
	opt.noise_figure = .5; % this is key!
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
		img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 0.3; % lungs 
	else
		test = 5
		keyboard
	end
	%m_frac = elem_select(fmdl, @(x,y,z) ...
	%	(x.^2)/0.6 + (y-0.25).^2/0.55 < 1);
	%img.elem_data = img.elem_data + m_frac*2;
	%keyboard
	vi = fwd_solve(img);
end