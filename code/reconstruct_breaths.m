function reconstruct_breaths(fname,seg_fname,CT_fname,ref_frame,pt_num)
	% Reconstruct all individual breaths in a file 
	% and the EA breath in a seperate figure
	[dd,auxdata,stim] = eidors_readdata(fname,"DRAEGER-EIT");
	pp.FR = stim(1).framerate;
	T = [0:1/pp.FR:(size(dd,2)/pp.FR)-1/pp.FR];
	pp.LPF = 0.8;
	pp = set_parameters(pp,dd);
	dd = preproc_data(dd,pp);
	% TODO select complete breaths to reconstruct...
	breaths = detect_breaths(dd);
	% Just do the first 5 breaths
	num_breaths = 5;
	seg_sel = [breaths(1).trgh(1):breaths(num_breaths).trgh(2)];
	clf;
	set(gcf,'renderer','painters');
	set(groot,'defaultAxesTickLabelInterpreter','latex');  
	set(groot,'defaulttextinterpreter','latex');
	set(groot,'defaultLegendInterpreter','latex');
	tiledlayout(5,5, 'Padding', 'none', 'TileSpacing', 'compact');
	nexttile([1 5])
	plot(T(seg_sel),sum(dd(:,seg_sel)))
	set(get(gca, 'XLabel'), 'String', 'time(s)');
	set(get(gca, 'YLabel'), 'String', '$\Delta$ Z');
	%xlim([3 21])
	axis tight
	ax = gca;
	ax.FontSize = 16; 
	hold on 
	for i=1:num_breaths
		xline((breaths(i).pk-1)/pp.FR,'r')
		xline((breaths(i).trgh(1)-1)/pp.FR,'b')
		xline((breaths(i).trgh(2)-1)/pp.FR,'b')
	end
	% Get the custom data to generate the plots
	seg_data = load(seg_fname);
	bounds = seg_data.segs.SRS00002.bounds;
	lung_masks = get_lung_masks(CT_fname,ref_frame);
	% Reconstruct the images 
	for j = 0:3 % use the 4 different types of models
		[imdl,ROIs,pp] = mk_imdl(pp,j,bounds,lung_masks,pt_num); 
		pp.ROIs = ROIs;

		for i=1:num_breaths
			img = inv_solve(imdl, dd(:,breaths(i).trgh(1)), dd(:,breaths(i).pk));
			nexttile 
			img.calc_colours.ref_level = 0; 
			show_slices(img)
			if j>1; view([90 90]); end
		end
	end
	% Show the breath reconstructions 
	set(gcf,'Position',[949          88        1397        1162])
	print(['../imgs/breath_imgs_PT0' num2str(pt_num)], '-dsvg');
	%TODO - make a 2x2 figure to show boundary differences 

	%TODO - make a 2x3 figure to show the lung shape differences 

	

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
 
%function pp = get_boundary(fmdl,pp);
	%fmdl2 = mdl2d_from3d(fmdl);
	%b = fmdl2.boundary; bs = size(b);
	%pp.boundaryx = reshape(fmdl2.nodes(b,1),bs)';
	%pp.boundaryy = reshape(fmdl2.nodes(b,2),bs)';

function [img,vh,vi] = set_elem_background(fmdl)
	img = mk_image(fmdl,1);
	vh = fwd_solve(img);
	%   img.elem_data(vertcat(fmdl.mat_idx{2:3}))= 2;
	if numel(fmdl.mat_idx)>1
		img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 2; % lungs 
	end
	m_frac = elem_select(fmdl, @(x,y,z) ...
		(x.^2)/0.6 + (y-0.25).^2/0.55 < 1);
	img.elem_data = img.elem_data + m_frac*2;
	%keyboard
	vi = fwd_solve(img);