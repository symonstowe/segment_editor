function output = thesis_figures(fig_num)
	fig_num = num2str(fig_num);
	data_dir = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data';
	% load the DICOM data
	switch fig_num
	case '1'
		% segmentation results image 
	case '2'
		% segmentation methods fig
		% show some of the intermedeate stages on the processs  
		% Load the data for reference subject 4
		dname = [data_dir, '/PTS4/SRS00002'];
		fileinfo = dir(fullfile(dname, '**', '*.DCM'));
		filenames = fullfile({fileinfo.folder}, {fileinfo.name});
		frame = length(filenames) - 107;
		for i=1:length(filenames)
			imgs(:,:,i) = mat2gray(dicomread(filenames{i}));
			%show_slices(imgs(:,:,i));
		end
		[e1,e2,e3,e4,e5,r1,r2,r3,r4,r5,r6,r7,r8,l1,l2,l3,l4] = step_seg_thorax(imgs,frame);
		% Steps to segment the external boundary
		f = figure(6); clf;
		set(gcf,'renderer','painters');
		set(groot,'defaultAxesTickLabelInterpreter','latex');  
		set(groot,'defaulttextinterpreter','latex');
		set(groot,'defaultLegendInterpreter','latex');
		tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact');
		% Print this as a figure (use subject 4?)
		xs = ["(A)","(B)","(C)","(D)","(E)"];
		for i=1:5
			nexttile 
			imshow(eval(['e',num2str(i)]))
			set(gca,'FontSize',25)
			set(get(gca, 'XLabel'), 'String', xs(i));
		end
		ext_plot = bwboundaries(e5);
		nexttile
		imshow(imgs(:,:,frame))
		hold on
		plot(ext_plot{1}(:,2),ext_plot{1}(:,1),'r','Linewidth',2)
		set(get(gca, 'XLabel'), 'String', '(F)');
		set(gca,'FontSize',25)
		xlim([0 512])
		ylim([0 512])
		set(gcf,'Position',[0          0        1768        1241])
		print('../imgs/boundary_seg_methods', '-dsvg');
		% Do the rib segmentation next 
		clf
		tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
		xs = ["(A)","(B)","(C)","(D)","(E)","(F)","(G)","(H)"];
		for i=1:8
			nexttile 
			imshow(eval(['r',num2str(i)]))
			set(gca,'FontSize',25)
			set(get(gca, 'XLabel'), 'String', xs(i));
		end
		set(gcf,'Position',[0          0        1968        1241])
		print('../imgs/chest_cavity_seg_methods', '-dsvg');
		% Do the lung segmentation methods 
		clf
		tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact');
		xs = ["(A)","(B)","(C)","(D)","(E)","(F)","(G)","(H)"];
		for i=1:4
			nexttile 
			imshow(eval(['l',num2str(i)]))
			set(gca,'FontSize',25)
			set(get(gca, 'XLabel'), 'String', xs(i));
		end
		set(gcf,'Position',[0          0        1241        1241])
		print('../imgs/lung_seg_methods', '-dsvg');
	case '3'
		% Show the 2D basic mesh 

	case '4'
		% Show the 2.5 D mesh 

		% A) with the 5 slices 

		% B) with 25 slices
	case '5'
		% reconstrusionts with with 3 meshes 
		% A) regular B) custom C) 2.5
		% This takes ages...
		for i=2:5 
			if i == 4
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432623_xueliu/2432623_xueliu_01.eit';
				CT_dir = [data_dir, '/PTS4/SRS00002'];
				ref_frame = 136-107;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS4.mat';
				reconstruct_breaths(fname,seg_data,CT_dir,ref_frame,i);
			elseif i == 2
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432616_xueliu/2432616_xueliu-peeptit_02_001.eit';
				CT_dir = [data_dir, '/PTS2/SRS00002'];
				ref_frame = 145-112;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS2.mat';
				reconstruct_breaths(fname,seg_data,CT_dir,ref_frame,i);
			elseif i == 3
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_23/2432000_xueliu__0-12pep/2432000_xueliu__peep_12_001.eit';
				CT_dir = [data_dir, '/PTS3/SRS00002'];
				ref_frame = 53-28;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS3.mat';
				reconstruct_breaths(fname,seg_data,CT_dir,ref_frame,i);
			elseif i == 5
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_29/2394977_xueliu-peep0/23394977_xueliu-peep0_01.eit';
				CT_dir = [data_dir, '/PTS5/SRS00002'];
				ref_frame = 73-42;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS5.mat';
				reconstruct_breaths(fname,seg_data,CT_dir,ref_frame,i);
			end
		end
	case '6' 
		% Reconstruct with the jpg slices 	
		jpg_dir = [data_dir '/JPG-EIT_data/1680127/'];
		% Identify a folder 
		img = imread([jpg_dir '1.JPG']);
		keyboard
	case '7'
		% Show the FEM for each of the 3 different methods 
		%fname = ['~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432623_xueliu/2432623_xueliu_01.eit'];
		for i=2:5 
			if i == 2
				ref_frame = 145-112;
			elseif i == 3
				ref_frame = 53-28;
			elseif i == 4
				ref_frame = 136-107;
			elseif i == 5
				ref_frame = 73-42;
			end
			CT_dir = [data_dir, '/PTS' num2str(i) '/SRS00002'];
			seg_data = ['~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS' num2str(i) '.mat'];

			% Get the custom data to generate the plots
			seg_data = load(seg_data);
			bounds = seg_data.segs.SRS00002.bounds;
			CT_fname = CT_dir;
			lung_masks = get_lung_masks(CT_fname,ref_frame);
			pp = [];
			clf;
			set(gcf,'renderer','opengl');
			set(groot,'defaultAxesTickLabelInterpreter','latex');  
			set(groot,'defaulttextinterpreter','latex');
			set(groot,'defaultLegendInterpreter','latex');
			tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
			title_let = ["(A)","(B)","(C)","(D)"];
			for j = 0:3 % use the 3 different types of models
				nexttile
				fmdl = mk_mdl(pp,j,bounds,lung_masks,i); 
				img = set_elem_background(fmdl);
				%pp.ROIs = ROIs;
				show_fem_enhanced(img,[0 1.04])
				if j <= 1 
					view([0 90])
				else
					view([90 90])
				end
				set(get(gca, 'Title'), 'String', title_let(j+1));
				set(gca,'FontSize',25)
				axis off
			end
			for j = 0:3 % use the 3 different types of models
				nexttile
				fmdl = mk_mdl(pp,j,bounds,lung_masks,i); 
				img = set_elem_background(fmdl);
				%pp.ROIs = ROIs;
				show_fem_enhanced(img,[0 1.04])
				if j <= 1 
					view([160 90+45])
				else
					view([90+160 90+45])
				end
				axis off
			end
			%if i == 4 % Only need the one figure 
				set(gcf,'Position',[465         500        1597         822])
				print(['../imgs/fem_models_PT0' num2str(i)], '-dsvg');
			%end
		end
	case '8'
% Show the FEM for each of the 3 different methods 
		%fname = ['~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432623_xueliu/2432623_xueliu_01.eit'];
		clf;
		set(gcf,'renderer','painters');
		set(groot,'defaultAxesTickLabelInterpreter','latex');  
		set(groot,'defaulttextinterpreter','latex');
		set(groot,'defaultLegendInterpreter','latex');
		% Do the basic meshes 
		tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact');
		for j =0:1
			nexttile
			fmdl = mk_mdl(pp,j,bounds,lung_masks,i); 
			img = set_elem_background(fmdl);
			%pp.ROIs = ROIs;
			show_fem_enhanced(img,[0 1.04])
			if j <= 1 
				view([0 90])
			else
				view([90 90])
			end
			set(get(gca, 'Title'), 'String', title_let(j+1));
			set(gca,'FontSize',25)
			axis off
		end
		for j = 0:2 % use the 3 different types of models
			nexttile
			fmdl = mk_mdl(pp,j,bounds,lung_masks,i); 
			img = set_elem_background(fmdl);
			%pp.ROIs = ROIs;
			show_fem_enhanced(img,[0 1.04])
			if j <= 1 
				view([160 90+45])
			else
				view([90+160 90+45])
			end
			axis off
		end
		set(gcf,'Position',[465         500        1397         822])
		print('../imgs/fem_models_generic', '-dsvg');
		
		for i=2:5 
			if i == 2
				ref_frame = 145-112;
			elseif i == 3
				ref_frame = 53-28;
			elseif i == 4
				ref_frame = 136-107;
			elseif i == 5
				ref_frame = 73-42;
			end
			CT_dir = [data_dir, '/PTS' num2str(i) '/SRS00002'];
			seg_data = ['~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS' num2str(i) '.mat'];

			% Get the custom data to generate the plots
			seg_data = load(seg_data);
			bounds = seg_data.segs.SRS00002.bounds;
			CT_fname = CT_dir;
			lung_masks = get_lung_masks(CT_fname,ref_frame);
			pp = [];
			clf;
			set(gcf,'renderer','painters');
			set(groot,'defaultAxesTickLabelInterpreter','latex');  
			set(groot,'defaulttextinterpreter','latex');
			set(groot,'defaultLegendInterpreter','latex');
			tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
			title_let = ["(A)","(B)","(C)","(D)","(E)","(F)"];
			for j = 2:2 % use the 3 different types of models
				nexttile(i)
				fmdl = mk_mdl(pp,j,bounds,lung_masks,i); 
				%pp.ROIs = ROIs;
				show_fem_enhanced(fmdl,[0 1.04])
				if j <= 1 
					view([0 90])
				else
					view([90 90])
				end
				set(get(gca, 'Title'), 'String', title_let(j+1));
				set(gca,'FontSize',25)
				axis off
			end
			for j = 2:2 % use the 3 different types of models
				nexttile(i+4)
				fmdl = mk_mdl(pp,j,bounds,lung_masks,i); 
				%pp.ROIs = ROIs;
				show_fem_enhanced(fmdl,[0 1.04])
				if j <= 1 
					view([160 90+45])
				else
					view([90+160 90+45])
				end
				axis off
			end
		end
		set(gcf,'Position',[465         500        1397         822])
		print('../imgs/fem_models_custom', '-dsvg');
	case '9'
		% Print the EA figure for each subject and then compute the metrics
		for i=2:5 
			if i == 4
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432623_xueliu/2432623_xueliu_01.eit';
				CT_dir = [data_dir, '/PTS4/SRS00002'];
				ref_frame = 136-107;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS4.mat';
			elseif i == 2
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432616_xueliu/2432616_xueliu-peeptit_02_001.eit';
				CT_dir = [data_dir, '/PTS2/SRS00002'];
				ref_frame = 145-112;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS2.mat';
			elseif i == 3
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_23/2432000_xueliu__0-12pep/2432000_xueliu__peep_12_001.eit';
				CT_dir = [data_dir, '/PTS3/SRS00002'];
				ref_frame = 53-28;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS3.mat';
			elseif i == 5
				fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_29/2394977_xueliu-peep0/23394977_xueliu-peep0_01.eit';
				CT_dir = [data_dir, '/PTS5/SRS00002'];
				ref_frame = 73-42;
				seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS5.mat';
			end
			[img_a,img_b,img_c,img_d] = ea_reconstruct_breaths(fname,seg_data,CT_dir,ref_frame,i);
			tiledlayout(2,5, 'Padding', 'none', 'TileSpacing', 'compact');
			ref_frame = 145-112;
			nexttile
			show_slices(img_a)
			set(get(gca, 'Title'), 'String', '(A)');
			set(gca,'FontSize',25)
			nexttile
			show_slices(img_b)
			set(get(gca, 'Title'), 'String', '(B)');
			set(gca,'FontSize',25)
			nexttile
			show_slices(img_c)
			set(get(gca, 'Title'), 'String', '(C)');
			set(gca,'FontSize',25)
			view([90 90])
			nexttile
			show_slices(img_d)
			set(get(gca, 'Title'), 'String', '(D)');
			set(gca,'FontSize',25)
			view([90 90])
			% Do the CT case as well  
		%keyboard
			[ct_img,c_of_m(i-1,1:2),ctD] = ct_c_of_m(CT_dir,ref_frame,seg_data);
			nexttile
			imshow(ct_img)
			set(get(gca, 'Title'), 'String', '(E)');
			set(gca,'FontSize',25)
			[c_of_m(i-1,3:4),D] = calc_c_of_m(img_a,0);
			nexttile
			b = barh(D,'stacked');
			xlim([-max(max(abs(D))) max(max(abs(D)))])
			xlim([-60 60])
			yticklabels({'1:10','11:20','21:30','31:40','41:50','51:60','61:70','71:80','81:90','91:100'})
			set(get(gca, 'YLabel'), 'String', 'Distance from model anterior');
			b(1).FaceColor =[49,130,189]/256; 
			b(2).FaceColor =[49,130,189]/256; 
			set(gca,'FontSize',25)
			[c_of_m(i-1,5:6),D] = calc_c_of_m(img_b,0);
			nexttile
			b = barh(D,'stacked');
			xlim([-max(max(abs(D))) max(max(abs(D)))])
			xlim([-60 60])
			b(1).FaceColor =[49,130,189]/256; 
			b(2).FaceColor =[49,130,189]/256; 
			set(gca,'FontSize',25)
			set(gca,'YTick', [])
			[c_of_m(i-1,7:8),D] = calc_c_of_m(img_c,1);
			nexttile
			b = barh(D,'stacked');
			xlim([-max(max(abs(D))) max(max(abs(D)))])
			xlim([-60 60])
			b(1).FaceColor =[49,130,189]/256; 
			b(2).FaceColor =[49,130,189]/256; 
			set(gca,'FontSize',25)
			set(gca,'YTick', [])
			set(get(gca, 'XLabel'), 'String', 'number of ventilated pixels');
			[c_of_m(i-1,9:10),D] = calc_c_of_m(img_d,1);
			nexttile
			b = barh(D,'stacked');
			xlim([-max(max(abs(D))) max(max(abs(D)))])
			xlim([-50 50])
			b(1).FaceColor =[49,130,189]/256; 
			b(2).FaceColor =[49,130,189]/256; 
			set(gca,'FontSize',25)
			set(gca,'YTick', [])
			nexttile
			set(gca,'FontSize',25)
			set(gca,'YTick', [])
			b = barh(ctD,'stacked');
			xlim([-60 60])
			b(1).FaceColor =[49,130,189]/256; 
			b(2).FaceColor =[49,130,189]/256; 
			set(gca,'YTick', [])
			set(gcf,'Position',[166         303        2180         947])
			set(gcf,'renderer','painters');
			set(gca,'FontSize',25)
			print(['../imgs/center_of_vent_PT0' num2str(i)], '-dsvg');
		end
		% Calcualte the error in the center of mass for each model and each EA breath
		diff_array = zeros(4,8);
		diff_array(:,[1,3,5,7]) =  c_of_m(:,[3,5,7,9]) -c_of_m(:,1);
		diff_array(:,[2,4,6,8]) =  c_of_m(:,[4,6,8,10])-c_of_m(:,2);
		dist_array = zeros(4);
		dist_array(:,1) = sqrt(diff_array(:,1).^2+diff_array(:,2).^2);
		dist_array(:,2) = sqrt(diff_array(:,3).^2+diff_array(:,4).^2);
		dist_array(:,3) = sqrt(diff_array(:,5).^2+diff_array(:,6).^2);
		dist_array(:,4) = sqrt(diff_array(:,7).^2+diff_array(:,8).^2);
		dist_array
		mean(dist_array)
		figure(7); clf;
		set(gcf,'renderer','painters');
		set(groot,'defaultAxesTickLabelInterpreter','latex');  
		set(groot,'defaulttextinterpreter','latex');
		set(groot,'defaultLegendInterpreter','latex');
		MeshType = repelem([{'Circular'}, {'Generic'}, {'Custom Extruded'}, {'Custom 3D lungs'}], [4 4 4 4])';  
		boxplot(dist_array,MeshType)
		set(gca,'FontSize',25)
		set(get(gca, 'YLabel'), 'String', 'Center of Mass error (pixels)');
		ax = gca;
    %ax.FontSize = 16; 
    	ax.TickLabelInterpreter = 'latex';
		keyboard
case '10'
	% Print the EA figure for each subject and then compute the metrics
	for i=2:5 
		if i == 4
			fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432623_xueliu/2432623_xueliu_01.eit';
			CT_dir = [data_dir, '/PTS4/SRS00002'];
			ref_frame = 136-107;
			seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS4.mat';
		elseif i == 2
			fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432616_xueliu/2432616_xueliu-peeptit_02_001.eit';
			CT_dir = [data_dir, '/PTS2/SRS00002'];
			ref_frame = 145-112;
			seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS2.mat';
		elseif i == 3
			fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_23/2432000_xueliu__0-12pep/2432000_xueliu__peep_12_001.eit';
			CT_dir = [data_dir, '/PTS3/SRS00002'];
			ref_frame = 53-28;
			seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS3.mat';
		elseif i == 5
			fname = '~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_29/2394977_xueliu-peep0/23394977_xueliu-peep0_01.eit';
			CT_dir = [data_dir, '/PTS5/SRS00002'];
			ref_frame = 73-42;
			seg_data = '~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS5.mat';
		end
		[img_a,img_b,img_c,img_d] = ea_reconstruct_breaths(fname,seg_data,CT_dir,ref_frame,i);
		tiledlayout(2,5, 'Padding', 'none', 'TileSpacing', 'compact');
		ref_frame = 145-112;
		nexttile
		show_slices(img_a)
		set(get(gca, 'Title'), 'String', '(A)');
		set(gca,'FontSize',25)
		nexttile
		show_slices(img_b)
		set(get(gca, 'Title'), 'String', '(B)');
		set(gca,'FontSize',25)
		nexttile
		show_slices(img_c)
		set(get(gca, 'Title'), 'String', '(C)');
		set(gca,'FontSize',25)
		view([90 90])
		nexttile
		show_slices(img_d)
		set(get(gca, 'Title'), 'String', '(D)');
		set(gca,'FontSize',25)
		view([90 90])
		for j=1:4
			if j==1
				img = img_a;
			elseif j==2
				img = img_b;
			elseif j==3
				img = img_c;
			elseif j==4
				img = img_d;
			end

			% Identify the lung region
			lung_pixels = img.elem_data(round([img.fwd_model.mat_idx{2}/2; ...
			img.fwd_model.mat_idx{3}]/2));
			%lung_pixels = imgr.elem_data(round([imdl.rec_model.mat_idx{2}; ...
			%imdl.rec_model.mat_idx{3}]/2));
			% using the difference  image and lung regions calculate the GI 
			GI(i-1,j) = abs(sum(abs(lung_pixels-median(lung_pixels)))/sum(lung_pixels));
		end
	end
keyboard
	figure(7); clf;
	set(gcf,'renderer','painters');
	set(groot,'defaultAxesTickLabelInterpreter','latex');  
	set(groot,'defaulttextinterpreter','latex');
	set(groot,'defaultLegendInterpreter','latex');
	MeshType = repelem([{'Circular'}, {'Generic'}, {'Custom Extruded'}, {'Custom 3D lungs'}], [4 4 4 4])';  
	boxplot(dist_array,MeshType)
	set(gca,'FontSize',25)
	set(get(gca, 'YLabel'), 'String', 'Center of Mass error (pixels)');
	ax = gca;
%ax.FontSize = 16; 
	ax.TickLabelInterpreter = 'latex';
	keyboard


end




end

%function [img,vh,vi] = set_elem_background(fmdl)
%	img = mk_image(fmdl,1);
%	%vh = fwd_solve(img);
%	%   img.elem_data(vertcat(fmdl.mat_idx{2:3}))= 2;
%	if numel(fmdl.mat_idx)>1
%		img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 0.2; % lungs 
%	else
%		test = 5
%		keyboard
%	end
%	%m_frac = elem_select(fmdl, @(x,y,z) ...
%	%	(x.^2)/0.6 + (y-0.25).^2/0.55 < 1);
%	%img.elem_data = img.elem_data + m_frac*2;
%	%keyboard
%	%vi = fwd_solve(img);
%end
