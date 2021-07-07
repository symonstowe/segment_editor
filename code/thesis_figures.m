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
		end
		[e1,e2,e3,e4,e5] = step_seg_thorax(imgs,frame);
		% Steps to segment the external boundary
		f = figure(1); clf;
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

	case '3'
		% Show the 2D basic mesh 

	case '4'
		% Show the 2.5 D mesh 

		% A) with the 5 slices 

		% B) with 25 slices
	case '5'
		% reconstrusionts with with 3 meshes 
		% A) regular B) custom C) 2.5
		fname = ['~/Dropbox/projects/2021/CT_2_MESH/EIT-CT/2020_03_24/2432623_xueliu/2432623_xueliu_01.eit'];
		CT_dir = [data_dir, '/PTS4/SRS00002'];
		ref_frame = 136-107;
		seg_data = ['~/Dropbox/projects/2021/CT_2_MESH/segment_editor/data/PTS4.mat'];
		reconstruct_breaths(fname,seg_data,CT_dir,ref_frame)
	case '6' 
		% Reconstruct with the jpg slices 	
		jpg_dir = [data_dir '/JPG-EIT_data/1680127/'];
		% Identify a folder 
		img = imread([jpg_dir '1.JPG']);
		keyboard


	end
end
