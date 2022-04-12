function [ct_img,vent_cent,D] = ct_c_of_m(ct_dir,frame,seg_fname)
	% Load all the CT images
	fileinfo = dir(fullfile(ct_dir, '**', '*.DCM'));
	filenames = fullfile({fileinfo.folder}, {fileinfo.name});
	for i=1:length(filenames)
		imgs(:,:,i) = mat2gray(dicomread(filenames{i}));
	end
	IMG = imgs;
	seg_data = load(seg_fname);
	bounds = seg_data.segs.SRS00002.bounds;
	[adjusted_slice, external, chest_cavity, lung_healthy, total_lung]=seg_thorax(IMG,frame);
	boundary_mask = poly2mask(bounds{3}.exterior(:,2),bounds{3}.exterior(:,1),512,512);
	bounds = boundary_mask;
	[~,y] = find(bounds);
	remove_cols = [1:min(y),max(y):length(bounds)];
	remove_rows = [1:floor(numel(remove_cols)/2),length(bounds)-ceil(numel(remove_cols)/2):length(bounds)];
	bounds(:,remove_cols) = [];
	bounds(remove_rows,:) = [];
	lung_healthy(:,remove_cols) = [];
	lung_healthy(remove_rows,:) = [];
	bound_array = imresize(bounds,[50 50]);
	vent_array = imresize(lung_healthy,[50 50]);
	% Center of mass of the boundary
	[I,J] = find(bound_array);
	max_bounds = max(I);
	min_bounds = min(I);
	bound_cent = [mean(J),mean(I)];
	% Center of mass of the lungs
	[I,J] = find(vent_array);
	vent_cent = [mean(J),mean(I)];
	ranges = floor(linspace(min_bounds,max_bounds,11));
	y= [];
%keyboard
	for i=1:10
		%vals = (i-1)*10+1:i*10;
		vals = ranges(i):ranges(i+1);
		img_temp = vent_array(vals,1:floor(bound_cent(2)));
		[~,R] = find(img_temp);
		if isempty(R)
			y(i,1) = 0;
		else
			y(i,1) = -length(R);% Sum on the left(actual right)
		end
		img_temp = vent_array(vals,ceil(bound_cent(2):end));
		[~,R] = find(img_temp);
		if isempty(R)
			y(i,2) = 0;
		else
			y(i,2) = length(R);% Sum on the left(actual right)
		end
	end
	D = flipud(y);
	% make the CT image	
	lungs = repmat(uint8(vent_array), 1, 1, 3);
	lungs(:,:,1) = lungs(:,:,1)*49;
	lungs(:,:,2) = lungs(:,:,2)*130;
	lungs(:,:,3) = lungs(:,:,3)*189;
	bounds = 256*repmat(uint8(bound_array), 1, 1, 3);
	[x,y] = find(vent_array);
	for i=1:numel(x)
		bounds(x(i),y(i),:) = lungs(x(i),y(i),:);
	end
	ct_img = bounds;
	%imshow(ct_img)
end