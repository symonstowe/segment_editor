function mask = get_lung_masks(CT_dir,ref_frame)
	% TODO could we have just used poly2mask?!?!?!?!?!?
	% Segment the lung frame for everywhere except the frames too close to the edge
	% create a meshgrid of the masks to return a point cloud for alphashape generation
	% I need to have all 5 slices corrected to do this properly...
	bot_lim =  10;
	top_lim =  10;
	fileinfo = dir(fullfile(CT_dir, '**', '*.DCM'));
	filenames = fullfile({fileinfo.folder}, {fileinfo.name});
	for i=1:length(filenames)
		imgs(:,:,i) = mat2gray(dicomread(filenames{i}));
	end
	c=0;
	for i=ref_frame-0:ref_frame+5
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
end