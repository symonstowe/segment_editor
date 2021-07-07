function [E1, E2, E3, E4, E5, chest_cavity, lung_healthy, total_lung]=step_seg_thorax(IMG,frame, NumberSlices,thresh_adjust)
	% IMG is a series of loaded CT images (512,512, number of frames)
	% Frame - a frame selected from the thorax
	% Number of slices (should be an odd number, used to define the rib region surrounding the slice)
	% E series is external boundary segmentation steps
	% E1 - adjusted slice 
	% E2 - ventilated lung tissue set to 0
	% E3 - Iobrcbr - dialated and reconstructed 
	% E4 - Binarized 
	% E5 - Body final
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
		E1t = J1;
		% Improve contrast level based on the image
		dist = imhist(J1); 
		dist(1:10) = 0;
		[~,locs] = findpeaks(-dist,'MinPeakProminence',100,'MinPeakDistance',15);
		locs = locs/255;
		level = 0.5; % Make sure it is somewhat close to the original code...
		[~,loc]=min(abs(locs-level));
		level = locs(loc); % This variable naming is weird - just finding threshold for lung tissue
		J1 = imadjust(J1,[level,1]); % increase contrast but better
		E2t = J1;
		%SEGMENT THE BODY BOUNDARY
		se = strel('disk', 20);
		Ie = imerode(J1, se);
		se = strel('disk', 10);
		Iobr = imreconstruct(Ie, J1);
		Iobrd = imdilate(Iobr, se);
		Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
		Iobrcbr = imcomplement(Iobrcbr); % This has a pretty defined heart region %%%%%% TODO!!
		E3t = Iobrcbr;
		bw = imbinarize(Iobrcbr, 0.5*max(max(Iobrcbr)));
		E4t = bw;
		bw = imfill(bw,'holes');
		bw = imclose(bw,strel('disk',2));
		body_final = imopen(bw,strel('disk', 5));
		E5t = body_final;
		if ((NumberSlices-1)/2-i+1) == 0 
		    % Save the external boundary
		    external = bwareafilt(body_final,1); % Largest object
		    adjusted_slice = J1;
		    lung_heart = Iobrcbr;
		    E1 = E1t; 
		    E2 = E2t; 
		    E3 = E3t; 
		    E4 = E4t; 
		    E5 = E5t;
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