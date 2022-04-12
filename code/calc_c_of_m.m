function [vent_cent,D] = calc_c_of_m(img,rot_flag)
	% Given an image calculage the center of mass of the ventilated region
	% the ventilated region is any region larger than 50% of the max ventilation signal
	img_arr = calc_slices(img);
	img_arr(img_arr>0)=0;
	img_min = min(min(img_arr));
	img_arr = img_arr/-img_min;
	img_arr(img_arr>-0.35)=0; % This may need a bit of changing! let's see
	if rot_flag == 1
		img_arr = rot90(rot90(rot90(img_arr)));
	end
	bound_array = abs(isnan(img_arr)-1);
	[I,J] = find(bound_array);
	max_bounds = max(I);
	min_bounds = min(I);
	bound_cent = [mean(J),mean(I)];
	% Compared to the boundary center of mass the vent center of mass is...
	vent_array = img_arr;
	vent_array(vent_array<0) = 1;
	vent_array(isnan(vent_array)) = 0;
	[I,J] = find(vent_array);
	vent_cent = [mean(J),mean(I)];
	img_arr(floor(bound_cent(2)),floor(bound_cent(1))) = 1;
	img_arr(floor(vent_cent(2)),floor(vent_cent(1))) = 1;
	%subplot(211)
	%show_slices(img_arr)
	% For each 5 slices - find amount on right and amount on left
	ranges = floor(linspace(min_bounds,max_bounds,11));
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
	%subplot(212)
	%barh(flipud(y),'stacked')
	%xlim([-max(max(abs(y))) max(max(abs(y)))])
	%keyboard
end