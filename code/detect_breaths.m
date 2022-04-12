function [breath,ea_breath,pk_loc] = detect_breaths(dd)
        breath_sig = sum(dd);
        breath_sig = rescale(breath_sig,0,1);
        prominence = 0.25; 
        pks = islocalmax(breath_sig,'MinProminence',prominence,'FlatSelection','first');
        p = prominence;
        while isempty(pks)
            p = p-0.05;
            pks = islocalmax(breath_sig,'MinProminence',p,'FlatSelection','first');
        end
        pks = find(pks);
        % Find the troughs using the signal min (end expiration)
        trghs_sig = -breath_sig;
        trghs = islocalmax(trghs_sig,'MinProminence',prominence,'FlatSelection','first');
        p = prominence;
        while isempty(pks)
            p = p-0.05;
            trghs = islocalmax(trghs_sig,'MinProminence',p,'FlatSelection','first');
        end
        trghs = find(trghs);
        %% Find the breath troughs using findchangepts
        %idx = findchangepts(breath_sig,'Statistic','linear','MinThreshold',0.3); % TODO needs to be adjusted base on Fs eventually
        %dx = diff(idx);
        %dy = diff(breath_sig(idx));
        %slope = dy./dx;
        %breathTroughs = idx(slope > 0.0001); % remove negative slopes
        % Remove the first trough because of the filtering
        trghs(1) = [];
        % discard peaks that don't have at least one trough before or after
        pks = pks(pks<max(trghs) & pks>min(trghs));
        breath = [];
%clf
%plot(breath_sig)
%hold on
%for i=1:numel(pks)
%    xline(pks(i))
%    xline(trghs(i))
%end
%keyboard
        if isempty(pks)
            keyboard
        end
        for i=1:numel(pks)
            breath(i).pk   = pks(i);%   = pks(pks>breathTroughs(i) & pks<breathTroughs(i+1));
            trghs_diff = trghs-pks(i);
            low_diffs  = trghs_diff;
            high_diffs = trghs_diff;
            low_diffs(low_diffs>0)    = inf;
            high_diffs(high_diffs<0)  = inf;
            [~,T1] = min(abs(low_diffs));
            [~,T2] = min(high_diffs);
            breath(i).trgh(:) = [trghs(T1), trghs(T2)];
            length_inhale(i) = breath(i).pk - breath(i).trgh(1);
            length_exhale(i) = breath(i).trgh(2) - breath(i).pk; 
        end
        % Remove high outliers in breath lengths - assume something weird is going on...
        while max(length_exhale)>1.5*mean(length_exhale) 
            [~,reject] = max(length_exhale);
            length_exhale(reject) = [];
            length_inhale(reject) = [];
            breath(reject) = [];
            pks(reject) = [];
        end
        % For ensemble average use the max length for inhale and exhale
        if size(dd,2)<pks(end)+max(length_exhale)
            pks(end) = [];
            breath(end) = [];
        end
        if size(dd,2)<pks(end)+max(length_exhale)
            pks(end) = [];
            breath(end) = [];
        end
        if size(dd,2)<pks(end)+max(length_exhale)
            pks(end) = [];
            breath(end) = [];
        end
        if pks(1)-max(length_inhale) < 0
            pks(1) = [];
            breath(1) = [];
        end
        for i=1:numel(pks)
            ea_dat(:,:,i) = dd(:,pks(i)-max(length_inhale):pks(i)+max(length_exhale));
        end
        ea_breath = mean(ea_dat,3);
        pk_loc = max(length_inhale)+1;
        breath_var = max(length_inhale)-min(length_inhale);