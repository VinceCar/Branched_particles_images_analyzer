%% First labelling and characterization of islands
label_islands=bwlabel(bw);
islands_index=(1:max(label_islands(:)))';

% report index, x,y , how many particles, how many rings, if it is
% contained in another patch, and time
islands_in_frame=[islands_index,zeros(max(islands_index),5),ii*ones(max(islands_index),1)];
support_matrix=regionprops(label_islands,'Centroid'); % momentary storage of center of masses
support_matrix=[support_matrix(:).Centroid];
islands_in_frame(:,2)=support_matrix(1:2:end);
islands_in_frame(:,3)=support_matrix(2:2:end);
islands_in_frame(:,6)=~ismember(islands_index,unique(imdilate(1-imfill(bw,'holes'),ones(3)).*label_islands));
