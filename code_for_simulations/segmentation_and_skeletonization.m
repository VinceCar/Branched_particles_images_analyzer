%% Attention to work on!!! Compare watersheds of enanchend_holes with the gradients (bw+2*(bwskel(imclose((watershed(-imadjust((imfilter(accumulator/max(accumulator(:)),hollow_object))),8)==0),minor_object)))+4*bwskel(imclose((watershed(imadjust(imgradient(imfilter(accumulator/max(accumulator(:)),hollow_object))),8)==0),minor_object)))


function [skel,bw]=segmentation_and_skeletonization(image_to_process,reference,dim_y_o,dim_x_o,dim_y,dim_x,diameter,thr_accumulator_reference)
support_matrix=image_to_process-min(image_to_process(:)); % disposable matrix used for this package, in order to not rewrite 
accumulator=ifft2(fft2(padarray(padarray(support_matrix, [ceil(size(reference)/2)],'symmetric','pre'),[floor(size(reference)/2)],'symmetric','post')).*fft2(padarray(reference, [(size(support_matrix))],0,'post')));
accumulator=accumulator(size(reference,1)+1:end,size(reference,2)+1:end);
thr=multithresh(accumulator); % initial otsu thresholding
accumulator=accumulator.*(accumulator>thr)+(accumulator<=thr).*imfilter(accumulator,ones(9)/81); % perform average filtering only on the background, needed to smooth the watershed
skel=double(((watershed(round(accumulator./(max(accumulator(:)+1e-6))*32),8)==0))); % perform a first skeletonization based on the watershed transform. note that the ridge lines will be part of both foreground and background since the background is not uniform!
thr=multithresh(accumulator(skel>0)); % more precise thresholding since now the populations are more balanced
bw=((accumulator.*(imfill(accumulator.*skel>thr,'holes')))); % first segmentation, it still has background in it, hopefully now in a bilanced sort
thr_accumulator=multithresh(bw(bw>0)); % final optimal threshold
if abs(thr_accumulator_reference/thr_accumulator-1)>0.10 && (thr_accumulator_reference>0) % if there is a reference, and we are too far away from it
    thr_accumulator=thr_accumulator_reference; % overwrite the threshold. we avoid doing this on the first place just because there may be some misalignments between images / different areas in the same folder
end
bw=accumulator>thr_accumulator; % actual binarization
bw=imgaussfilt(bw+0,3)>0.8; % smooth the mask
bw=bwareaopen(bw,120,8)+0; % remove small stuff... 120 is hardcoded and emperically found
%% preparation for filling of scars and removal of full area
small_object=strel('disk',round(diameter*2),8); % a circle slightly smaller than a particle
big_object=strel('disk',round(diameter*4),8); % a circle slightly bigger than a particle
half_difference_in_size=(size(big_object.Neighborhood)-size(small_object.Neighborhood))/2; % variable to relate the sizes of big and small disks
hollow_object=big_object.Neighborhood./sum(big_object.Neighborhood(:))-[zeros(half_difference_in_size(1),size(big_object.Neighborhood,2));zeros(size(small_object.Neighborhood,2),half_difference_in_size(2)),small_object.Neighborhood./sum(small_object.Neighborhood(:)),zeros(size(small_object.Neighborhood,2),half_difference_in_size(2));zeros(half_difference_in_size(1),size(big_object.Neighborhood,2))];
support_matrix=bw>0;
%support_matrix=1-(imfilter(1-(imfilter(support_matrix,small_object.Neighborhood+0)>0),small_object.Neighborhood+0)>0);
support_matrix=imclose(support_matrix,small_object)+0;
%% removal of full area
accumulator=accumulator./(max(accumulator(:)+1e-6));
accumulator=accumulator-prctile(accumulator(:),1);
accumulator(accumulator<0)=0;
thr_up=(prctile(accumulator(:),99)+1e-6);
accumulator(accumulator>thr_up)=thr_up;
accumulator=accumulator*1/thr_up;
accumulator(size(accumulator,1)+size(hollow_object,1),size(accumulator,2)+size(hollow_object,2))=0; % 0 pad for fourier transform
hollow_object(size(accumulator,1),size(accumulator,2))=0; % 0 pad for fourier transform
enhanced_holes=ifft2(fft2(accumulator).*fft2(hollow_object)); % enhancing the holes through normalization, contrast adjustment and filtering using an hollow kernel
enhanced_holes=enhanced_holes(floor(size(big_object.Neighborhood,1)/2)+1:end,floor(size(big_object.Neighborhood,2)/2)+1:end);
enhanced_holes=enhanced_holes(1:(end-floor(size(big_object.Neighborhood,1)/2)-1),1:(end-floor(size(big_object.Neighborhood,2)/2)-1));

enhanced_holes=(imgaussfilt(enhanced_holes,diameter,'Padding','symmetric')); % additional gaussian filtering with a kernel that has a deviation equal to 1/3 of the diameter of a particle (in enhanced image)
enhanced_holes=enhanced_holes-min(enhanced_holes(:)); % put min to 0
enhanced_holes=round(enhanced_holes./(max(enhanced_holes(:))+1e-6).*50); % discretize the gray levels in the image. 50 was an empirical value for optimal skeletonization using watershed
preliminary_skeleton=(watershed(-enhanced_holes,8)==0); % candidate areas are holes that contain only background
finding_the_candidates_to_empty=bwlabel(imopen(imfill(preliminary_skeleton.*support_matrix,'holes'),ones(5)).*(1-preliminary_skeleton),4); % smoothing and identification of holes to make
region_of_interests=regionprops(finding_the_candidates_to_empty,(finding_the_candidates_to_empty.*(bw)>0)+0,'Area','MaxIntensity','MinIntensity','WeightedCentroid'); % calculate useful infos about the holes
index_to_empty=find(([region_of_interests(:).MinIntensity]==1).*([region_of_interests(:).Area]>diameter)); % region to empty satisfy this criteria.

to_remove=zeros(size(support_matrix)); % intialization of a matrix that is going to contain what to remove form the original segmentation to make the holes
if ~isempty(index_to_empty)
    support_matrix_2=[region_of_interests(index_to_empty').WeightedCentroid];
    to_remove(sub2ind(size(support_matrix),round(support_matrix_2(:,2:2:end)'),round(support_matrix_2(:,1:2:end)')))=1; % filling it up
    to_remove=(imfilter(to_remove+0,small_object.Neighborhood+0)>0)+0; % dilation (for size reason)
end
support_matrix=support_matrix-to_remove; % subtract from the support matrix
bw=bw.*support_matrix; % apply modification to the black and white

%% filling scars
minor_object=strel('disk',round(diameter*3/4),8); % object used during filling of the scars

additional_pieces=bwskel((watershed((imgaussfilt(bw+0,diameter,'Padding','symmetric')))==0).*(1-bw)>0); % first 4-connectivity skeleton first passage smooth the new bw, then apply watershed, then force 4 connectivity and consider only the lines that are included in the  mask
endp=bwmorph(additional_pieces,'endpoints'); % find endpoints
branchp=bwmorph(additional_pieces,'branchpoints'); % find branchpoints
distances=(bwdistgeodesic(additional_pieces>0,endp)); % geodesic distance transforms from the endpoints
distances(isnan(distances))=0; % put NaN value to 0, they are background/not connected to endpoints anyway
distances(distances==Inf)=0; % same as above for Inf
distances=distances.*(1-branchp); % remove branchpoints
additional_pieces=bwlabel(additional_pieces.*(1-branchp)); % identify different subsegments
distances=regionprops(additional_pieces>0,distances,'MaxIntensity'); % register the maximum distance in each segment
index_to_keep=find([distances(:).MaxIntensity]<=round(diameter/2/0.33)); % keep segments that are only half the theoretical size of a particle
additional_pieces=ismember(additional_pieces,index_to_keep); % keeping
bw=imdilate(additional_pieces,ones(9)).*imfill(bw,'holes')+bw>0; % include the new segments in the original bw, but only if there was background in a hole
support_matrix=imfill(bw,'holes'); % support matrix to perform a smother filling
bw=bw+(imdilate(additional_pieces,ones(11)).*(1-support_matrix)); % final smother integration of the new segments in the bw

%% make it smoother
bw=imresize(bw+0,round([dim_y/2,dim_x/2])); % we will make it smooth through a resizing... 0.3 is an arbitrary threshold selected because we are reducing the images by a factor of 0.3
bw=(imresize(bw.*(bw>0.25),[dim_y,dim_x])>0.5)+0;

%% Skeletonization

% % support_matrix=[zeros(size(bw)), flipud(bw),zeros(size(bw));fliplr(bw),bw,fliplr(bw);zeros(size(bw)), flipud(bw),zeros(size(bw))]; % to avoid finite size artefacts we tile the image, mirroring it along the edges... but on smaller versions
% support_matrix=[zeros(size(bw(1:10,1:10))), flipud(bw(1:10,:)),zeros(size(bw(1:10,1:10)));fliplr(bw(:,1:10)),bw,fliplr(bw(:,end-9:end));zeros(size(bw(1:10,1:10))), flipud(bw(end-9:end,:)),zeros(size(bw(end-9:end,1:10)))]; % to avoid finite size artefacts we tile the image, mirroring it along the edges... but on smaller versions
% 
% skel=bwmorph(support_matrix,'skel',Inf)+bwskel(support_matrix>0)>0; % combine 8-connectivity and 4-connectivity (trick used to have sharper angles)
% %skel=skel((size(skel,1)/3+1):(2*size(skel,1)/3),(size(skel,2)/3+1):(2*size(skel,2)/3)); % take only the skelton in the central frame
% skel=skel(11:end-10,11:end-10);
bw(1,:)=0;
bw(:,1)=0;
bw(end,:)=0;
bw(:,end)=0;
support_matrix=bw;
skel=bwmorph(support_matrix,'skel',Inf)+bwskel(support_matrix>0)>0; % combine 8-connectivity and 4-connectivity (trick used to have sharper angles)

support_matrix=bwlabel(imfill(skel,'holes')-skel,4)+0; % fill the spaces between the two skeletons
to_not_use=unique((1-bw).*support_matrix); % avoid using the ones that are part of the background or holes through this trick. label them and then if they are not entirely on the white delete them
to_not_use(to_not_use==0)=[];
support_matrix=(support_matrix>0)-ismember(support_matrix,to_not_use);
skel=bwmorph(support_matrix+skel>0,'skel',Inf); % skeletonize and correct for minor errors
skel=bwskel(skel);
% skel=skel.*(imfilter(skel,[0 1 0; 1 2 1; 0 1 0])<5);

%% Pruning
% this is a sequential pruning
skel=(bwmorph(skel,'spur',31));
skel=bwskel(skel);
end
