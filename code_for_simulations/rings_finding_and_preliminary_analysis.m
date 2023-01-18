% detection of the rings and first classification in convex or concave, and
% nearby concave
shrinked_skel=bwareaopen(bwmorph(skel,'shrink',Inf),2);
labelled_rings=imfill(shrinked_skel,'holes'); % filling the skeleton, we store it in this variable because we do not need it for long
branchmap=bwmorph(skel,'branchpoints').*imdilate(labelled_rings-skel,ones(3)); % here for speed and memory reasons... it is useful afterwards
labelled_rings=double(bwlabel(labelled_rings-shrinked_skel,4)); % rings are labelled, mainly needed to avoid the 8 connectivity forced by regionprops
rings_index=unique(labelled_rings(labelled_rings>0)); % useful to keep the identity of the holes within a variable
labelled_rings=labelled_rings.*(imopen((labelled_rings>0)-imdilate(shrinked_skel,ones(3))>0,ones(3))); % rings are labelled, mainly needed to avoid the 8 connectivity forced by regionprops
perimiters_holes=bwmorph((bwperim(labelled_rings>0,4)),'skel',Inf);
region_property=regionprops(labelled_rings,'Centroid','Area','ConvexArea'); % calculate the occupied area and the occupied convex area

rings_index=rings_index(:); % to make sure rings_index is a column vector
index_conc=[region_property(:).ConvexArea]./[region_property(:).Area]>1.1; % 1.1 is an empirically found threshold, if the convex area diverge too much from the real area, then the polygon is concave

support_matrix=labelled_rings; % this is a support matrix that we use to perform our calculations
support_matrix(imfill(imdilate(skel,ones(3))>0,'holes')-(labelled_rings>0)>0)=Inf; % to performe erosion
matrix_sorted=imerode(support_matrix,ones(11)).*(labelled_rings>0); % the name is not indicative of the function, we just call it this way to discard (exploiting garbage collector)
support_matrix=labelled_rings; % this is a support matrix that we use % to perform our calculations
support_matrix(imfill(imdilate(skel,ones(3))>0,'holes')-(labelled_rings>0)>0)=-Inf; % to performe erosion
support_matrix=imdilate(support_matrix,ones(11)).*(labelled_rings>0);
support_matrix=support_matrix.*(support_matrix~=labelled_rings);
matrix_sorted=matrix_sorted.*(matrix_sorted~=labelled_rings);
matrix_sorted=(matrix_sorted.*(matrix_sorted>support_matrix)+support_matrix.*(support_matrix>matrix_sorted)).*perimiters_holes; % the name is not indicative of the function, we just call it this way to discard (exploiting garbage collector)% now each polygon has been labelled with the confining polygons
matrix_sorted=(matrix_sorted+1).*perimiters_holes; % add 1 so that you do not need to go crazy afterwards
support_matrix=matrix_sorted;
support_matrix(matrix_sorted==0)=-Inf;
vertex_matrix=imdilate(support_matrix,ones(3))>support_matrix;
vertex_matrix(perimiters_holes==0)=0;
support_matrix=matrix_sorted;
support_matrix(matrix_sorted==0)=Inf;
vertex_matrix=vertex_matrix+(imerode(support_matrix,ones(3))<support_matrix);

vertex_matrix(perimiters_holes==0)=0;
support_matrix=bwlabel(perimiters_holes.*imdilate(branchmap,ones(17)));
support_matrix=support_matrix.*(1-ismember(support_matrix,support_matrix(vertex_matrix>0)))>0;
vertex_matrix=vertex_matrix+support_matrix;
vertex_matrix=bwmorph(vertex_matrix,'shrink',Inf); % where to store the verteces
matrix_sorted=matrix_sorted-1;
index_nearby_concave=false(size(index_conc)); 

index_nearby_concave(labelled_rings(logical(ismember(matrix_sorted(:),[0;find(index_conc)']))))=true;
index_nearby_concave(index_conc)=false;

%% preparing what is needed to recognize the polygons
convex_perimiter_rings=bwperim(bwconvhull(ismember(labelled_rings,find(~index_conc)).*(labelled_rings>0),'objects')).*perimiters_holes+vertex_matrix>0; % this is a box reporting candidate points to become verteces. it already performs a reduction  based on the concept that real verteces would be part of the original match

convex_perimiter_rings=convex_perimiter_rings.*labelled_rings;
%% applying the Ramer-Douglas-Pucker algorithm, starting from the convex polygons and branches if present
index_to_use=rings_index(~index_conc);
index_to_use([region_property(~index_conc).ConvexArea]==0)=[];
if ~isempty(index_to_use)
    [vertex_matrix]=apply_RDP(vertex_matrix+0,convex_perimiter_rings,index_to_use,diameter);
    support_matrix=bwlabel(perimiters_holes.*imdilate(vertex_matrix>0,ones(17)));
    support_matrix=support_matrix.*(1-ismember(support_matrix,support_matrix(vertex_matrix>0)))>0;
    vertex_matrix=(vertex_matrix>0)+support_matrix;
    vertex_matrix=bwmorph(vertex_matrix,'shrink',Inf); % where to store the verteces
end
% note even if the variable is called convex_perimeter now we are dealing
% with concave perimiter

convex_perimiter_rings=bwlabel(bwperim(ismember(labelled_rings,find(index_conc)).*perimiters_holes+vertex_matrix).*(ismember(labelled_rings,find(index_conc)))>0); % this is a box reporting candidate points to become verteces. it already performs a reduction  based on the concept that real verteces would be part of the original match
% the label will be the different segments we can divide them into! to
% avoid repetitions, let's add the maximum id + 1 so that we can create a
% new index
convex_perimiter_rings(convex_perimiter_rings>0)=convex_perimiter_rings(convex_perimiter_rings>0)+max(rings_index);
convex_perimiter_rings(convex_perimiter_rings==0)=-1;
convex_perimiter_rings(convex_perimiter_rings>0)=convex_perimiter_rings(convex_perimiter_rings>0)+matrix_sorted(convex_perimiter_rings>0); % assign an identity to the different segments
rings_index_2=unique(convex_perimiter_rings);
rings_index_2(rings_index_2==-1)=[];
if ~isempty(rings_index_2)
    [vertex_matrix]=apply_RDP_to_segments(vertex_matrix+0,convex_perimiter_rings,rings_index_2,diameter);
    vertex_matrix=bwmorph(vertex_matrix,'shrink',Inf); % by default we will have 2 points, let's reduce to 1
end
if ~isempty(rings_index)
	type_of_rings=(histcounts(vertex_matrix(labelled_rings>0).*labelled_rings(labelled_rings>0),(min(rings_index)-0.5):(max(rings_index)+0.5)))'; % identify polygons based on verteces
end
% now we will use convex perimiter to locate where the segments are
branchmap=bwmorph(skel,'branchpoints');
convex_perimiter_rings=bwlabel(bwareaopen(skel-imopen(imfill(skel,'holes'),ones(3))-branchmap>0,2));
rings_index_2=unique(convex_perimiter_rings);
rings_index_2(rings_index_2==0)=[];
if ~isempty(rings_index_2)
    [vertex_matrix]=apply_RDP_to_segments(vertex_matrix+0,convex_perimiter_rings,rings_index_2,diameter);
    vertex_matrix=vertex_matrix.*(1-imdilate(branchmap+bwmorph(skel,'endpoints'),ones(11)))+branchmap+bwmorph(skel,'endpoints');
    vertex_matrix=bwmorph(imdilate(vertex_matrix,ones(13)).*skel>0,'shrink',Inf);
end
if ~isempty(rings_index)
    type_of_rings(index_conc)=-type_of_rings(index_conc);
end

%% Update of lablled_rings in order to make it correspond to holes nearby the skeleton so that later it is easy to track what particles are part of the holes
labelled_holes=bwlabel(imfill(skel,'holes')-skel,4);
labelled_rings=double(labelled_holes);
region_property=regionprops(labelled_rings,'Centroid','Area','ConvexArea'); % calculate the occupied area and the occupied convex area
support_matrix=[region_property(:).Centroid];
%% find out which level from the boarder the holes are
level_to_store=zeros(length(type_of_rings),1); % store a level for each ring
level=0; % the starting level is the background (already identified) or non convex polygon
matrix_level=labelled_rings; % matrix that we will update with the levels
matrix_level(ismember(matrix_level,find(type_of_rings<0)))=0; % make concave polygons to level 0
matrix_sorted(ismember(matrix_sorted,find(type_of_rings<0)))=0;
original_max_level=max(matrix_level(:)); % variable useful for a technique we apply in the following
max_matrix_level=original_max_level;
while level<max(matrix_level(:))-1 % to avoid infinite loop we will do this research until the last original label is equal to the actual level to look into
    level=level+1; % actual level we are interested in
    max_matrix_level=max_matrix_level+1; % we do a circular shift if a label that we still need to indentify is equal to the level of interest
    matrix_level(matrix_level==level)=max_matrix_level; % do the subs
    matrix_sorted(matrix_sorted==level)=max_matrix_level; % do the subs
    on_boarder=find(logical(ismember(rings_index,unique(labelled_rings(matrix_sorted==level-1))).*(level_to_store==0)));
    on_boarder(on_boarder==level)=max_matrix_level;
    matrix_level(ismember(matrix_level,on_boarder))=level;
    matrix_sorted(ismember(matrix_sorted,on_boarder))=level;
    
    on_boarder=mod(on_boarder,original_max_level);
    on_boarder(on_boarder==0)=original_max_level;
    level_to_store(on_boarder)=level;
end
if level>0
    %image_by_grains=double(watershed(-imdilate(matrix_level,ones(3)))); v1
    % v2
%     matrix_level=imdilate(matrix_level,ones(3));
%     matrix_level(bwmorph(matrix_level>0,'remove')>0)=-Inf;
%     image_by_grains=double(watershed(-imerode(matrix_level,ones(3))));
    % v3
    matrix_level=imclose(matrix_level,ones(3));
    matrix_level(matrix_level>0)=matrix_level(matrix_level>0)+1;
    matrix_level=matrix_level+imdilate(bw,ones(3)).*(matrix_level==0);
    matrix_level(bwmorph(matrix_level>0,'remove')>0)=-Inf;
    image_by_grains=double(watershed(-imerode(matrix_level,ones(3))));
    image_by_grains(imclose(ismember(labelled_holes,find(index_conc)),ones(3)))=0;
    image_by_grains((labelled_holes>0)+bw==0)=0;
end
% the structure of ring is x,y, type, if it is in the boarder, which patch
% and time
if ~isempty(rings_index)
    ring_in_frame=[support_matrix(1:2:end)',support_matrix(2:2:end)',type_of_rings, level_to_store,zeros(length(type_of_rings),1) ,ii*ones(length(type_of_rings),1)];
end
function [vertex_matrix]=apply_RDP(vertex_matrix,labelled_reduced_rings,vector_ids,diameter)
prop=regionprops(labelled_reduced_rings,'PixelIdxList');
sizes=size(labelled_reduced_rings);
[branch_pos_Y,branch_pos_X]=find(vertex_matrix);

for pp=vector_ids' 
    [y,x]=ind2sub(sizes,[prop(pp).PixelIdxList]); % find the verteces
    [theta,rho]=cart2pol(x-mean(x),y-mean(y)); % transform in polar coordinates

    index_reduced=sortrows([(1:length(theta))',theta,rho],2); % sort based on angle
    x=x(index_reduced(1:end,1)); % apply sort
    y=y(index_reduced(1:end,1)); % apply sort
    [~,pos]=max(index_reduced(:,3)); % put the furthest point from the center of mass as first
    candidate_verteces=circshift([y,x],-pos); % initial candidates
    y=candidate_verteces(:,1);
    x=candidate_verteces(:,2);
    % to order we can start from the furthest point frmo the center of
    % mass, (as first) and the last point is the furthest point from
    % there
    furthest=find(vertex_matrix(sub2ind(size(labelled_reduced_rings),y,x))); % identify branchpoints that are nearby the vertex they would be good start points for the RDP algorithm
    pos=furthest;
    if isempty(furthest) % if there aren't branches, the next best option is to find the furthest vertex
        cm=mean(candidate_verteces); % center of mass
        [~,furthest]=max(sum((candidate_verteces-cm).^2,2)); % the furthest vertex
        branch_ids=1;
    end
    if size(unique(furthest),1)==1
        candidate_verteces=circshift(candidate_verteces,-furthest(1)+1);  % reordering, making the first element of furthest first
        candidate_verteces=[candidate_verteces;candidate_verteces(1,:)]; % closing the polygon.
        verteces=reducepoly2(candidate_verteces,12); % 9 is an empirical identify threshold, interestingly would correspond to 1 nm distance in real units
    else % make a subdivision based on where the branches are... for each branch take the single point that is closest to the branch itself
        [~,pos_max]=max((x-x(furthest)').^2+(y-y(furthest)').^2); % find the furthest points to those points
        to_remove=false(size(candidate_verteces,1),1);
        furthest=sort(unique([furthest;pos_max'])); % the verteces to start from
        iterations=length(furthest); % the number of segments to subdivide our perimeter
        verteces=false(size(candidate_verteces,1),1); % prepare a flag to avoid changing size
        for tt=1:iterations
            if tt<iterations
                vertex_of_interest=candidate_verteces(furthest(tt):furthest(tt+1),:);
                vertex_of_interest(to_remove(furthest(tt):furthest(tt+1)),:)=[];
            else
                vertex_of_interest=candidate_verteces(([furthest(tt):size(candidate_verteces,1),1:furthest(1)]),:); % circularity
                vertex_of_interest(to_remove([furthest(tt):end,1:furthest(1)]),:)=[];
            end
            vertex_of_interest=reducepoly2(vertex_of_interest,9);
            verteces=verteces+sum((candidate_verteces(:,1)-vertex_of_interest(:,1)').^2+(candidate_verteces(:,2)-vertex_of_interest(:,2)').^2==0,2);
        end
        candidate_verteces=candidate_verteces(verteces>0,:);
        to_remove=false(size(candidate_verteces,1),1);
%         to_remove(furthest)=true;
%         to_remove(furthest(pos))=false;
        furthest=find((vertex_matrix(sub2ind(size(labelled_reduced_rings),candidate_verteces(:,1),candidate_verteces(:,2))))); % the verteces to start from
        iterations=length(furthest); % the number of segments to subdivide our perimeter
        verteces=false(size(candidate_verteces,1),1); % prepare a flag to avoid changing size
        for tt=1:iterations
            if tt<iterations
                vertex_of_interest=candidate_verteces(furthest(tt):furthest(tt+1),:);
                vertex_of_interest(to_remove(furthest(tt):furthest(tt+1)),:)=[];
            else
                vertex_of_interest=candidate_verteces(([furthest(tt):size(candidate_verteces,1),1:furthest(1)]),:); % circularity
                vertex_of_interest(to_remove([furthest(tt):end,1:furthest(1)]),:)=[];
            end
            vertex_of_interest=reducepoly2(vertex_of_interest,12);
            verteces=verteces+sum((candidate_verteces(:,1)-vertex_of_interest(:,1)').^2+(candidate_verteces(:,2)-vertex_of_interest(:,2)').^2==0,2);
        end
        verteces=candidate_verteces(verteces>0,:);
    end
    vertex_matrix(sub2ind(size(vertex_matrix),verteces(:,1),verteces(:,2)))=pp;
    if isempty(verteces)==0
        index_new=sum((verteces(:,1)-branch_pos_Y').^2+(verteces(:,2)-branch_pos_X').^2>0,2)==length(branch_pos_Y);
        index_new=logical(index_new.*(sum(tril(verteces(:,1)==verteces(:,1)').*tril(verteces(:,2)==verteces(:,2)'),2)==1));
        branch_pos_Y=([branch_pos_Y;verteces(index_new,1)]);
        branch_pos_X=([branch_pos_X;verteces(index_new,2)]);
    end
end
end


function [vertex_matrix]=apply_RDP_to_segments(vertex_matrix,labelled_reduced_rings,vector_ids,diameter)
prop=regionprops(labelled_reduced_rings,'PixelIdxList');
sizes=size(labelled_reduced_rings);
for pp=vector_ids'
%    pp % for debugging purposes
    [y,x]=ind2sub(sizes,[prop(pp).PixelIdxList]);
    [branch_pos_Y,branch_pos_X]=find(bwmorph(labelled_reduced_rings(min(y):max(y),min(x):max(x))==pp,'endpoints'));
    % find endpoints
    if isempty(branch_pos_Y)
        labelled_reduced_rings_2=zeros(sizes);
        labelled_reduced_rings_2(prop(pp).PixelIdxList)=1;
        [vertex_matrix]=apply_RDP(vertex_matrix,labelled_reduced_rings_2,1,diameter);
    else
    distances=bwdistgeodesic(labelled_reduced_rings(min(y):max(y),min(x):max(x))==pp,branch_pos_X(1),branch_pos_Y(1));
    distances=distances(sub2ind(size(distances),y-min(y)+1,x-min(x)+1));
    index_reduced=sortrows([(1:length(x))',distances(:)],2); % sort based on distance
    x=x(index_reduced(1:end,1)); % apply sort
    y=y(index_reduced(1:end,1)); % apply sort
    candidate_verteces=[y,x]; % initial candidates
    y=candidate_verteces(:,1);
    x=candidate_verteces(:,2);
    % to order we can start from the furthest point frmo the center of
    % mass, (as first) and the last point is the furthest point from
    % there
    furthest=find(vertex_matrix(sub2ind(size(labelled_reduced_rings),y,x))); % identify branchpoints that are nearby the vertex they would be good start points for the RDP algorithm
    pos=furthest;
    if isempty(furthest) % if there aren't branches, the next best option is to find the furthest vertex
        cm=mean(candidate_verteces); % center of mass
        [~,furthest]=max(sum((candidate_verteces-cm).^2,2)); % the furthest vertex
        branch_ids=1;
    end
    if size(unique(furthest),1)==1
        candidate_verteces=circshift(candidate_verteces,-furthest(1));  % reordering, making the first element of furthest first
        candidate_verteces=[candidate_verteces;candidate_verteces(1,:)]; % closing the polygon.
        verteces=reducepoly2(candidate_verteces,12); % 9 is an empirical identify threshold, interestingly would correspond to 1 nm distance in real units
    else % make a subdivision based on where the branches are... for each branch take the single point that is closest to the branch itself
        [~,pos_max]=max((x-x(furthest)').^2+(y-y(furthest)').^2); % find the furthest points to those points
        to_remove=false(size(candidate_verteces,1),1);
        furthest=sort(unique([furthest;pos_max'])); % the verteces to start from
        iterations=length(furthest); % the number of segments to subdivide our perimeter
        verteces=false(size(candidate_verteces,1),1); % prepare a flag to avoid changing size
        for tt=1:iterations
            if tt<iterations
                vertex_of_interest=candidate_verteces(furthest(tt):furthest(tt+1),:);
                vertex_of_interest(to_remove(furthest(tt):furthest(tt+1)),:)=[];
            else
                vertex_of_interest=candidate_verteces(([furthest(tt):size(candidate_verteces,1),1:furthest(1)]),:); % circularity
                vertex_of_interest(to_remove([furthest(tt):end,1:furthest(1)]),:)=[];
            end
            vertex_of_interest=reducepoly2(vertex_of_interest,9);
            verteces=verteces+sum((candidate_verteces(:,1)-vertex_of_interest(:,1)').^2+(candidate_verteces(:,2)-vertex_of_interest(:,2)').^2==0,2);
        end
        candidate_verteces=candidate_verteces(verteces>0,:);
        to_remove=false(size(candidate_verteces,1),1);
%         to_remove(furthest)=true;
%         to_remove(furthest(pos))=false;
        furthest=find((vertex_matrix(sub2ind(size(labelled_reduced_rings),candidate_verteces(:,1),candidate_verteces(:,2))))); % the verteces to start from
        iterations=length(furthest); % the number of segments to subdivide our perimeter
        verteces=false(size(candidate_verteces,1),1); % prepare a flag to avoid changing size
        for tt=1:iterations
            if tt<iterations
                vertex_of_interest=candidate_verteces(furthest(tt):furthest(tt+1),:);
                vertex_of_interest(to_remove(furthest(tt):furthest(tt+1)),:)=[];
            else
                vertex_of_interest=candidate_verteces(([furthest(tt):size(candidate_verteces,1),1:furthest(1)]),:); % circularity
                vertex_of_interest(to_remove([furthest(tt):end,1:furthest(1)]),:)=[];
            end
            vertex_of_interest=reducepoly2(vertex_of_interest,12);
            verteces=verteces+sum((candidate_verteces(:,1)-vertex_of_interest(:,1)').^2+(candidate_verteces(:,2)-vertex_of_interest(:,2)').^2==0,2);
        end
        verteces=candidate_verteces(verteces>0,:);
    end
    vertex_matrix(sub2ind(size(vertex_matrix),verteces(:,1),verteces(:,2)))=pp;
    end
end
end
