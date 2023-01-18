[y,x]=find(vertex_matrix); % coordinates of where the particles are
island_id_per_particle=label_islands((sub2ind([dim_y,dim_x],y,x))); % find out which island each particle belongs too.
label_fill_islands=imfill(label_islands,'holes');

if ~isempty(island_id_per_particle)

labelled_vertex_matrix=vertex_matrix+0; % convert to double
index_particle=(1:length(x))'; % a unique id [in frame] for the particles
labelled_vertex_matrix(sub2ind([dim_y,dim_x],y,x))=index_particle; % create a label matrix of each vertex
%labelled_vertex_matrix=imdilate(labelled_vertex_matrix,ones(3)); % make them one pixel bigger
% labelled_vertex_matrix=labelled_vertex_matrix.*(bwmorph((labelled_vertex_matrix>0)-(labelled_rings==0)>0,'shrink',Inf)); % we have one pixel representing each vertex per hole! redundant, but needed at this point
labelled_vertex_matrix=imdilate(labelled_vertex_matrix,ones(5)); % make them one pixel bigger
labelled_vertex_matrix=labelled_vertex_matrix.*bwmorph((labelled_vertex_matrix>0)-(imerode(labelled_rings,ones(3))==0)>0,'shrink',Inf);
index=find(labelled_vertex_matrix); % find the indeces (in matrix) of where the verteces are.
particle_id=labelled_vertex_matrix(index); % each vertex what particle represent. Nb: from one particle we may either have 1, 2 or 3 points
ring_particle=labelled_rings(index); % find out to what rings the particles belong to (either nothing, 1, 2 or 3)
if ~isempty(rings_index)
    % this is to find out to what islands the rings are part of, to do so we use regionprops
    labelled_centers=bwmorph(labelled_rings>0,'shrink',Inf).*double(labelled_rings);
    support_matrix=zeros(size(ring_in_frame,1),1);
    support_matrix(labelled_centers(labelled_centers>0))=double(label_fill_islands(labelled_centers>0));
	ring_in_frame(:,5)=support_matrix; % add the info on what patch contains the ring
end
% compleating the informations of patches with number of particles and
% number of rings (column 5)
    islands_in_frame(:,4)=histcounts(island_id_per_particle,[islands_index-0.5;max(islands_index)+0.5]); % count how many particles each island has

if ~isempty(rings_index)
     islands_in_frame(:,5)=histcounts(ring_in_frame(:,5),[islands_index-0.5;max(islands_index)+0.5]); % count how many rings each island has
     islands_in_frame(:,6)=histcounts(label_fill_islands(bwmorph(imfill(imerode(image_by_grains.*(matrix_level>0)>0,ones(3)),'holes').*(label_fill_islands>0),"shrink",Inf)>0),[islands_index-0.5;max(islands_index)+0.5]); % find out how many subparts are the patches made of 
     
end
% finding which islands the particles are connected to
missing_particle_amount=holes_per_particle-histcounts(particle_id,[(1:length(x))'-0.5;length(x)+0.5]); % find out if a particle is on a boarder! at maxium a particle can be in contact with "holes_per_particle", we then count how many holes it is in actually in contact with and then subtract
for kk=find(missing_particle_amount<0) % if you have more holes per particle print a warning and force compatibility (really weird need to look into)
    id_to_remove=find(particle_id==kk);
    particle_id(id_to_remove(4:end))=[]; % force to save a maximum of 3 holes, by keeping a maximum of 3 verteces per id
    ring_particle(id_to_remove(4:end))=[];
    fprintf(['Warning! In frame ',num2str(ii),' with particle ', num2str(kk),'. Too many rings in contact with a single particle!\n'])
end

ring_particle=[ring_particle; zeros(length(x)*3-length(ring_particle),1)]; % 0 pad how many rings per particle we have, so that we are able to have always three entries per particle
index_particle=index_particle(missing_particle_amount>0); % find the index of the particles that are on the boarders (at least one ring in contact has label 0)
missing_particle_amount=missing_particle_amount(missing_particle_amount>0); % report it only for those particles that actualy are not in contact with all the possible holes
support_matrix=repelem(index_particle,missing_particle_amount); % create replicas of the indeces so that everything has 3 replicas
support_matrix=sortrows([[particle_id;support_matrix(:)],ring_particle],1);
particle_in_frame=[x,y,reshape(support_matrix(:,2),[],holes_per_particle),ii*ones(length(x),1)];

% transferring infos to record in time
particle_infos_list=[particle_infos_list;particle_in_frame];
if ~isempty(rings_index)
	ring=[ring;ring_in_frame];
end
islands_infos_list=[islands_infos_list;islands_in_frame];
end