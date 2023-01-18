%% parts I am scared to remove
%%% Opening the files
% if isfile([name_folder_r,'data.mat'])
%     load([name_folder_r,'data.mat']);
% end
 reference=creating_filter(diameter); % create a filter based on the expected size of the particle... we do not need it anymore
 to_build_model=3*diameter+(1-mod(3*diameter,2));
%% variable intialization specific to the folders of analysis

delete([name_folder_r,'._*']); % necessary if you opened the images in Mac... delete data you cannot access anyway

list_of_pictures=dir(fullfile(name_folder_r,data_format)); % find the full list of files
start_of_pictures=1; % Discard pictures from 1 to start_of_pictures
end_of_pictures=length(list_of_pictures); % nth final picture to consider
list_of_pictures=list_of_pictures(start_of_pictures:end_of_pictures); % list of files to open
number_of_pictures=length(list_of_pictures); % total number of pictures to analyse

%% where to collect data
if end_of_pictures>0
    total_area=[total_area,zeros(1,number_of_pictures)]; % allocate the necessary space for storing the infos of the subfolder
    regions_number=[regions_number;zeros(end_of_pictures-start_of_pictures+1,1)]; % allocate the necessary space for storing the infos of the subfolder
    polygons=[polygons;sparse(end_of_pictures-start_of_pictures+1,maximum_polygon)]; % allocate the necessary space for storing the infos of the subfolder
    polygons_on_the_boarders=[polygons_on_the_boarders;sparse(end_of_pictures-start_of_pictures+1,maximum_polygon)]; % allocate the necessary space for storing the infos of the subfolder
    total_number_of_particle=[total_number_of_particle,zeros(1,number_of_pictures)]; % allocate the necessary space for storing the infos of the subfolder
    ideal_number_of_rings=[ideal_number_of_rings,zeros(1,number_of_pictures)]; % allocate the necessary space for storing the infos of the 
    % reference_picture=number_of_pictures;
    ring=[]; % additional variables
    connectivity=[]; % additional variables, is it used?
    time=[time,int32(0:number_of_pictures-1)]; % intialize time axes
end
retrace=[]; % where the retrace is stored
trace=[]; % where the trace is stored
frame_axis=1:sense:number_of_pictures; % indeces of the frames to analyse
if sense<0
    frame_axis=fliplr(frame_axis); % go from higher to lower index
end
for ii=frame_axis
    particle_in_frame=[]; % the particles in a single frame
    open_afm_image % open image and do flattening
    type_of_rings=0;
    %% recover information
    % time
    C=cell2mat(textscan(list_of_pictures(ii).date,'%*d-%*c%*c%*c-%*d %d:%d:%d')); % specific of the format used in the HS-AFM
    if time_series==true
     time(ii+reference_ii)=C(1).*3600+C(2).*60+C(3); % retrive absolute time in seconds (within the day)
     if ii+reference_ii==1
         time_starting=time(ii+reference_ii); % when the subfolder started the measurments
         delay_x=0; % to use to correlate patches when the field of view changes, displacement in X
         delay_y=0; % to use to correlate patches when the field of view changes, displacement in Y
     end
    time(ii+reference_ii)=time(ii+reference_ii)-time_starting; % TIME in seconds from the beginning of the experiment
    end
    % useful info in SPACE dimensions
    particle_size=ceil(diameter/sqrt(W*L));
    total_area(ii+reference_ii)=non_segmented_image.xreal.*non_segmented_image.yreal/to_m^2; % area in SPACE UNITS

    %% segmentation and skeletonization

    [skel,bw]=segmentation_and_skeletonization(image_to_process,reference,dim_y_o,dim_x_o,dim_y,dim_x,diameter,thr_accumulator_reference);
    %% Identifying different islands and analysing them for the first time
    first_island_analysis
    %% Time recording and checking for jumps (a jump will have an higher dt between the images)
%     if ii>1
%        if time(ii)-time(ii-1)>(dt*1.5)
%        crr=xcorr2(double(target),double(previous>0));
%        [~,pos_lag]=max(crr(:));
%        [indy,indx]=ind2sub(size(crr),pos_lag);
%        [dim_yb,dim_xb]=size(target);
%        Ydel=-dim_yb:dim_yb;
%        Xdel=-dim_xb:dim_xb;
%        delay_x=Xdel(indx);
%        delay_y=Ydel(indy);
%        else
%         delay_x=0;
%         delay_y=0;
%        end
%     end

    %% Plotting original image
    support_matrix=imresize((non_segmented_image.data-min(non_segmented_image.data(:)))./max((non_segmented_image.data(:)-min(non_segmented_image.data(:))))*255,[dim_y_o,dim_x_o]);
    path_to_image=[name_folder_r,list_of_pictures(ii).name(1:end-4)];
    tt=printing_image_simple(support_matrix,path_to_image,'_original',mappa);
    close(tt)
    clear tt
    %% Plotting segmented image
    support_matrix=imresize((non_segmented_image.data-min(non_segmented_image.data(:)))./max((non_segmented_image.data(:)-min(non_segmented_image.data(:))))*255.*bw,[dim_y_o,dim_x_o])>0.2;
    tt=printing_image_simple(support_matrix,path_to_image,'_segmented',[0 0 0; 1 1 1]);
    close(tt)
    clear tt
    %% Finding voids and first Polygon classification (concave/concentric). it is going to be useful for misscaractherized particles
    rings_finding_and_preliminary_analysis
%     finding_the_voids_hugh
    %% Looking for particles and finalising analysis of rings and patches

    final_analysis
    %% Plotting image where the data will be represented
    total_number_of_particle(ii+reference_ii)=size(particle_in_frame,1);
    density=total_number_of_particle(1:ii+reference_ii)./total_area(1:ii+reference_ii);

    tt=preparing_image(imresize((non_segmented_image.data-min(non_segmented_image.data(:)))./max((non_segmented_image.data(:)-min(non_segmented_image.data(:))))*255,[dim_y_o,dim_x_o]),mappa);
    % copy the figure in a new one as subplot taken from https://ch.mathworks.com/matlabcentral/answers/92538-how-can-i-copy-an-existing-figure-onto-another-figure-as-a-subplot-using-matlab-7-10-r2010a
    vv=figure('visible','off');
    ax1=subplot(1,2,1);
    ax2=subplot(1,2,1);
    axes_to_be_copied = findobj(tt,'type','axes');
    % Identify the children of this axes 
    chilred_to_be_copied = get(axes_to_be_copied,'children'); 
    % Identify orientation of the axes 
    [az,el] = view; 
    % Copy the children of the axes 
    copyobj(chilred_to_be_copied,ax1); 
    % Set the limits and orientation of the subplot as the original figure 
    set(ax1,'Xlim',get(axes_to_be_copied,'XLim')) 
    set(ax1,'Ylim',get(axes_to_be_copied,'YLim')) 
    set(ax1,'Zlim',get(axes_to_be_copied,'ZLim')) 
    view(ax1,[az,el]) 
    % One may set other properties such as labels, ticks etc. using the same 
    % idea 
    % Close the figure 
    close(tt);     
    figure(vv)
    set(vv,'visible','off')
    set(vv,'Position',[10,10,1024,512])
    subplot(1,2,1)
    axis(ax1,'square')
    ax1.XTick=[];
    ax1.YTick=[];
    colormap(mappa)
    hold on
    %% Plotting the centers
    if isempty(particle_in_frame)==0
        scatter(round(particle_in_frame(:,1)/dim_x*dim_x_o),round(particle_in_frame(:,2)/dim_y*dim_y_o),10,[84,39,136]/255,'filled','p');
    end
    hold on
    subplot(1,2,2)
    plot(time(1:ii+reference_ii),density(1:ii+reference_ii))
    xlabel('Time [s]')
    ylabel('Particle density [1/nm^2]')
    hold on
    plot(time(ii+reference_ii),density(ii+reference_ii),'o','LineWidth',2)
    axis([0 double(max(time(ii+reference_ii),dt)) 0 max([1.8e-3*total_area(reference_ii+ii)/64000;max(density(1:(ii+reference_ii)))])])
    curtick = get(gca, 'XTick');
    set(gca, 'XTick', unique(round(curtick)));
    
    print(vv,[name_folder_r,list_of_pictures(ii).name(1:end-4),'_particle'],'-dpng','-r0');

    %% Plotting the concave region 
    % could have pulled in the previous loop, but for code transferability
    % I keep outside
    subplot(1,2,1)
%     for ss=list_region(logical(conc))'
%         [row,col]=find(lab==ss);
%         plot(col.*resized_f,row.*resized_f,'.','Color',colors_concave)
%     end

    
if isempty(particle_in_frame)==0    
    %% plotting the defects
    subplot(1,2,1)
    labelled_rings=imerode(labelled_rings.*bwmorph(imerode(labelled_rings>0,ones(3)),'thin',1/2*diameter),ones(3));
    if ~isempty(rings_index)
        for zz=rings_index(logical((type_of_rings>=minimum_polygon) .* (type_of_rings<=maximum_polygon) ))'
            [row,col]=find(labelled_rings==zz);
            plot(round(col/dim_x*dim_x_o),round(row/dim_y*dim_y_o),'.','Color',colors_map(type_of_rings(zz),:)) % Color classification
        end
    end  
    polygons(ii+reference_ii,minimum_polygon:maximum_polygon)=histcounts(type_of_rings,[(minimum_polygon-0.5):(maximum_polygon+0.5)]);
    %% adding legends and bar graphs
    ax2=subplot(1,2,2);
    delete(ax2)
    ax2=subplot(1,2,2);
    ff=bar(poly_x_axes,full(polygons(ii+reference_ii,4:10))/total_area(ii+reference_ii));
    ff.FaceColor='flat';
    for kk=1:7
        ff.CData(kk,:)=colors_map(kk+3,:);
    end
    title('Density of polygons')
   
    ylabel('Density [nm⁻²]')
    hsp2 = get(gca, 'Position');                   % Get 'Position' for (2,1,2)
%     set(gca, 'Position', [hsp2(1:3)  hsp1(4)]) 
    curtick = get(gca, 'XTick');
       %Create list of positions for the centres
    count_particles_found=length(particle_in_frame(:,1));
end
%     vv=subplot(1,2,2);
%     hold off
%     density=full(polygons(1:ii+reference_ii,4:10))./total_area(1:ii+reference_ii)';
%  
% 
%     for kk=1:7
%         h(kk)=plot(time(1:ii+reference_ii),density(1:ii+reference_ii,kk),'Color',colors_map(kk+3,:));
%         
%         hold on
%         xlim(vv,[0 max(time(ii+reference_ii),dt)])
%         ylim(vv,[0 1.1*max([max(density),1./(total_area(1:ii+reference_ii))])])
%     end
%     h_legend=legend({'Quadrilateral','Pentagon','Hexagon','Heptagon','Octagon','Ennagon','Decagon'},'AutoUpdate','off');
%     asa=plot(time(ii+reference_ii),density(ii+reference_ii,kk),'.','Color',colors_map(kk+3,:),'LineWidth',3);
%     title('Polygons counted in the frame')
%     ylabel('Density [nm⁻²]')
%     xlabel('Time [s]')
% %     hsp2 = get(gca, 'Position');                   % Get 'Position' for (2,1,2)
% %     set(gca, 'Position', [hsp2(1:3)  hsp1(4)]) 
%     curtick = get(gca, 'XTick');

%     set(gca, 'XTick', unique(round(curtick)));
%     set(h_legend, 'location', 'northeastoutside')

    print(vv,[name_folder_r,list_of_pictures(ii).name(1:end-4),'_defects'],'-dpng','-r0');
    
    % Save image of patch color coded
%    prepare the color coded image
    
%     vv=figure('visible','off');
%     photo_1=imshow(imresize(color_coded/(color_value+1),[dim_x_o,dim_y_o],'nearest'));
%     set(vv,'Position',[10,10,512,512])
%     axis off
%     colormap(vv,map_color_patches)    
%     print(vv,[name_folder,list_of_pictures(ii).name(1:end-4),'_color_coded'],'-dpng','-r0'); 

%% Closing all
    map_islands=parula(max(image_by_grains(:)));
    map_islands(2:end,:)=map_islands(1+randperm(size(map_islands,1)-1),:);
    tt=printing_image_simple(imresize(image_by_grains,[dim_y_o,dim_x_o],'nearest'),path_to_image,'_grains',map_islands);
    clear tt

    close all;
    fclose('all');
    clear vv;
    if mod(ii,50)==0
       save('final.mat') 
    elseif ii==number_of_pictures
       save('final.mat') 
    end
end
 save([name_folder_r,'final','.mat']);

%% Make a video 
reference_ii=reference_ii+number_of_pictures;
%  videosm
 


function tt=printing_image_simple(matrix,path_to_image,suffix,mappa)
    tt=preparing_image(matrix,mappa);
    print(tt,[path_to_image,suffix],'-dpng','-r600');
end

function tt=preparing_image(matrix,mappa)
tt=figure('visible','off');
    set(tt,'Position',[10,10,512,512])
     photo_1=image(matrix);
     photo_1.CDataMapping='direct';
     axis off
     colormap(photo_1.Parent,mappa)
end
