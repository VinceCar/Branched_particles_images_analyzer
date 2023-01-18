close all
clear all
%% This is a script to plot the final steady states. It has to be used interactivtely
path_main_folder='/Volumes/One Touch/AFM/last_analysis/HS_AFM_definitive/grains_analysis/';
%monomers=[{'sGCGC'};{'sGCTA'};{'sTATA'};{'sTATA'};{'lGCGC'};{'lGCTA'};{'lTATA'};];
%subfolders=[{'/D'};{'/C'};{'_1/A'};{'_3/2021-07-15_exp37ct105a'};{'/2021-04-27_exp08112C'};{'/2021-04-27_exp04118'};{'/2021-04-27_exp0581'};];
monomers=[{'sGCGC'};{'sGCTA'};{'sTATA'};{'lGCGC'};{'lGCTA'};{'lTATA'};];
subfolders=[{'/D'};{'/C'};{'_1/A'};{'/2021-04-27_exp08112C'};{'/2021-04-27_exp11118'};{'/2021-04-27_exp0581'};];
color_monomers=[0.752,0.243,0.137;1,0.406,0;0.999,0.576,0;0.787,0.344,0.256;1,0.478,0.153;1,0.678,0.249];

rows_subplot=2;
cols_subplot=3;
list_of_subplot_where_to_put_x_label=[2;5];
list_of_subplot_where_to_put_y_label=[1;4];

monomers=[{'sGCGC'};{'lGCGC'}];
subfolders=[{'/D'};{'/2021-04-27_exp08112C'}];
color_monomers=[0.752,0.243,0.137;0.787,0.344,0.256];

rows_subplot=2;
cols_subplot=3;
list_of_subplot_where_to_put_x_label=[1;2];
list_of_subplot_where_to_put_y_label=[1;2];


% monomers=[{'sGCTA'};{'lGCTA'};{'lcGCTA'}];
% subfolders=[{'/C'};{'/2021-04-27_exp04118'};{'/long core3_20hz'}];
frames_to_discard=[{[]};{[]};{[]};{[]};{[]};{[]};{[]}];
file_name='final.mat';
statistics_rings=zeros(length(monomers),3);
statistics_hexagon_percentage=zeros(length(monomers),2);

number_of_pics_to_use=15;
collecting_MNDs=zeros(size(monomers,1),number_of_pics_to_use);
collecting_hexagons=collecting_MNDs;
xf=0;   % Screen position of the figure
yf=0;   % Screen position of the figure
width=1024; % Width in pixels
height=512; % height in pixels
Fontsize_xaxis=12;
Fontsize_yaxis=12;
Fontsize_label_x=12;
Fontsize_label_y=12;
for monomer_number=1:length(monomers)
    counter_plotted_images=0;
    Smoothin_param=3.99956076419875e-07;
    load([path_main_folder,monomers{monomer_number},subfolders{monomer_number},'/',file_name],'total_area','ring','particle_infos_list','frame_axis','ii','minimum_polygon','maximum_polygon','islands_infos_list','diameter','list_of_pictures');
    number_of_pictures_analysed=find(frame_axis==ii);
    %% Recompute the time cause the way we detected it before was based on the creation of the file (but we copied it ... so it was messed up)
    time=zeros(1,number_of_pictures_analysed);
    for tt=frame_axis(1:number_of_pictures_analysed)
        C=cell2mat(textscan(list_of_pictures(tt).name,'%*d-%*d-%*d_%d-%d-%d__%*d.gwy')); % specific of the format used in the HS-AFM
        time(tt)=double(C(1).*3600+C(2).*60+C(3));
    end
    time=time-min(time);

    %% Remove data relative to a previous measurment that has not been separated in another folder
    if strcmp(monomers{monomer_number},'sGCTA') && strcmp(subfolders{monomer_number},'/C')
        frames_to_remove=find(time<306);
        total_area(frames_to_remove)=[];
        ring(ismember(ring(:,end),frames_to_remove),:)=[];
        particle_infos_list(ismember(particle_infos_list(:,end),frames_to_remove),:)=[];
        frame_axis(ismember(frame_axis,frames_to_remove))=[];
        islands_infos_list(ismember(islands_infos_list(:,end),frames_to_remove),:)=[];
        time(frames_to_remove)=[];
        time=time-time(1);
        min_frame_ax=min(frame_axis);
        islands_infos_list(:,end)=islands_infos_list(:,end)-min_frame_ax+1;
        particle_infos_list(:,end)=particle_infos_list(:,end)-min_frame_ax+1;       
        ring(:,end)=ring(:,end)-min_frame_ax+1;       
        frame_axis=frame_axis-min_frame_ax+1;        
    end
    if strcmp(monomers{monomer_number},'sTATA') && strcmp(subfolders{monomer_number},'_1/A')
        frames_to_remove=find(time<330);
        total_area(frames_to_remove)=[];
        ring(ismember(ring(:,end),frames_to_remove),:)=[];
        particle_infos_list(ismember(particle_infos_list(:,end),frames_to_remove),:)=[];
        frame_axis(ismember(frame_axis,frames_to_remove))=[];
        islands_infos_list(ismember(islands_infos_list(:,end),frames_to_remove),:)=[];
        time(frames_to_remove)=[];
        time=time-time(1);
        min_frame_ax=min(frame_axis);
        islands_infos_list(:,end)=islands_infos_list(:,end)-min_frame_ax+1;
        particle_infos_list(:,end)=particle_infos_list(:,end)-min_frame_ax+1;       
        ring(:,end)=ring(:,end)-min_frame_ax+1;       
        frame_axis=frame_axis-min_frame_ax+1;        
    end
    number_of_pictures_analysed=length(time);
    %% Islands to consider based on minimu size
    island_to_consider=islands_infos_list(:,4)>3;
    %% Basi quantites such as number of rings, particles and hexagons per frame
    axis_of_interest_for_histogram=[frame_axis-0.5,frame_axis(end)+0.5];
    total_number_of_rings_per_frame=histcounts(ring(:,end),axis_of_interest_for_histogram);
%    total_number_of_particle_per_frame=histcounts(particle_infos_list(:,end),axis_of_interest_for_histogram);
    total_number_of_particle_per_frame=histcounts(repelem(islands_infos_list(island_to_consider,end),islands_infos_list(island_to_consider,4)),axis_of_interest_for_histogram);
    total_number_of_hexagons=histcounts(ring(ring(:,3)==6,end),axis_of_interest_for_histogram);
    dens_particle=total_number_of_particle_per_frame./total_area*diameter^2;
    dens_rings=total_number_of_rings_per_frame./total_area*diameter^2;
    gradient_data=0;
    dens_particle_monotone=makedatamonotone(time,dens_particle);
    while sum(gradient_data<=0)>0
        fitted_dens_particle=(feval(smoothing_spline(time',dens_particle_monotone,Smoothin_param),time))';
        Smoothin_param=Smoothin_param/8;
        gradient_data=gradient(fitted_dens_particle);
    end
   dens_hexagons=total_number_of_hexagons./total_number_of_rings_per_frame;

    %% Plot the density of particles in function of time
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(time,dens_particle,[],color_monomers(monomer_number,:),'filled')
    hold on
    plot(time,fitted_dens_particle,'LineWidth',2,'Color','k')
    % title('Number of particles in function of particle density')
    title(monomers{monomer_number},'Fontsize',14)
    legend('Measured data','Monotone fitted spline')
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Time [s]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    ylim([0 0.25]);
    xlim([0 1500])

    %% Compute and plot the density of islands in function of time
    islands_in_function_of_time=histcounts(islands_infos_list(island_to_consider,end),axis_of_interest_for_histogram);
    dens_islands=islands_in_function_of_time./total_area*diameter^2;
    smoothed_dens_islands=(feval(smoothing_spline(time,dens_islands,Smoothin_param),time))';
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(time,dens_islands,[],color_monomers(monomer_number,:),'filled')
    hold on
    plot(time,smoothed_dens_islands,'LineWidth',2,'Color','k')
    % title('Number of particles in function of particle density')
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Time [s]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Islands density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    ylim([0 0.02]);
    xlim([0 1500])

    %% Identify outlayers forcing continuity on the density of particle and islands
    identify_outlayers=logical((abs(fitted_dens_particle-dens_particle)>0.1)+(abs(smoothed_dens_islands-dens_islands)>0.003));
    frames_to_remove=find(identify_outlayers);
    frames_to_keep=~identify_outlayers;

    %% Plot the density of rings in function of particle density
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(fitted_dens_particle(frames_to_keep),dens_rings(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    [P,S]=polyfit(fitted_dens_particle+rand(1,number_of_pictures_analysed)*min(fitted_dens_particle)/1e2,dens_rings,1);
    % title(['Rings density in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Rings density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1])
    ylim([0 2.5e-1])
    yfit = polyval(P,dens_particle);          % Estimated  Regression Line
    hold on
    plot([0 sort(dens_particle) 2.5e-1],polyval(P,[0 sort(dens_particle) 2.5e-1]),'LineWidth',2,'Color','k')
    SStot = sum((dens_rings-mean(dens_rings)).^2);                    % Total Sum-Of-Squares
    SSres = sum((dens_rings-yfit).^2);                       % Residual Sum-Of-Squares
    Rsq = 1-SSres/SStot;                            % R^2
    xlim([0 2.5e-1]);
    ylim([0 2e-1])
    xl = xlim;
    yl = ylim;
    xt = 0.05 * (xl(2)-xl(1)) + xl(1);
    yt = 0.90 * (yl(2)-yl(1)) + yl(1);
    caption = sprintf('y = %f * x + %f\n R^2 = %f', P(1), P(2), Rsq);
    text(xt, yt, caption, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
    statstics_rings(ii,:)=[P(2), P(1),Rsq];
    
    %% Plot the density of hexagon in function of particle density
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(dens_particle(frames_to_keep),dens_hexagons(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    [P,S]=polyfit(dens_particle+rand(1,number_of_pictures_analysed)*min(dens_particle)/1e2,dens_hexagons,0);
    % title(['Fraction of hexagons in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Fraction','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 1])

    statistics_hexagon_percentage(ii,:)=[P(1),Rsq];   
    

    %% Ideally you will want to do a vector that records the position of the maximum island in function of time
    % to do that you can use the time labelling and using particle number
    % as intensity masks in regionprops 
    frames_with_islands=ismember(frame_axis,unique(islands_infos_list((islands_infos_list(:,5)>0),end)));
    maximum_size_in_function_of_time=zeros(number_of_pictures_analysed,1);
    support_matrix=regionprops(islands_infos_list(islands_infos_list(:,4)>0,end),islands_infos_list(islands_infos_list(:,4)>0,4),'MaxIntensity','Area');
    maximum_size_in_function_of_time(frames_with_islands)=[support_matrix(:).MaxIntensity];
    density_islands=zeros(number_of_pictures_analysed,1);
    density_islands(frames_with_islands)=[support_matrix(frames_with_islands).Area]'./([total_area(frames_with_islands)]')*diameter^2;
    % one option that is memory intensive is to compare each element in
    % time and particles
    index_major_island=find(sum((islands_infos_list(:,end)==(frames_with_islands)).*(islands_infos_list(:,4)==(maximum_size_in_function_of_time(frames_with_islands)')),2));

    %% plot the maximum size in function of time
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(dens_particle(frames_to_keep),maximum_size_in_function_of_time(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    % title(['Size of maximum cluster in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
       ylabel('Number of monomers','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 400]);
    %% Compute the solidity of the biggest patch
    solidity=islands_infos_list(:,5)./(1+islands_infos_list(:,4)/2-sqrt(3*islands_infos_list(:,4)/2));
    if length(unique(islands_infos_list(index_major_island,end)))~=length(index_major_island)
        hisorep=histcounts(islands_infos_list(index_major_island,end),axis_of_interest_for_histogram);
        repeated_frames=find(hisorep>1);
        solidity_big_patches=solidity(index_major_island);
        matrix_of_reference=sortrows([index_major_island,islands_infos_list(index_major_island,end),solidity_big_patches],[2,-3]);
        check_for_ascending_order=[true(1);diff(matrix_of_reference(:,2))>0]; % take only the first one per each patch
        index_major_island=[matrix_of_reference(check_for_ascending_order)];
    end

    solidity_major_island=zeros(number_of_pictures_analysed,1);
    solidity_major_island(frames_with_islands)=solidity(index_major_island);

    %% Plot the solidity in function of particle density
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    %title(['MND of the maximum cluster in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    scatter(dens_particle(frames_to_keep),solidity_major_island(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('MND','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 1.2]);

    %% Computing the rings in maximum patch
    rings_in_maximum_patch_in_function_of_time=zeros(number_of_pictures_analysed,1);
    rings_in_maximum_patch_in_function_of_time(frames_with_islands)=islands_infos_list(index_major_island,5);
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(dens_particle(frames_to_keep),rings_in_maximum_patch_in_function_of_time(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    %title(['Rings in maximum patch in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Count','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 100]);
    %% Computing the fraction of rings that are on the boarder
    % this rely on the same approach as before, we need to find which rings
    % are part of the biggest island at the given time. we have the id of
    % the island for each time frame there are islands. Remember that cannot be a ring that is not an island...
    id_rings_in_maximum_patch=find(sum((ring(:,5)==(islands_infos_list(index_major_island,1)')).*(ring(:,end)==(islands_infos_list(index_major_island,end)')),2));
    % now for each time frame we have the rings in the patch... if we
    % multiply the flag of the boarder contained in column 4 with the time
    % frame and we divide by the histcounts of how many rings we have per
    % time frame in the maximum patch, we can get the fraction...
    ring(ring(:,4)==0,4)=1; %% correct for bug in program
    rings_on_boarder_maximum_patch=histcounts((ring(id_rings_in_maximum_patch,4)==1).*ring(id_rings_in_maximum_patch,end),axis_of_interest_for_histogram)';
    fraction_of_rings_on_boarder_maximum_patch=rings_on_boarder_maximum_patch./rings_in_maximum_patch_in_function_of_time;
    fraction_of_rings_on_boarder_maximum_patch(isnan(fraction_of_rings_on_boarder_maximum_patch))=0;
    if sum(isinf(fraction_of_rings_on_boarder_maximum_patch))
        fprintf('Error!')
    end
    %% Plot the fraction of rings on boarder for the maximum patch
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(dens_particle(frames_to_keep),fraction_of_rings_on_boarder_maximum_patch(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    %title(['Fraction of rings in maximum patch that are on boarder in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Fraction','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 1]);

    %% Same in function of time
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(time(frames_to_keep),fraction_of_rings_on_boarder_maximum_patch(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    %title(['Fraction of rings in maximum patch that are on boarder in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    title(monomers{monomer_number},'Fontsize',14)
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Time [s]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Fraction','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 1500]);
    ylim([0 1]);
    %% Computing the fraction of particles that are on the boarder for the maximum patch
    % here we need a bit of rearrangment of particle_infos_list cause I
    % saved the label of rings in a frame and not of the total counted, we
    % need therefore to add the cumsum of the total rings counted in the
    % previous instant, if the particle is not on the boarder
    % this one can be particularly memory intensive if we are not
    % careful... an easy option is first to devide in particles that are
    % not on the boarders. they must be part of rings, so you can then
    % check the particles that are part of specific rings at given times.
    % in case the memory runs out, run a for cycle
    internal_particles_index=find(prod(particle_infos_list(:,3:5)>0,2));
    internal_particles_subset=particle_infos_list(internal_particles_index,:);
    quantity_to_add=[0;cumsum(total_number_of_rings_per_frame(1:end-1))'];
    % now identify which patch are they part of by comparing times and
    % rings at given time
%     index_internal_particles_maximum_patch=find(sum(((internal_particles_subset(:,3)+quantity_to_add(internal_particles_subset(:,end)))==((id_rings_in_maximum_patch'))).*(internal_particles_subset(:,end)==((ring(id_rings_in_maximum_patch,end))')),2));
%     internal_particles_per_frame_maximum_patch=histcounts(internal_particles_subset(index_internal_particles_maximum_patch,end),axis_of_interest_for_histogram);
%     fraction_of_particles_on_boarder_maximum_patch=internal_particles_per_frame_maximum_patch'./maximum_size_in_function_of_time;
%     fraction_of_particles_on_boarder_maximum_patch(isnan(fraction_of_particles_on_boarder_maximum_patch))=0;
%     % fraction of particles that are on a boarder in function of time
%     figure
%     scatter(dens_particle,fraction_of_rings_on_boarder_maximum_patch,'LineWidth',2)
%     title(['Fraction of particles in maximum patch that are on boarder in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
%     xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
%     ylabel('Fraction','Fontsize',12,'FontWeight','bold')  
    % Plot the grains per islands
    %% average size of islands in function of time
    weight=2; % use 1 for aritmetic averages and 2 for weighted averages
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    weighted_sum_of_sizes=histcounts(repelem(islands_infos_list(island_to_consider,end),islands_infos_list(island_to_consider,4).^weight),axis_of_interest_for_histogram);
    if weight==2
        weighted_averages_size_islands=weighted_sum_of_sizes./total_number_of_particle_per_frame;
    else
        weighted_averages_size_islands=weighted_sum_of_sizes./islands_in_function_of_time;
    end
    subplot(rows_subplot,cols_subplot,monomer_number)
    plot(time(frames_to_keep),weighted_averages_size_islands(frames_to_keep),'LineWidth',2,'Color',color_monomers(monomer_number,:))
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Time [s]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Average size of islands [monomers]','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 1500])
    ylim([0 80])
    %% Do it but in function of particle density
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(dens_particle(frames_to_keep),weighted_averages_size_islands(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Average size of islands [monomers]','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 400])
    %% average size of grains in function of time?
    islands_infos_list(islands_infos_list(:,6)==0,6)=1; % the minimum level is 0, this happens for misscoding on the grain part. on the bright side all the others number should be correct
    weight=2; % use 1 for aritmetic averages and 2 for weighted averages

    grains_in_function_of_time=histcounts(repelem(islands_infos_list(island_to_consider,end),islands_infos_list(island_to_consider,6)),axis_of_interest_for_histogram);
    weighted_sum_of_grains=histcounts(repelem(islands_infos_list(island_to_consider,end),islands_infos_list(island_to_consider,6).^weight),axis_of_interest_for_histogram);
    if weight==2
        weighted_averages_grains_islands=weighted_sum_of_grains./grains_in_function_of_time;
    else
        weighted_averages_grains_islands=weighted_sum_of_grains./islands_in_function_of_time;
    end
    figure(counter_plotted_images+1)
    set(gcf,'Position',[xf,yf,width,height]);
    counter_plotted_images=counter_plotted_images+1;
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(dens_particle(frames_to_keep),weighted_averages_grains_islands(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Particle density [n.u.]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Average size of grains','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 2.5e-1]);
    ylim([0 10])
    %% average size of islands in terms of grains in function of time?
    figure(counter_plotted_images+1)
    counter_plotted_images=counter_plotted_images+1;
    set(gcf,'Position',[xf,yf,width,height]);
    subplot(rows_subplot,cols_subplot,monomer_number)
    scatter(weighted_averages_size_islands(frames_to_keep),weighted_averages_grains_islands(frames_to_keep),[],color_monomers(monomer_number,:),'filled')
    if sum(monomer_number==list_of_subplot_where_to_put_x_label)
        xlabel('Islands size [monomers]','Fontsize',12,'FontWeight','bold')
    end
    if sum(monomer_number==list_of_subplot_where_to_put_y_label)
        ylabel('Average size of islands [grains]','Fontsize',12,'FontWeight','bold')
    end
    xlim([0 200])
    ylim([0 15])
end
min_x_pos=Inf;
max_x_pos=-Inf;
min_y_pos=Inf;
max_y_pos=-Inf;
for ii=1:length(monomers)
    for jj=1:counter_plotted_images
        figure(jj)
        subplot(rows_subplot,cols_subplot,ii)
        ax=gca;
        ax.XAxis.FontSize = Fontsize_xaxis;
        ax.YAxis.FontSize = Fontsize_yaxis;
        ax.XLabel.FontSize = Fontsize_label_x;
        ax.YLabel.FontSize = Fontsize_label_y;
        ax.XAxisLocation='origin';
        ax.YAxisLocation='origin';

    end
end
