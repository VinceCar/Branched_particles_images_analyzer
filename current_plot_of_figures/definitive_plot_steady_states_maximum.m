%% This is a script to plot the final steady states. It has to be used interactivtely
path_main_folder='/Volumes/One Touch/AFM/last_analysis/';
monomers=[{'sGCGC'};{'sGCTA'};{'sTATA'};{'lGCGC'};{'lGCTA'};{'lTATA'};];
% monomers=[{'lGCTA'};{'lGCTA_1_2'};{'lGCTA_1_3'}]
%monomers=[{'sGCTA'};{'sGCTA_1_2'};{'sGCTA_1_3'}]

file_name='final.mat';
statistics_rings=zeros(length(monomers),3);
statistics_hexagon_percentage=zeros(length(monomers),2);

number_of_pics_to_use=15;
collecting_MNDs=zeros(size(monomers,1),number_of_pics_to_use);
collecting_hexagons=collecting_MNDs;
for monomer_number=1:length(monomers)
    load([path_main_folder,monomers{monomer_number},'/',file_name],'total_area','ring','particle_infos_list','frame_axis','ii','minimum_polygon','maximum_polygon','islands_infos_list','diameter');
    number_of_pictures_analysed=find(frame_axis==ii);
    axis_of_interest_for_histogram=[frame_axis(1:number_of_pictures_analysed)-0.5,frame_axis(number_of_pictures_analysed)+0.5];
    total_number_of_rings_per_frame=histcounts(ring(:,end),axis_of_interest_for_histogram);
    total_number_of_particle_per_frame=histcounts(particle_infos_list(:,end),axis_of_interest_for_histogram);
    total_number_of_hexagons=histcounts(ring(ring(:,3)==6,end),axis_of_interest_for_histogram);
    dens_particle=total_number_of_particle_per_frame./total_area(1:number_of_pictures_analysed)*diameter^2;
    dens_rings=total_number_of_rings_per_frame./total_area(1:number_of_pictures_analysed)*diameter^2;
    figure()
    scatter(dens_particle,dens_rings,'LineWidth',2)
    [P,S]=polyfit(dens_particle+rand(1,number_of_pictures_analysed)*min(dens_particle)/1e2,dens_rings,1);
    title(['Rings density in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('Rings density [n.u.]','Fontsize',12)
    yfit = polyval(P,dens_particle);          % Estimated  Regression Line
    hold on
    plot([0 sort(dens_particle) 5e-1],polyval(P,[0 sort(dens_particle) 5e-1]),'LineWidth',2)
    SStot = sum((dens_rings-mean(dens_rings)).^2);                    % Total Sum-Of-Squares
    SSres = sum((dens_rings-yfit).^2);                       % Residual Sum-Of-Squares
    Rsq = 1-SSres/SStot;                            % R^2
    xlim([0 4e-1]);
    ylim([0 2e-1])
    xl = xlim;
    yl = ylim;
    xt = 0.05 * (xl(2)-xl(1)) + xl(1);
    yt = 0.90 * (yl(2)-yl(1)) + yl(1);
    caption = sprintf('y = %f * x + %f\n R^2 = %f', P(1), P(2), Rsq);
    text(xt, yt, caption, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
    statstics_rings(ii,:)=[P(2), P(1),Rsq];
    
    dens_hexagons=total_number_of_hexagons./total_number_of_rings_per_frame;
    figure()
    scatter(dens_particle,dens_hexagons,'LineWidth',2)
    [P,S]=polyfit(dens_particle+rand(1,number_of_pictures_analysed)*min(dens_particle)/1e2,dens_hexagons,0);
    title(['Fraction of hexagons in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('Fraction','Fontsize',12)

    yfit = polyval(P,dens_particle);          % Estimated  Regression Line
    hold on
    plot([0 sort(dens_particle) 5e-1],polyval(P,[0 sort(dens_particle) 5e-1]),'LineWidth',2)
    SStot = sum((dens_hexagons-mean(dens_hexagons)).^2);                    % Total Sum-Of-Squares
    SSres = sum((dens_hexagons-yfit).^2);                       % Residual Sum-Of-Squares
    Rsq = 1-SSres/SStot;                            % R^2
    xlim([0 5e-1]);
    ylim([0 1])

    xl = xlim;
    yl = ylim;
    xt = 0.05 * (xl(2)-xl(1)) + xl(1);
    yt = 0.90 * (yl(2)-yl(1)) + yl(1);
    
    caption = sprintf('y = %f\n std = %f', P(1), std(dens_hexagons));
    text(xt, yt, caption, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
    statistics_hexagon_percentage(ii,:)=[P(1),Rsq];   
    %% Ideally you will want to do a vector that records the position of the maximum island in function of time
    % to do that you can use the time labelling and using particle number
    % as intensity masks in regionprops 
    frames_with_islands=unique(islands_infos_list(:,end));
    maximum_size_in_function_of_time=zeros(number_of_pictures_analysed,1);
    support_matrix=regionprops(islands_infos_list(islands_infos_list(:,4)>0,end),islands_infos_list(islands_infos_list(:,4)>0,4),'MaxIntensity','Area');
    maximum_size_in_function_of_time(frames_with_islands)=[support_matrix(:).MaxIntensity];
    density_islands=zeros(number_of_pictures_analysed,1);
    density_islands(frames_with_islands)=[support_matrix(:).Area]'./([total_area(1:number_of_pictures_analysed)]')*diameter^2;
    % one option that is memory intensive is to compare each element in
    % time and particles
    index_major_island=find(sum((islands_infos_list(:,end)==(frames_with_islands')).*(islands_infos_list(:,4)==(maximum_size_in_function_of_time(frames_with_islands)')),2));

    %% plot the maximum size in function of time
    figure()
    scatter(dens_particle,maximum_size_in_function_of_time,'LineWidth',2)
    title(['Size of maximum cluster in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('Number of monomers','Fontsize',12)

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
    figure()
    scatter(dens_particle,solidity_major_island,'LineWidth',2)
    title(['MND of the maximum cluster in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('MND','Fontsize',12)
  
    [P,S]=polyfit(dens_particle+rand(1,number_of_pictures_analysed)*min(dens_particle)/1e2,solidity_major_island,0);

    xlim([0 5e-1]);
    ylim([0 1])
%     yfit = polyval(P,dens_particle);          % Estimated  Regression Line
%     hold on
%     plot([0 sort(dens_particle) 5e-1],polyval(P,[0 sort(dens_particle) 5e-1]),'LineWidth',2)
%     SStot = sum((dens_hexagons-mean(dens_hexagons)).^2);                    % Total Sum-Of-Squares
%     SSres = sum((dens_hexagons-yfit).^2);                       % Residual Sum-Of-Squares
%     Rsq = 1-SSres/SStot;                            % R^2
    xlim([0 5e-1]);
    ylim([0 1])

    xl = xlim;
    yl = ylim;
    xt = 0.05 * (xl(2)-xl(1)) + xl(1);
    yt = 0.90 * (yl(2)-yl(1)) + yl(1);
%     
%     caption = sprintf('y = %f\n std = %f', P(1), std(dens_hexagons));
%     text(xt, yt, caption, 'FontSize', 14, 'Color', 'k', 'FontWeight', 'bold');
    figure(5)
    hold on
    scatter(maximum_size_in_function_of_time,solidity_major_island,'LineWidth',2)
    title('MND in function of maximum island size','FontSize',12);
    xlabel('Size [monomers]','Fontsize',12)
    ylabel('MND','FontSize',12)
collecting_MNDs(monomer_number,:)=solidity_major_island(end-number_of_pics_to_use+1:end);
collecting_hexagons(monomer_number,:)=dens_hexagons(end-number_of_pics_to_use+1:end);


  %% Computing the rings in maximum patch
    rings_in_maximum_patch_in_function_of_time=zeros(number_of_pictures_analysed,1);
    rings_in_maximum_patch_in_function_of_time(frames_with_islands)=islands_infos_list(index_major_island,5);

    figure
    scatter(dens_particle,rings_in_maximum_patch_in_function_of_time,'LineWidth',2)
    title(['Rings in maximum patch in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('Count','Fontsize',12)  
    %% Computing the fraction of rings that are on the boarder
    % this rely on the same approach as before, we need to find which rings
    % are part of the biggest island at the given time. we have the id of
    % the island for each time frame there are islands. Remember that cannot be a ring that is not an island...
    id_rings_in_maximum_patch=find(sum((ring(:,5)==(islands_infos_list(index_major_island,1)')).*(ring(:,end)==(islands_infos_list(index_major_island,end)')),2));
    % now for each time frame we have the rings in the patch... if we
    % multiply the flag of the boarder contained in column 4 with the time
    % frame and we divide by the histcounts of how many rings we have per
    % time frame in the maximum patch, we can get the fraction...

    rings_on_boarder_maximum_patch=histcounts(ring(id_rings_in_maximum_patch,4).*ring(id_rings_in_maximum_patch,end),axis_of_interest_for_histogram)';
    fraction_of_rings_on_boarder_maximum_patch=rings_on_boarder_maximum_patch./rings_in_maximum_patch_in_function_of_time;
    fraction_of_rings_on_boarder_maximum_patch(isnan(fraction_of_rings_on_boarder_maximum_patch))=0;
    if sum(isinf(fraction_of_rings_on_boarder_maximum_patch))
        fprintf('Error!')
    end
    figure
    scatter(dens_particle,fraction_of_rings_on_boarder_maximum_patch,'LineWidth',2)
    title(['Fraction of rings in maximum patch that are on boarder in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('Fraction','Fontsize',12)  
    
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
    index_internal_particles_maximum_patch=find(sum(((internal_particles_subset(:,3)+quantity_to_add(internal_particles_subset(:,end)))==((id_rings_in_maximum_patch'))).*(internal_particles_subset(:,end)==((ring(id_rings_in_maximum_patch,end))')),2));
    internal_particles_per_frame_maximum_patch=histcounts(internal_particles_subset(index_internal_particles_maximum_patch,end),axis_of_interest_for_histogram);
    fraction_of_particles_on_boarder_maximum_patch=internal_particles_per_frame_maximum_patch'./maximum_size_in_function_of_time;
    fraction_of_particles_on_boarder_maximum_patch(isnan(fraction_of_particles_on_boarder_maximum_patch))=0;
    % fraction of particles that are on a boarder in function of time
    figure
    scatter(dens_particle,fraction_of_rings_on_boarder_maximum_patch,'LineWidth',2)
    title(['Fraction of particles in maximum patch that are on boarder in function of particle density (',monomers{monomer_number},')'],'Fontsize',14)
    xlabel('Particle density [n.u.]','Fontsize',12)
    ylabel('Fraction','Fontsize',12)  

end
figure()
boxplot(collecting_MNDs','Labels',monomers)
ylim([0 1])
    title(['MND of largest island in function of design'],'Fontsize',14)
    xlabel('Design','Fontsize',12)
    ylabel('MND','Fontsize',12)

    figure()
boxplot(collecting_hexagons','Labels',monomers)
ylim([0 1])
    title(['Hexagons fraction in function of design'],'Fontsize',14)
    xlabel('Design','Fontsize',12)
    ylabel('Fraction','Fontsize',12)
