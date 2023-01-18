close all
clear all
%% Common variables
%% Set unit of measure, number of dimensions
SPACE='nm';
TIME='s';
dimensions=2;

to_m=1e-9; % from SPACE to meter
to_s=1; % from TIME to s

dt=256/100; % ideal time step in TIME between each picture, but not really used afterwards
fs=1/dt; % Hz sampling frequency, again not really used but always nice to check
diameter=12; % diameter within which only a particle can exist in SPACE units  12 for short 16 for long 

expected_DNA_width=7; % the DNA two arms width is approximatevely the geometric mean between the two arms diameter and the radius tip
particle_occupied_area=3*(diameter/2-expected_DNA_width*sqrt(3)/2/3)*expected_DNA_width+expected_DNA_width.^2*sqrt(3)/4; % eastimed occupied area

holes_per_particle=3;
%% variables to look for
minimum_polygon=4; % minimum polygon expected  (or of interest)
maximum_polygon=12; % maximum polygon expected  (or of interest)
additional_variables=5; % patch size, flag for internal patch, total number of rings, identity of the patch and  and flag for coalesced, to store together in a patch
polygons=[]; % where to store infos about the polygons in each frame
polygons_on_the_boarders=[]; % where to store infos about the polygons in each frame that are on the boarder of an island
total_number_of_particle=[]; % storing the total number of particles in function of time... just useful to have it ready already
ideal_number_of_rings=[]; % storing the ideal number of rings in function of time .... just useful to have ir ready already
total_area=[]; % storing infos about the total area in function of time... in case there is a change of the field of view
regions_number=[];% storing infos about the number of islands in function of time ... it is redundant but convenient to have
particle_infos_list=[]; % where to store infos about single particles in function of time
islands_infos_list=[]; % where to save infos about the single islands in function of time
solidity_measure=[]; % where to store the solidity in function of time ... not really used
time=[];

%% Setting up where the files are
experimental=true; % if we are dealing with simulations

% name_folder_original_r='/Volumes/One Touch/AFM/last_analysis/HS_AFM_definitive/lTATA/2021-04-27_exp0581/';
% name_folder_original_r='/Volumes/One Touch/AFM/last_analysis/HS_AFM_definitive/sGCGC/D/';
% name_folder_original_r='/Volumes/One Touch/AFM/last_analysis/HS_AFM_definitive/sGCTA/C/';
% name_folder_original_r='/Volumes/One Touch/AFM/last_analysis/HS_AFM_definitive/lGCTA/2021-04-27_exp11118/'
name_folder_original_r='/Volumes/One Touch/AFM/last_analysis/sGCGC/'
files=dir(name_folder_original_r); % look the content of the folder
dirFlags=[files.isdir]; % identify subdirectories
subFolders=files(dirFlags); % store the subdirectories
data_format='*.gwy'; % format of the pictures

subFolders(2)=[]; % needed to avoid going one level back from the main folder
start_folder=1;
end_folder=length(subFolders); % the maximum foler you may go for
%% Parameters for image treatment, for filtering

reference_ii=0;
time_series=true; % if we are dealing with images to compare in time, therefore concatenate in the time matrix


channel_fwd=0; % Where the forward height table is contained in the gwy file
% channel_bwd=3;% Where the backward height table is contained in the gwy file
channel_bwd=6;% Where the backward height table is contained in the gwy file
use_retrace=0; % wether to use the backward trace
use_line_corr=1; %wether to activate line correction
vertical_correction=1; %wether to use median alignment along columns
periodic_noise=false; % wether to activate periodic noise correction

recount=true;

resizing=0;

expx=4; % maximum grade of polynomial for polynomial fitting of the plane... x direction
expy=4; % maximum grade of polynomial for polynomial fitting of the plane... x direction

expected_populations=2; % how many popoulations we expect in the image: DNA and background in our case

tarres=0.33; % target resolution in nm, magnify images to this resolution
ws = warning('off','all');  % Turn off warning
% thr_accumulator_reference=0; % intialization of a reference threshold ... most probably not used

% thr_accumulator_reference=4.980475710800370e5; 
% thr_accumulator_reference=4.2e5;
% thr_accumulator_reference=9.731910511775464e5; 
% thr_accumulator_reference=7.6e5; % intialization of a reference threshold ... most probably not used
thr_accumulator_reference=0;
%% Colors and axes for histograms
colors_concave=[216,179,101]/255;
colors_map=lines(maximum_polygon); % automatic colors for the different polygons out there
colors_map(4:10,:)=[79, 91,102; 8, 81, 113; 49, 136, 189; 150, 163, 172; 232, 142, 133; 223, 78, 113; 213, 32, 75;]/255; % specific colors for quadrilaters to heptagon
poly_x_axes=categorical({'Quadrilateral','Pentagon','Hexagon','Heptagon','Octagon','Ennagon','Decagon'}); % intializing the axis for histograms
poly_x_axes=reordercats(poly_x_axes,{'Quadrilateral','Pentagon','Hexagon','Heptagon','Octagon','Ennagon','Decagon'}); % ordering the axis
load('gwyddion_color_map.mat'); % load the colormap

color_value=2000; % value to intialize the colors for the islands... not really thought this through
map_color_patches= colormap(jet(color_value)); % prepare the pallet of colors
map_color_patches=[0,0,0; map_color_patches(randperm(color_value),:);1,1,1]; % randomly order it 
offset_base=0.2; % offset for aligning the graphs
%% open the files and analysing
started=false; % do not change... it is foundamental to intialize some variables we cannot intialize now.
new_identity=1; % do not change ... it is intialized for the part of islands identification
identity=1; % do not change ... it is intialized for the part of islands identification
sense=1; % how to read the folder, positive goes forward, negative goes backwards in time
for sf=start_folder:end_folder %cycle over each of the folders we are interested in
    name_folder_r=[name_folder_original_r,subFolders(sf).name,'/']; %target folder
        image_analysis_experiments % the script for the image analysis
end

% Sampling_time=number_of_pictures*dt; % total time spent for the image sampling in TIME unit
% 
% %% Parameters to reconstruct the trajectories (track in traj_construction)
% % Define parameters for track
% param.mem =0; % this is the number of time steps that a particle can be 'lost' and then recovered again
% param.good=0; %Discard every trajectory less than param.good
% param.quiet=0; %Show information
% param.dim=dimensions; % set dimension: 2 for 2D
% 
%                                 
% %% Constructing the trajectories
% traj_construction
% 
% %% Using MSD analyser the output is analysis
% MSD_analysis
