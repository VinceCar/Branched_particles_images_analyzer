%% Open afm image
% acquiring the data

%some times the channel may not be were we expect it to be but later
%(especially for one subset of images) in case it is not there, look at
%the channel afterwards
non_segmented_image=0;
non_segmented_image=looking_for_channel(non_segmented_image,[name_folder_r,list_of_pictures(ii).name],channel_fwd);
non_segmented_image.data(1,:)=[]; % remove first line because it is often affected by artefacts

%% Dimensional informations extracted from the gwy file
[dim_y_o, dim_x_o]=size(non_segmented_image.data); % Getting size of the image
if use_line_corr
    % first line correction, align based on the median of differences
    % between lines
    non_segmented_image.data=non_segmented_image.data-polynomial_fit_1(true(dim_y_o-1,dim_x_o),cumsum([zeros(1,size(non_segmented_image.data,2));median(diff(non_segmented_image.data),2).*ones(size(non_segmented_image.data,1)-1,size(non_segmented_image.data,2))]),0,1); 
end
L=non_segmented_image.xreal/to_m/dim_x_o; % spatial resolution in x SPA
% CE/pixel
W=non_segmented_image.yreal/to_m/(dim_y_o+1); % spatial resolution in y, take into account you removed a pixel!

non_segmented_image.data=non_segmented_image.data/to_m; % conver to m

%% Starting a first processing of the image, flattening the background
image_to_process=normalize_matrix(non_segmented_image.data); %  set 0 the minimum and normalize
image_to_process=polynomial_fit_1(ones(size(image_to_process))>0,image_to_process,1,1); % apply a plane removal by fitting all the points
image_to_process=normalize_matrix(image_to_process); % renomarilze
mask_b=(image_to_process<adaptthresh(image_to_process)); % prepare a preliminary mask through adaptthresh. 1 is pixels that are probably background
image_to_process=non_segmented_image.data; % return to work the original data
if expx>=0 || expy>=0
    image_to_process=polynomial_fit_1(mask_b,image_to_process,expx,expy); % do a polynomial removal with the exponents given in the main. use for the fitting only the point under the mask
    if vertical_correction
        % remove the median of the backgrounds
        non_segmented_image.data=image_to_process-median(image_to_process);
        % now repeat the fitting using image_to_process
        image_to_process=normalize_matrix(non_segmented_image.data); %  set 0 the minimum and normalize
        image_to_process=polynomial_fit_1(ones(size(image_to_process))>0,image_to_process,1,1); % apply a plane removal by fitting all the points
        image_to_process=normalize_matrix(image_to_process); % renomarilze
        mask_b=(image_to_process<adaptthresh(image_to_process)); % prepare a preliminary mask through adaptthresh. 1 is pixels that are probably background
        image_to_process=non_segmented_image.data; % return to work the original data
        prct_of_interest=sum(mask_b)/size(mask_b,1)*100/2; % find the percentage of the background and divide by 2 now that we have a more reliable mask
        vector_of_medians=diag(prctile(image_to_process,prct_of_interest))'; % each median of the background per column
        image_to_process=image_to_process-vector_of_medians;

        image_to_process=polynomial_fit_1(mask_b,image_to_process,expx,expy); % do a polynomial removal with the exponents given in the main. use for the fitting only the point under the mask
    end
end

non_segmented_image.data=image_to_process; % store it as image with background 
%% resizing to desired resolution
dim_y=dim_y_o*ceil(W/tarres); % pixels in final y axis
dim_x=dim_x_o*ceil(L/tarres); % pixels in final x axis
non_segmented_image.data=imresize(non_segmented_image.data,[dim_y,dim_x],'Method','lanczos3'); % actual resizing    
L=L*dim_x_o/dim_x; % new pixel resolution
W=W*dim_y_o/dim_y; % new pixel resolution


%% furhter flattening of background subtracting the median of differences of the background... consider that the expected foreground is wide we want to remove fast spikes
image_to_process=normalize_matrix(non_segmented_image.data);
image_to_process=median_line_correction_of_background(image_to_process,non_segmented_image.data,expy,'differences');
non_segmented_image.data=image_to_process; % storing the variabel in the non segmented

%% %% furhter flattening of background subtracting the median of the background... consider that the expected foreground is wide, we want to remove fast spikes
image_to_process=normalize_matrix(non_segmented_image.data);
image_to_process=median_line_correction_of_background(image_to_process,(non_segmented_image.data),expy,'none');

%% Further smooth the image
image_to_process=imgaussfilt(image_to_process,2); % preliminary gaussian filtering to remove discontinuities

% cap 2 percent below and 2 percent on top
low_limit=prctile(image_to_process(:),2);
image_to_process(image_to_process<low_limit)=low_limit;
%% support functions
function image_to_process=normalize_matrix(input_matrix)
    % just to put minimum to 0 and rescale maximum to 1
    image_to_process=input_matrix-min(input_matrix(:)); % minimum to 0
    image_to_process=input_matrix./max(input_matrix(:)); % normalize to max
end

function original_image=looking_for_channel(original_image,path_to_file,intial_channel)
    %some times the channel may not be were we expect it to be but later
    %(especially for one subset of images) in case it is not there, look at
    %the channel afterwards
    while isstruct(original_image)==false
        original_image=readgwychannel(path_to_file,intial_channel);
        intial_channel=intial_channel+1;
        fclose all;
    end
end

function image_to_process=median_line_correction_of_background(input_image,matrix_to_calculate_the_median_about,exp_y,type_of_correction)
    % This function perform median line correction, based on the inputs you
    % can do 2 subtypes of median line corrections
    % INPUTS:
    % input_image, matrix bounded between 0 and 1, it contains the image we
    % want to process
    % matrix_to_calculate_the_median_about, as per name it contains the
    % matrix that we want to use to compute the median to subtract, use the
    % not normalized values of input_image 
    % type_of_correction, if equal to difference perform median of
    % differences line correction otherwise median differences
    % OUTPUTS:
    % image_to_process final result
    image_to_process=matrix_to_calculate_the_median_about;
    yy=(1:size(matrix_to_calculate_the_median_about,1))'-1; % support variable to fit the line to add obtained by vec_removal to a polynomial of order given by exp_y 
    mask_b=(input_image<adaptthresh(input_image)); % segmentation of background
    mask_b=mask_b-bwmorph(mask_b,'remove'); % be sure of having background by discarding the pixels at the interface
    if strcmp(type_of_correction,'differences')
        matrix_to_calculate_the_median_about=diff(matrix_to_calculate_the_median_about); % matrix of differences line by line
        mask_b=mask_b(2:end,:); % mask of the difference matrix
    end
    matrix_to_calculate_the_median_about(logical(1-mask_b))=max(matrix_to_calculate_the_median_about(:))+1; % avoid foreground to be taken into account
    vec_removal=(sum(mask_b,2))*100/2/size(matrix_to_calculate_the_median_about,2); % find what percentage of the pixels per line is occupied by background. divide by 2 to obtain the location of the median of the background
    % what is done afterwards is simply to avoid a for over the values of
    % vec_removal 
    vec_removal_2=unique(sort(vec_removal)); % make it ordered
    support=prctile(matrix_to_calculate_the_median_about,vec_removal_2,2); % compute the given percentiles per all the lines
    support(:,vec_removal_2==0)=0; % makes sure that the 0 percentile starts at 0!
    vec_removal=[(sum(support.*(vec_removal_2'==vec_removal),2))];
    if strcmp(type_of_correction,'differences')
        vec_removal=[0;cumsum(sum(support.*(vec_removal_2'==vec_removal),2))];  % in each line find where the percentile is the one we are interested in and consider that the difference between the lines needs to be integrated to have the flatten image
    end
    % now vec_removal is a column vector reporting the median of background
    % per each raw of matrix_to_calculate_the_median about
    support=fliplr(yy.^([1,2:exp_y])); % support variable to do the fitting
    coeff=support\vec_removal; % compute the coefficients
    vec_removal=vec_removal-sum(coeff'.*support,2); % subtract the fitted line 
    image_to_process=image_to_process-vec_removal;
end