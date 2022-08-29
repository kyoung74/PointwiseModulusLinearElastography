%% This program calculates the pointwise modulus (PWM) for each cell for each cell region (lamellipodia or perinuclear)
clc;clear; %clear workspace and command window

%% Load information about each forcemap from a master spreadsheet
%Make sure the current directory is the one with the master spreadsheet
sample_info = readtable(''); %Populate this with spreadsheet directory 
sample_names = string(sample_info.SampleName); %name of sample file
sample_xdim = double(sample_info.xDimLength_um_); %vertical range of forcemap
sample_ydim = double(sample_info.yDimLength_um_); %horizontal range of forcemap
sample_exclude = double(sample_info.ExcludedPoints); %1 = some points were excluded, 0 = no points excluded
sample_indenter_radius = double(sample_info.IndenterRadius);
sample_align_R = double(sample_info.AligningRMatrix);

%% Iterate through each sample
for aa = 1:length(sample_names) + 1


    sample_name = sample_names(aa);
    
    dir_forcemap_csv = ''; %POPULATE with directory for file containing the forcemap data
    dir_code_results = sprintf('/%s',sample_name); %POPULATE with directory for folder you want results in
    dir_FvI_txts = sprintf('/%s/%s',sample_name,sample_name); %POPULATE with directory of the Force vs Indentation CSV file
    dir_sample_info = ''; %POPULATE with the directory the sample info CSV file is in
    
    cd(dir_forcemap_csv)
    csv_data = readtable(sprintf('%s.csv',sample_name)); %load the foremap CSV file
    
    %check if column 'converged' consists of either 0/1 or TRUE/FALSE, then
    %convert the data from a table into a matrix of doubles
    converged = string(csv_data.converged);
    if ismember(converged,'TRUE') %if convergence is indicated by "TRUE" or "FALSE", convert to 1 or 0
        for zz = 1:length(converged)
            if converged(zz) == 'TRUE'
                converged(zz) = 1;
            elseif converged(zz) == 'FALSE'
                converged(zz) = 0;
            end
        end
        dub_data = [table2array(csv_data(:,1:6)) str2double(converged) table2array(csv_data(:,8:9))];
    else %if convergence is already indicated by 1 or 0
        dub_data = table2array(csv_data);
    end
    
    % calculate height of the cell
    maxcp = max(dub_data(:,3)); %get the maximum height of the cell
    og_height = maxcp - dub_data(:,3); %the original height without any filtering
    
    mask1 = dub_data(:,7) == 0; %locations of nonconverged points
    mask2 = og_height < 2.5e-6; %any height less than 2.5E-6 is considered glass
    
    height = og_height;
    height(mask1) = 0; %truncate nonconverged points to zero
    height(mask2) = 0; %truncate height of glass to zero
    
    %Exclude points that don't meet criteria
    if sample_exclude(aa) == 1 
    
        cd(dir_sample_info)
        excluded_pts = xlsread('Excluded Points',sprintf(sample_name)); %load the excel file detailing coordinates of excluded points
        [num_pts, ~] = size(excluded_pts);
        exclude = nan(num_pts,1);
        
        %find coordinates for excluded points
        for qq = 1:num_pts
            exclude(qq) = find(dub_data(:,8) == excluded_pts(qq,1) & dub_data(:,9) == excluded_pts(qq,2) & dub_data(:,7) == 1);
        end
    
        height(exclude) = 0; %height of points to exclude is considered 0
    end
    
    % arrange txt file in Line, Point, and Height string matrix
    A = [dub_data(:,8) dub_data(:,9) height];
    
    % get number of rows and columns 
    rows = max(A(:,1)) + 1; %number of rows
    cols = max(A(:,2)) + 1; %number of columns
    
    og_height_mat = nan(rows,cols);
    height_mat = nan(rows,cols);
    
    % For troubleshooting: check to see if indices in csv file are in order
    for i = 1:length(dub_data)-1
        if dub_data(i+1,1) == dub_data(i,1) + 1
        else
            break
        end
    end
    
    %For troubleshooting: check to see which line has an unequal number of points
    %This suggests that you forgot to exclude a point that didn't converge
    %and helps you locate it.
    unequal_lines = false;
    for i = 0:rows-1
        num_of_lines = size(find(dub_data(:,8) == i));
        if num_of_lines < rows
            unequal_lines = true;
        end
    end
    
    height_mat(1) = 0;
    for bb = 1:rows*cols
        if dub_data(bb,7) == 0
            height_mat(dub_data(bb,8)+1,dub_data(bb,9)+1) = 0;
        else
            og_height_mat(dub_data(bb,8)+1,dub_data(bb,9)+1) = og_height(bb);
            height_mat(dub_data(bb,8)+1,dub_data(bb,9)+1) = height(bb);
        end
    end
    

    figure(9)
    imagesc(og_height_mat)
    set(gca, 'YDir','reverse')
    
    figure(10)
    imagesc(height_mat)
    
    
    % value of ticker increments (size of afm measurement (in m) / number of points measured)
    val_xaxis = sample_xdim(aa)*10^-6 / (cols - 1);
    
    % final vector of all values
    A_fin = A;
    
    % Multipy the lines and points by their actual x and y value (value in
    % meters) and add to vector
    first_col_updated = val_xaxis .* A_fin(:,2);
    % convert Line indices (start at 0) to Row indices (start at 1)
    ind_line = A_fin(:,1) + 1;
    % convert Point indices (start at 0) to Column indices (start at 1)
    ind_point = A_fin(:,2) + 1;
    % add line indices to end of matrix
    A_fin = [A_fin(:,1) first_col_updated A_fin(:,3) ind_line];
    % Multipy the lines and points by their actual x and y value (value in
    % meters) and add to vector
    val_yaxis = sample_ydim(aa)*10^-6 / (rows - 1); %45 is length of scan, y axis in microns, 15 is # of boxes
    second_col_updated = val_yaxis .* A_fin(:,1);
    % add indices to end of matrix
    A_fin = [second_col_updated A_fin(:,2) A_fin(:,3:end) ind_point]; 
    
    
    %A_fin columns are [col [m], row [m], height [m], col [index] row [index]
    
    %% Calculate Radii of curvature
    
    %start with finding radii of curvature in the x direction
    Radius_x = []; % empty Rx matrix - to be filled in with the Rx value for each point on the cell
    
    start_sel = 1; % select first value of col
    end_sel = cols; % select last value of col
    
    for t = 1: rows %measure radius one row at a time
        radius_x_complete = [];
    
        first_line = A_fin(start_sel:end_sel,2:end); %load data for the current row (row t) ss - I think this should be called 'first_row', and should only load relevant data
    
        % de-serpentine the pattern
        if first_line(1,4) == cols %if the first col in row t is = to # of cols, then row is out of order
            first_line = flip(A_fin(start_sel:end_sel,2:end)); %flip row to be in column order
        end
    
    % %   find radii of curvature
        for j = (2):(cols - 1) %don't iterate through edgepoints because radius cannot be measured there
            syms x y
         % %   if both the height and the start location are 0 for the three consecutive values automatically make radius = 0
            if first_line(j,1) == 0 & first_line(j,2) == 0 & first_line(j+1,1) == 0 & first_line(j+1,2) == 0 & first_line(j-1,1) == 0 & first_line(j-1,2) == 0
                radius_x = 0;
                radius_x_complete = [radius_x_complete radius_x];
            elseif  (first_line(j,1) == 0 & first_line(j+1,1)==0 & first_line(j-1,2)==0) | (first_line(j,2) == 0 & first_line(j+1,2)==0 & first_line(j-1,2)==0)
                     radius_x = 0;
                radius_x_complete = [radius_x_complete radius_x];
            else
    % %     find eq of circle fit to 3 points
            eq1 = (first_line(j,1) - x)^2 + (first_line(j,2)- y)^2;
            eq2 = (first_line(j+1,1) - x)^2 + (first_line(j+1,2)- y)^2;
            eq3 = (first_line(j-1,1) - x)^2 + (first_line(j-1,2)- y)^2;
            [solx,soly] = vpasolve(eq1 == eq2, eq1 == eq3, [x y]);
            if length(solx) == 0 | length(soly) == 0
                radius_x = 0;
                radius_x_complete = [radius_x_complete radius_x];
                continue
            end
            diff = (first_line(j,1) - solx)^2 +(first_line(j,2)- soly)^2;
            radius_x = sqrt(diff);
        %     place each radius value in apprpriate point in the line
            radius_x_complete = [radius_x_complete radius_x];
            end
        end
          %  add line to overall Rx matrix
        Radius_x = [Radius_x; radius_x_complete];
    
    %     check to see if you are done evaluating all data points
        if end_sel == length(A_fin) %length of sorted_matrix for cols
            start_sel = 1;
            end_sel = cols; %cols
        else
            start_sel = end_sel + 1;
            end_sel = start_sel + cols - 1; %cols - 1
        end
    end
    
    % sort matrix based on point
    sorted_matrix = [];
    for g = 1:cols %cols
    %     start from point = 1
      ind = find(A_fin(:,5) == g);
    %   get all rows where point = 1
      sort_m = A_fin(ind,:);
    %   Create matrix sorted on points
      sorted_matrix = [sorted_matrix; sort_m];
    end
    
    % create Radius in the y direction matrix
    Radius_y = [];
    % select first value of row
    start_sel = 1;
    % select last value of row
    end_sel = rows; %rows
    for t = 1: cols
        radius_y_complete = [];
        % %  select first line in order to measure the radius
        first_line = sorted_matrix(start_sel:end_sel,:); %this sould be called first point, or first column
    % %   find radius
            for h = (2):(rows - 1) %rows - 1
                if t == 1
                    disp('stop')
                end
                syms x y
                 % %   if both the height and the start location are 0 for the three consecutive values automatically make radius = 0
                if first_line(h,1) == 0 & first_line(h,3) == 0 & first_line(h+1,1) == 0 & first_line(h+1,3) == 0 & first_line(h-1,1) == 0 & first_line(h-1,3) == 0
                    radius_y = 0;
                    radius_y_complete = [radius_y_complete; radius_y];
                elseif (first_line(h,1) == 0 & first_line(h+1,1) == 0 & first_line(h-1,1) == 0) | (first_line(h,3) == 0 & first_line(h+1,3) == 0 & first_line(h-1,3) == 0)
                    radius_y = 0;
                    radius_y_complete = [radius_y_complete; radius_y]; %i found an error where this was [radius_y radius_y_complete]. This made the Radius_y vector look weird0
                else
    % %     find eq
                eq1 = (first_line(h,1) - x)^2 + (first_line(h,3)- y)^2;
                eq2 = (first_line(h+1,1) - x)^2 + (first_line(h+1,3)- y)^2;
                eq3 = (first_line(h-1,1) - x)^2 + (first_line(h-1,3)- y)^2;
                [solx,soly] = vpasolve(eq1 == eq2, eq1 == eq3, [x y]);
                diff = (first_line(h,1) - solx)^2 +(first_line(h,3)- soly)^2;
                radius_y = sqrt(diff);
    %     place each radius value in apprpriate point in the line
                radius_y_complete = [radius_y_complete; radius_y]; %adds to the beginning
                end
            end
          %  add line to overall Rx matrix
        Radius_y = [Radius_y radius_y_complete];
    %     check to see if you are done evaluating all data points
        if end_sel == length(A_fin) %sorted matrix for cols
        start_sel = 1;
        end_sel = rows; %rows
        else
        start_sel = end_sel + 1;
        end_sel = start_sel + rows - 1; %rows - 1
        end
    end
    
    % fix aspect ratio of radii of x and y (i.e zoom in)
    Radius_x = Radius_x(2:end-1,:); 
    Radius_y = Radius_y(:,2:end-1);
    % get the reciprocal
    One_ov_Radius_x = 1 ./ Radius_x;
    One_ov_Radius_y = 1 ./ Radius_y;
    % the sample surface curvature calculation
    One_ov_Radius_s = (One_ov_Radius_x + One_ov_Radius_y) ./ 2;
    % get indenter radius and then take reciprocal
    indenter_radius = sample_indenter_radius(aa);
    one_ov_indenter = ( 1 ./ indenter_radius);
    % get final radius of sample
    one_ov_R = one_ov_indenter + One_ov_Radius_s;
    R = double((1 ./ one_ov_R));
    
    % create matrix based on hertzian 2 sphere model
    [X,Y] = meshgrid(1:size(R,2), 1:size(R,1));
    [X2,Y2] = meshgrid(1:0.01:size(R,2), 1:0.01:size(R,1));
    outData = interp2(X, Y, R, X2, Y2, 'linear');
    figure(1)
    imagesc(outData);
    set(gca, 'XTick', linspace(1,size(X2,2),size(X,2)));
    set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
    set(gca, 'XTickLabel', 1:size(X,2));
    set(gca, 'YTickLabel', 1:size(X,1));
    ax =gca;
    ax.YDir = 'normal';
    h = colorbar;
    colormap(pink(512));
    axis image
    xlabel('?m')
    ylabel('?m')
    ylabel(h, 'height (?m)')
    
    %%
    
    figure(2)
    [C2,h2]= imcontour(flipdim(outData,1));
    xlabel('?m')
    ylabel('?m')
    
    %%
    
    figure(3)
    surfl(outData)
    colormap(pink(512))
    shading interp
    hold on
    imagesc(outData)
    xlabel('?m')
    ylabel('?m')
    zlabel('?m')
    %%
    figure(4)
    
    [new_r,new_c] = size(A);
    Mat = [];
    new_column =  1;
    for g = 1:new_r
        Line_no = A(g,new_column) + 1;
        Point_no = A(g,new_column + 1) + 1;
        height = A(g, new_column + 2);
        Mat(Line_no,Point_no) = height;
    end
    Mat = Mat * 1e6;
    [X,Y] = meshgrid(1:size(Mat,2), 1:size(Mat));
    [X2,Y2] = meshgrid(1:0.01:size(Mat,2), 1:0.01:size(Mat));
    outData2 = interp2(X, Y, Mat, X2, Y2, 'linear');
    maskstat = find(outData2 > 1.05e-5);
    imagesc(outData2);
    set(gca, 'XTick', linspace(1,size(X2,2),size(X,2))); 
    set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
    set(gca, 'XTickLabel', 1:size(X,2));
    set(gca, 'YTickLabel', 1:size(X,1));
    ax =gca;
    ax.YDir = 'normal';
    h = colorbar;
    colormap(pink(512));
    axis image
    ylabel('Line')
    ylabel(h, 'point')
    cd(dir_code_results)
    saveas(gcf,'heightcolormap.png')
    
    figure(5)
    contour(outData2,'ShowText','on')
    axis image
    cd(dir_code_results)
    saveas(gcf,'contourmap.png')
    %%
    
    figure(6)
    surfl(outData2)
    colormap(pink(512))
    shading interp
    saveas(gcf, 'celltopography.png')
    
    %%
    cd(dir_FvI_txts)
    files = dir('*.txt'); %load all .txt file names
    filmat =[];
    NaN_spots = sum(dub_data(:,7) == 0); %count the amount of NaN spots (spots where force curve didn't converge)
    l = length(files); %length should be #rows*#cols
    Mat_1 = [0];
    R_vec = [];
    length_R = size(R,1);
    
    %Rearrange the R matrix so that it aligns with the PWM calculation. Each
    %matrix is aligned differently, so that's why there is an if statement with
    %different alignment procedures depending on the sample
    align_R = sample_align_R(aa);
    
    if align_R == 1
        for qq = 1:length_R 
            R_v = R(qq,:)';
            R_vec = [R_vec; R_v];
        end
    elseif align_R == 2
        for qq = length_R:-1:1
            R_v = R(qq,:)';
    %         R_v = flip(R_v);
            R_vec = [R_vec; R_v];
        end
    elseif align_R == 3
        for qq = 1:length_R 
            R_v = R(qq,:)';
            R_vec = [R_vec; R_v];
        end
    elseif align_R == 4
        for qq = 1:length_R
            R_v = R(qq,:)';
    %         R_v = flip(R_v);
            R_vec = [R_vec; R_v];
        end
    else
        for qq = 1:length_R 
            R_v = R(qq,:)';
            R_v = flip(R_v);
            R_vec = [R_vec; R_v];
        end
    end

    %collect all the.txt filenames with Force vs Indentation data into one vector 
    for ii=1:l - NaN_spots
        names = files(ii).name ;
        filmat = [filmat; {names}];
    end
    
    Filmat = []; % to be populated with .txt filenames excluding those on the border
    for jj=1:l - NaN_spots %change depending on the size
        A = filmat{jj};
    
        %skip points on the edge of the forcemap (where R could not be calculated)
        if contains(A,'Point0000')|contains(A,sprintf('Point00%d',cols - 1)) %ss-edit based on the size
            continue
        elseif contains(A,sprintf('Line00%d',rows - 1))|contains(A,'Line0000') %ss-edit based on the size
            continue
        else
          Filmat = [Filmat;{A}];
        end
    end
    
    %% Calculate the PWM for every point on the cell

    %Commented code in this section is for trobleshooting, especially for
    %when pwm calculation results in "NaN"

    %check to see if the current point is NaN
    % NaN_Ind = find(dub_data(:,7) == 0);
    % NaN_Coord = [dub_data(NaN_Ind,8) dub_data(NaN_Ind,9)];
    % 
    % NaN_Track = 0;
    for rv = 1:length(Filmat) %iterate through .txt filenames
        A = Filmat{rv};
    % 
    %     line = str2double(A(7:8));
    %     point = str2double(A(16:17));
    %     
    %     for NN = 1:length(NaN_Coord)
    %         if line-1 == NaN_Coord(NN,1) && point-1 == NaN_Coord(NN,2)
    %             NaN_Track = NaN_Track + 1;
    % 
    %             prev_NaN = 1;
    %             while prev_NaN > 0 
    %                 if line-prev_NaN-1 == NaN_Coord(NN,1) && point-prev_NaN-1 == NaN_Coord(NN,2)
    %                     NaN_Track = NaN_Track + 1;
    %                     prev_NaN = prev_NaN + 1;
    %                 else
    %                     prev_NaN = 0;
    %                 end
    %             end
    %         end
    %     end
    
    %     if A == 'Line0014Point0011.txt'
    %         disp(' ')
    %     end

    %     if statement for troubleshooting individual force v indent. curves
        if rv == 165 %force curve you want to stop at
            disp('put a stopper on this line') %put a stopper on this line to pause code
        end

    
        res = importdata(A);
        C = res.data;
        fun_call = @pwmreader;
         O_P = fun_call(C,R_vec(rv)); %run pwmreader function (see at the end of code) with inputs C (force v indentation data) and R (calculated radii of curvature)
    %      O_P = fun_call(C,R_vec(rv + NaN_Track));
        out = str2double(regexp(A,'[\d.]+','match')); 
        row = out(1);
        column = out(2);
        Mat_1(row,column) = O_P; %insert calculated PWM into the appropriate spot in the PWM matrix
    end
    
    Mask = (Mat_1 < 0 );
    Mat_1(Mask) = 0; %any spots with negative PWM values are considered 0
    Mat_1 = Mat_1 .* 1e6; %unit conversion
    
    cd(dir_code_results)
    csvwrite(sprintf('%s_PWM_map.csv',sample_name),Mat_1) %save the PWM map as a CSV file
    
    figure(7)
    subplot(2,2,1)
    find(double(Radius_x > 0))
    imagesc(double(Radius_x))
    title('Rx')
    subplot(2,2,2)
    imagesc(double(Radius_y))
    title('Ry')
    subplot(2,2,3)
    imagesc(R)
    title('R')
    subplot(2,2,4)
    imagesc(Mat_1)
    title('PWM')
    saveas(gcf,'ROC.png')
    
    [X,Y] = meshgrid(1:size(Mat_1,2), 1:size(Mat_1,1));
    [X2,Y2] = meshgrid(1:0.01:size(Mat_1,2), 1:0.01:size(Mat_1,1));
    outData3 = interp2(X, Y, Mat_1, X2, Y2, 'linear');
    figure(8)
    imagesc(outData3);
    set(gca, 'XTick', linspace(1,size(X2,2),size(X,2)));
    set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
    set(gca, 'XTickLabel', 1:size(X,2));
    set(gca, 'YTickLabel', 1:size(X,1));
    ax =gca;
    ax.YDir = 'normal';
    h = colorbar;
    colormap(jet(512));
    axis image
    ylabel('Line')
    ylabel(h, 'Pointwise Modulus')
    cd(dir_code_results)
    saveas(gcf,'pwmcolormap.png')
    
    %% Determine points in the lamelopodial and perinuclear regions
    figure(11)
    subplot(1,2,1)
    contour(outData2,linspace(0,max(max(outData2)),10),'ShowText','on') %create a contour plot
    subplot(1,2,2)
    imagesc(outData2)
    cd(dir_code_results)
    saveas(gcf,'contour_and_PWM.png') 
    
    %see if perinuclear bound has already been set. If not, ask user to set it
    %based on the contour lines
    if sample_info.ContourLineBoundary(aa) > 0
        perinuc_bound = sample_info.ContourLineBoundary(aa);
    else
        disp(linspace(0,max(A_fin(:,3)),10))
        perinuc_bound = input('Which line # is the outermost line of the two lines that have the greatest distance between them?');
    end
    
    
    div_height = linspace(0,max(A_fin(:,3)),10); %vector of heights of contour lines 
    
    %find perinuclear points (height > contour boundary)
    nuc_points = find(A_fin(:,3) >= div_height(perinuc_bound) & A_fin(:,4) > 0 & A_fin(:,4) < rows-1 & A_fin(:,5) > 0 & A_fin(:,5) < cols-1);
    nuc_rows = A_fin(nuc_points,4);
    nuc_cols = A_fin(nuc_points,5);
    
    %find lamelipodial points (height < contour boundary)
    lamel_points = find(A_fin(:,3) < div_height(perinuc_bound) & A_fin(:,3) > 0 & A_fin(:,4) > 0 & A_fin(:,4) < rows-1 & A_fin(:,5) > 0 & A_fin(:,5) < cols-1);
    lamel_rows = A_fin(lamel_points,4);
    lamel_cols = A_fin(lamel_points,5);
    
    lamel_or_nuc = nan(size(Mat_1)); %matrix showing which points were labeled perinuclear (0) or lameipodial (1)
    
    %For visualization of sorting of points into lamellipodial or
    %perinuclear

    %If a point is in the nuclear region, visualize it as a '0' in a matrix
    counter = 1;
    for dd = 1:length(nuc_points)
        lamel_or_nuc(nuc_rows(dd),nuc_cols(dd)) = 0;
         if Mat_1(nuc_rows(dd),nuc_cols(dd)) > 0
            pwm_nuc(counter) = Mat_1(nuc_rows(dd),nuc_cols(dd));
            nuc_rows_filt(counter) = nuc_rows(counter);
            nuc_cols_filt(counter) = nuc_cols(counter);
            counter = counter + 1;
        end
    end
    
    %If a point is in the lamel. region, visualize it as a '1' in a matrix
    counter = 1;
    for cc = 1:length(lamel_points)
        lamel_or_nuc(lamel_rows(cc),lamel_cols(cc)) = 1;
        if Mat_1(lamel_rows(cc),lamel_cols(cc)) > 0
            pwm_lamel(counter) = Mat_1(lamel_rows(cc),lamel_cols(cc));
            lamel_rows_filt(counter) = lamel_rows(counter);
            lamel_cols_filt(counter) = lamel_cols(counter);
            counter = counter + 1;
        end
    end
   
    %% Calculate values of interest from PWM data and save
    avg_pwc_nuc = mean(pwm_nuc); %average PWM of nuclear region
    std__pwc_nuc = std(pwm_nuc); %standard deviation of PWMs of nuclear region
    
    avg_pwc_lamel = mean(pwm_lamel); %average PWM of lamellipodial region
    std_pwc_lamel = std(pwm_lamel); %standard deviation of PWMs of lamellipodial region

    %load calculated values into the sample_info spreadsheet
    sample_info.PerinuclearPWM(aa) = avg_pwc_nuc;
    sample_info.PeriuclcearPWMSTDev(aa) = std__pwc_nuc;
    sample_info.LamelopodialPWM(aa) = avg_pwc_lamel;
    sample_info.LamelopodialPWMSTDev(aa) = std_pwc_lamel;
    sample_info.ContourLineBoundary(aa) = perinuc_bound;
    
    cd(dir_sample_info)
    writetable(sample_info,'Force Map Avg PWMs.xlsx') %save updated spreadsheet
    
    %load individual PWM values from perinuclar and lamellipodial regions
    %with their locations into a spreadsheet
    cd(dir_code_results)
    ind_nuc_pwms = table(pwm_nuc',nuc_rows_filt',nuc_cols_filt','VariableNames',{'PWM','Row numer','Column number'});
    ind_lamel_pwms = table(pwm_lamel',lamel_rows_filt',lamel_cols_filt','VariableNames',{'PWM','Row numer','Column number'});
    writetable(ind_nuc_pwms,'Individual Perinuclear PWM Values.xlsx')
    writetable(ind_lamel_pwms,'Individual Lamellipodia PWM Values.xlsx')

end

%% pwmreader function - calculates PWM at each point of cell
function [out] = pwmreader(in1,in2) %in1: force/indentation data, in2: radius of curvature
v=0.5;
R = in2;
in1 = rmmissing(in1); %remove missing rows or columns
ind= in1(1:end,1); %indentation data
force= in1(1:end,2); %force data
% the following codes determine the point-wise Young's modulus
pntindex=find(ind>0 & force>0); %indexes where tip is in contact with cell
if isempty(pntindex) %if tip is never in contact, force and indentaion are considered 0
    i=0;
    f=0;
else
i=ind(pntindex); %i is indentation data only while in contact [m]
f=force(pntindex); %f is force data only while in contact [N]
end
% force and indendation have to be positive
phi=(4/3).*sqrt(R.*(i.^3)); %denominator of Hertz eqn.
E_a=(f./phi); %divide numberator by denominator of Hertz eqn.
% E_a is the effective/reduced stiffness(Pa)or apparent Young's modulus
E=zeros(size(E_a));
% initialize E
for k=1:length(pntindex)
E(k)=(1-v^2)/(1/(E_a(k))-(1-0.17^2)/74000000000); %ss - I couldn't find this equation in the paper
end
E = E/1000000; %[MPa]
i = i*1e6; %[um]
if length(i) == 1
    b  = 0;
else
min_1 = diff(i);

% find the difference between values in the X-axis
first_der_E = diff(E)./min_1;
% take the first derivative (slope) of the difference between the y-values
% divided by the difference of x-values
len = length(min_1);
% find the length of the x-values difference vector
min_1(len) = [];
%remove the last value of the vector
second_deriv_E = diff(first_der_E)./ min_1;
%take the second derivative of the y-values
absval_E = abs(first_der_E);
%find the absolute value of first derivative numbers to remove the errors
%with the signage
absval_E_2 = abs(second_deriv_E);
%find the absolute value of second derivative numbers to remove the errors
%with the signage

faf = find(absval_E < 0.1); %0.1 kPa/nm (I think)

% find the indices of the numbers where the absolute value of the first
% derivative approachs zero

faf_2 = find(absval_E_2 < 1); 

% find the indices of the numbers where the absolute value of the second
% derivative approachs zero
chec = ismember(faf_2,faf);
% find the indices where both the first and second derivative approach zero
faf_2= faf_2(chec);
% create a vector of region where 1st and 2nd derivative approach zero
% (linear region)

E_2 = E(faf);
% find the y-values that correspond to the region
i_2 = i(faf_2);
% find the x-values that correspond to the region
% use properties of fused quartz(silica) as the indenter, poisson's ration: 0.17, Etip=74.00GPa
% the poisson's ration of the sample is assumed to be 0.33.
% hold on
% plot(i,E)% plot indentation vs Young's modulus
% plot(i_2,E_2)
% xlabel('Indentation/um');
% ylabel('Elastic Modulus/MPa');
% title('Point wise elastic modulus');
b = mean(E_2);
end

if b > 0.01 
    b = 0;
end

out = b;
% find the mean of the linear region

end