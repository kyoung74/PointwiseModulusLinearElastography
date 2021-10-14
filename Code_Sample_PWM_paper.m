% open txt file
fileid = fopen('height.txt');
% extract txt
Str   = fileread('height.txt');
%arrange txt file in Line, Point, and Height string matrix
A = str2num(Str);
% Find number of colums and rows in matrix (should be 3 columns and an odd number of rows)
[r,c] = size(A)
% get number of rows needed for image matrix
count_box_rows= sqrt(256)
% value of ticker increments (size of afm measurement (in m) / number of points measured)
val_xaxis = 40e-6 ./ (count_box_rows -1)
s = r.*c;
% create a new vector of double values (to increment)
A_updated=[];
% final vector of all values
A_fin =[];
% change every string value into double value
for i = 1:c
    for k = 2:r
   x = str2num(A(k,i));
   A_updated = [A_updated; x];
    end
    A_fin = [A_fin A_updated];
    A_updated=[];
end
fclose(fileid);

%%
% Multipy the lines and points by their actual x and y value (value in
% meters) and add to vector
first_col_updated = val_xaxis .* A_fin(:,1)
% get Line indices and add to vector (increment by 1 for matlab)
ind_line = A_fin(:,1) + 1
% get point indices and add to vector (increment by 1 for matlab)
ind_point = A_fin(:,2) + 1
% add line indices to end of matrix
A_fin = [first_col_updated A_fin(:,2:end) ind_line]
% Multipy the lines and points by their actual x and y value (value in
% meters) and add to vector
second_col_updated = val_xaxis .* A_fin(:,2)
% add indices to end of matrix
A_fin = [A_fin(:,1) second_col_updated A_fin(:,3:end) ind_point]
% find length of rows and columns of new matrix of doubles
[new_r,new_c] = size(A_fin);
% create the Rx vector
Radius_x = []
% select first value of row
start_sel = 1
% select last value of row
end_sel = 16
for t = 1: count_box_rows
    radius_x_complete = []
    % %  select first line in order to measure the radius
    first_line = A_fin(start_sel:end_sel,2:end)
    % % de-serpentine the pattern
if first_line(1,4) == 1
     first_line = A_fin(start_sel:end_sel,2:end)
elseif first_line(1,4) == 16
    first_line = flip(A_fin(start_sel:end_sel,2:end))
end
% %   find radius
        for j = (2):(15)
    syms x y
     % %   if both the height and the start location are 0 for the three consecutive values automatically make radius = 0
    if first_line(j,1) == 0 & first_line(j,2) == 0 & first_line(j+1,1) == 0 & first_line(j+1,2) == 0 & first_line(j-1,1) == 0 & first_line(j-1,2) == 0
        radius_x = 0
        radius_x_complete = [radius_x_complete radius_x]
    elseif  (first_line(j,1) == 0 & first_line(j+1,1)==0 & first_line(j-1,2)==0) | (first_line(j,2) == 0 & first_line(j+1,2)==0 & first_line(j-1,2)==0)
             radius_x = 0
        radius_x_complete = [radius_x_complete radius_x]
    else
% %     find eq
    eq1 = (first_line(j,1) - x)^2 + (first_line(j,2)- y)^2;
    eq2 = (first_line(j+1,1) - x)^2 + (first_line(j+1,2)- y)^2;
    eq3 = (first_line(j-1,1) - x)^2 + (first_line(j-1,2)- y)^2;
    [solx,soly] = vpasolve(eq1 == eq2, eq1 == eq3, [x y]);
    diff = (first_line(j,1) - solx)^2 +(first_line(j,2)- soly)^2;
    radius_x = sqrt(diff)
%     place each radius value in apprpriate point in the line
    radius_x_complete = [radius_x_complete radius_x]
    end
        end
      %  add line to overall Rx matrix
    Radius_x = [Radius_x; radius_x_complete]
%     check to see if you are done evaluating all data points
    if end_sel == 256
    start_sel = 1
    end_sel = 16
    else
    start_sel = end_sel + 1
    end_sel = start_sel + 15
    end
end
% sort matrix based on point
sorted_matrix = []
for g = 1:16
%     start from point = 1
  ind = find(A_fin(:,5) == g)
%   get all rows where point = 1
  sort_m = A_fin(ind,:)
%   Create matrix sorted on points
  sorted_matrix = [sorted_matrix; sort_m]
end
% create Radius in the y direction matrix
Radius_y = []
% select first value of row
start_sel = 1
% select last value of row
end_sel = 16
for t = 1: count_box_rows
    radius_y_complete = []
    % %  select first line in order to measure the radius
    first_line = sorted_matrix(start_sel:end_sel,:)
% %   find radius
        for h = (2):(15)
    syms x y
     % %   if both the height and the start location are 0 for the three consecutive values automatically make radius = 0
    if first_line(h,1) == 0 & first_line(h,3) == 0 & first_line(h+1,1) == 0 & first_line(h+1,3) == 0 & first_line(h-1,1) == 0 & first_line(h-1,3) == 0
        radius_y = 0
        radius_y_complete = [radius_y_complete; radius_y]
    elseif (first_line(h,1) == 0 & first_line(h+1,1) == 0 & first_line(h-1,1) == 0) | (first_line(h,3) == 0 & first_line(h+1,3) == 0 & first_line(h-1,3) == 0)
        radius_y = 0
        radius_y_complete = [radius_y_complete; radius_y]
    else
% %     find eq
    eq1 = (first_line(h,1) - x)^2 + (first_line(h,3)- y)^2;
    eq2 = (first_line(h+1,1) - x)^2 + (first_line(h+1,3)- y)^2;
    eq3 = (first_line(h-1,1) - x)^2 + (first_line(h-1,3)- y)^2;
    [solx,soly] = vpasolve(eq1 == eq2, eq1 == eq3, [x y]);
    diff = (first_line(h,1) - solx)^2 +(first_line(h,3)- soly)^2;
    radius_y = sqrt(diff)
%     place each radius value in apprpriate point in the line
    radius_y_complete = [radius_y; radius_y_complete]
    end
        end
      %  add line to overall Rx matrix
    Radius_y = [Radius_y radius_y_complete]
%     check to see if you are done evaluating all data points
    if end_sel == 256
    start_sel = 1
    end_sel = 16
    else
    start_sel = end_sel + 1
    end_sel = start_sel + 15
    end
end
% fix aspect ratio of radii of x and y (i.e zoom in)
Radius_x = Radius_x(2:end-1,:)
Radius_y = Radius_y(:,2:end-1)
% get the reciprocal
One_ov_Radius_x = 1 ./ Radius_x
One_ov_Radius_y = 1 ./ Radius_y
% the sample surface curvature calculation
One_ov_Radius_s = (One_ov_Radius_x + One_ov_Radius_y) ./ 2
% get indenter radius and then take reciprocal
one_ov_indenter = ( 1 ./ 3.6600e-06)
% get final radius of sample
one_ov_R = one_ov_indenter + One_ov_Radius_s
R = double((1 ./ one_ov_R))
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
h = colorbar
colormap(pink(512));
axis image
xlabel('?m')
ylabel('?m')
ylabel(h, 'height (?m)')

%%

figure(2)
[C2,h2]= imcontour(flipdim(outData,1))
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
fileid = fopen('height.txt')
Str   = fileread('height.txt')
A = str2num(Str)
[r,c] = size(A)
s = r.*c
A_updated=[]
A_fin =[]
for i = 1:c
    for k = 2:r
   x = str2num(A(k,i))
   A_updated = [A_updated; x]
    end
    A_fin = [A_fin A_updated]
    A_updated=[]
end
fclose(fileid)

[new_r,new_c] = size(A_fin)
Mat = []
new_column =  1
for g = 1:new_r
 Line_no = A_fin(g,new_column) + 1
 Point_no = A_fin(g,new_column + 1) + 1
 height = A_fin(g, new_column + 2)
 Mat(Line_no,Point_no) = height
end
[X,Y] = meshgrid(1:size(Mat,2), 1:size(Mat));
[X2,Y2] = meshgrid(1:0.01:size(Mat,2), 1:0.01:size(Mat));
outData2 = interp2(X, Y, Mat, X2, Y2, 'linear');
imagesc(outData2);
set(gca, 'XTick', linspace(1,size(X2,2),size(X,2)));
set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
set(gca, 'XTickLabel', 1:size(X,2));
set(gca, 'YTickLabel', 1:size(X,1));
ax =gca;
ax.YDir = 'normal';
h = colorbar
colormap(pink(512));
axis image
ylabel('Line')
ylabel(h, 'point')

figure(5)
contour(outData2)
axis image
%%

figure(6)
surfl(outData2)
colormap(pink(512))
shading interp
%%
files = dir('*.txt')
filmat =[];
l = length(files);
Mat_1 = [0];
R_vec = []
length_R = length(R)
for q = length_R:-1:1
    R_v = R(q,:)'
    R_v = flip(R_v)
    R_vec = [R_vec; R_v]
end

for i=1:l
    names = files(i).name ;
    filmat = [filmat; {names}];
end
Filmat = []
for k=1:256
    A = filmat{k};
    if contains(A,'Point0000')|contains(A,'Point0015')
        continue
    elseif contains(A,'Line0015')|contains(A,'Line0000')
        continue
    else
      Filmat = [Filmat;{A}]
    end
end
  for rv = 1:length(R_vec)
    A = Filmat{rv}
    res = importdata(A);
    C = res.data;
    fun_call = @pwmreader;
     O_P = fun_call(C,R_vec(rv));
    out = str2double(regexp(A,'[\d.]+','match'));
    row = out(1);
    column = out(2);
    Mat_1(row,column) = O_P;
     end
Mask = (Mat_1 < 0 )
Mat_1(Mask) = 0
csvwrite('CAm01_new.csv',Mat_1)
[X,Y] = meshgrid(1:size(Mat_1,2), 1:size(Mat_1,1));
[X2,Y2] = meshgrid(1:0.01:size(Mat_1,2), 1:0.01:size(Mat_1,1));
outData3 = interp2(X, Y, Mat_1, X2, Y2, 'linear');
figure(7)
imagesc(outData3);
set(gca, 'XTick', linspace(1,size(X2,2),size(X,2)));
set(gca, 'YTick', linspace(1,size(X2,1),size(X,1)));
set(gca, 'XTickLabel', 1:size(X,2));
set(gca, 'YTickLabel', 1:size(X,1));
ax =gca;
ax.YDir = 'normal';
h = colorbar
colormap(jet(512));
axis image
ylabel('Line')
ylabel(h, 'Pointwise Modulus')
%%
function [out] = pwmreader(in1,in2)
v=0.5;
R = in2;
in1 = rmmissing(in1);
ind= in1(1:end,1);
force= in1(1:end,2);
% the following codes determine the point-wise Young's modulus
pntindex=find(ind>0 & force>0);
if isempty(pntindex)
    i=0
    f=0
else
i=ind(pntindex);
f=force(pntindex);
end
% force and indendation have to be positive
phi=(4/3).*sqrt(R.*(i.^3))
E_a=(f./phi);
% E_a is the effective/reduced stiffness(Pa)or apparent Young's modulus
E=zeros(size(E_a));
% initialize E
for k=1:length(pntindex)
E(k)=(1-v^2)/(1/(E_a(k))-(1-0.17^2)/74000000000);
end
E = E/1000000;
i = i*1e6;
if length(i) == 1
    b  = 0
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
faf = find(absval_E < 0.001);
% find the indices of the numbers where the absolute value of the first
% derivative approachs zero
faf_2 = find(absval_E_2 < 0.5);
% find the indices of the numbers where the absolute value of the second
% derivative approachs zero
chec = ismember(faf_2,faf);
% find the indices where both the first and second derivative approach zero
faf_2= faf_2(chec);
% create a vector of region where 1st and 2nd derivative approach zero
% (linear region)
E_2 = E(faf_2);
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
if b > 0.0016
    b = 0;
else
    b = b;
end

out = b;
% find the mean of the linear region


end
