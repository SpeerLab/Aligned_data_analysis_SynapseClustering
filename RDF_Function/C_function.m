clear;
clc;

expfolder = 'X:\Chenghang\4_Color\Raw\12.21.2020_P8EA_B_V2\';
voxel = [15.5,15.5,70];

analysis_folder = [expfolder 'analysis\'];
Result_folder = [analysis_folder 'Result\'];
VGluT2_folder = [Result_folder '3_Vglut2\'];
V_Syn_folder = [Result_folder '5_V_Syn\'];
Img_folder = [analysis_folder 'elastic_align\storm_merged\'];

files = [dir([Img_folder '*.tif']) dir([Img_folder '*.png'])];
infos = imfinfo([Img_folder files(1,1).name]);
num_images = numel(files);

load([V_Syn_folder 'R_paired_3.mat']);
%
Weighted_centroid_1 = zeros(numel(statsRwater_ssss),3);
Weighted_centroid_2 = zeros(numel(statsRwater_sssn),3);
for i = 1:numel(statsRwater_ssss)
    Weighted_centroid_1(i,:) = statsRwater_ssss(i).WeightedCentroid;
end
for i = 1:numel(statsRwater_sssn)
    Weighted_centroid_2(i,:) = statsRwater_sssn(i).WeightedCentroid;
end
disp('Done loading data. ');
%%
dist_list_11 = [];
for i = 1:size(Weighted_centroid_1,1)
    disp(i);
    dist = [];
    for j = 1:size(Weighted_centroid_1)
        dist_temp = sqrt(((Weighted_centroid_1(i,1) - Weighted_centroid_1(j,1))*0.0155)^2 + ...
            ((Weighted_centroid_1(i,2) - Weighted_centroid_1(j,2))*0.0155)^2 + ...
            ((Weighted_centroid_1(i,3) - Weighted_centroid_1(j,3))*0.07)^2);
        dist = cat(2,dist,dist_temp);
    end
    dist = sort(dist);
    dist_list_11 = cat(1,dist_list_11,dist);
end
dist_list_11 = tril(dist_list_11);
%
dist_list_22 = [];
for i = 1:size(Weighted_centroid_2,1)
    disp(i);
    dist = [];
    for j = 1:size(Weighted_centroid_2)
        dist_temp = sqrt(((Weighted_centroid_2(i,1) - Weighted_centroid_2(j,1))*0.0155)^2 + ...
            ((Weighted_centroid_2(i,2) - Weighted_centroid_2(j,2))*0.0155)^2 + ...
            ((Weighted_centroid_2(i,3) - Weighted_centroid_2(j,3))*0.07)^2);
        dist = cat(2,dist,dist_temp);
    end
    dist = sort(dist);
    dist_list_22 = cat(1,dist_list_22,dist);
end
dist_list_22 = tril(dist_list_22);
%
dist_list_12 = [];
for i = 1:size(Weighted_centroid_1,1)
    disp(i);
    dist = [];
    for j = 1:size(Weighted_centroid_2)
        dist_temp = sqrt(((Weighted_centroid_1(i,1) - Weighted_centroid_2(j,1))*0.0155)^2 + ...
            ((Weighted_centroid_1(i,2) - Weighted_centroid_2(j,2))*0.0155)^2 + ...
            ((Weighted_centroid_1(i,3) - Weighted_centroid_2(j,3))*0.07)^2);
        dist = cat(2,dist,dist_temp);
    end
    dist = sort(dist);
    dist_list_12 = cat(1,dist_list_12,dist);
end
%%
bin_num = 100;
lim = 2;
test = dist_list_11(:);
test(test == 0) = [];
test = test(test <= lim);
[x,y] = hist(test,bin_num);
%x = x./(4*pi*y.^2*lim/bin_num);
y2 = y + y(2) - y(1);
x = x./(4/3*pi*y2.^3 - 4/3*pi*y.^3);
x = x/x(bin_num);
%x = cat(2,0,x);
%y = cat(2,0,y);
figure;plot(y,x);
xlim([0,lim]);
hold on
plot([0,lim],[1,1],'r');

%

test = dist_list_22(:);
test(test == 0) = [];
test = test(test <= lim);
[x,y] = hist(test,bin_num);
%x = x./(4*pi*y.^2*lim/bin_num);
y2 = y + y(2) - y(1);
x = x./(4/3*pi*y2.^3 - 4/3*pi*y.^3);
x = x/x(bin_num);
%x = cat(2,0,x);
%y = cat(2,0,y);
figure;plot(y,x);
xlim([0,lim]);
hold on
plot([0,lim],[1,1],'r');

%

test = dist_list_12(:);
test(test == 0) = [];
test = test(test <= lim);
[x,y] = hist(test,bin_num);
%x = x./(4*pi*y.^2*lim/bin_num);
y2 = y + y(2) - y(1);
x = x./(4/3*pi*y2.^3 - 4/3*pi*y.^3);
x = x/x(bin_num);
%x = cat(2,0,x);
%y = cat(2,0,y);
figure;plot(y,x);
xlim([0,lim]);
hold on
plot([0,lim],[1,1],'r');
%%
%Random package: 
Weighted_centroid_1 = zeros(numel(statsRwater_ssss),3);
Weighted_centroid_2 = zeros(numel(statsRwater_sssn),3);
for i = 1:numel(statsRwater_ssss)
    temp = zeros(1,3);
    temp(1) = rand(1) * infos.Width;
    temp(2) = rand(1) * infos.Height;
    temp(3) = rand(1) * num_images;
    Weighted_centroid_1(i,:) = temp;
end
for i = 1:numel(statsRwater_sssn)
    temp = zeros(1,3);
    temp(1) = rand(1) * infos.Width;
    temp(2) = rand(1) * infos.Height;
    temp(3) = rand(1) * num_images;
    Weighted_centroid_2(i,:) = temp;
end
disp('Done loading data. ');