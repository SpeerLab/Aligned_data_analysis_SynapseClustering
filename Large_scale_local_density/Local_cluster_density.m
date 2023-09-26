clear;
clc;

expfolder = 'X:\Chenghang\4_Color\Raw\12.21.2020_P8EA_B_V2\';
voxel = [15.5,15.5,70];
Search_radius = 12000; %nm.
Search_depth = 24; %Sections
Center_section_id = 25;

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
% disp('Done loading data. ');
%
%Exclude ranomized centered from the soma region.
Soma_folder = [Result_folder '1_soma\'];
S_filter = zeros(infos.Height,infos.Width,num_images,'uint8');
for i = 1:num_images
    S_filter(:,:,i) = imread([Soma_folder,'F_',sprintf('%03d',i),'.tif']);
end
%
Weighted_centroid_1(round(Weighted_centroid_1(:,3))<1,3) = 1;
Weighted_centroid_1(round(Weighted_centroid_1(:,3))>num_images,3) = num_images;
Weighted_centroid_1_new = [];
for i = 1:size(Weighted_centroid_1,1)
    %disp(i);
    temp = [floor(Weighted_centroid_1(i,2)),floor(Weighted_centroid_1(i,1)),floor(Weighted_centroid_1(i,3))];
    temp(temp<=0) =1;
    if S_filter(temp(1),temp(2),temp(3)) > 0
        Weighted_centroid_1_new = cat(1,Weighted_centroid_1_new,Weighted_centroid_1(i,:));
    end
end
Weighted_centroid_1 = Weighted_centroid_1_new;

Weighted_centroid_2(round(Weighted_centroid_2(:,3))<1,3) = 1;
Weighted_centroid_2(round(Weighted_centroid_2(:,3))>num_images,3) = num_images;
Weighted_centroid_2_new = [];
for i = 1:size(Weighted_centroid_2,1)
    %disp(i);
    temp = [floor(Weighted_centroid_2(i,2)),floor(Weighted_centroid_2(i,1)),floor(Weighted_centroid_2(i,3))];
    temp(temp<=0) =1;
    if S_filter(temp(1),temp(2),temp(3)) > 0
        Weighted_centroid_2_new = cat(1,Weighted_centroid_2_new,Weighted_centroid_2(i,:));
    end
end
Weighted_centroid_2 = Weighted_centroid_2_new;
%disp('Done loading data. ');
%
X_list = ceil(Search_radius/voxel(1)) + 1 : floor((infos.Width*15.5-Search_radius)/voxel(1))-1;
X_list = X_list(1:10:numel(X_list));
Y_list = ceil(Search_radius/voxel(1)) + 1 : floor((infos.Height*15.5-Search_radius)/voxel(1))-1;
Y_list = Y_list(1:10:numel(Y_list));
Local_density_1 = zeros(numel(Y_list),numel(X_list));
Local_density_2 = zeros(numel(Y_list),numel(X_list));
Local_volume = pi * Search_radius^2 * (Search_depth * 2 +1) * voxel(3);

for i = 1:numel(X_list)
    for j = 1:numel(Y_list)
        Search_radius_pixel = Search_radius/voxel(1);
        center = [X_list(i),Y_list(j),Center_section_id];
        num_temp_1 = find_local_density(Weighted_centroid_1,center,Search_radius_pixel,Search_depth);
        num_temp_2 = find_local_density(Weighted_centroid_2,center,Search_radius_pixel,Search_depth);
        den_1 = num_temp_1 / Local_volume * 1000000000;
        den_2 = num_temp_2 / Local_volume * 1000000000;
        Local_density_1(j,i) = den_1;
        Local_density_2(j,i) = den_2;
    end
end
%
r = corr2(Local_density_1,Local_density_2);
disp(r);
%%
figure;imagesc(Local_density_1);
colorbar;
figure;imagesc(Local_density_2);
colorbar;
%%
a = Local_density_1(:);
b = Local_density_2(:);
figure;scatter(a,b,'.');
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
%%
%Exclude ranomized centered from the soma region.
Soma_folder = [Result_folder '1_soma\'];
S_filter = zeros(infos.Height,infos.Width,num_images,'uint8');
for i = 1:num_images
    S_filter(:,:,i) = imread([Soma_folder,'F_',sprintf('%03d',i),'.tif']);
end
%
Weighted_centroid_1(round(Weighted_centroid_1(:,3))<1,3) = 1;
Weighted_centroid_1(round(Weighted_centroid_1(:,3))>num_images,3) = num_images;
Weighted_centroid_1_new = [];
for i = 1:size(Weighted_centroid_1,1)
    disp(i);
    temp = [round(Weighted_centroid_1(i,2)),round(Weighted_centroid_1(i,1)),round(Weighted_centroid_1(i,3))];
    if S_filter(temp(2),temp(1),temp(3)) > 0
        Weighted_centroid_1_new = cat(1,Weighted_centroid_1_new,Weighted_centroid_1(i,:));
    end
end
Weighted_centroid_1 = Weighted_centroid_1_new;

Weighted_centroid_2(round(Weighted_centroid_2(:,3))<1,3) = 1;
Weighted_centroid_2(round(Weighted_centroid_2(:,3))>num_images,3) = num_images;
Weighted_centroid_2_new = [];
for i = 1:size(Weighted_centroid_2,1)
    disp(i);
    temp = [round(Weighted_centroid_2(i,2)),round(Weighted_centroid_2(i,1)),round(Weighted_centroid_2(i,3))];
    if S_filter(temp(2),temp(1),temp(3)) > 0
        Weighted_centroid_2_new = cat(1,Weighted_centroid_2_new,Weighted_centroid_2(i,:));
    end
end
Weighted_centroid_2 = Weighted_centroid_2_new;
%%
%randomize of clusters.
statsG_temp = statsRwater_sssn;
Weighted_centroid_rand = [];
new_G = zeros(infos.Height,infos.Width,num_images,'uint8');
for i = 1:numel(statsG_temp)
    disp(i);
    while 1
        x_coor = round(rand(1) * infos.Width);
        z_coor = round(rand(1) * num_images);
        if x_coor<1||x_coor>infos.Width...
                ||z_coor<1||z_coor>num_images
            continue;
        end
        pick_list = find((new_G(:,x_coor,z_coor) == 0));
        y_coor = randsample(pick_list,1);
        PixelList = statsG_temp(i).PixelList;
        PixelList = PixelList - statsG_temp(i).WeightedCentroid...
            + [x_coor_round,y_coor_round,z_coor_round];
        for j = 1:size(PixelList,1)
            try
                x = PixelList(j,2);
                y = PixelList(j,1);
                z = PixelList(j,3);
                new_G(x,y,z) = uint8(1);
            end
        end
        Weighted_centroid_rand = cat(1,Weighted_centroid_rand,[x_coor,y_coor,z_coor]);
        break;
    end
end
%%
Weighted_centroid_2 = Weighted_centroid_rand;
%%
function [l_num] = find_local_density(Weighted_Centroid,center,radius,depth)
%radius in the unit of pixels.
l_num = 0;
for i = 1:size(Weighted_Centroid,1)
    if (Weighted_Centroid(i,3) >= (center(3) - depth)) && (Weighted_Centroid(i,3) <= (center(3) + depth))
        dist_temp = sqrt((Weighted_Centroid(i,1)-center(1))^2 + (Weighted_Centroid(i,2)-center(2))^2);
        if dist_temp <= radius
            l_num = l_num+1;
        end
    end
end
end
