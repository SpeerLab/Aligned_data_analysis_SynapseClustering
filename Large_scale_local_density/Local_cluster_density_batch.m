%%
clear all;clc
%
data_path = 'X:\Chenghang\4_Color_Continue\';

pathname = strings(18,1);
pathname(1) = ['X:\Chenghang\Backup_Raw_Data\1.2.2021_P2EA_B\']; %#ok<*NBRAK>
pathname(2) = ['X:\Chenghang\Backup_Raw_Data\1.4.2021_P2EB_B\'];
pathname(3) = ['X:\Chenghang\4_Color\Raw\1.6.2021_P2EC_B\'];
pathname(4) = ['X:\Chenghang\Backup_Raw_Data\7.29.2020_P4EB\'];
pathname(5) = ['X:\Chenghang\Backup_Raw_Data\9.25.2020_P4EC_B\'];
pathname(6) = ['X:\Chenghang\Backup_Raw_Data\12.5.2020_P4ED_B\'];
pathname(7) = ['X:\Chenghang\Backup_Raw_Data\12.21.2020_P8EA_B\'];
pathname(8) = ['X:\Chenghang\4_Color\Raw\12.23.2020_P8EB_B\'];
pathname(9) = ['X:\Chenghang\4_Color\Raw\1.12.2021_P8EC_B\'];
pathname(10) = ['X:\Chenghang\Backup_Raw_Data\9.29.2020_B2P2A_B\'];
pathname(11) = ['X:\Chenghang\4_Color\Raw\12.13.2020_B2P2B_B\'];
pathname(12) = ['X:\Chenghang\Backup_Raw_Data\12.18.2020_B2P2C_B\'];
pathname(13) = ['X:\Chenghang\Backup_Raw_Data\10.3.2020_B2P4A_B\'];
pathname(14) = ['X:\Chenghang\Backup_Raw_Data\10.27.2020_B2P4B_B\'];
pathname(15) = ['X:\Chenghang\Backup_Raw_Data\12.8.2020_B2P4C_B\'];
pathname(16) = ['X:\Chenghang\Backup_Raw_Data\12.12.2020_B2P8A_B\'];
pathname(17) = ['X:\Chenghang\4_Color\Raw\1.13.2021_B2P8B_B\'];
pathname(18) = ['X:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B\'];
%%
r = zeros(18,1);
linear_r = zeros(18,1);
for file_ID = 1:18
    disp(file_ID);
    expfolder = char(pathname(file_ID));
    voxel = [15.5,15.5,70];
    Search_radius = 4000; %nm.
    Search_depth = 24; %Sections
    Center_section_id = [25];

    for center_id = 1
        analysis_folder = [expfolder 'analysis\'];
        Result_folder = [analysis_folder 'Result\'];
        CTB_folder = [Result_folder '4_CTB\'];
        Img_folder = [analysis_folder 'elastic_align\storm_merged\'];

        files = [dir([Img_folder '*.tif']) dir([Img_folder '*.png'])];
        infos = imfinfo([Img_folder files(1,1).name]);
        num_images = numel(files);

        load([CTB_folder 'R_paired_VC.mat']);
        %
        %         Weighted_centroid_1 = zeros(numel(statsRwater_ssss),3);
        %         Weighted_centroid_2 = zeros(numel(statsRwater_sssn),3);
        %         for i = 1:numel(statsRwater_ssss)
        %             Weighted_centroid_1(i,:) = statsRwater_ssss(i).WeightedCentroid;
        %         end
        %         for i = 1:numel(statsRwater_sssn)
        %             Weighted_centroid_2(i,:) = statsRwater_sssn(i).WeightedCentroid;
        %         end

        Weighted_centroid_1 = zeros(numel(statsRwater_sssss),3);
        Weighted_centroid_2 = zeros(numel(statsRwater_ssssn),3);
        for i = 1:numel(statsRwater_sssss)
            temp = zeros(1,3);
            temp(1) = rand(1) * infos.Width;
            temp(2) = rand(1) * infos.Height;
            temp(3) = rand(1) * num_images;
            Weighted_centroid_1(i,:) = temp;
        end
        for i = 1:numel(statsRwater_ssssn)
            temp = zeros(1,3);
            temp(1) = rand(1) * infos.Width;
            temp(2) = rand(1) * infos.Height;
            temp(3) = rand(1) * num_images;
            Weighted_centroid_2(i,:) = temp;
        end
        %disp('Done loading data. ');
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
            temp = [round(Weighted_centroid_1(i,2)),round(Weighted_centroid_1(i,1)),round(Weighted_centroid_1(i,3))];
            temp(temp<=0) =1;
            if ~S_filter(temp(1),temp(2),temp(3)) > 0
                Weighted_centroid_1_new = cat(1,Weighted_centroid_1_new,Weighted_centroid_1(i,:));
            end
        end
        Weighted_centroid_1 = Weighted_centroid_1_new;

        Weighted_centroid_2(round(Weighted_centroid_2(:,3))<1,3) = 1;
        Weighted_centroid_2(round(Weighted_centroid_2(:,3))>num_images,3) = num_images;
        Weighted_centroid_2_new = [];
        for i = 1:size(Weighted_centroid_2,1)
            temp = [round(Weighted_centroid_2(i,2)),round(Weighted_centroid_2(i,1)),round(Weighted_centroid_2(i,3))];
            temp(temp<=0) =1;
            if ~S_filter(temp(1),temp(2),temp(3)) > 0
                Weighted_centroid_2_new = cat(1,Weighted_centroid_2_new,Weighted_centroid_2(i,:));
            end
        end
        Weighted_centroid_2 = Weighted_centroid_2_new;
        %
        %
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
                center = [X_list(i),Y_list(j),Center_section_id(center_id)];
                num_temp_1 = find_local_density(Weighted_centroid_1,center,Search_radius_pixel,Search_depth);
                num_temp_2 = find_local_density(Weighted_centroid_2,center,Search_radius_pixel,Search_depth);
                den_1 = num_temp_1 / Local_volume * 1000000000;
                den_2 = num_temp_2 / Local_volume * 1000000000;
                Local_density_1(j,i) = den_1;
                Local_density_2(j,i) = den_2;
            end
        end
        figure;subplot(1,2,1);imagesc(Local_density_1);colorbar;
        subplot(1,2,2);imagesc(Local_density_2);colorbar;
        set(gcf,'position',[0,0,1200,400]);
        saveas(gcf,[data_path sprintf('%03d',file_ID) '_heatmap.png']);
        saveas(gcf,[data_path sprintf('%03d',file_ID) '_heatmap.fig']);
        close;
        a = Local_density_1(:);
        b = Local_density_2(:);
        figure;scatter(a,b,'.');
        set(gcf,'position',[0,0,600,500]);
        saveas(gcf,[data_path sprintf('%03d',file_ID) '_linearfit.png']);
        saveas(gcf,[data_path sprintf('%03d',file_ID) '_linearfit.fig']);
        close;
        %
        R = corrcoef(cat(2,a,b));
        linear_r(file_ID) = R(2,1)^2;
        r(file_ID,center_id) = corr2(Local_density_1,Local_density_2);
    end
end
%%
% mean_r = mean(r,2);
% disp(mean_r);
mean_r = r;
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
    if S_filter(temp(1),temp(2),temp(3)) > 0
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
    if S_filter(temp(1),temp(2),temp(3)) > 0
        Weighted_centroid_2_new = cat(1,Weighted_centroid_2_new,Weighted_centroid_2(i,:));
    end
end
Weighted_centroid_2 = Weighted_centroid_2_new;
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
