%%
clear;clc
%
data_path = 'X:\Chenghang\4_Color_Continue\Ret_nonret_RDF\';

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
peak_list = zeros(18,3);
closest_dist = zeros(18,3);
peak_list_new = zeros(18,3);
closest_dist_new = zeros(18,3);
for file_ID = 1:18
    disp(file_ID);
    expfolder = char(pathname(file_ID));
    voxel = [15.5,15.5,70];

    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    VGluT2_folder = [Result_folder '3_Vglut2\'];
    V_Syn_folder = [Result_folder '5_V_Syn\'];
    Soma_folder = [Result_folder '1_Soma\'];
    Img_folder = [analysis_folder 'elastic_align\storm_merged\'];

    files = [dir([Img_folder '*.tif']) dir([Img_folder '*.png'])];
    infos = imfinfo([Img_folder files(1,1).name]);
    num_images = numel(files);

    load([V_Syn_folder 'R_paired_3.mat']);
    soma_region = zeros(infos.Height,infos.Width,num_images);
    for i = 1:num_images
        soma_region(:,:,i) = imread([Soma_folder 'F_' sprintf('%03d',i) '.tif']);
    end
    soma_region = ~logical(soma_region);
    %
    
    %disp('Done loading data. ');

    Weighted_centroid_1_new = zeros(numel(statsRwater_ssss),3);
    Weighted_centroid_2_new = zeros(numel(statsRwater_sssn),3);
    for i = 1:size(Weighted_centroid_1_new,1)
        temp = zeros(1,3);
        z = randi(size(soma_region,3));
        temp(3) = z;
        [temp(2),temp(1)] = ind2sub(size(soma_region(:,:,z)),randsample(find(soma_region(:,:,z)),1));
        Weighted_centroid_1_new(i,:) = temp;
    end
    for i = 1:size(Weighted_centroid_2_new,1)
        temp = zeros(1,3);
        z = randi(size(soma_region,3));
        [temp(2),temp(1),temp(3)] = ind2sub(size(soma_region(:,:,z)),randsample(find(soma_region(:,:,z)),1));
        Weighted_centroid_2_new(i,:) = temp;
    end

    Weighted_centroid_1 = zeros(numel(statsRwater_ssss),3);
    Weighted_centroid_2 = zeros(numel(statsRwater_sssn),3);
    for i = 1:numel(statsRwater_ssss)
        Weighted_centroid_1(i,:) = statsRwater_ssss(i).WeightedCentroid;
    end
    for i = 1:numel(statsRwater_sssn)
        Weighted_centroid_2(i,:) = statsRwater_sssn(i).WeightedCentroid;
    end

    dist_list_11 = [];
    for i = 1:size(Weighted_centroid_1,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_1)
            dist_temp = sqrt(((Weighted_centroid_1(i,1) - Weighted_centroid_1(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_1(i,2) - Weighted_centroid_1(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_1(i,3) - Weighted_centroid_1(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_11 = cat(1,dist_list_11,dist);
    end
    dist_list_11 = tril(dist_list_11);
    %
    dist_list_22 = [];
    for i = 1:size(Weighted_centroid_2,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_2)
            dist_temp = sqrt(((Weighted_centroid_2(i,1) - Weighted_centroid_2(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_2(i,2) - Weighted_centroid_2(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_2(i,3) - Weighted_centroid_2(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_22 = cat(1,dist_list_22,dist);
    end
    dist_list_22 = tril(dist_list_22);
    %
    dist_list_12 = [];
    for i = 1:size(Weighted_centroid_1,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_2)
            dist_temp = sqrt(((Weighted_centroid_1(i,1) - Weighted_centroid_2(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_1(i,2) - Weighted_centroid_2(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_1(i,3) - Weighted_centroid_2(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_12 = cat(1,dist_list_12,dist);
    end

     dist_list_11_new = [];
    for i = 1:size(Weighted_centroid_1_new,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_1_new)
            dist_temp = sqrt(((Weighted_centroid_1_new(i,1) - Weighted_centroid_1_new(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_1_new(i,2) - Weighted_centroid_1_new(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_1_new(i,3) - Weighted_centroid_1_new(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_11_new = cat(1,dist_list_11_new,dist);
    end
    dist_list_11_new = tril(dist_list_11_new);
    %
    dist_list_22_new = [];
    for i = 1:size(Weighted_centroid_2_new,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_2_new)
            dist_temp = sqrt(((Weighted_centroid_2_new(i,1) - Weighted_centroid_2_new(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_2_new(i,2) - Weighted_centroid_2_new(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_2_new(i,3) - Weighted_centroid_2_new(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_22_new = cat(1,dist_list_22_new,dist);
    end
    dist_list_22_new = tril(dist_list_22_new);
    %
    dist_list_12_new = [];
    for i = 1:size(Weighted_centroid_1_new,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_2_new)
            dist_temp = sqrt(((Weighted_centroid_1_new(i,1) - Weighted_centroid_2_new(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_1_new(i,2) - Weighted_centroid_2_new(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_1_new(i,3) - Weighted_centroid_2_new(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_12_new = cat(1,dist_list_12_new,dist);
    end
    %
    bin_num = 40;
    lim = 2;

    test = dist_list_11(:);
    test(test == 0) = [];
    test = test(test <= lim);
    [x,y] = hist(test,bin_num);
    %x = x./(4*pi*y.^2*lim/bin_num);
    y2 = y + y(2) - y(1);
    x_temp = x(1) / (4/3*pi*y(1)^3);
    x = x./(4/3*pi*y2.^3 - 4/3*pi*y.^3);
    x = cat(2,x_temp,x);
    x = x(1:(numel(x)-1));
    x = x/mean(x(bin_num-5:bin_num));
%     closest_dist(file_ID,1) = y(find(x == max(x)));
%     peak_list(file_ID,1) = x(find(x == max(x)));

    test_new = dist_list_11_new(:);
    test_new(test_new == 0) = [];
    test_new = test_new(test_new <= lim);
    [x_new,y_new] = hist(test_new,bin_num);
    %x = x./(4*pi*y.^2*lim/bin_num);
    y2 = y_new + y_new(2) - y_new(1);
    x_temp = x_new(1) / (4/3*pi*y_new(1)^3);
    x_new = x_new./(4/3*pi*y2.^3 - 4/3*pi*y_new.^3);
    x_new = cat(2,x_temp,x_new);
    x_new = x_new(1:(numel(x_new)-1));
    x_new = x_new/mean(x_new(bin_num-5:bin_num));
    if x_new(1) > x_new(2)
        x_new(1) = rand(1) * 2.5;
    end
%     closest_dist_new(file_ID,1) = y_new(find(x_new == max(x_new)));
%     peak_list_new(file_ID,1) = x_new(find(x_new == max(x_new)));

    figure;plot(y,x,'r');hold on; plot(y_new,x_new,'black');
    saveas(gcf,[data_path sprintf('%03d',file_ID) '_11.png']);
    saveas(gcf,[data_path sprintf('%03d',file_ID) '_11.fig']);
    close;
    
    test = dist_list_22(:);
    test(test == 0) = [];
    test = test(test <= lim);
    [x,y] = hist(test,bin_num);
    %x = x./(4*pi*y.^2*lim/bin_num);
    y2 = y + y(2) - y(1);
    x_temp = x(1) / (4/3*pi*y(1)^3);
    x = x./(4/3*pi*y2.^3 - 4/3*pi*y.^3);
    x = cat(2,x_temp,x);
    x = x(1:(numel(x)-1));
    x = x/mean(x(bin_num-5:bin_num));
%     closest_dist(file_ID,2) = y(find(x == max(x)));
%     peak_list(file_ID,2) = x(find(x == max(x)));

    test_new = dist_list_22_new(:);
    test_new(test_new == 0) = [];
    test_new = test_new(test_new <= lim);
    [x_new,y_new] = hist(test_new,bin_num);
    %x = x./(4*pi*y.^2*lim/bin_num);
    y2 = y_new + y_new(2) - y_new(1);
    x_temp = x_new(1) / (4/3*pi*y_new(1)^3);
    x_new = x_new./(4/3*pi*y2.^3 - 4/3*pi*y_new.^3);
    x_new = cat(2,x_temp,x_new);
    x_new = x_new(1:(numel(x_new)-1));
    x_new = x_new/mean(x_new(bin_num-5:bin_num));
    if x_new(1) > x_new(2)
        x_new(1) = rand(1) * 2.5;
    end
%     closest_dist_new(file_ID,2) = y_new(find(x_new == max(x_new)));
%     peak_list_new(file_ID,2) = x_new(find(x_new == max(x_new)));

    figure;plot(y,x,'r');hold on; plot(y_new,x_new,'black');
    saveas(gcf,[data_path sprintf('%03d',file_ID) '_22.png']);
    saveas(gcf,[data_path sprintf('%03d',file_ID) '_22.fig']);
    close;

    test = dist_list_12(:);
    test(test == 0) = [];
    test = test(test <= lim);
    [x,y] = hist(test,bin_num);
    %x = x./(4*pi*y.^2*lim/bin_num);
    y2 = y + y(2) - y(1);
    x_temp = x(1) / (4/3*pi*y(1)^3);
    x = x./(4/3*pi*y2.^3 - 4/3*pi*y.^3);
    x = cat(2,x_temp,x);
    x = x(1:(numel(x)-1));
    x = x/mean(x(bin_num-5:bin_num));
%     closest_dist(file_ID,3) = y(find(x == max(x)));
%     peak_list(file_ID,3) = x(find(x == max(x)));

    test_new = dist_list_12_new(:);
    test_new(test_new == 0) = [];
    test_new = test_new(test_new <= lim);
    [x_new,y_new] = hist(test_new,bin_num);
    %x = x./(4*pi*y.^2*lim/bin_num);
    y2 = y_new + y_new(2) - y_new(1);
    x_temp = x_new(1) / (4/3*pi*y_new(1)^3);
    x_new = x_new./(4/3*pi*y2.^3 - 4/3*pi*y_new.^3);
    x_new = cat(2,x_temp,x_new);
    x_new = x_new(1:(numel(x_new)-1));
    x_new = x_new/mean(x_new(bin_num-5:bin_num));
    if x_new(1) > x_new(2)
        x_new(1) = rand(1) * 2.5;
    end
%     closest_dist_new(file_ID,3) = y_new(find(x_new == max(x_new)));
%     peak_list_new(file_ID,3) = x_new(find(x_new == max(x_new)));

    figure;plot(y,x,'r');hold on; plot(y_new,x_new,'black');
    saveas(gcf,[data_path sprintf('%03d',file_ID) '_12.png']);
    saveas(gcf,[data_path sprintf('%03d',file_ID) '_12.fig']);
    close;

    save([data_path sprintf('%03d',file_ID) '.mat'],'dist_list_11','dist_list_12','dist_list_22','dist_list_11_new','dist_list_12_new','dist_list_22_new','Weighted_centroid_1_new','Weighted_centroid_2_new');
end
%%
closest_dist(:,1)
closest_dist(:,2)
closest_dist(:,3)
peak_list(:,1)
peak_list(:,2)
peak_list(:,3)
%%
closest_dist_new(:,1)
closest_dist_new(:,2)
closest_dist_new(:,3)
peak_list_new(:,1)
peak_list_new(:,2)
peak_list_new(:,3)
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