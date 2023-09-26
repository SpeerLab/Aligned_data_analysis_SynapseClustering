%%
clear;clc
data_path = 'X:\Chenghang\4_Color\Complex_Syn\';
out_path = 'X:\Chenghang\4_Color\Complex_Syn\Comp_C_Func\';

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
closest_dist = zeros(18,3);
peak_list = zeros(18,3);
closest_dist_new = zeros(18,3);
peak_list_new = zeros(18,3);
for file_ID = 1:18
    disp(file_ID);
    expfolder = char(pathname(file_ID));
    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    VGluT2_folder = [Result_folder '3_Vglut2\'];
    V_Syn_folder = [Result_folder '5_V_Syn\'];
    Img_folder = [analysis_folder 'elastic_align\storm_merged\'];

    files = [dir([Img_folder '*.tif']) dir([Img_folder '*.png'])];
    infos = imfinfo([Img_folder files(1,1).name]);
    num_images = numel(files);

    load([V_Syn_folder 'R_paired_3.mat']);
    voxel = [15.5,15.5,70];
    load([data_path 'Vawter' sprintf('%03d',file_ID) '.mat']); %Load statsVwater_ss with B_ID.

    Centroid_normal = [];
    Centroid_complex = [];
    No_Matching_V_count = 0;
    for i = 1:numel(statsVwater_ss)
        B_ID = sort(statsVwater_ss(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
%             for j = 1:numel(B_ID)
%                 Centroid_complex = cat(1,Centroid_complex,statsRwater_ssss(B_ID(j)).WeightedCentroid);
%             end
            Centroid_complex = cat(1,Centroid_complex,statsRwater_ssss(B_ID(1)).WeightedCentroid);
        else
            try
                Centroid_normal = cat(1,Centroid_normal,statsRwater_ssss(B_ID(1)).WeightedCentroid);
            catch
                No_Matching_V_count = No_Matching_V_count+1;
            end
        end
    end
    disp(["The number of VGluT2 clusters far from Bassoon " num2str(No_Matching_V_count)]);
    Centroid_normal_rand = zeros(numel(statsRwater_ssss),3);
    Centroid_complex_rand = zeros(numel(statsRwater_sssn),3);
    for i = 1:numel(statsRwater_ssss)
        temp = zeros(1,3);
        temp(1) = rand(1) * infos.Width;
        temp(2) = rand(1) * infos.Height;
        temp(3) = rand(1) * num_images;
        Centroid_normal_rand(i,:) = temp;
    end
    for i = 1:numel(statsRwater_sssn)
        temp = zeros(1,3);
        temp(1) = rand(1) * infos.Width;
        temp(2) = rand(1) * infos.Height;
        temp(3) = rand(1) * num_images;
        Centroid_complex_rand(i,:) = temp;
    end

    dist_list_11 = [];
    for i = 1:size(Centroid_normal,1)
        %disp(i);
        dist = [];
        for j = 1:size(Centroid_normal)
            dist_temp = sqrt(((Centroid_normal(i,1) - Centroid_normal(j,1))*0.0155)^2 + ...
                ((Centroid_normal(i,2) - Centroid_normal(j,2))*0.0155)^2 + ...
                ((Centroid_normal(i,3) - Centroid_normal(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_11 = cat(1,dist_list_11,dist);
    end
    dist_list_11 = tril(dist_list_11);
    %
    dist_list_22 = [];
    for i = 1:size(Centroid_complex,1)
        %disp(i);
        dist = [];
        for j = 1:size(Centroid_complex)
            dist_temp = sqrt(((Centroid_complex(i,1) - Centroid_complex(j,1))*0.0155)^2 + ...
                ((Centroid_complex(i,2) - Centroid_complex(j,2))*0.0155)^2 + ...
                ((Centroid_complex(i,3) - Centroid_complex(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_22 = cat(1,dist_list_22,dist);
    end
    dist_list_22 = tril(dist_list_22);
    %
    dist_list_12 = [];
    for i = 1:size(Centroid_normal,1)
        %disp(i);
        dist = [];
        for j = 1:size(Centroid_complex)
            dist_temp = sqrt(((Centroid_normal(i,1) - Centroid_complex(j,1))*0.0155)^2 + ...
                ((Centroid_normal(i,2) - Centroid_complex(j,2))*0.0155)^2 + ...
                ((Centroid_normal(i,3) - Centroid_complex(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_12 = cat(1,dist_list_12,dist);
    end

    dist_list_11_rand = [];
    for i = 1:size(Centroid_normal_rand,1)
        %disp(i);
        dist = [];
        for j = 1:size(Centroid_normal_rand)
            dist_temp = sqrt(((Centroid_normal_rand(i,1) - Centroid_normal_rand(j,1))*0.0155)^2 + ...
                ((Centroid_normal_rand(i,2) - Centroid_normal_rand(j,2))*0.0155)^2 + ...
                ((Centroid_normal_rand(i,3) - Centroid_normal_rand(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_11_rand = cat(1,dist_list_11_rand,dist);
    end
    dist_list_11_rand = tril(dist_list_11_rand);
    %
    dist_list_22_rand = [];
    for i = 1:size(Centroid_complex_rand,1)
        %disp(i);
        dist = [];
        for j = 1:size(Centroid_complex_rand)
            dist_temp = sqrt(((Centroid_complex_rand(i,1) - Centroid_complex_rand(j,1))*0.0155)^2 + ...
                ((Centroid_complex_rand(i,2) - Centroid_complex_rand(j,2))*0.0155)^2 + ...
                ((Centroid_complex_rand(i,3) - Centroid_complex_rand(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_22_rand = cat(1,dist_list_22_rand,dist);
    end
    dist_list_22_rand = tril(dist_list_22_rand);
    %
    dist_list_12_rand = [];
    for i = 1:size(Centroid_normal_rand,1)
        %disp(i);
        dist = [];
        for j = 1:size(Centroid_complex_rand)
            dist_temp = sqrt(((Centroid_normal_rand(i,1) - Centroid_complex_rand(j,1))*0.0155)^2 + ...
                ((Centroid_normal_rand(i,2) - Centroid_complex_rand(j,2))*0.0155)^2 + ...
                ((Centroid_normal_rand(i,3) - Centroid_complex_rand(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        dist_list_12_rand = cat(1,dist_list_12_rand,dist);
    end
    %
    bin_num = 20;
    lim = 1.5;

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
    closest_dist(file_ID,1) = y(find(x == max(x)));
    peak_list(file_ID,1) = x(find(x == max(x)));

    test_new = dist_list_11_rand(:);
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
    closest_dist_new(file_ID,1) = y_new(find(x_new == max(x_new)));
    peak_list_new(file_ID,1) = x_new(find(x_new == max(x_new)));

    figure;plot(y,x,'r');hold on; plot(y_new,x_new,'black');
    saveas(gcf,[out_path sprintf('%03d',file_ID) '_11.png']);
    saveas(gcf,[out_path sprintf('%03d',file_ID) '_11.fig']);
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
    closest_dist(file_ID,2) = y(find(x == max(x)));
    peak_list(file_ID,2) = x(find(x == max(x)));

    test_new = dist_list_22_rand(:);
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
    closest_dist_new(file_ID,2) = y_new(find(x_new == max(x_new)));
    peak_list_new(file_ID,2) = x_new(find(x_new == max(x_new)));

    figure;plot(y,x,'r');hold on; plot(y_new,x_new,'black');
    saveas(gcf,[out_path sprintf('%03d',file_ID) '_22.png']);
    saveas(gcf,[out_path sprintf('%03d',file_ID) '_22.fig']);
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
    closest_dist(file_ID,3) = y(find(x == max(x)));
    peak_list(file_ID,3) = x(find(x == max(x)));

    test_new = dist_list_12_rand(:);
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
    closest_dist_new(file_ID,3) = y_new(find(x_new == max(x_new)));
    peak_list_new(file_ID,3) = x_new(find(x_new == max(x_new)));

    figure;plot(y,x,'r');hold on; plot(y_new,x_new,'black');
    saveas(gcf,[out_path sprintf('%03d',file_ID) '_12.png']);
    saveas(gcf,[out_path sprintf('%03d',file_ID) '_12.fig']);
    close;
    %
end
save([out_path 'Complex_Syn_batch.mat'],"closest_dist","peak_list","closest_dist_new","peak_list_new");
%%
closest_dist(:,1)
closest_dist(:,2)
closest_dist(:,3)
peak_list(:,1)
peak_list(:,2)
peak_list(:,3)
%%
figure;bar(peak_list(:,1));
figure;bar(peak_list(:,2));
figure;bar(peak_list(:,3));
%%
figure;plot(y,x)
%%
%Random package:
Centroid_normal = zeros(numel(statsRwater_ssss),3);
Centroid_complex = zeros(numel(statsRwater_sssn),3);
for i = 1:numel(statsRwater_ssss)
    temp = zeros(1,3);
    temp(1) = rand(1) * infos.Width;
    temp(2) = rand(1) * infos.Height;
    temp(3) = rand(1) * num_images;
    Centroid_normal(i,:) = temp;
end
for i = 1:numel(statsRwater_sssn)
    temp = zeros(1,3);
    temp(1) = rand(1) * infos.Width;
    temp(2) = rand(1) * infos.Height;
    temp(3) = rand(1) * num_images;
    Centroid_complex(i,:) = temp;
end
disp('Done loading data. ');