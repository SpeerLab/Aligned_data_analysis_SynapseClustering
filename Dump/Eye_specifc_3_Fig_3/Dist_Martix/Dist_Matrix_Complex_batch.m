%%
clear;clc
data_path = 'X:\Chenghang\4_Color\Complex_Syn_2\';
out_path = 'X:\Chenghang\4_Color\Complex_Syn_3\';
Centroid_path = 'X:\Chenghang\4_Color\Complex_Syn_3\Weighted_centroid\Comp_pos_simp_neg\';

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

for file_ID = 1:18
    disp(file_ID);
    expfolder = char(pathname(file_ID));
    voxel = [15.5,15.5,70];

    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    Soma_folder = [Result_folder '1_soma\'];
    CTB_folder = [Result_folder '4_CTB\'];
    Img_folder = [analysis_folder 'elastic_align\storm_merged\'];

    files = [dir([Img_folder '*.tif']) dir([Img_folder '*.png'])];
    infos = imfinfo([Img_folder files(1,1).name]);
    num_images = numel(files);

    load([data_path 'Vwater_ss_' sprintf('%03d',file_ID) '.mat']);
    load([data_path 'Vwater_sn_' sprintf('%03d',file_ID) '.mat']);
    %statsVwater_ss
%     soma_region = zeros(infos.Height,infos.Width,num_images);
%     for i = 1:num_images
%         soma_region(:,:,i) = imread([Soma_folder 'F_' sprintf('%03d',i) '.tif']);
%     end
%     soma_region = ~logical(soma_region);

    Centroid_simp_ss = [];
    Centroid_complex_ss = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = sort(statsVwater_ss(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            Centroid_complex_ss = cat(1,Centroid_complex_ss,statsVwater_ss(i).WeightedCentroid);
        else
            Centroid_simp_ss = cat(1,Centroid_simp_ss,statsVwater_ss(i).WeightedCentroid);
        end
    end

    Centroid_simp_sn = [];
    Centroid_complex_sn = [];
    for i = 1:numel(statsVwater_sn)
        B_ID = sort(statsVwater_sn(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            Centroid_complex_sn = cat(1,Centroid_complex_sn,statsVwater_sn(i).WeightedCentroid);
        else
            Centroid_simp_sn = cat(1,Centroid_simp_sn,statsVwater_sn(i).WeightedCentroid);
        end
    end
    
    Weighted_centroid_1 = Centroid_complex_ss;
    Weighted_centroid_2 = Centroid_simp_sn;

%     Weighted_centroid_1_new = zeros(size(Weighted_centroid_1,1),3);
%     Weighted_centroid_2_new = zeros(size(Weighted_centroid_2,1),3);
%     for i = 1:size(Weighted_centroid_1,1)
%         temp = zeros(1,3);
%         z = randi(size(soma_region,3));
%         temp(3) = z;
%         [temp(2),temp(1)] = ind2sub(size(soma_region(:,:,z)),randsample(find(soma_region(:,:,z)),1));
%         Weighted_centroid_1_new(i,:) = temp;
%     end
%     for i = 1:size(Weighted_centroid_2,1)
%         temp = zeros(1,3);
%         z = randi(size(soma_region,3));
%         [temp(2),temp(1),temp(3)] = ind2sub(size(soma_region(:,:,z)),randsample(find(soma_region(:,:,z)),1));
%         Weighted_centroid_2_new(i,:) = temp;
%     end
    load([Centroid_path sprintf('%03d',file_ID) '.mat']);

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

    save([out_path sprintf('%03d',file_ID) '.mat'],'dist_list_12','dist_list_12_new');
end
%%