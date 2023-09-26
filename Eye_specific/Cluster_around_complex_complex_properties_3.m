%%
%Similar as property 2 but for same eye comp-comp distance vs. simple
%synaps
clear;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn\Comp_Simp_RDF\Comp_neg_comp_pos\';
V_path = 'X:\Chenghang\4_Color\Complex_Syn\CTB_Specific\';
out_path = 'X:\Chenghang\4_Color\Complex_Syn\Comp_Simp_RDF\Comp_neg_comp_pos_property_2\';

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
ave_near = zeros(18,1);
ave_near_rand = zeros(18,1);

for file_ID = 1:18
    expfolder = char(pathname(file_ID));
    disp(expfolder)
    voxel = [15.5,15.5,70];
    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    
    load([V_path 'Vwater_ss_' sprintf('%03d',file_ID) '.mat']);
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
    
    Weighted_centroid_1 = Centroid_complex_ss;
    Weighted_centroid_2 = Centroid_complex_ss;

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

    dist_list_12_sorted = sort(dist_list_12,2);
    dist_list_12_sorted = dist_list_12_sorted(:,1:min(size(dist_list_12_sorted,2),40));
    dist_list_12_sorted_mean = mean(dist_list_12_sorted);
    dist_list_12_sorted_se = std(dist_list_12_sorted) / sqrt(size(dist_list_12_sorted,1));

    dist_list_12_simp = dist_list_12_sorted;
    num_near = zeros(size(dist_list_12_sorted,1),1);
    for i = 1:size(dist_list_12_sorted,1)
        num_near(i) = numel(find(dist_list_12_sorted(i,:) < 1.0));
    end

    load([data_path sprintf('%03d',file_ID) '.mat']);
    %
    near_dist = dist_list_12_sorted(:,1);
    near_dist_rand = dist_list_12_sorted_new(:,1);

    writematrix(num_near',[out_path 'num_near.csv'],'WriteMode','append');
    writematrix(near_dist',[out_path 'near_dist.csv'],'WriteMode','append');
    writematrix(near_dist_rand',[out_path 'near_dist_rand.csv'],'WriteMode','append');

    %figure;
    %boxplot(near_dist,num_near);
    %figure;
    %boxplot(near_dist_rand,num_near);


    %ave_near(file_ID) = near_dist;
    %ave_near_rand(file_ID) = near_dist_rand;
end
%%
disp(ave_near);
disp(ave_near_rand);