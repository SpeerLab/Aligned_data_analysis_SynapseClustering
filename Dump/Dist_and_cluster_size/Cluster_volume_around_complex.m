%%
clear all;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn\';

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
    expfolder = char(pathname(file_ID));
    disp(expfolder)
    voxel = [15.5,15.5,70];
    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    %V_Syn_folder = [Result_folder '5_V_Syn\'];
    VGluT2_folder = [Result_folder '3_Vglut2\'];

    %load([V_Syn_folder 'R_paired_3.mat']);
    load([data_path 'Vawter' sprintf('%03d',file_ID) '.mat']);
    %
    comp_centroids= [];
    simp_centroids = [];
    simp_volume = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(i).B_ID;
        B_ID = sort(B_ID);
        B_ID(1) = [];
        if numel(B_ID) > 1
            comp_centroids = cat(1,comp_centroids,statsVwater_ss(i).WeightedCentroid);
        else
            simp_centroids = cat(1,simp_centroids,statsVwater_ss(i).WeightedCentroid);
            simp_volume = cat(1,simp_volume,statsVwater_ss(i).Volume1_0);
        end
    end
    
    simp_check = zeros(size(simp_centroids,1),1);
    for i = 1:size(comp_centroids,1)
        for j = 1:size(simp_centroids,1)
            dist_temp = sqrt(((comp_centroids(i,1) - simp_centroids(j,1))*0.0155)^2 + ...
                ((comp_centroids(i,2) - simp_centroids(j,2))*0.0155)^2 + ...
                ((comp_centroids(i,3) - simp_centroids(j,3))*0.07)^2);
            if dist_temp < 1.0
                simp_check(j) = simp_check(j) + 1;
            end
        end
    end

    Volume_list = simp_volume(logical(simp_check));
    writematrix(Volume_list',[data_path 'volume_simp_near.csv'],'WriteMode','append');
end