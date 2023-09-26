%%
clear;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn\Comp_Simp_RDF\Comp_neg_comp_pos\';
outpath = 'X:\Chenghang\4_Color\Complex_Syn\Comp_Simp_RDF\Comp_neg_comp_pos\';

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

    load([data_path sprintf('%03d',file_ID) '.mat']);
    %
    near_dist = mean(dist_list_12_sorted(:,1));
    near_dist_rand = mean(dist_list_12_sorted_new(:,1));
    ave_near(file_ID) = near_dist;
    ave_near_rand(file_ID) = near_dist_rand;

    %num_near = numel(find(dist_list_12_sorted(:)<=4.0));
    %num_near_rand = numel(find(dist_list_12_sorted_new(:)<=4.0));
    %ave_near(file_ID) = num_near;
    %ave_near_rand(file_ID) = num_near_rand;

end
%%
disp(ave_near);
disp(ave_near_rand);