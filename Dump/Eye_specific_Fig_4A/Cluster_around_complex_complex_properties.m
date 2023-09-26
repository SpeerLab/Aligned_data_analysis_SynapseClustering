%%
clear;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn_1\Comp_Simp_RDF\Comp_neg_comp_pos\';
V_path = 'X:\Chenghang\4_Color\Complex_Syn_2\';
outpath = 'C:\Users\Chenghang\Desktop\My_m_Codes\4_color_continue\Eye_specific_Fig_4A\';

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

    load([V_path 'Vwater_sn_' sprintf('%03d',file_ID) '.mat']);
    Volume_comp_sn = [];
    for i = 1:numel(statsVwater_sn)
        B_ID = sort(statsVwater_sn(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            Volume_comp_sn = cat(1,Volume_comp_sn,numel(B_ID));
        else
            %Centroid_simp_sn = cat(1,Centroid_simp_sn,statsVwater_sn(i).WeightedCentroid);
        end
    end

    load([data_path sprintf('%03d',file_ID) '.mat']);
    %
    
    near_dist = dist_list_12_sorted(:,1);
    writematrix(near_dist',[outpath 'Near_dist.csv'],'WriteMode','append');
    writematrix(Volume_comp_sn',[outpath 'B_ID.csv'],'WriteMode','append');
%     near_dist_rand = dist_list_12_sorted_new(:,1);
%     figure;
%     boxplot(near_dist,Volume_comp_sn);
%     figure;
%     scatter(near_dist_rand,Volume_comp_sn,'.');

    %ave_near(file_ID) = near_dist;
    %ave_near_rand(file_ID) = near_dist_rand;

end
%%
disp(ave_near);
disp(ave_near_rand);