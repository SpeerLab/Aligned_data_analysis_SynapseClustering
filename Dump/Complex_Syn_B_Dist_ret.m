%%
clear;clc
%
output_path = 'X:\Chenghang\4_Color_Continue\';
load_path = 'X:\Chenghang\4_Color\Complex_Syn\';

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
for i = 1:18
    disp(i);
    load([load_path 'Vawter' sprintf('%03d',i) '.mat']);
    %statsVwater_ss
    base_path = char(pathname(i));
    load([base_path 'analysis\Result\5_V_Syn\R_paired_3.mat']);
    
    min_dist_list = [];
    for j = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(j).B_ID;
        if numel(B_ID) > 2
            Centroids_list = [];
            for k = 1:numel(B_ID) %Delete 0 element
                if B_ID(k) ~= 0
                    temp = statsRwater_ssss(B_ID(k)).WeightedCentroid;
                    %disp(temp);
                    Centroids_list = cat(1,Centroids_list,temp);
                end
            end
            dist_all = [];
            for k = 1:size(Centroids_list,1)
                dist_list = [];
                for m = 1:size(Centroids_list,1)
                    dist = sqrt(((Centroids_list(k,1)-Centroids_list(m,1))*0.0155)^2 + ...
                        ((Centroids_list(k,2)-Centroids_list(m,2))*0.0155)^2 + ...
                        ((Centroids_list(k,3)-Centroids_list(m,3))*0.070)^2);
                    dist_list = cat(1,dist_list,dist);
                end
                dist_all = cat(2,dist_all,dist_list);
            end
            dist_all = dist_all(:);
            dist_all= dist_all(find(dist_all));
            min_dist = min(dist_all);
            min_dist_list = cat(1,min_dist_list,min_dist);
        end
    end
    writematrix(min_dist_list',[output_path 'mind_dist_list.csv'],'WriteMode','append');
end