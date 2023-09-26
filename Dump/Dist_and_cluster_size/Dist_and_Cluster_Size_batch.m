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
volume_mean_all = [];
for file_ID = 1:18
    expfolder = char(pathname(file_ID));
    disp(expfolder)
    voxel = [15.5,15.5,70];
    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    V_Syn_folder = [Result_folder '5_V_Syn\'];
    VGluT2_folder = [Result_folder '3_Vglut2\'];

    %load([V_Syn_folder 'R_paired_3.mat']);
    load([data_path 'Vawter' sprintf('%03d',file_ID) '.mat']);
    %
    num_paired_idx = [];
    VGluT2_volume = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(i).B_ID;
        B_ID = sort(B_ID);
        B_ID(1) = [];
        if numel(B_ID) > 0
            num_paired_idx = cat(1,num_paired_idx,numel(B_ID));
            VGluT2_volume = cat(1,VGluT2_volume,statsVwater_ss(i).Volume1_0);
        end
    end
    
    Volume_mean = zeros(1,4);
    for i = 1:4
        Volume_mean(i) = mean(VGluT2_volume(num_paired_idx == i));
    end
    volume_mean_all = cat(1,volume_mean_all,Volume_mean);
end
%%
volume_y_WTP2 = volume_mean_all(1:3,:);
volume_y_WTP4 = volume_mean_all(4:6,:);
volume_y_WTP8 = volume_mean_all(7:9,:);
volume_y_B2P2 = volume_mean_all(10:12,:);
volume_y_B2P4 = volume_mean_all(13:15,:);
volume_y_B2P8 = volume_mean_all(16:18,:);
figure;
x = [1,2,3,4];
subplot(2,3,1);bar(x,mean(volume_y_WTP2));ylim([0,15000]);
subplot(2,3,2);bar(x,mean(volume_y_WTP4));ylim([0,15000]);
subplot(2,3,3);bar(x,mean(volume_y_WTP8));ylim([0,15000]);
subplot(2,3,4);bar(x,mean(volume_y_B2P2));ylim([0,15000]);
subplot(2,3,5);bar(x,mean(volume_y_B2P4));ylim([0,15000]);
subplot(2,3,6);bar(x,mean(volume_y_B2P8));ylim([0,15000]);
% subplot(2,3,1);bar(x,mean(volume_y_WTP2));ylim([-5,0.0]);
% subplot(2,3,2);bar(x,mean(volume_y_WTP4));ylim([-5,0.0]);
% subplot(2,3,3);bar(x,mean(volume_y_WTP8));ylim([-5,0.0]);
% subplot(2,3,4);bar(x,mean(volume_y_B2P2));ylim([-5,0.0]);
% subplot(2,3,5);bar(x,mean(volume_y_B2P4));ylim([-5,0.0]);
% subplot(2,3,6);bar(x,mean(volume_y_B2P8));ylim([-5,0.0]);
%%
save([data_path 'volume_mean_all.mat'],"volume_mean_all");
%%
writematrix(volume_mean_all,[data_path 'volume_mean_all.csv'],'WriteMode','append');