%%
clear;clc
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
num_comp_B_list = zeros(18,1);
num_simp_B_list = zeros(18,1);
num_comp_V_list = zeros(18,1);
num_simp_V_list = zeros(18,1);
for file_ID = 1:18
    disp(file_ID);
    
    expfolder = char(pathname(file_ID));
    analysis_folder = [expfolder 'analysis\'];
    Result_folder = [analysis_folder 'Result\'];
    V_Syn_folder = [Result_folder '5_V_Syn\'];
    load([V_Syn_folder 'R_paired_3.mat']);
    load([data_path 'Vawter' sprintf('%03d',file_ID) '.mat']);
    %Complex or simple Bassoon clusters
    comp = [];
    simp = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(i).B_ID;
        B_ID = sort(B_ID);
        B_ID(1) = [];
        if numel(B_ID) > 1
            for j = 1:numel(B_ID)
                comp = cat(1,comp,statsRwater_ssss(B_ID(j)));
            end
        else 
            if numel(B_ID) == 1
                simp = cat(1,simp,statsRwater_ssss(B_ID(1)));
            end
        end
    end
    
    %Complex or simple VGluT2 clusters
    comp_V = [];
    simp_V = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(i).B_ID;
        B_ID = sort(B_ID);
        B_ID(1) = [];
        if numel(B_ID) > 1
                comp_V = cat(1,comp_V,statsVwater_ss(i));
        else 
            if numel(B_ID) == 1
                simp_V = cat(1,simp_V,statsVwater_ss(i));
            end
        end
    end

    num_comp_B_list(file_ID) = numel(comp);
    num_simp_B_list(file_ID) = numel(simp);
    num_simp_V_list(file_ID) = numel(simp_V);
    num_comp_V_list(file_ID) = numel(comp) / numel(comp_V);
end
%save([data_path 'complex_denstiy.mat'],'num_comp_B_list','num_simp_B_list','num_simp_V_list','num_comp_V_list')