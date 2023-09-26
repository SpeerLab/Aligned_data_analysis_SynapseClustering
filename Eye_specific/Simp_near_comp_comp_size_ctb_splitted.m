%%
clear;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn\CTB_Specific\';
outpath = 'X:\Chenghang\4_Color\Complex_Syn\Simp_near_comp_size_ctb_splitted\';

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
    CTB_folder = [Result_folder '4_CTB\'];
    Img_folder = [analysis_folder 'elastic_align\storm_merged\'];

    files = [dir([Img_folder '*.tif']) dir([Img_folder '*.png'])];
    infos = imfinfo([Img_folder files(1,1).name]);
    num_images = numel(files);

    load([data_path 'Vwater_ss_' sprintf('%03d',file_ID) '.mat']);
    load([data_path 'Vwater_sn_' sprintf('%03d',file_ID) '.mat']);
    %statsVwater_ss
    comp_ss = [];
    simp_ss = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(i).B_ID;
        B_ID = sort(B_ID);
        B_ID(1) = [];
        if numel(B_ID) > 1
            for j = 1:numel(B_ID)
                comp_ss = cat(1,comp_ss,statsVwater_ss(i));
            end
        else 
            if numel(B_ID) == 1
                simp_ss = cat(1,simp_ss,statsVwater_ss(i));
            end
        end
    end
    comp_sn = [];
    simp_sn = [];
    for i = 1:numel(statsVwater_sn)
        B_ID = statsVwater_sn(i).B_ID;
        B_ID = sort(B_ID);
        B_ID(1) = [];
        if numel(B_ID) > 1
            for j = 1:numel(B_ID)
                comp_sn = cat(1,comp_sn,statsVwater_sn(i));
            end
        else 
            if numel(B_ID) == 1
                simp_sn = cat(1,simp_sn,statsVwater_sn(i));
            end
        end
    end
    
    clear statsRwater_ssss statsRwater_sssn
    statsRwater_ssss = simp_ss;
    statsRwater_sssn = comp_ss;


    Weighted_centroid_1 = zeros(numel(statsRwater_ssss),3);
    Weighted_centroid_2 = zeros(numel(statsRwater_sssn),3);
    for i = 1:numel(statsRwater_ssss)
        Weighted_centroid_1(i,:) = statsRwater_ssss(i).WeightedCentroid;
    end
    for i = 1:numel(statsRwater_sssn)
        Weighted_centroid_2(i,:) = statsRwater_sssn(i).WeightedCentroid;
    end
    %test
    thre = 2;
   
    sel_id = zeros(size(Weighted_centroid_1,1),1);
    for i = 1:size(Weighted_centroid_1,1)
        %disp(i);
        dist = [];
        for j = 1:size(Weighted_centroid_2)
            dist_temp = sqrt(((Weighted_centroid_1(i,1) - Weighted_centroid_2(j,1))*0.0155)^2 + ...
                ((Weighted_centroid_1(i,2) - Weighted_centroid_2(j,2))*0.0155)^2 + ...
                ((Weighted_centroid_1(i,3) - Weighted_centroid_2(j,3))*0.07)^2);
            dist = cat(2,dist,dist_temp);
        end
        if numel(find(dist<thre)) >0
            sel_id(i) = 1;
        end
    end

    size_list = [statsRwater_ssss.Volume1_0];
    sel = size_list(find(sel_id));
    nonsel = size_list(find(~logical(sel_id)));

    writematrix(sel,[outpath 'CTBPos_near.csv'],'WriteMode','append');
    writematrix(nonsel,[outpath 'CTBPos_far.csv'],'WriteMode','append');
end
%%