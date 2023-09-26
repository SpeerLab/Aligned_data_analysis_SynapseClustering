%%
clear;clc
%
data_path = 'C:\Users\Chenghang\Desktop\Data\';

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
ratio_complex_list = [];
Num_paired_600 = [];
Num_paired = [];
for file_ID = 1:18
    disp(file_ID);
    base_path = char(pathname(file_ID));
    analysis_folder = [base_path 'analysis\'];
    path = [analysis_folder  'elastic_align\'];
    mergedpath = [path 'storm_merged\'];

    mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
    num_images = numel(mergedfiles);
    info = imfinfo([mergedpath mergedfiles(1,1).name]);
    %disp(num_images)
    %set voxel size in nm
    voxel=[15.5, 15.5, 70];

    Result_folder = [analysis_folder 'Result\'];
    VGluT2_folder = [Result_folder '3_Vglut2\'];
    V_Syn_folder = [Result_folder '5_V_Syn\'];
    %
    load([VGluT2_folder 'V_paired.mat']);
    load([V_Syn_folder 'R_paired_3.mat']);

    clear BP BG bg2
    %disp('allocating arrays')
    BP = zeros(info.Height, info.Width, num_images,'uint16');
    BP2 = zeros(info.Height, info.Width, num_images,'uint16');
    for i = 1:numel(statsVwater_ss)
        %disp(int2str(i))
        BP(statsVwater_ss(i).PixelIdxList)=uint16(statsVwater_ss(i).PixelValues);
    end

    %Stopped here! 1.30.2023.
    for i = 1:numel(statsRwater_ssss)
        %disp(int2str(i))
        BP2(statsRwater_ssss(i).PixelIdxList)=uint16(i);
    end
    %
    %disp('loading data')
    %
    %if not present, load in stats list

    %make categories for new variables in parfor loop and slice BP for use in
    %parfor loop
    for i=numel(statsVwater_ss):-1:1
        %disp(i)
        statsVwater_ss(i).B_ID = [];

        %try preproccessing the BP varialbe in a for loop
        minpix = min(statsVwater_ss(i).PixelList);  maxpix = max(statsVwater_ss(i).PixelList);
        min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
        max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
        if min1 < 1; min1=1; end
        if min2 < 1; min2=1; end
        if min3 < 1; min3=1; end
        if max1 > info.Width; max1=info.Width; end
        if max2 > info.Height; max2=info.Height; end
        if max3 > num_images; max3=num_images; end
        BPp(i).mat = BP(min2:max2,min1:max1,min3:max3);
        BP2p(i).mat = BP2(min2:max2,min1:max1,min3:max3);
    end
    %

    parfor jj=1:numel(statsVwater_ss)
        %
        %disp(jj)
        minpix = min(statsVwater_ss(jj).PixelList);  maxpix = max(statsVwater_ss(jj).PixelList);
        min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
        max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
        if min1 < 1; min1=1; end
        if min2 < 1; min2=1; end
        if min3 < 1; min3=1; end
        if max1 > info.Width; max1=info.Width; end
        if max2 > info.Height; max2=info.Height; end
        if max3 > num_images; max3=num_images; end

        size1 = max1-min1 + 1; size2 = max2-min2 + 1; size3 = max3-min3 + 1;
        curr2 = false(size2, size1, size3);

        for j=1: numel(statsVwater_ss(jj).PixelList(:,1))
            curr2(statsVwater_ss(jj).PixelList(j,2)-min2+1, ...
                statsVwater_ss(jj).PixelList(j,1)-min1+1, ...
                statsVwater_ss(jj).PixelList(j,3)-min3+1)= 1;
        end
        curr1a = BPp(jj).mat;
        curr1b = BP2p(jj).mat;
        Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
        %

        size2= 140;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m140 = [x,y,z];
        RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
        for kk=1:numel(PixelList_m140(:,1))
            RPixelValues_m140(kk) = curr1b(PixelList_m140(kk,2),...
                PixelList_m140(kk,1), PixelList_m140(kk,3));
        end

        statsVwater_ss(jj).B_ID = unique([RPixelValues_m140]);
    end
    %
    %B_ID should have at least 2 elements. Complex synapse are defined as
    %VGluT2 terminals that have more than 2 B)ID elements.
    num_paired_list = zeros(numel(statsVwater_ss),1);
    for i = 1:numel(statsVwater_ss)
        num_paired = numel(statsVwater_ss(i).B_ID);
        num_paired_list(i) = num_paired - 1;
    end
    %
    ratio_complex = numel(find(num_paired_list>1)) / numel(statsVwater_ss);
    ratio_complex_list = cat(1,ratio_complex_list,ratio_complex);
    %
    %Distance of Bassoon clusters that are paired under one
    Paired_dist = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = statsVwater_ss(i).B_ID;
        if numel(B_ID) > 2
            Centroids_list = [];
            for j = 1:numel(B_ID) %Delete 0 element
                if B_ID(j) ~= 0
                    temp = statsRwater_ssss(B_ID(j)).WeightedCentroid;
                    %disp(temp);
                    Centroids_list = cat(1,Centroids_list,temp);
                end
            end
            dist_all = [];
            for j = 1:size(Centroids_list,1)
                dist_list = [];
                for k = 1:size(Centroids_list,1)
                    dist = sqrt(((Centroids_list(j,1)-Centroids_list(k,1))*0.0155)^2 + ...
                        ((Centroids_list(j,2)-Centroids_list(k,2))*0.0155)^2 + ...
                        ((Centroids_list(j,3)-Centroids_list(k,3))*0.070)^2);
                    dist_list = cat(1,dist_list,dist);
                end
                dist_all = cat(2,dist_all,dist_list);
            end
            dist_all = tril(dist_all);
            Paired_dist = cat(1,Paired_dist,dist_all(:));
        end
    end
    %
    Paired_dist_nozero = [];
    for i = 1:numel(Paired_dist)
        if Paired_dist(i) > 0
            Paired_dist_nozero = cat(1,Paired_dist_nozero,Paired_dist(i));
        end
    end
    %
    Num_paired_600 = cat(1,Num_paired_600,numel(find(Paired_dist_nozero<=0.6)));
    Num_paired = cat(1,Num_paired,numel(Paired_dist_nozero));
end
%
save(['X:\Chenghang\4_Color\Complex_Syn\Complex_Syn.mat'],'Num_paired','Num_paired_600','ratio_complex_list');