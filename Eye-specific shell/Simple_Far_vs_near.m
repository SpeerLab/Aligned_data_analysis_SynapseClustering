%%
clear;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn_2\';
out_path = 'X:\Chenghang\4_Color\Complex_Syn_4_shell\';

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
for cur_path = 1:18
    %cur_path = 1;
    disp(cur_path);
    mergedpath = [char(pathname(cur_path)) 'analysis\elastic_align\conv_488\'];
    mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
    num_images = numel(mergedfiles);
    info = imfinfo([mergedpath mergedfiles(1,1).name]);
    voxel=[15.5, 15.5, 70];
    Syn_path = [char(pathname(cur_path)) 'analysis\Result\4_CTB\'];
    %
    %
    % Pos_near
    load([Syn_path 'R_paired_VC.mat']);
    load([data_path,'Vwater_ss_' sprintf('%03d',cur_path) '.mat']);
    %statsRwater_sssss and statsRwater_ssssn
    %statsVwater_ss and stastVwater_sn
    BID_simp = [];
    BID_comp = [];
    simp_ss = [];
    comp_ss = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = sort(statsVwater_ss(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            BID_comp = cat(1,BID_comp,statsVwater_ss(i).B_ID);
            comp_ss = cat(1,comp_ss,statsVwater_ss(i));
        else
            BID_simp = cat(1,BID_simp,statsVwater_ss(i).B_ID);
            simp_ss = cat(1,simp_ss,statsVwater_ss(i));
        end
    end
    BID_comp = unique(BID_comp);
    BID_simp = unique(BID_simp);
    BID_simp = BID_simp>0;
    statsRwater_sssss = statsRwater_sssss(BID_simp);
    Weighted_centroid_1 = [];
    Weighted_centroid_2 = [];
    for i = 1:numel(statsRwater_sssss)
        Weighted_centroid_1 = cat(1,Weighted_centroid_1,statsRwater_sssss(i).WeightedCentroid);
    end
    for i = 1:numel(comp_ss)
        Weighted_centroid_2 = cat(1,Weighted_centroid_2,comp_ss(i).WeightedCentroid);
    end
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
    min_dist_list_12 = zeros(1,size(dist_list_12,1));
    for i = 1:size(dist_list_12)
        temp_list = dist_list_12(i,:);
        temp_list = sort(temp_list);
        min_dist_list_12(i) = temp_list(1);
    end
    threshold = 1.5;
    statsRwater_sssss = statsRwater_sssss(min_dist_list_12 < threshold);
    clear BP BG bg2
    disp('allocating arrays')
    BP = zeros(info.Height, info.Width, num_images,'uint8');
    BP2 = zeros(info.Height, info.Width, num_images,'uint8');
    for i = 1:numel(statsRwater_sssss)
        %disp(int2str(i))
        BP(statsRwater_sssss(i).PixelIdxList)=statsRwater_sssss(i).PixelValues;
    end
    for i = 1:numel(statsVwater_ss)
        %disp(int2str(i))
        BP2(statsVwater_ss(i).PixelIdxList)=statsVwater_ss(i).PixelValues;
    end
    disp('loading data')
    for i=numel(statsRwater_sssss):-1:1
        %disp(i)
        %Area means all pixels while volume contains only positive pixels. 
        statsRwater_sssss(i).tints_p16 = [];
        statsRwater_sssss(i).volume_p16 = [];
        statsRwater_sssss(i).area_p16 = [];
        statsRwater_sssss(i).WeightedCentroid_p16 = [];
        statsRwater_sssss(i).tints_p32 = [];
        statsRwater_sssss(i).volume_p32 = [];
        statsRwater_sssss(i).area_p32 = [];
        statsRwater_sssss(i).WeightedCentroid_p32 = [];
        statsRwater_sssss(i).tints_p48 = [];
        statsRwater_sssss(i).volume_p48 = [];
        statsRwater_sssss(i).area_p48 = [];
        statsRwater_sssss(i).WeightedCentroid_p48 = [];
        statsRwater_sssss(i).tints_p64 = [];
        statsRwater_sssss(i).volume_p64 = [];
        statsRwater_sssss(i).area_p64 = [];
        statsRwater_sssss(i).WeightedCentroid_p64 = [];
        statsRwater_sssss(i).tints_p80 = [];
        statsRwater_sssss(i).volume_p80 = [];
        statsRwater_sssss(i).area_p80 = [];
        statsRwater_sssss(i).WeightedCentroid_p80 = [];
        minpix = min(statsRwater_sssss(i).PixelList);  maxpix = max(statsRwater_sssss(i).PixelList);
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
    parfor jj=1:numel(statsRwater_sssss)
        minpix = min(statsRwater_sssss(jj).PixelList);  maxpix = max(statsRwater_sssss(jj).PixelList);
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
        for j=1: numel(statsRwater_sssss(jj).PixelList(:,1))
            curr2(statsRwater_sssss(jj).PixelList(j,2)-min2+1, ...
                statsRwater_sssss(jj).PixelList(j,1)-min1+1, ...
                statsRwater_sssss(jj).PixelList(j,3)-min3+1)= 1;
        end
        curr1a = BPp(jj).mat;
        curr1b = BP2p(jj).mat;
        Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
        size2= 16;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m16 = [x,y,z];
        RPixelValues_m16 = zeros(numel(PixelList_m16(:,1)),1);
        for kk=1:numel(PixelList_m16(:,1))
            RPixelValues_m16(kk) = curr1b(PixelList_m16(kk,2),...
                PixelList_m16(kk,1), PixelList_m16(kk,3));
        end
        statsRwater_sssss(jj).tints_p16 = sum([RPixelValues_m16]);
        statsRwater_sssss(jj).area_p16 = numel([RPixelValues_m16]);
        statsRwater_sssss(jj).volume_p16 = numel(find(RPixelValues_m16));
        statsRwater_sssss(jj).WeightedCentroid_p16(1) = sum([PixelList_m16(:,1)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_sssss(jj).WeightedCentroid_p16(2) = sum([PixelList_m16(:,2)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_sssss(jj).WeightedCentroid_p16(3) =  sum([PixelList_m16(:,3)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        size2= 32;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m32 = [x,y,z];
        RPixelValues_m32 = zeros(numel(PixelList_m32(:,1)),1);
        for kk=1:numel(PixelList_m32(:,1))
            RPixelValues_m32(kk) = curr1b(PixelList_m32(kk,2),...
                PixelList_m32(kk,1), PixelList_m32(kk,3));
        end
        statsRwater_sssss(jj).tints_p32 = sum([RPixelValues_m32]);
        statsRwater_sssss(jj).area_p32 = numel([RPixelValues_m32]);
        statsRwater_sssss(jj).volume_p32 = numel(find(RPixelValues_m32));
        statsRwater_sssss(jj).WeightedCentroid_p32(1) = sum([PixelList_m32(:,1)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_sssss(jj).WeightedCentroid_p32(2) = sum([PixelList_m32(:,2)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_sssss(jj).WeightedCentroid_p32(3) =  sum([PixelList_m32(:,3)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        size2= 48;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m48 = [x,y,z];
        RPixelValues_m48 = zeros(numel(PixelList_m48(:,1)),1);
        for kk=1:numel(PixelList_m48(:,1))
            RPixelValues_m48(kk) = curr1b(PixelList_m48(kk,2),...
                PixelList_m48(kk,1), PixelList_m48(kk,3));
        end
        statsRwater_sssss(jj).tints_p48 = sum([RPixelValues_m48]);
        statsRwater_sssss(jj).area_p48 = numel([RPixelValues_m48]);
        statsRwater_sssss(jj).volume_p48 = numel(find(RPixelValues_m48));
        statsRwater_sssss(jj).WeightedCentroid_p48(1) = sum([PixelList_m48(:,1)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_sssss(jj).WeightedCentroid_p48(2) = sum([PixelList_m48(:,2)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_sssss(jj).WeightedCentroid_p48(3) =  sum([PixelList_m48(:,3)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        size2= 64;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m64 = [x,y,z];
        RPixelValues_m64 = zeros(numel(PixelList_m64(:,1)),1);
        for kk=1:numel(PixelList_m64(:,1))
            RPixelValues_m64(kk) = curr1b(PixelList_m64(kk,2),...
                PixelList_m64(kk,1), PixelList_m64(kk,3));
        end
        statsRwater_sssss(jj).tints_p64 = sum([RPixelValues_m64]);
        statsRwater_sssss(jj).area_p64 = numel([RPixelValues_m64]);
        statsRwater_sssss(jj).volume_p64 = numel(find(RPixelValues_m64));
        statsRwater_sssss(jj).WeightedCentroid_p64(1) = sum([PixelList_m64(:,1)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_sssss(jj).WeightedCentroid_p64(2) = sum([PixelList_m64(:,2)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_sssss(jj).WeightedCentroid_p64(3) =  sum([PixelList_m64(:,3)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        size2= 80;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m80 = [x,y,z];
        RPixelValues_m80 = zeros(numel(PixelList_m80(:,1)),1);
        for kk=1:numel(PixelList_m80(:,1))
            RPixelValues_m80(kk) = curr1b(PixelList_m80(kk,2),...
                PixelList_m80(kk,1), PixelList_m80(kk,3));
        end
        statsRwater_sssss(jj).tints_p80 = sum([RPixelValues_m80]);
        statsRwater_sssss(jj).area_p80 = numel([RPixelValues_m80]);
        statsRwater_sssss(jj).volume_p80 = numel(find(RPixelValues_m80));
        statsRwater_sssss(jj).WeightedCentroid_p80(1) = sum([PixelList_m80(:,1)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_sssss(jj).WeightedCentroid_p80(2) = sum([PixelList_m80(:,2)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_sssss(jj).WeightedCentroid_p80(3) =  sum([PixelList_m80(:,3)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
    end
    area = [statsRwater_sssss.area_p16];
    tints = [statsRwater_sssss.tints_p16];
    volume = [statsRwater_sssss.volume_p16];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p16.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p16.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p16.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p16.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p32];
    tints = [statsRwater_sssss.tints_p32];
    volume = [statsRwater_sssss.volume_p32];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p32.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p32.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p32.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p32.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p48];
    tints = [statsRwater_sssss.tints_p48];
    volume = [statsRwater_sssss.volume_p48];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p48.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p48.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p48.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p48.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p64];
    tints = [statsRwater_sssss.tints_p64];
    volume = [statsRwater_sssss.volume_p64];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p64.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p64.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p64.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p64.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p80];
    tints = [statsRwater_sssss.tints_p80];
    volume = [statsRwater_sssss.volume_p80];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p80.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p80.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p80.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p80.csv'],'WriteMode','append');

    %
    %
    % Pos_far
    load([Syn_path 'R_paired_VC.mat']);
    load([data_path,'Vwater_ss_' sprintf('%03d',cur_path) '.mat']);
    %statsRwater_sssss and statsRwater_ssssn
    %statsVwater_ss and stastVwater_sn
    BID_simp = [];
    BID_comp = [];
    simp_ss = [];
    comp_ss = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = sort(statsVwater_ss(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            BID_comp = cat(1,BID_comp,statsVwater_ss(i).B_ID);
            comp_ss = cat(1,comp_ss,statsVwater_ss(i));
        else
            BID_simp = cat(1,BID_simp,statsVwater_ss(i).B_ID);
            simp_ss = cat(1,simp_ss,statsVwater_ss(i));
        end
    end
    BID_comp = unique(BID_comp);
    BID_simp = unique(BID_simp);
    BID_simp = BID_simp>0;
    statsRwater_sssss = statsRwater_sssss(BID_simp);
    Weighted_centroid_1 = [];
    Weighted_centroid_2 = [];
    for i = 1:numel(statsRwater_sssss)
        Weighted_centroid_1 = cat(1,Weighted_centroid_1,statsRwater_sssss(i).WeightedCentroid);
    end
    for i = 1:numel(comp_ss)
        Weighted_centroid_2 = cat(1,Weighted_centroid_2,comp_ss(i).WeightedCentroid);
    end
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
    min_dist_list_12 = zeros(1,size(dist_list_12,1));
    for i = 1:size(dist_list_12)
        temp_list = dist_list_12(i,:);
        temp_list = sort(temp_list);
        min_dist_list_12(i) = temp_list(1);
    end
    threshold = 1.5;
    statsRwater_sssss = statsRwater_sssss(min_dist_list_12 >= threshold);
    clear BP BG bg2
    disp('allocating arrays')
    BP = zeros(info.Height, info.Width, num_images,'uint8');
    BP2 = zeros(info.Height, info.Width, num_images,'uint8');
    for i = 1:numel(statsRwater_sssss)
        %disp(int2str(i))
        BP(statsRwater_sssss(i).PixelIdxList)=statsRwater_sssss(i).PixelValues;
    end
    for i = 1:numel(statsVwater_ss)
        %disp(int2str(i))
        BP2(statsVwater_ss(i).PixelIdxList)=statsVwater_ss(i).PixelValues;
    end
    disp('loading data')
    for i=numel(statsRwater_sssss):-1:1
        %disp(i)
        %Area means all pixels while volume contains only positive pixels. 
        statsRwater_sssss(i).tints_p16 = [];
        statsRwater_sssss(i).volume_p16 = [];
        statsRwater_sssss(i).area_p16 = [];
        statsRwater_sssss(i).WeightedCentroid_p16 = [];
        statsRwater_sssss(i).tints_p32 = [];
        statsRwater_sssss(i).volume_p32 = [];
        statsRwater_sssss(i).area_p32 = [];
        statsRwater_sssss(i).WeightedCentroid_p32 = [];
        statsRwater_sssss(i).tints_p48 = [];
        statsRwater_sssss(i).volume_p48 = [];
        statsRwater_sssss(i).area_p48 = [];
        statsRwater_sssss(i).WeightedCentroid_p48 = [];
        statsRwater_sssss(i).tints_p64 = [];
        statsRwater_sssss(i).volume_p64 = [];
        statsRwater_sssss(i).area_p64 = [];
        statsRwater_sssss(i).WeightedCentroid_p64 = [];
        statsRwater_sssss(i).tints_p80 = [];
        statsRwater_sssss(i).volume_p80 = [];
        statsRwater_sssss(i).area_p80 = [];
        statsRwater_sssss(i).WeightedCentroid_p80 = [];
        minpix = min(statsRwater_sssss(i).PixelList);  maxpix = max(statsRwater_sssss(i).PixelList);
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
    parfor jj=1:numel(statsRwater_sssss)
        minpix = min(statsRwater_sssss(jj).PixelList);  maxpix = max(statsRwater_sssss(jj).PixelList);
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
        for j=1: numel(statsRwater_sssss(jj).PixelList(:,1))
            curr2(statsRwater_sssss(jj).PixelList(j,2)-min2+1, ...
                statsRwater_sssss(jj).PixelList(j,1)-min1+1, ...
                statsRwater_sssss(jj).PixelList(j,3)-min3+1)= 1;
        end
        curr1a = BPp(jj).mat;
        curr1b = BP2p(jj).mat;
        Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
        size2= 16;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m16 = [x,y,z];
        RPixelValues_m16 = zeros(numel(PixelList_m16(:,1)),1);
        for kk=1:numel(PixelList_m16(:,1))
            RPixelValues_m16(kk) = curr1b(PixelList_m16(kk,2),...
                PixelList_m16(kk,1), PixelList_m16(kk,3));
        end
        statsRwater_sssss(jj).tints_p16 = sum([RPixelValues_m16]);
        statsRwater_sssss(jj).area_p16 = numel([RPixelValues_m16]);
        statsRwater_sssss(jj).volume_p16 = numel(find(RPixelValues_m16));
        statsRwater_sssss(jj).WeightedCentroid_p16(1) = sum([PixelList_m16(:,1)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_sssss(jj).WeightedCentroid_p16(2) = sum([PixelList_m16(:,2)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_sssss(jj).WeightedCentroid_p16(3) =  sum([PixelList_m16(:,3)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        size2= 32;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m32 = [x,y,z];
        RPixelValues_m32 = zeros(numel(PixelList_m32(:,1)),1);
        for kk=1:numel(PixelList_m32(:,1))
            RPixelValues_m32(kk) = curr1b(PixelList_m32(kk,2),...
                PixelList_m32(kk,1), PixelList_m32(kk,3));
        end
        statsRwater_sssss(jj).tints_p32 = sum([RPixelValues_m32]);
        statsRwater_sssss(jj).area_p32 = numel([RPixelValues_m32]);
        statsRwater_sssss(jj).volume_p32 = numel(find(RPixelValues_m32));
        statsRwater_sssss(jj).WeightedCentroid_p32(1) = sum([PixelList_m32(:,1)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_sssss(jj).WeightedCentroid_p32(2) = sum([PixelList_m32(:,2)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_sssss(jj).WeightedCentroid_p32(3) =  sum([PixelList_m32(:,3)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        size2= 48;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m48 = [x,y,z];
        RPixelValues_m48 = zeros(numel(PixelList_m48(:,1)),1);
        for kk=1:numel(PixelList_m48(:,1))
            RPixelValues_m48(kk) = curr1b(PixelList_m48(kk,2),...
                PixelList_m48(kk,1), PixelList_m48(kk,3));
        end
        statsRwater_sssss(jj).tints_p48 = sum([RPixelValues_m48]);
        statsRwater_sssss(jj).area_p48 = numel([RPixelValues_m48]);
        statsRwater_sssss(jj).volume_p48 = numel(find(RPixelValues_m48));
        statsRwater_sssss(jj).WeightedCentroid_p48(1) = sum([PixelList_m48(:,1)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_sssss(jj).WeightedCentroid_p48(2) = sum([PixelList_m48(:,2)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_sssss(jj).WeightedCentroid_p48(3) =  sum([PixelList_m48(:,3)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        size2= 64;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m64 = [x,y,z];
        RPixelValues_m64 = zeros(numel(PixelList_m64(:,1)),1);
        for kk=1:numel(PixelList_m64(:,1))
            RPixelValues_m64(kk) = curr1b(PixelList_m64(kk,2),...
                PixelList_m64(kk,1), PixelList_m64(kk,3));
        end
        statsRwater_sssss(jj).tints_p64 = sum([RPixelValues_m64]);
        statsRwater_sssss(jj).area_p64 = numel([RPixelValues_m64]);
        statsRwater_sssss(jj).volume_p64 = numel(find(RPixelValues_m64));
        statsRwater_sssss(jj).WeightedCentroid_p64(1) = sum([PixelList_m64(:,1)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_sssss(jj).WeightedCentroid_p64(2) = sum([PixelList_m64(:,2)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_sssss(jj).WeightedCentroid_p64(3) =  sum([PixelList_m64(:,3)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        size2= 80;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m80 = [x,y,z];
        RPixelValues_m80 = zeros(numel(PixelList_m80(:,1)),1);
        for kk=1:numel(PixelList_m80(:,1))
            RPixelValues_m80(kk) = curr1b(PixelList_m80(kk,2),...
                PixelList_m80(kk,1), PixelList_m80(kk,3));
        end
        statsRwater_sssss(jj).tints_p80 = sum([RPixelValues_m80]);
        statsRwater_sssss(jj).area_p80 = numel([RPixelValues_m80]);
        statsRwater_sssss(jj).volume_p80 = numel(find(RPixelValues_m80));
        statsRwater_sssss(jj).WeightedCentroid_p80(1) = sum([PixelList_m80(:,1)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_sssss(jj).WeightedCentroid_p80(2) = sum([PixelList_m80(:,2)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_sssss(jj).WeightedCentroid_p80(3) =  sum([PixelList_m80(:,3)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
    end
    area = [statsRwater_sssss.area_p16];
    tints = [statsRwater_sssss.tints_p16];
    volume = [statsRwater_sssss.volume_p16];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p16.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p16.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p16.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p16.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p32];
    tints = [statsRwater_sssss.tints_p32];
    volume = [statsRwater_sssss.volume_p32];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p32.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p32.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p32.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p32.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p48];
    tints = [statsRwater_sssss.tints_p48];
    volume = [statsRwater_sssss.volume_p48];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p48.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p48.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p48.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p48.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p64];
    tints = [statsRwater_sssss.tints_p64];
    volume = [statsRwater_sssss.volume_p64];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p64.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p64.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p64.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p64.csv'],'WriteMode','append');
    area = [statsRwater_sssss.area_p80];
    tints = [statsRwater_sssss.tints_p80];
    volume = [statsRwater_sssss.volume_p80];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p80.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p80.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p80.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_sssss)],[out_path 'Null_num_p80.csv'],'WriteMode','append');
    %
    %
    % Neg_near
    load([Syn_path 'R_paired_VC.mat']);
    load([data_path,'Vwater_sn_' sprintf('%03d',cur_path) '.mat']);
    %statsRwater_ssssn and statsRwater_ssssn
    %statsVwater_sn and stastVwater_sn
    BID_simp = [];
    BID_comp = [];
    simp_sn = [];
    comp_sn = [];
    for i = 1:numel(statsVwater_sn)
        B_ID = sort(statsVwater_sn(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            BID_comp = cat(1,BID_comp,statsVwater_sn(i).B_ID);
            comp_sn = cat(1,comp_sn,statsVwater_sn(i));
        else
            BID_simp = cat(1,BID_simp,statsVwater_sn(i).B_ID);
            simp_sn = cat(1,simp_sn,statsVwater_sn(i));
        end
    end
    BID_comp = unique(BID_comp);
    BID_simp = unique(BID_simp);
    BID_simp = BID_simp>0;
    statsRwater_ssssn = statsRwater_ssssn(BID_simp);
    Weighted_centroid_1 = [];
    Weighted_centroid_2 = [];
    for i = 1:numel(statsRwater_ssssn)
        Weighted_centroid_1 = cat(1,Weighted_centroid_1,statsRwater_ssssn(i).WeightedCentroid);
    end
    for i = 1:numel(comp_sn)
        Weighted_centroid_2 = cat(1,Weighted_centroid_2,comp_sn(i).WeightedCentroid);
    end
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
    min_dist_list_12 = zeros(1,size(dist_list_12,1));
    for i = 1:size(dist_list_12)
        temp_list = dist_list_12(i,:);
        temp_list = sort(temp_list);
        min_dist_list_12(i) = temp_list(1);
    end
    threshold = 1.5;
    statsRwater_ssssn = statsRwater_ssssn(min_dist_list_12 < threshold);
    clear BP BG bg2
    disp('allocating arrays')
    BP = zeros(info.Height, info.Width, num_images,'uint8');
    BP2 = zeros(info.Height, info.Width, num_images,'uint8');
    for i = 1:numel(statsRwater_ssssn)
        %disp(int2str(i))
        BP(statsRwater_ssssn(i).PixelIdxList)=statsRwater_ssssn(i).PixelValues;
    end
    for i = 1:numel(statsVwater_sn)
        %disp(int2str(i))
        BP2(statsVwater_sn(i).PixelIdxList)=statsVwater_sn(i).PixelValues;
    end
    disp('loading data')
    for i=numel(statsRwater_ssssn):-1:1
        %disp(i)
        %Area means all pixels while volume contains only positive pixels. 
        statsRwater_ssssn(i).tints_p16 = [];
        statsRwater_ssssn(i).volume_p16 = [];
        statsRwater_ssssn(i).area_p16 = [];
        statsRwater_ssssn(i).WeightedCentroid_p16 = [];
        statsRwater_ssssn(i).tints_p32 = [];
        statsRwater_ssssn(i).volume_p32 = [];
        statsRwater_ssssn(i).area_p32 = [];
        statsRwater_ssssn(i).WeightedCentroid_p32 = [];
        statsRwater_ssssn(i).tints_p48 = [];
        statsRwater_ssssn(i).volume_p48 = [];
        statsRwater_ssssn(i).area_p48 = [];
        statsRwater_ssssn(i).WeightedCentroid_p48 = [];
        statsRwater_ssssn(i).tints_p64 = [];
        statsRwater_ssssn(i).volume_p64 = [];
        statsRwater_ssssn(i).area_p64 = [];
        statsRwater_ssssn(i).WeightedCentroid_p64 = [];
        statsRwater_ssssn(i).tints_p80 = [];
        statsRwater_ssssn(i).volume_p80 = [];
        statsRwater_ssssn(i).area_p80 = [];
        statsRwater_ssssn(i).WeightedCentroid_p80 = [];
        minpix = min(statsRwater_ssssn(i).PixelList);  maxpix = max(statsRwater_ssssn(i).PixelList);
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
    parfor jj=1:numel(statsRwater_ssssn)
        minpix = min(statsRwater_ssssn(jj).PixelList);  maxpix = max(statsRwater_ssssn(jj).PixelList);
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
        for j=1: numel(statsRwater_ssssn(jj).PixelList(:,1))
            curr2(statsRwater_ssssn(jj).PixelList(j,2)-min2+1, ...
                statsRwater_ssssn(jj).PixelList(j,1)-min1+1, ...
                statsRwater_ssssn(jj).PixelList(j,3)-min3+1)= 1;
        end
        curr1a = BPp(jj).mat;
        curr1b = BP2p(jj).mat;
        Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
        size2= 16;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m16 = [x,y,z];
        RPixelValues_m16 = zeros(numel(PixelList_m16(:,1)),1);
        for kk=1:numel(PixelList_m16(:,1))
            RPixelValues_m16(kk) = curr1b(PixelList_m16(kk,2),...
                PixelList_m16(kk,1), PixelList_m16(kk,3));
        end
        statsRwater_ssssn(jj).tints_p16 = sum([RPixelValues_m16]);
        statsRwater_ssssn(jj).area_p16 = numel([RPixelValues_m16]);
        statsRwater_ssssn(jj).volume_p16 = numel(find(RPixelValues_m16));
        statsRwater_ssssn(jj).WeightedCentroid_p16(1) = sum([PixelList_m16(:,1)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_ssssn(jj).WeightedCentroid_p16(2) = sum([PixelList_m16(:,2)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_ssssn(jj).WeightedCentroid_p16(3) =  sum([PixelList_m16(:,3)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        size2= 32;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m32 = [x,y,z];
        RPixelValues_m32 = zeros(numel(PixelList_m32(:,1)),1);
        for kk=1:numel(PixelList_m32(:,1))
            RPixelValues_m32(kk) = curr1b(PixelList_m32(kk,2),...
                PixelList_m32(kk,1), PixelList_m32(kk,3));
        end
        statsRwater_ssssn(jj).tints_p32 = sum([RPixelValues_m32]);
        statsRwater_ssssn(jj).area_p32 = numel([RPixelValues_m32]);
        statsRwater_ssssn(jj).volume_p32 = numel(find(RPixelValues_m32));
        statsRwater_ssssn(jj).WeightedCentroid_p32(1) = sum([PixelList_m32(:,1)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_ssssn(jj).WeightedCentroid_p32(2) = sum([PixelList_m32(:,2)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_ssssn(jj).WeightedCentroid_p32(3) =  sum([PixelList_m32(:,3)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        size2= 48;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m48 = [x,y,z];
        RPixelValues_m48 = zeros(numel(PixelList_m48(:,1)),1);
        for kk=1:numel(PixelList_m48(:,1))
            RPixelValues_m48(kk) = curr1b(PixelList_m48(kk,2),...
                PixelList_m48(kk,1), PixelList_m48(kk,3));
        end
        statsRwater_ssssn(jj).tints_p48 = sum([RPixelValues_m48]);
        statsRwater_ssssn(jj).area_p48 = numel([RPixelValues_m48]);
        statsRwater_ssssn(jj).volume_p48 = numel(find(RPixelValues_m48));
        statsRwater_ssssn(jj).WeightedCentroid_p48(1) = sum([PixelList_m48(:,1)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_ssssn(jj).WeightedCentroid_p48(2) = sum([PixelList_m48(:,2)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_ssssn(jj).WeightedCentroid_p48(3) =  sum([PixelList_m48(:,3)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        size2= 64;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m64 = [x,y,z];
        RPixelValues_m64 = zeros(numel(PixelList_m64(:,1)),1);
        for kk=1:numel(PixelList_m64(:,1))
            RPixelValues_m64(kk) = curr1b(PixelList_m64(kk,2),...
                PixelList_m64(kk,1), PixelList_m64(kk,3));
        end
        statsRwater_ssssn(jj).tints_p64 = sum([RPixelValues_m64]);
        statsRwater_ssssn(jj).area_p64 = numel([RPixelValues_m64]);
        statsRwater_ssssn(jj).volume_p64 = numel(find(RPixelValues_m64));
        statsRwater_ssssn(jj).WeightedCentroid_p64(1) = sum([PixelList_m64(:,1)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_ssssn(jj).WeightedCentroid_p64(2) = sum([PixelList_m64(:,2)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_ssssn(jj).WeightedCentroid_p64(3) =  sum([PixelList_m64(:,3)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        size2= 80;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m80 = [x,y,z];
        RPixelValues_m80 = zeros(numel(PixelList_m80(:,1)),1);
        for kk=1:numel(PixelList_m80(:,1))
            RPixelValues_m80(kk) = curr1b(PixelList_m80(kk,2),...
                PixelList_m80(kk,1), PixelList_m80(kk,3));
        end
        statsRwater_ssssn(jj).tints_p80 = sum([RPixelValues_m80]);
        statsRwater_ssssn(jj).area_p80 = numel([RPixelValues_m80]);
        statsRwater_ssssn(jj).volume_p80 = numel(find(RPixelValues_m80));
        statsRwater_ssssn(jj).WeightedCentroid_p80(1) = sum([PixelList_m80(:,1)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_ssssn(jj).WeightedCentroid_p80(2) = sum([PixelList_m80(:,2)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_ssssn(jj).WeightedCentroid_p80(3) =  sum([PixelList_m80(:,3)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
    end
    area = [statsRwater_ssssn.area_p16];
    tints = [statsRwater_ssssn.tints_p16];
    volume = [statsRwater_ssssn.volume_p16];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p16.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p16.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p16.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p16.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p32];
    tints = [statsRwater_ssssn.tints_p32];
    volume = [statsRwater_ssssn.volume_p32];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p32.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p32.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p32.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p32.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p48];
    tints = [statsRwater_ssssn.tints_p48];
    volume = [statsRwater_ssssn.volume_p48];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p48.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p48.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p48.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p48.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p64];
    tints = [statsRwater_ssssn.tints_p64];
    volume = [statsRwater_ssssn.volume_p64];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p64.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p64.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p64.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p64.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p80];
    tints = [statsRwater_ssssn.tints_p80];
    volume = [statsRwater_ssssn.volume_p80];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p80.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p80.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p80.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p80.csv'],'WriteMode','append');

    %
    %
    % Neg_far
    load([Syn_path 'R_paired_VC.mat']);
    load([data_path,'Vwater_sn_' sprintf('%03d',cur_path) '.mat']);
    %statsRwater_ssssn and statsRwater_ssssn
    %statsVwater_sn and stastVwater_sn
    BID_simp = [];
    BID_comp = [];
    simp_sn = [];
    comp_sn = [];
    for i = 1:numel(statsVwater_sn)
        B_ID = sort(statsVwater_sn(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            BID_comp = cat(1,BID_comp,statsVwater_sn(i).B_ID);
            comp_sn = cat(1,comp_sn,statsVwater_sn(i));
        else
            BID_simp = cat(1,BID_simp,statsVwater_sn(i).B_ID);
            simp_sn = cat(1,simp_sn,statsVwater_sn(i));
        end
    end
    BID_comp = unique(BID_comp);
    BID_simp = unique(BID_simp);
    BID_simp = BID_simp>0;
    statsRwater_ssssn = statsRwater_ssssn(BID_simp);
    Weighted_centroid_1 = [];
    Weighted_centroid_2 = [];
    for i = 1:numel(statsRwater_ssssn)
        Weighted_centroid_1 = cat(1,Weighted_centroid_1,statsRwater_ssssn(i).WeightedCentroid);
    end
    for i = 1:numel(comp_sn)
        Weighted_centroid_2 = cat(1,Weighted_centroid_2,comp_sn(i).WeightedCentroid);
    end
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
    min_dist_list_12 = zeros(1,size(dist_list_12,1));
    for i = 1:size(dist_list_12)
        temp_list = dist_list_12(i,:);
        temp_list = sort(temp_list);
        min_dist_list_12(i) = temp_list(1);
    end
    threshold = 1.5;
    statsRwater_ssssn = statsRwater_ssssn(min_dist_list_12 >= threshold);
    clear BP BG bg2
    disp('allocating arrays')
    BP = zeros(info.Height, info.Width, num_images,'uint8');
    BP2 = zeros(info.Height, info.Width, num_images,'uint8');
    for i = 1:numel(statsRwater_ssssn)
        %disp(int2str(i))
        BP(statsRwater_ssssn(i).PixelIdxList)=statsRwater_ssssn(i).PixelValues;
    end
    for i = 1:numel(statsVwater_sn)
        %disp(int2str(i))
        BP2(statsVwater_sn(i).PixelIdxList)=statsVwater_sn(i).PixelValues;
    end
    disp('loading data')
    for i=numel(statsRwater_ssssn):-1:1
        %disp(i)
        %Area means all pixels while volume contains only positive pixels. 
        statsRwater_ssssn(i).tints_p16 = [];
        statsRwater_ssssn(i).volume_p16 = [];
        statsRwater_ssssn(i).area_p16 = [];
        statsRwater_ssssn(i).WeightedCentroid_p16 = [];
        statsRwater_ssssn(i).tints_p32 = [];
        statsRwater_ssssn(i).volume_p32 = [];
        statsRwater_ssssn(i).area_p32 = [];
        statsRwater_ssssn(i).WeightedCentroid_p32 = [];
        statsRwater_ssssn(i).tints_p48 = [];
        statsRwater_ssssn(i).volume_p48 = [];
        statsRwater_ssssn(i).area_p48 = [];
        statsRwater_ssssn(i).WeightedCentroid_p48 = [];
        statsRwater_ssssn(i).tints_p64 = [];
        statsRwater_ssssn(i).volume_p64 = [];
        statsRwater_ssssn(i).area_p64 = [];
        statsRwater_ssssn(i).WeightedCentroid_p64 = [];
        statsRwater_ssssn(i).tints_p80 = [];
        statsRwater_ssssn(i).volume_p80 = [];
        statsRwater_ssssn(i).area_p80 = [];
        statsRwater_ssssn(i).WeightedCentroid_p80 = [];
        minpix = min(statsRwater_ssssn(i).PixelList);  maxpix = max(statsRwater_ssssn(i).PixelList);
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
    parfor jj=1:numel(statsRwater_ssssn)
        minpix = min(statsRwater_ssssn(jj).PixelList);  maxpix = max(statsRwater_ssssn(jj).PixelList);
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
        for j=1: numel(statsRwater_ssssn(jj).PixelList(:,1))
            curr2(statsRwater_ssssn(jj).PixelList(j,2)-min2+1, ...
                statsRwater_ssssn(jj).PixelList(j,1)-min1+1, ...
                statsRwater_ssssn(jj).PixelList(j,3)-min3+1)= 1;
        end
        curr1a = BPp(jj).mat;
        curr1b = BP2p(jj).mat;
        Dg = bwdistsc(curr2,[voxel(1),voxel(2),voxel(3)]);
        size2= 16;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m16 = [x,y,z];
        RPixelValues_m16 = zeros(numel(PixelList_m16(:,1)),1);
        for kk=1:numel(PixelList_m16(:,1))
            RPixelValues_m16(kk) = curr1b(PixelList_m16(kk,2),...
                PixelList_m16(kk,1), PixelList_m16(kk,3));
        end
        statsRwater_ssssn(jj).tints_p16 = sum([RPixelValues_m16]);
        statsRwater_ssssn(jj).area_p16 = numel([RPixelValues_m16]);
        statsRwater_ssssn(jj).volume_p16 = numel(find(RPixelValues_m16));
        statsRwater_ssssn(jj).WeightedCentroid_p16(1) = sum([PixelList_m16(:,1)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_ssssn(jj).WeightedCentroid_p16(2) = sum([PixelList_m16(:,2)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        statsRwater_ssssn(jj).WeightedCentroid_p16(3) =  sum([PixelList_m16(:,3)].*...
            double([RPixelValues_m16]))/(sum([RPixelValues_m16]));
        size2= 32;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m32 = [x,y,z];
        RPixelValues_m32 = zeros(numel(PixelList_m32(:,1)),1);
        for kk=1:numel(PixelList_m32(:,1))
            RPixelValues_m32(kk) = curr1b(PixelList_m32(kk,2),...
                PixelList_m32(kk,1), PixelList_m32(kk,3));
        end
        statsRwater_ssssn(jj).tints_p32 = sum([RPixelValues_m32]);
        statsRwater_ssssn(jj).area_p32 = numel([RPixelValues_m32]);
        statsRwater_ssssn(jj).volume_p32 = numel(find(RPixelValues_m32));
        statsRwater_ssssn(jj).WeightedCentroid_p32(1) = sum([PixelList_m32(:,1)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_ssssn(jj).WeightedCentroid_p32(2) = sum([PixelList_m32(:,2)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        statsRwater_ssssn(jj).WeightedCentroid_p32(3) =  sum([PixelList_m32(:,3)].*...
            double([RPixelValues_m32]))/(sum([RPixelValues_m32]));
        size2= 48;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m48 = [x,y,z];
        RPixelValues_m48 = zeros(numel(PixelList_m48(:,1)),1);
        for kk=1:numel(PixelList_m48(:,1))
            RPixelValues_m48(kk) = curr1b(PixelList_m48(kk,2),...
                PixelList_m48(kk,1), PixelList_m48(kk,3));
        end
        statsRwater_ssssn(jj).tints_p48 = sum([RPixelValues_m48]);
        statsRwater_ssssn(jj).area_p48 = numel([RPixelValues_m48]);
        statsRwater_ssssn(jj).volume_p48 = numel(find(RPixelValues_m48));
        statsRwater_ssssn(jj).WeightedCentroid_p48(1) = sum([PixelList_m48(:,1)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_ssssn(jj).WeightedCentroid_p48(2) = sum([PixelList_m48(:,2)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        statsRwater_ssssn(jj).WeightedCentroid_p48(3) =  sum([PixelList_m48(:,3)].*...
            double([RPixelValues_m48]))/(sum([RPixelValues_m48]));
        size2= 64;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m64 = [x,y,z];
        RPixelValues_m64 = zeros(numel(PixelList_m64(:,1)),1);
        for kk=1:numel(PixelList_m64(:,1))
            RPixelValues_m64(kk) = curr1b(PixelList_m64(kk,2),...
                PixelList_m64(kk,1), PixelList_m64(kk,3));
        end
        statsRwater_ssssn(jj).tints_p64 = sum([RPixelValues_m64]);
        statsRwater_ssssn(jj).area_p64 = numel([RPixelValues_m64]);
        statsRwater_ssssn(jj).volume_p64 = numel(find(RPixelValues_m64));
        statsRwater_ssssn(jj).WeightedCentroid_p64(1) = sum([PixelList_m64(:,1)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_ssssn(jj).WeightedCentroid_p64(2) = sum([PixelList_m64(:,2)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        statsRwater_ssssn(jj).WeightedCentroid_p64(3) =  sum([PixelList_m64(:,3)].*...
            double([RPixelValues_m64]))/(sum([RPixelValues_m64]));
        size2= 80;
        curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
        [y,x,z] = ind2sub(size(curr4b),find(curr4b));
        PixelList_m80 = [x,y,z];
        RPixelValues_m80 = zeros(numel(PixelList_m80(:,1)),1);
        for kk=1:numel(PixelList_m80(:,1))
            RPixelValues_m80(kk) = curr1b(PixelList_m80(kk,2),...
                PixelList_m80(kk,1), PixelList_m80(kk,3));
        end
        statsRwater_ssssn(jj).tints_p80 = sum([RPixelValues_m80]);
        statsRwater_ssssn(jj).area_p80 = numel([RPixelValues_m80]);
        statsRwater_ssssn(jj).volume_p80 = numel(find(RPixelValues_m80));
        statsRwater_ssssn(jj).WeightedCentroid_p80(1) = sum([PixelList_m80(:,1)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_ssssn(jj).WeightedCentroid_p80(2) = sum([PixelList_m80(:,2)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
        statsRwater_ssssn(jj).WeightedCentroid_p80(3) =  sum([PixelList_m80(:,3)].*...
            double([RPixelValues_m80]))/(sum([RPixelValues_m80]));
    end
    area = [statsRwater_ssssn.area_p16];
    tints = [statsRwater_ssssn.tints_p16];
    volume = [statsRwater_ssssn.volume_p16];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p16.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p16.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p16.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p16.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p32];
    tints = [statsRwater_ssssn.tints_p32];
    volume = [statsRwater_ssssn.volume_p32];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p32.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p32.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p32.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p32.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p48];
    tints = [statsRwater_ssssn.tints_p48];
    volume = [statsRwater_ssssn.volume_p48];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p48.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p48.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p48.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p48.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p64];
    tints = [statsRwater_ssssn.tints_p64];
    volume = [statsRwater_ssssn.volume_p64];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p64.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p64.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p64.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p64.csv'],'WriteMode','append');
    area = [statsRwater_ssssn.area_p80];
    tints = [statsRwater_ssssn.tints_p80];
    volume = [statsRwater_ssssn.volume_p80];
    num_null = numel(find(volume == 0));
    area = area(volume>0);
    tints = tints(volume>0);
    volume = volume(volume>0);
    writematrix(area,[out_path 'Area_p80.csv'],'WriteMode','append');
    writematrix(tints,[out_path 'Tints_p80.csv'],'WriteMode','append');
    writematrix(volume,[out_path 'Volume_p80.csv'],'WriteMode','append');
    writematrix([num_null,numel(statsRwater_ssssn)],[out_path 'Null_num_p80.csv'],'WriteMode','append');
end