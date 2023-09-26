clear;clc;
expfolder = 'X:\Chenghang\Backup_Raw_Data\1.2.2021_P2EA_B\';
analysis_folder = [expfolder 'analysis\'];
Result_folder = [analysis_folder 'Result\'];
V_Syn_folder = [Result_folder '5_V_Syn\'];
VGluT2_folder = [Result_folder '3_Vglut2\'];

load([V_Syn_folder 'R_paired_3.mat']);
load([VGluT2_folder 'V_paired.mat']);
%
Centroid_R = zeros(numel(statsRwater_ssss),3);
for i = 1:numel(statsRwater_ssss)
    Centroid_R(i,:) = statsRwater_ssss(i).WeightedCentroid;
end
Centroid_V = zeros(numel(statsVwater_ss),3);
for i = 1:numel(statsVwater_ss)
    Centroid_V(i,:) = statsVwater_ss(i).WeightedCentroid;
end
%%
Dist_num = [];
for i = 1:numel(statsRwater_ssss)
    cur_cent = Centroid_R(i,:);
    cur_dist = [];
    for j = 1:numel(statsRwater_ssss)
        dist_temp = sqrt(((cur_cent(1)-Centroid_R(j,1))*0.0155)^2 +...
            ((cur_cent(2)-Centroid_R(j,2))*0.0155)^2+...
            ((cur_cent(3)-Centroid_R(j,3))*0.07)^2);
        cur_dist = cat(1,cur_dist,dist_temp);
    end
    Dist_num = cat(1,Dist_num,numel(find(cur_dist<0.6)));
end
Dist_num = Dist_num - 1;
%%
figure;hist(Dist_num,max(Dist_num),20);
%%
Closest_V_Volume = [];
for i = 1:numel(statsRwater_ssss)
    cur_cent = Centroid_R(i,:);
    cur_dist = [];
    for j = 1:numel(statsVwater_ss)
        dist_temp = sqrt(((cur_cent(1)-Centroid_V(j,1))*0.0155)^2 +...
            ((cur_cent(2)-Centroid_V(j,2))*0.0155)^2+...
            ((cur_cent(3)-Centroid_V(j,3))*0.07)^2);
        cur_dist = cat(1,cur_dist,dist_temp);
    end
    volume_temp = statsVwater_ss(find(cur_dist == min(cur_dist))).TintsG;
    Closest_V_Volume = cat(1,Closest_V_Volume,volume_temp);
end
%%
figure;scatter(Dist_num,log10(Closest_V_Volume*0.0155*0.0155*0.07),'.');
%%
dist_list = unique(Dist_num);
dist_list = sort(dist_list);
dist_x = [];
volume_y = [];
%for i = 1:numel(dist_list)
for i = 1:4
    Volume_list = Closest_V_Volume(Dist_num == dist_list(i));
    dist_x = cat(1,dist_x,dist_list(i));
    volume_y = cat(1,volume_y,mean(Volume_list));
end
%figure;bar(dist_x,log10(volume_y*0.0155*0.0155*0.07));
figure;bar(dist_x,volume_y);