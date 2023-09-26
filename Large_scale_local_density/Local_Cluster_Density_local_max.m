%%
clear;clc
output_path = 'X:\Chenghang\4_Color_Continue\';

All_max_list_ret = [];
Num_maximum_list_ret = zeros(18,1);
for file_id = 1:18
    disp(file_id);
    load([output_path 'Density_Large_Small_Ratio_Retinal\' sprintf('%03d',file_id) '.mat']);
    load([output_path 'Density_martix\' sprintf('%03d',file_id) '.mat']);

    Local_density_1(~Local_density_1_large_logic) = 0;
    Local_density_1_filt = imgaussfilt(Local_density_1,25);
    BW = imregionalmax(Local_density_1_filt);
    [row_list,column_list] = ind2sub(size(BW),find(BW));
    max_list = cat(2,row_list,column_list);
    max_list = cat(2,max_list,file_id .* ones(size(max_list,1),1));
    Num_maximum = size(max_list,1);
    
    All_max_list_ret = cat(1,All_max_list_ret,max_list);
    Num_maximum_list_ret(file_id) = Num_maximum;
end
save([output_path 'Max_ret_large.mat'],'All_max_list_ret',"Num_maximum_list_ret"');

%
output_path = 'X:\Chenghang\4_Color_Continue\';

All_max_list_nonret = [];
Num_maximum_list_nonret = zeros(18,1);
for file_id = 1:18
    disp(file_id);
    load([output_path 'Density_Large_Small_Ratio_Nonretinal\' sprintf('%03d',file_id) '.mat']);
    load([output_path 'Density_martix\' sprintf('%03d',file_id) '.mat']);

    Local_density_2(~Local_density_2_large_logic) = 0;
    Local_density_2_filt = imgaussfilt(Local_density_2,25);
    BW = imregionalmax(Local_density_2_filt);
    [row_list,column_list] = ind2sub(size(BW),find(BW));
    max_list = cat(2,row_list,column_list);
    max_list = cat(2,max_list,file_id .* ones(size(max_list,1),1));
    Num_maximum = size(max_list,1);
    
    All_max_list_nonret = cat(1,All_max_list_nonret,max_list);
    Num_maximum_list_nonret(file_id) = Num_maximum;
end
save([output_path 'Max_nonret_large.mat'],'All_max_list_nonret',"Num_maximum_list_nonret");
%%
output_path = 'X:\Chenghang\4_Color_Continue\Density_Large_Small_Local_Hopspot\';
load([output_path 'Max_nonret_large.mat']);
load([output_path 'Max_ret_large.mat']);

%%

for i = 1:18
    ret_temp = All_max_list_ret(All_max_list_ret(:,3)==i,1:2);
    nonret_temp = All_max_list_nonret(All_max_list_nonret(:,3)==i,1:2);
    figure;scatter(ret_temp(:,1),ret_temp(:,2),'r.');
    hold on;scatter(nonret_temp(:,1),nonret_temp(:,2),'b.');
end
%%
for i = 1:18
    ret_temp = All_max_list_ret(All_max_list_ret(:,3)==i,1:2);
    nonret_temp = All_max_list_nonret(All_max_list_nonret(:,3)==i,1:2);
end