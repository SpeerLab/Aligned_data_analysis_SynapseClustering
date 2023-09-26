%%
clear;clc
%
data_path = 'X:\Chenghang\4_Color\Complex_Syn\Dist_matrix_comp_simp\';

%%
dist_list_11_all = [];
dist_list_22_all = [];
dist_list_12_all = [];
dist_list_11_new_all = [];
dist_list_22_new_all = [];
dist_list_12_new_all = [];

for file_ID = 1:18
    disp(file_ID);
    load([data_path sprintf('%03d',file_ID) '.mat'],'dist_list_11_sorted','dist_list_12_sorted','dist_list_22_sorted','dist_list_11_sorted_new','dist_list_12_sorted_new','dist_list_22_sorted_new');

    dist_list_11_sorted_mean = mean(dist_list_11_sorted);
    dist_list_22_sorted_mean = mean(dist_list_22_sorted);
    dist_list_12_sorted_mean = mean(dist_list_12_sorted);
    dist_list_11_sorted_new_mean = mean(dist_list_11_sorted_new);
    dist_list_22_sorted_new_mean = mean(dist_list_22_sorted_new);
    dist_list_12_sorted_new_mean = mean(dist_list_12_sorted_new);

    dist_list_11_all = cat(1,dist_list_11_all,dist_list_11_sorted_mean);
    dist_list_22_all = cat(1,dist_list_22_all,dist_list_22_sorted_mean);
    dist_list_12_all = cat(1,dist_list_12_all,dist_list_12_sorted_mean);
    dist_list_11_new_all = cat(1,dist_list_11_new_all,dist_list_11_sorted_new_mean);
    dist_list_22_new_all = cat(1,dist_list_22_new_all,dist_list_22_sorted_new_mean);
    dist_list_12_new_all = cat(1,dist_list_12_new_all,dist_list_12_sorted_new_mean);
end
dist_list_11_all_norm = dist_list_11_all./dist_list_11_new_all;
dist_list_22_all_norm = dist_list_22_all./dist_list_22_new_all;
dist_list_12_all_norm = dist_list_12_all./dist_list_12_new_all;
%%
WTP2_Mean = mean(dist_list_11_all_norm(1:3,:));
WTP2_SE = std(dist_list_11_all_norm(1:3,:)) / sqrt(3);
WTP4_Mean = mean(dist_list_11_all_norm(4:6,:));
WTP4_SE = std(dist_list_11_all_norm(4:6,:)) / sqrt(3);
WTP8_Mean = mean(dist_list_11_all_norm(7:9,:));
WTP8_SE = std(dist_list_11_all_norm(7:9,:)) / sqrt(3);
B2P2_Mean = mean(dist_list_11_all_norm(10:12,:));
B2P2_SE = std(dist_list_11_all_norm(10:12,:)) / sqrt(3);
B2P4_Mean = mean(dist_list_11_all_norm(13:15,:));
B2P4_SE = std(dist_list_11_all_norm(13:15,:)) / sqrt(3);
B2P8_Mean = mean(dist_list_11_all_norm(16:18,:));
B2P8_SE = std(dist_list_11_all_norm(16:18,:)) / sqrt(3);
%
x=1:1:20;
figure;errorbar(x,WTP2_Mean,WTP2_SE,'r');
hold on;errorbar(x,WTP4_Mean,WTP4_SE,'g');
hold on;errorbar(x,WTP8_Mean,WTP8_SE,'b');
ylim([0.0,1]);
xlim([0,10]);
saveas(gcf,[data_path 'WT_11.png']);
saveas(gcf,[data_path 'WT_11.fig']);
saveas(gcf,[data_path 'WT_11.eps']);

close;

figure;errorbar(x,B2P2_Mean,B2P2_SE,'r');
hold on;errorbar(x,B2P4_Mean,B2P4_SE,'g');
hold on;errorbar(x,B2P8_Mean,B2P8_SE,'b');
ylim([0.0,1]);
xlim([0,10]);
saveas(gcf,[data_path 'B2_11.fig']);
saveas(gcf,[data_path 'B2_11.png']);
saveas(gcf,[data_path 'B2_11.eps']);
close;
%%
WTP2_Mean = mean(dist_list_22_all_norm(1:3,:));
WTP2_SE = std(dist_list_22_all_norm(1:3,:)) / sqrt(3);
WTP4_Mean = mean(dist_list_22_all_norm(4:6,:));
WTP4_SE = std(dist_list_22_all_norm(4:6,:)) / sqrt(3);
WTP8_Mean = mean(dist_list_22_all_norm(7:9,:));
WTP8_SE = std(dist_list_22_all_norm(7:9,:)) / sqrt(3);
B2P2_Mean = mean(dist_list_22_all_norm(10:12,:));
B2P2_SE = std(dist_list_22_all_norm(10:12,:)) / sqrt(3);
B2P4_Mean = mean(dist_list_22_all_norm(13:15,:));
B2P4_SE = std(dist_list_22_all_norm(13:15,:)) / sqrt(3);
B2P8_Mean = mean(dist_list_22_all_norm(16:18,:));
B2P8_SE = std(dist_list_22_all_norm(16:18,:)) / sqrt(3);
%
x=1:1:20;
figure;errorbar(x,WTP2_Mean,WTP2_SE,'r');
hold on;errorbar(x,WTP4_Mean,WTP4_SE,'g');
hold on;errorbar(x,WTP8_Mean,WTP8_SE,'b');
ylim([0.6,1]);
xlim([0,10]);
saveas(gcf,[data_path 'WT_22.png']);
saveas(gcf,[data_path 'WT_22.fig']);
saveas(gcf,[data_path 'WT_22.eps']);

close;

figure;errorbar(x,B2P2_Mean,B2P2_SE,'r');
hold on;errorbar(x,B2P4_Mean,B2P4_SE,'g');
hold on;errorbar(x,B2P8_Mean,B2P8_SE,'b');
ylim([0.6,1]);
xlim([0,10]);
saveas(gcf,[data_path 'B2_22.fig']);
saveas(gcf,[data_path 'B2_22.png']);
saveas(gcf,[data_path 'B2_22.eps']);
close;
%%
WTP2_Mean = mean(dist_list_12_all_norm(1:3,:));
WTP2_SE = std(dist_list_12_all_norm(1:3,:)) / sqrt(3);
WTP4_Mean = mean(dist_list_12_all_norm(4:6,:));
WTP4_SE = std(dist_list_12_all_norm(4:6,:)) / sqrt(3);
WTP8_Mean = mean(dist_list_12_all_norm(7:9,:));
WTP8_SE = std(dist_list_12_all_norm(7:9,:)) / sqrt(3);
B2P2_Mean = mean(dist_list_12_all_norm(10:12,:));
B2P2_SE = std(dist_list_12_all_norm(10:12,:)) / sqrt(3);
B2P4_Mean = mean(dist_list_12_all_norm(13:15,:));
B2P4_SE = std(dist_list_12_all_norm(13:15,:)) / sqrt(3);
B2P8_Mean = mean(dist_list_12_all_norm(16:18,:));
B2P8_SE = std(dist_list_12_all_norm(16:18,:)) / sqrt(3);
%
x=1:1:20;
figure;errorbar(x,WTP2_Mean,WTP2_SE,'r');
hold on;errorbar(x,WTP4_Mean,WTP4_SE,'g');
hold on;errorbar(x,WTP8_Mean,WTP8_SE,'b');
ylim([0.6,1.2]);
xlim([0,10]);
saveas(gcf,[data_path 'WT_12.png']);
saveas(gcf,[data_path 'WT_12.fig']);
saveas(gcf,[data_path 'WT_12.eps']);

close;

figure;errorbar(x,B2P2_Mean,B2P2_SE,'r');
hold on;errorbar(x,B2P4_Mean,B2P4_SE,'g');
hold on;errorbar(x,B2P8_Mean,B2P8_SE,'b');
ylim([0.6,1.2]);
xlim([0,10]);
saveas(gcf,[data_path 'B2_12.fig']);
saveas(gcf,[data_path 'B2_12.png']);
saveas(gcf,[data_path 'B2_12.eps']);
close;
%%
writematrix(dist_list_22_all_norm,[data_path '22.csv'],'WriteMode','append');