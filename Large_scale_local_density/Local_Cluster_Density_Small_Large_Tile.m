%%
clear;clc
input_path = 'X:\Chenghang\4_Color_Continue\Density_martix\';
output_path = 'X:\Chenghang\4_Color_Continue\';

high_thre_list = zeros(18,1);
low_thre_list = zeros(18,1);
Ratio_large_list = zeros(18,1);
Ratio_small_list = zeros(18,1);
for file_id = 1:18
    disp(file_id);
    load([input_path sprintf('%03d',file_id) '.mat']);
    figure;
    a = histogram(Local_density_1(:));a.NumBins = 20;a.FaceColor = 'r';a.FaceAlpha= 0.5;hold on;
    b = histogram(Local_density_1_new(:));b.BinWidth = a.BinWidth;b.BinLimits(1) = a.BinLimits(1);b.FaceColor = 'g';b.FaceAlpha= 0.5;
    fit_x = b.BinEdges;
    fit_y = b.Values;
    binsize = (fit_x(2) - fit_x(1))/2;
    fit_x = fit_x + binsize;
    fit_x(end) = [];
    f = fit(fit_x',fit_y','gauss1');
    plot(f,fit_x,fit_y)
    high_thre = f.b1 + 2*f.c1;
    low_thre = f.b1 - 2*f.c1;
    close;

    high_thre_list(file_id) = high_thre;
    low_thre_list(file_id) = low_thre;
    
    Density_size = numel(Local_density_1(:));
    Num_large = numel(find(Local_density_1(:) > high_thre));
    Num_small = numel(find(Local_density_1(:) < low_thre));
    Ratio_large = Num_large / Density_size;
    Ratio_small = Num_small / Density_size;

    Ratio_large_list(file_id) = Ratio_large;
    Ratio_small_list(file_id) = Ratio_small;

    Local_density_1_large_logic = Local_density_1 > high_thre;
    Local_density_1_small_logic = Local_density_1 < low_thre;
    save([output_path sprintf('%03d',file_id) '.mat'],'Local_density_1_large_logic','Local_density_1_small_logic');
end
save([output_path 'ratiolist.mat'],'high_thre_list','low_thre_list','Ratio_large_list','Ratio_small_list');

%%
clear;clc
input_path = 'X:\Chenghang\4_Color_Continue\Density_martix\';
output_path = 'X:\Chenghang\4_Color_Continue\';

high_thre_list = zeros(18,1);
low_thre_list = zeros(18,1);
Ratio_large_list = zeros(18,1);
Ratio_small_list = zeros(18,1);
for file_id = 1:18
    load([input_path sprintf('%03d',file_id) '.mat']);
    figure;
    a = histogram(Local_density_2(:));a.NumBins = 20;a.FaceColor = 'r';a.FaceAlpha= 0.5;hold on;
    b = histogram(Local_density_2_new(:));b.BinWidth = a.BinWidth;b.BinLimits(1) = a.BinLimits(1);b.FaceColor = 'g';b.FaceAlpha= 0.5;
    fit_x = b.BinEdges;
    fit_y = b.Values;
    binsize = (fit_x(2) - fit_x(1))/2;
    fit_x = fit_x + binsize;
    fit_x(end) = [];
    f = fit(fit_x',fit_y','gauss1');
    plot(f,fit_x,fit_y)
    high_thre = f.b1 + 2*f.c1;
    low_thre = f.b1 - 2*f.c1;
    close;

    high_thre_list(file_id) = high_thre;
    low_thre_list(file_id) = low_thre;
    
    Density_size = numel(Local_density_2(:));
    Num_large = numel(find(Local_density_2(:) > high_thre));
    Num_small = numel(find(Local_density_2(:) < low_thre));
    Ratio_large = Num_large / Density_size;
    Ratio_small = Num_small / Density_size;

    Ratio_large_list(file_id) = Ratio_large;
    Ratio_small_list(file_id) = Ratio_small;

    Local_density_2_large_logic = Local_density_2 > high_thre;
    Local_density_2_small_logic = Local_density_2 < low_thre;
    save([output_path sprintf('%03d',file_id) '.mat'],'Local_density_2_large_logic','Local_density_2_small_logic');
end
save([output_path 'ratiolist.mat'],'high_thre_list','low_thre_list','Ratio_large_list','Ratio_small_list');