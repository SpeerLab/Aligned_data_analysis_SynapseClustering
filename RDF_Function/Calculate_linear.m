inpath = 'X:\Chenghang\4_Color\Complex_Syn\CTB_Specific\';
for i = 1:18
    disp(i);
    load([inpath 'Vwater_ss_' sprintf('%03d',i) '.mat']);
    load([inpath 'Vwater_sn_' sprintf('%03d',i) '.mat']);

    comp_ss = [];
    simp_ss = [];
    for i = 1:numel(statsVwater_ss)
        B_ID = sort(statsVwater_ss(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            comp_ss = cat(1,comp_ss,statsVwater_ss(i));
        else
            simp_ss = cat(1,simp_ss,statsVwater_ss(i));
        end
    end
    comp_sn = [];
    simp_sn = [];
    for i = 1:numel(statsVwater_sn)
        B_ID = sort(statsVwater_sn(i).B_ID);
        B_ID = B_ID(B_ID>0);
        if numel(B_ID) > 1
            comp_sn = cat(1,comp_sn,statsVwater_sn(i));
        else
            simp_sn = cat(1,simp_sn,statsVwater_sn(i));
        end
    end


    max_lin_ss = zeros(1,numel(comp_ss));
    for j =1:numel(comp_ss)
        PixList = comp_ss(j).PixelList;
        minpix = min(PixList);
        maxpix = max(PixList);
        box_size = maxpix - minpix + [1,1,1];
        box_size = box_size.*[15.5,15.5,70];
        box_size = sqrt(box_size(1)^2+box_size(2)^2+box_size(3)^2);
        %box_size = max(box_size);
        max_lin_ss(j) = box_size;
    end
    max_lin_sn = zeros(1,numel(comp_sn));
    for j =1:numel(comp_sn)
        PixList = comp_sn(j).PixelList;
        minpix = min(PixList);
        maxpix = max(PixList);
        box_size = maxpix - minpix + [1,1,1];
        box_size = box_size.*[15.5,15.5,70];
        box_size = sqrt(box_size(1)^2+box_size(2)^2+box_size(3)^2);
        %box_size = max(box_size);
        max_lin_sn(j) = box_size;
    end
    writematrix(max_lin_ss,[inpath 'max_lin_ss.csv'],'WriteMode','append');
    writematrix(max_lin_sn,[inpath 'max_lin_sn.csv'],'WriteMode','append');
end