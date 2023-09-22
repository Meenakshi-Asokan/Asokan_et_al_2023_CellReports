%Meenakshi M Asokan's code for CAmy cell counts and density in Figure 3 from Asokan et al., Cell Rep 2023
%August 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;
root = fullfile(currentFolder,'/Data/Figure3');
data_filename = fullfile(root,'Amy_Aud_Areas.csv');
T2 = readtable(data_filename);
%Add to path the folder Utils including all subfolders
addpath(genpath(fullfile(currentFolder,'Utils')));
%%
%Initialize variables
%5 mice, 10 slices overall
mouse_nums = [2,3,6,9,10];
slices_mice = cell(1,length(mouse_nums));
slices_mice{1} = 2;
slices_mice{2} = [1,2,3];
slices_mice{3} = 1;
slices_mice{4} = [2,3,4];
slices_mice{5} = [1,3];
field_names = {'AuP','AuV','TeA'};
%%
%1) Use ImageJ to mark the cells, add to the ROI manager, and save the coordinates
%2) Mark the pial surface as a line and run the macro - Macro_DistanceToPialSurface.ijm (function written by Michael Cammer, Microscopy core NYU Langone Medical Center)
%to compute the distance of the cells from the pial surface
%3) Save results for each slice, each cell's distance from pia - following
%code will read these results from each subfolder

ars = T2.Area;
neuron_counts = struct;
k = 1;
for mouse = 1:length(mouse_nums)
    mouse_num = mouse_nums(mouse);
    slices = slices_mice{mouse};
    dir_results = fullfile(root,sprintf('/AAA%02d/Results',mouse_num));
    dist_to_pia = cell(1,length(field_names));
    area_rois = cell(1,length(field_names));

    for slice = 1:length(slices)
        for field = 1:length(field_names)
            T = readtable(fullfile(dir_results,sprintf('AAA%02d_tilescan%d_%s_results.txt',mouse_num,slices(slice),field_names{field})));
            dist_to_pia{field} = T.Var4;%in microns
            x = intersect(intersect(find(T2.Mouse == mouse_num), find(T2.Slice == slices(slice))),find(strcmp(T2.Field,field_names{field})));
            area_rois{field} = (ars(x))/1e6;%mm square
        end
        neuron_counts(k).mouse_num = mouse_num;
        neuron_counts(k).slice_num = slices(slice);
        neuron_counts(k).AuP_dist = dist_to_pia{1};
        neuron_counts(k).AuP_area = area_rois{1};
        neuron_counts(k).AuV_dist = dist_to_pia{2};
        neuron_counts(k).AuV_area = area_rois{2};
        neuron_counts(k).TeA_dist = dist_to_pia{3};
        neuron_counts(k).TeA_area = area_rois{3};
        k = k+1;
%         slice
    end
    mouse
end
  
%%
%pool counts across all slices and mice and then plot histograms
dist_to_pia = cell(1,3);
for slice = 1:length(neuron_counts)
    dist_to_pia{1,1} = [dist_to_pia{1,1} neuron_counts(slice).AuP_dist'];
    dist_to_pia{1,2} = [dist_to_pia{1,2} neuron_counts(slice).AuV_dist'];
    dist_to_pia{1,3} = [dist_to_pia{1,3} neuron_counts(slice).TeA_dist'];
end

figure();
for field = 1:3
    subplot(3,1,field);
    data_histplot = dist_to_pia{field};
    colr = CT(field,:);
    bw = 10;
    myfunc_histogram(data_histplot,colr,bw);
    xlim([0 1200]);
    hold on
% [values, edges] = histcounts(dist_to_pia{field}, 'Normalization', 'probability','BinWidth',BW);
% centers = (edges(1:end-1)+edges(2:end))/2;
% plot(centers, values, 'k-')
% vline(0,'k--');
box off
if (field == 3)
    xlabel('Distance from pia (microns)');
end
if (field ==2)
    ylabel('Percentage of cells (/100)');
end
title(field_names{field});
set(gca,'TickDir','out');
set(gca,'fontsize',12);
end
%%
%raincloud plots for each slice separately
close all
for slice = 1:length(neuron_counts)
figure();
for field = 1:3
%     subplot(3,1,field);
    if (field ==1)
        dist = neuron_counts(slice).AuP_dist;
    elseif (field ==2)
        dist = neuron_counts(slice).AuV_dist;
    else
        dist = neuron_counts(slice).TeA_dist;
    end
h1 = raincloud_plot(dist, 'box_on', 0, 'color', CT(field,:), 'cloud_edge_col', CT(field,:), 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .5, 'dot_dodge_amount', (field)*.25,...
     'box_col_match', 0,'bandwidth', 0.5, 'density_type', 'ks');
h1{1}.EdgeAlpha = 0;
h1{2}.SizeData = 30;% handles 1-6 are the cloud area, scatterpoints, and boxplot elements respectively
h1{2}.MarkerFaceAlpha = 0.5;

view([90 -90]);
box off
end
legend('AuP','','AuV','','TeA','');
set(gca,'Xdir','reverse');
xlim([0 1200]);
xlabel('Depth from Pia (microns)');
set(gca,'fontsize',12);
end

%%
%cell density box plot
density = cell(1,3);
for slice = 1:length(neuron_counts)
    num_cells = length(neuron_counts(slice).AuP_dist);
    area = neuron_counts(slice).AuP_area;
    density{1} = [density{1} num_cells/area];
    
    num_cells = length(neuron_counts(slice).AuV_dist);
    area = neuron_counts(slice).AuV_area;
    density{2} = [density{2} num_cells/area];
    
    num_cells = length(neuron_counts(slice).TeA_dist);
    area = neuron_counts(slice).TeA_area;
    density{3} = [density{3} num_cells/area];
end

figure();
data = density;
scatter_width = 0.1;
marker_size = 30; marker_color = [0.5 0.5 0.5]; marker_alpha = 0.25;
myfunc_scatterplot(data,scatter_width,marker_size,marker_color,marker_alpha);
box_colors = 'k';
box_width = 0.2;
myfunc_boxplot_customized(data,box_colors,box_width);
hold on
hline(0,'k--');
% ylim([-1 1]);
% xlim([0 3]);
box off
ylabel('Density of amygdala projecting cortical neurons (Cells/mm^2)');
xticklabels(field_names);
% title('Post-Conditioning');
set(gca,'fontsize',12);
%%
%layer assignment?
% 0 - 100 - L1
% 100 - 300 - L2/3
% 300 - 400 - L4
% 400 - 800 - L5
% 800 - 1000 - L6a
% 1000 - 1200 - L6b + wm
%%
%STATS
%STEP1 : NORMALITY, HOMOSCEDASTICITY?
%%lillietest
ltest = [];
for i = 1:length(data)
    ltest(i) = lillietest(data{i})%0 => data are normal
    ktest(i) = kstest((data{i} - mean(data{i}))./std(data{i}))%0 => data are normal
end

%%
%data are normal!
%completely randomized (unbalanced) one way anova
y = [];
for i = 1:length(data)
    y = [y data{i}];
end
group = cell(1,length(y));
g = cell(1,3);
g{1} = 'AuP';g{2} = 'AuV';g{3} = 'TeA';
k = 1;
for i = 1:length(data)
    for j = 1:length(data{i})
        group{k} = g{i};
        k = k+1;
    end
end
[p,tbl,stats] = anova1(y,group);
%%
%post-hoc pairwise comparisons with bonferroni correction
tbl = multcompare(stats,'CType','bonferroni')
pvalues_bonf = tbl(:,6);

%post-hoc pairwise comparisons with bonf-holm correction
pvalues = [];
coh_d = [];
k = 1;
for i = 1:length(data)-1
    for j = i+1:length(data)
        [h,pvalues(k)] = ttest2(data{i},data{j});
        coh_d(k,1) = computeCohen_d(data{i},data{j}, 'independent');
        k = k+1;
        
    end
end
[corrected_p, h]=bonf_holm(pvalues);
pvalues_bonfholm = corrected_p';