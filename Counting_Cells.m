clearvars
close all
field_names = {'AuP','AuV','TeA'};
root = 'F:\Meurika histology\AAA_counting';

mouse_nums = [2,3,6,9,10];
slices_mice = cell(1,length(mouse_nums));
slices_mice{1} = 2;
slices_mice{2} = [1,2,3];
slices_mice{3} = 1;
slices_mice{4} = [2,3,4];
slices_mice{5} = [1,3];

%%
T2 = readtable(fullfile(root,'Amy_Aud_Areas.csv'));
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
close all
CT = cbrewer('qual','Dark2',3);
for slice = 1:length(neuron_counts)
figure();
for field = 1:3
    subplot(3,1,field);
    if (field ==1)
        dist = neuron_counts(slice).AuP_dist;
    elseif (field ==2)
        dist = neuron_counts(slice).AuV_dist;
    else
        dist = neuron_counts(slice).TeA_dist;
    end

% histfit(dist,50,'kernel')
h1 = histogram(dist,length(dist));
% h1.DisplayStyle = 'stairs';
h1.FaceColor = CT(field,:);
h1.Normalization = 'probability';
BW = 15;
h1.BinWidth = BW;  
xlim([0 1200]);
hold on
[values, edges] = histcounts(dist, 'Normalization', 'probability','BinWidth',BW);
centers = (edges(1:end-1)+edges(2:end))/2;
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
% histfit(dist,50,'kernel')
h1 = histogram(dist_to_pia{field},length(dist_to_pia{field}));
% h1.DisplayStyle = 'stairs';
h1.FaceColor = CT(field,:);
h1.Normalization = 'probability';
BW = 10;
h1.BinWidth = BW;  
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
%centers of prob distributions - mean across mice
% close all
figure();
hist_values = cell(1,3);
for slice = 1:length(neuron_counts)
% figure();
for field = 1:3
%     subplot(3,1,field);
    if (field ==1)
        dist = neuron_counts(slice).AuP_dist;
    elseif (field ==2)
        dist = neuron_counts(slice).AuV_dist;
    else
        dist = neuron_counts(slice).TeA_dist;
    end
    BW = 30;
    edges = 0:BW:1200;
    values = histcounts(dist,'BinEdges',edges, 'Normalization', 'probability');
    centers = (edges(1:end-1)+edges(2:end))/2;
%     plot(centers, values, 'k-');
    hist_values{field}(slice,:) = values;
%     xlim([0 1200]);
%     hold on
end
end

sm = 5;
for field = 1:3
    subplot(3,1,field);
    data = smoothdata(hist_values{field},2,'gaussian',sm);
    for i = 1:size(data,1)
        plot(centers, data(i,:),'color',CT(field,:),'LineWidth',0.1);
        hold on
    end
    line_sem_plot(centers,data,CT(field,:),1.5,0.5);
    xlim([0 1200]);
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
figure();
sm = 5;
for field = 1:3
    data = smoothdata(hist_values{field},2,'gaussian',sm);
    line_sem_plot(centers,data,CT(field,:),1.5,0.5);
    xlim([0 1200]);
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
%raincloud plots
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
%cell density
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
for i = 1:length(data)
a = i-0.1;
b = i+0.1;
rand_x = (b-a).*rand(1,length(data{i})) + a;
scatter(rand_x,data{i},30,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
hold on
end
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