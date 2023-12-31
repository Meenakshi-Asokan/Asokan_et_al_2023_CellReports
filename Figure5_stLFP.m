%Meenakshi M Asokan's code to plot stLFP 
%Related to Figure 5 from Asokan et al., Cell Rep 2023
%September 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;

data_filename = fullfile(currentFolder,'/Data/Figure5/units_stlfp_pre.mat');
data = load(data_filename);
units_stlfp_pre = data.units_stlfp_pre;

data_filename = fullfile(currentFolder,'/Data/Figure5/units_stlfp_post.mat');
data = load(data_filename);
units_stlfp_post = data.units_stlfp_post;

data_filename = fullfile(currentFolder,'/Data/Figure3/cell_types_names.mat');
data = load(data_filename);
cell_types = data.cell_type;

CT_sound = cell(1,2);
CT_sound{1}=cbrewer('seq', 'RdPu', 7);
CT_sound{2}=cbrewer('seq', 'Greens', 7);
%%
close all
% epoch = 1;%Pre sound epoch
epoch = 2;%During sound epoch
% epoch = 3;%Post sound epoch
win = 50;
for day = 1:2
    if day ==1
        units_type = units_stlfp_pre;
    else
        units_type = units_stlfp_post;
    end
    %merging together all RS HO-AC units into type2
    units_type{2} = cat(2,units_type{2},units_type{1});

%for each unit stLFP has been computed using linear deconvolution for lfp from each channel in the other region and then averaged across all channels)
for type = 2 %2 - RS HO-AC units and LA lfp; 4 - RS LA units and HO-AC lfp
mat_stLFP = cell(1,2);
k = 1;
for clu = 1:length(units_type{type})
    stlfp_type_main = units_type{type}(clu).sig_stLFP_4{epoch};
    stlfp_type = stlfp_type_main;
    if ~isempty(stlfp_type)
        for sound = 1:2
            sta = stlfp_type(sound,:);
            mat_stLFP{sound}(k,:) = sta;
        end
        k = k+1;
    end
end
%%display after sorting - template
%sorting based on latency
sorted_mat_stLFP = cell(1,2);
where_min = [];
for i = 1:size(mat_stLFP{1},1)
    where_min(i) = find(mat_stLFP{1}(i,:) == min(mat_stLFP{1}(i,1000+(-win:win))));
end
[out,idx] = sort(where_min);

%sort
for sound = 1:2
    for i = 1:size(mat_stLFP{sound},1)
        sorted_mat_stLFP{sound}(i,:) = mat_stLFP{sound}(idx(i),:);
    end
end
figure();
for sound = 1:2
subplot(1,3,sound);
CT=cbrewer('div','RdBu',20,'pchip');
colormap(flipud(CT));
imagesc(sorted_mat_stLFP{sound});
hold on
vline(1000,'k');
caxis([-4 4]);
box off
xticks(0:200:2000);
xticklabels(-1000:200:1000);
xlabel('Time wrt spike(ms)');
ylabel('Units');
set(gcf, 'Color', 'w');
set(gca,'fontsize',12);
xlim([500 1500]);
end
subplot(1,3,3);
for sound = 1:2
x = 1:size(sorted_mat_stLFP{sound},2);
y = sorted_mat_stLFP{sound};
color = CT_sound{sound}(5,:);
lw = 1.5;
fa = 0.2;
line_sem_plot(x, y, color,lw,fa);
hold on
vline(1000,'k');
xticks(0:200:2000);
xticklabels(-1000:200:1000);
xlim([500 1500]);
ylim([-3.5 1.5]);
xlabel('Time wrt spike(ms)');
ylabel('stLFP');
set(gca, 'TickDir', 'in')
set(gcf, 'Color', 'w');
set(gca,'fontsize',12);
end
end

end
%%

%assymetry indices
close all
win = 50;
data = cell(2,length(cell_types));%days,unit_types
sta_amps = cell(2,2,length(cell_types));%days,sounds,units_type

figure();
for day = 1:2
    if day ==1
        units_type = units_stlfp_pre;
    else
        units_type = units_stlfp_post;
    end
    %merging together all RS HO-AC units into type2
    units_type{2} = cat(2,units_type{2},units_type{1});

type = 2;%HO-AC spikes to LA LFP
% type = 4;%LA spikes to HO-AC LFP
    sta_amp = [];
    k = 1;
    for clu = 1:length(units_type{type})
        stlfp_type_main = units_type{type}(clu).sig_stLFP_4{epoch};%for epoch 2, doesn't quite matter whetehr _,2,3,4

        stlfp_type = stlfp_type_main;
        if ~isempty(stlfp_type)
        for sound = 1:2
            sta_amp(k,sound) = max(abs(stlfp_type(sound,1000+(-win:win))));
        end
        k = k+1;
        end

    end
measure = sta_amp;

for sound = 1:2
sta_amps{day,sound,type} = measure(:,sound)';
end
data{day,type} = ((measure(:,1) - measure(:,2))./(measure(:,1) + measure(:,2)))';

a = day-0.25;
b = day+0.25;
rand_x = (b-a).*rand(1,length(data{day,type})) + a;
scatter(rand_x,data{day,type},30,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
hold on
end
box_colors = 'k';
box_width = 0.7;
myfunc_boxplot_customized(data(:,type),box_colors,box_width)
hold on
hline(0,'k--');
box off

xticklabels({'Hab','Recall'});
ylabel('CS+ vs CS- stLFP assymetry index');
set(gca,'fontsize',12);
hold on


%cs+ cs- comparisons of amplitudes
figure();
for day = 1:2%1 - Hab; 2 - Recall
    subplot(1,2,day);
    for sound = 1:2
        a = sound-0.25;
        b = sound+0.25;
        rand_x = (b-a).*rand(1,length(sta_amps{day,sound,type})) + a;
        scatter(rand_x,sta_amps{day,sound,type},30,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
        hold on
    end
    box_colors = 'k';
    box_width = 0.7;
    myfunc_boxplot_customized(sta_amps(day,:,type),box_colors,box_width)
    hold on
    box off
    ylim([0 7]);
    xticklabels({'CS+','CS-'});
    ylabel('stLFP amplitude');
    if day ==1
        title('Hab');
    else
        title('Recall');
    end
    set(gca,'fontsize',12);
end
%%
%Stats
data_filename = fullfile(currentFolder,'/Data/Figure5/sta_amps.mat');
% data_filename = fullfile(currentFolder,'/Data/Figure5/sta_amps_pseudocond.mat');
data_temp = load(data_filename);
sta_amps = data_temp.sta_amps;

type = 2;
data = cell(1,2);
for day = 1:2
    for sound = 1:2
        data{day}(sound,:) = sta_amps{day,sound,type};
    end
end
%%
%2 x 2 mixed model anova: within subject factor/repeated measure - sound,
%between subject factor - session (because units in diff session are diff)
data_mat = cat(2,data{1},data{2})';
between_factor = cat(1,ones(size(data{1},2),1),2*ones(size(data{2},2),1));
[tbl,rm] = simple_mixed_anova(data_mat,between_factor,{'Sound'},{'Session'})

%%
%assym indices
data = cell(1,2);

figure();
for day = 1:2
    measure = [];
    for sound = 1:2
    measure(:,sound) = sta_amps{day,sound,type};
    end
    data{1,day} = ((measure(:,1) - measure(:,2))./(measure(:,1) + measure(:,2)))';
    % figure();
    % for i = 1:size(data,2)
    a = day-0.25;
    b = day+0.25;
    rand_x = (b-a).*rand(1,length(data{1,day})) + a;
    scatter(rand_x,data{1,day},30,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
    hold on
end

box_colors = 'k';
box_width = 0.7;
myfunc_boxplot_customized(data,box_colors,box_width)
hold on
hline(0,'k--');
ylim([-0.85 0.85]);
box off
xticklabels({'Hab','Recall'});
ylabel('CS+ vs CS- stLFP assymetry index');
set(gca,'fontsize',12);

%%
for r_num = 1:length(data)
    values1 = data{r_num};
    ktest(r_num) = kstest((values1 - mean(values1))./std(values1));
end
%ktest 0 => data are normal

%one sample t test against population mean of 0
stats_mod = cell(1,length(data));
h = [];
p = [];
t = [];
d = [];
for r_num = 1:length(data)
    [h(r_num),p(r_num),ci,stats_mod{r_num}] = ttest(data{r_num},0);
    t(r_num) = stats_mod{r_num}.tstat;
    d(r_num) = computeCohen_d(data{r_num},0);
end
h
p
d

%unpaired 2 sample ttest
stats_mod = cell(1,length(data));
h = [];p = [];d = [];
[h,p,ci,stats] = ttest2(data{2},data{1})
p
t = stats.tstat
d = computeCohen_d(data{2},data{1})
