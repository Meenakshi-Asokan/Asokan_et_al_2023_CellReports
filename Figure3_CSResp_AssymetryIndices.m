%Meenakshi M Asokan's code to compute CS-evoked responses and CS+ vs CS-
%assymetry indices
%Related to Figure 3 from Asokan et al., Cell Rep 2023
%September 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;
data_filename = fullfile(currentFolder,'/Data/Figure3/cell_types_names.mat');
data = load(data_filename);
cell_types = data.cell_type;

data_filename = fullfile(currentFolder,'/Data/Figure3/units_all_types_pre.mat');
data = load(data_filename);
units_pre = data.pre_units_type;
data_filename = fullfile(currentFolder,'/Data/Figure3/units_all_types_post.mat');
data = load(data_filename);
units_post = data.post_units_type;

%Add to path the folder Utils including all subfolders
addpath(genpath(fullfile(currentFolder,'Utils')));
%%
%Compute CS resp as AUC and then compute assymetry indices

mouse_numbers = [9,10,11,12,13,14,15,17];

cs_plus = [1,2,1,2,2,1,2,1];
cs_minus = [2,1,2,1,1,2,1,2];
spont_win = 1000;
figure();
data = cell(length(cell_types),2);

for type = 1:length(cell_types)
% data = cell(1,2);
for prepost = 1:2
    if prepost ==1
    units = units_pre{type};
    else
    units = units_post{type};        
    end
trialwise_spike_counts = [];
z_psth = [];
fr_psth = [];
for clu = 1:length(units)
    m_num = find(mouse_numbers == units(clu).mouse_num);
    cs_p = cs_plus(m_num);
    cs_m = cs_minus(m_num);
    cs = [cs_p cs_m];
    for sound = 1:2  
        raster = units(clu).sound_raster{1,cs(sound)};
        psth = mean(raster,1)*1e3;
        mean_spont = mean(psth(1:spont_win));
        std_spont = std(psth(1:spont_win));
        z_psth(clu,sound,:) = (psth - mean_spont)/std_spont;
    end
    clu
end

num_reps = 5;
spont_win = 1000;
sound_win = 5000;
sm = 50;
area_clus = [];
max_resp = [];
c_amy = [];
for clu = 1:length(units)
    for sound = 1:2  
    z = squeeze(z_psth(clu,sound,:));
    z_psth_sm = smoothdata(z,'gaussian',sm);
    temp_sound = [];
    for rep = 1:num_reps
        temp_sound(rep,:) = z_psth_sm(1000*rep-250+(1:1000));
    end
    burst_resp_mean = mean(temp_sound,1);
    %Only the positive area under the curve
    pos_data = z_psth_sm;
    neg = find(pos_data<0);
    pos_data(neg) =0;
    cum_area = cumtrapz(pos_data);
    %AUC across all bursts
   area_clus(clu,sound) = cum_area(spont_win+sound_win) - cum_area(spont_win);   
    end
    clu
end

%assym indices
cs_area_assym = [];
for clu = 1:length(units)
cs_area_assym(clu) = (area_clus(clu,1) - area_clus(clu,2))/(area_clus(clu,1) + area_clus(clu,2));
end

data{type,prepost} = cs_area_assym;

subplot(1,length(cell_types),type);
%Plot data points as scatter
a = prepost-0.25;
b = prepost+0.25;
rand_x = (b-a).*rand(1,length(data{type,prepost})) + a;
scatter(rand_x,data{type,prepost},20,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on
end
%box plots
box_colors = [0.5 0.5 0.5];
box_width = 0.5;
myfunc_boxplot_customized(data(type,:),box_colors,box_width)
hold on
hline(0,'k--');
ylim([-1 1]);
box off
xticklabels({'Hab','Recall'});
ylabel('CS+ vs CS- response assymetry index');
title(cell_types{type});
set(gca,'fontsize',12);
end
%%
%Visualize sound rasters, laser rasters, waveform shapes together

units = units_post{1};%CAmy units post - change pre vs post and cell type to visualize other units

close all
clear CT
CT{1}=cbrewer('seq', 'RdPu', 7);
CT{2}=cbrewer('seq', 'Greens', 7);
spont_win = 1000;

num_reps = 5;
spont_win_micros = 100*10;

%First plot rasters & PSTH for CS-evoked activity (1s baseline, 5s CS, 1s post)
for clu = 1:length(units)
    m_num = find(mouse_numbers == units(clu).mouse_num);
    cs_p = cs_plus(m_num);
    cs_m = cs_minus(m_num);
    cs = [cs_p cs_m];
    figure();
    for sound = 1:2  
        raster = units(clu).sound_raster{1,cs(sound)};
        psth = mean(raster,1)*1e3;
        sm = 100;
        psth_sm = smoothdata(psth,'gaussian',sm);
        mean_spont = mean(psth(1:spont_win));
        std_spont = std(psth(1:spont_win));
        z_psth = (psth - mean_spont)/std_spont;
        z_psth_sm = smoothdata(z_psth,'gaussian',sm);
        subplot(3,3,sound*3-2);
        L = logical(raster);
        LineFormat = struct();
        LineFormat.Color = CT{sound}(5,:);
        LineFormat.LineWidth = 1;
        LineFormat.LineStyle = '-';   
        plotSpikeRaster(L,'PlotType','vertline','LineFormat',LineFormat);
%         MarkerFormat.Color = CT{sound}(5,:);
%         plotSpikeRaster(L,'PlotType','scatter','MarkerFormat',MarkerFormat);
        xticks(0:1000:11000);
        xticklabels(-1:1:10);
        xlim([0 7000]);
        box off 
        y_lowbound = 1-0.5;
        y_upbound = size(raster,1)+0.5;
        ylim([y_lowbound y_upbound]);
        for ii = 1:num_reps
            x_on = (0+ii)*1000;
            x_off = (0.50+ii)*1000;
            x_fill = [x_on, x_on, x_off, x_off];
            y_fill = [y_lowbound, y_upbound, y_upbound, y_lowbound];
            patch(x_fill, y_fill, 'c','FaceAlpha',.1,'EdgeAlpha',0);
        end
        if (sound ==1)
        title(sprintf('mouse num = %d, cluster = %d',units(clu).mouse_num, units(clu).cluster));
        end
        set(gca,'fontsize',12);

        subplot(3,3,7);
        hold on
%         plot(1:length(z_psth_sm),z_psth_sm,colr{sound});
        plot(1:length(psth_sm),psth_sm,'color',CT{sound}(5,:),'linewidth',1);
        hold on
        xticks(0:1000:11000);
        xticklabels(-1:1:10);
            xlim([0 7000]);
        xlabel('Time wrt. sound onset (s)');
        ylabel('Firing rate (sp/s)');
        set(gca,'fontsize',12);
        
        %Zoom in around the first burst (0.5 s before and after)
        subplot(3,3,1+sound*3-2);
        L = logical(raster);
        LineFormat = struct();
        LineFormat.Color = CT{sound}(5,:);
        LineFormat.LineWidth = 1;
        LineFormat.LineStyle = '-';
        plotSpikeRaster(L,'PlotType','vertline','LineFormat',LineFormat);
        xticks(0:100:11000);
        xticklabels(-1000:100:10000);
        xlim([500 1500]);
        box off 
        y_lowbound = 1-0.5;
        y_upbound = size(raster,1)+0.5;
        ylim([y_lowbound y_upbound]);
        for ii = 1:num_reps
            x_on = (0+ii)*1000;
            x_off = (0.50+ii)*1000;
            x_fill = [x_on, x_on, x_off, x_off];
            y_fill = [y_lowbound, y_upbound, y_upbound, y_lowbound];
            patch(x_fill, y_fill, 'c','FaceAlpha',.1,'EdgeAlpha',0);
        end
        set(gca,'fontsize',12);
        subplot(3,3,8);
        hold on
%         plot(1:length(z_psth_sm),z_psth_sm,colr{sound});
        plot(1:length(psth_sm),psth_sm,'color',CT{sound}(5,:),'linewidth',1);
        hold on
        hold on
        xticks(0:100:11000);
        xticklabels(-1000:100:10000);
        xlim([500 1500]);
        xlabel('Time wrt. sound onset (ms)');
        ylabel('Firing rate (sp/s)');
        set(gca,'fontsize',12);
    end
    
    %Laser response, fsl, fsj
    subplot(3,3,3);
    raster = units(clu).laser_burstwise_rasters;
    fsl = [];
    for ii = 1:size(raster,1)
        first_spike = find(raster(ii,spont_win_micros+(1:500))==1);%50 ms response window after laser onset
        if ~(isempty(first_spike))
            fsl = [fsl min(first_spike)];
        end
    end
    latency = mean(fsl)/10;
    jitter = std(fsl)/10;       
    L = logical(raster);
    MarkerFormat.Color = 'r';
    plotSpikeRaster(L,'PlotType','scatter','MarkerFormat',MarkerFormat);
    box off
    xticks(0:500:2000);
    xticklabels([-1000:500:1000]);
    xlabel('Time wrt laser onset (0.1 ms)');
    ylabel('Repetitions');
    title(sprintf('FSL = %0.2f ms, FSJ = %0.2f ms',latency, jitter));
    set(gca,'fontsize',12);
    psth = mean(raster,1)*1e4;
    sm = 5;
    psth_sm = smoothdata(psth,'gaussian',sm);
    subplot(3,3,6);
    plot(1:length(psth_sm),psth_sm,'r');
    box off
    xlabel('Time wrt. laser onset (0.1 ms)');
    ylabel('Firing rate (sp/s)');
    set(gca,'fontsize',12);
    
    %Spike waveform shape
    subplot(3,3,9);
    wf_shape = units(clu).wf_shape;
    plot(1:length(wf_shape),wf_shape,'color','k','linewidth',1.5);
    xlim([1,81]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    xticks(1:10:81);
    xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
    box off
    title(sprintf('Peak - trough interval = %0.2f ms',units(clu).peak_trough_interval));
    set(gca,'fontsize',12);
    xlabel('Time wrt. AP trough (ms)');
    ylabel('AP amplitude (z)');
end
set(gcf,'color','w');
%%
%Stats
data_filename = fullfile(currentFolder,'/Data/Figure3/su_cs_spiking_data_areas.mat');%DTC group
% data_filename = fullfile(currentFolder,'/Data/Figure3/su_cs_spiking_data_pseudocond_areas.mat');%pseudoconditioning group
data_temp = load(data_filename);
clear data
data = data_temp.data;

%Checking for normality
ltest = [];
ktest = [];
mean_val = [];
stderr = [];
for r_num = 1:length(cell_types)
stats_data = data(r_num,:);
for prepost = 1:2
values = stats_data{prepost};
mean_val(prepost,r_num) = mean(values);
stderr(prepost,r_num) = std(values)/sqrt(length(values));
ktest(r_num,prepost) = kstest((values - mean(values))./std(values));
end
end
%k = 0 =>normal as per kstest

%2 sample t-test (unpaired)
h = [];
p = [];
t = [];
d = [];
stats = cell(1,length(cell_types));
for r_num = 1:length(cell_types)
stats_data = data(r_num,:);
%2 sample t-test
values1 = stats_data{1};
values2 = stats_data{2};
[h(r_num),p(r_num),ci,stats{r_num}] = ttest2(values2,values1);
t(r_num) = stats{r_num}.tstat;
d(r_num) = computeCohen_d(values2, values1, 'independent');
end
h
p
d