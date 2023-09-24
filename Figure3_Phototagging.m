%Meenakshi M Asokan's code to identify photo-tagged cells
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

units_all = cell(1,length(cell_types));
for type = 1:length(cell_types)
    units_all{type} = cat(2,units_pre{type},units_post{type});
end

%Add to path the folder Utils including all subfolders
addpath(genpath(fullfile(currentFolder,'Utils')));
%%
%Plot the mean spike waveform shape of each unit within each cell type (CAmy, RS-HOAC, FS-HOAC,
%RS-LA, FS-LA)

wf = cell(1,length(units_all));
for type = 1:length(units_all)
    units = units_all{type};
    for clu = 1:length(units)
        wf{type}(clu,:) = units(clu).wf_shape;
    end   
end

figure();
clear CT
CT=cbrewer('qual', 'Set1', 5);
CT1 = flipud(CT(1:3,:));
subplot(2,2,1);
for type = [2,1,3]
    wfs = wf{type};
    for i = 1:size(wfs,1)
        plot(1:size(wfs,2),wfs(i,:),'color',CT1(type,:),'linewidth',0.5);
        hold on
    end
end
xlim([1,81]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xticks(1:10:81);
xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
box off
title('ACtx');
set(gca,'fontsize',12);
% xlabel('Time wrt. AP trough (ms)');
ylabel('AP amplitude (z)');

subplot(2,2,3);
for type = 1:3
    wfs = wf{type};
line_sem_plot(1:size(wfs,2),wfs,CT1(type,:),1.5,0.5);   
end
legend(sprintf('n = %d',size(wf{1,1},1)),'CAmy',sprintf('n = %d',size(wf{1,2},1)),'RS',sprintf('n = %d',size(wf{1,3},1)),'FS');
xlim([1,81]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xticks(1:10:81);
xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
box off
% title('ACtx');
set(gca,'fontsize',12);
xlabel('Time wrt. AP trough (ms)');
ylabel('AP amplitude (z)');

subplot(2,2,2);
for type = 4:5
    wfs = wf{type};
    for i = 1:size(wfs,1)
        plot(1:size(wfs,2),wfs(i,:),'color',CT(type,:),'linewidth',0.5);
        hold on
    end
%     line_sem_plot(1:size(wfs,2),wfs,CT(type-2,:),1.5,0.5);   
end
xlim([1,81]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xticks(1:10:81);
xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
box off
title('LA');
set(gca,'fontsize',12);
% xlabel('Time wrt. AP trough (ms)');
% ylabel('AP amplitude (z)');

subplot(2,2,4);
for type = 4:5
    wfs = wf{type};
line_sem_plot(1:size(wfs,2),wfs,CT(type,:),1.5,0.5);   
end
legend(sprintf('n = %d',size(wf{1,4},1)),'RS',sprintf('n = %d',size(wf{1,5},1)),'FS');
xlim([1,81]);
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xticks(1:10:81);
xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
box off
set(gca,'fontsize',12);
xlabel('Time wrt. AP trough (ms)');
% ylabel('AP amplitude (z)');
%%
%Compute FSL, FSJ and plot as a scatter
figure();
CT0=cbrewer('qual', 'Set1', 5);
CT = cat(1,flipud(CT0(1:3,:)),CT0(4:5,:));
for type = 1:length(units_all)
    units = units_all{type};
    resp_clus = find([units.laser_resp] == 1);


%Latency vs Jitter on a log scale
latency = [];
jitter = [];
spont_win = 100;
resp_win = 100;
spont_win_micros = spont_win*10;
resp_win_micros = resp_win*10;

for i = 1:length(resp_clus)
clu = resp_clus(i);
laser_resp_rasters = units(clu).laser_burstwise_rasters;
raster = laser_resp_rasters(:,spont_win_micros+(1:resp_win_micros));
    fsl = [];
    for ii = 1:size(raster,1)
        first_spike = find(raster(ii,1:500)==1);%50 ms response window for fsl
        if ~(isempty(first_spike))
            fsl = [fsl min(first_spike)];
        end
    end
    latency(i) = mean(fsl)/10;
    jitter(i) = std(fsl)/10;
    units_all{type}(clu).laser_fsl_fsj = [latency(i) jitter(i)];
end
scatter(jitter,latency,30,CT(type,:),'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0);
hold on
end
%Setting latency and jitter thresholds to define a unit as CAmy
latency_thresh = 5;
jit_thresh = 0.75;
set(gca,'xscale','log')
hold on
set(gca,'yscale','log')
hold on
xlim([0.1 20]);
hold on
ylim([1 30]);
hold on
hline(latency_thresh,'k--');
hold on
vline(jit_thresh,'k--');
% ylim([0 30]);
axis square
xlabel('FSL Jitter (ms)');
ylabel('FSL (ms)');
set(gcf,'color','w');
box off
set(gca,'fontsize',12);
%%
%Plot laser response rasters in ms resolution 
close all
for type = 1%CAmy cells (change the type to plot the laser responses for other cell types)
    units = units_all{type};
    for clu = 1:length(units)
        raster = units(clu).laser_burstwise_rasters;
        raster_compress = zeros(size(raster,1),size(raster,2)/10);%saved in 0.1 ms resolution, therefore need to compress
        for i = 1:size(raster,2)/10
            for j = 1:size(raster,1)
            raster_compress(j,i) = double(sum(raster(j,(i-1)*10+(1:10)))>0);
            end
        end
        figure();
        L = logical(raster_compress);
        MarkerFormat.Color = CT(type,:);
        plotSpikeRaster(L,'PlotType','scatter','MarkerFormat',MarkerFormat);
        box off
        xticks(0:20:200);
%         xticklabels(-spont_win:100:resp_win);
        xticklabels([-100:20:100]);
        xlabel('Time wrt laser onset (ms)');
        ylabel('Repetitions');
        title(sprintf('AAP %d, Sess %d, Clu %d, FSL = %0.2f ms, FSJ = %0.2f ms',units(clu).mouse_num,units(clu).sess_num,units(clu).cluster,units(clu).laser_fsl_fsj(1),units(clu).laser_fsl_fsj(2)));
        set(gca,'fontsize',12);
    end
end