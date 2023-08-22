%Meenakshi M Asokan's code for the neurograms in Figure 2 from Asokan et al., Cell Rep 2023
%August 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;
data_filename = fullfile(currentFolder,'/Data/Figure2/Habituation_FM_psth_allunits.mat');
data = load(data_filename);
psth_all = data.psth_all;
cell_type = cell(1,2);
cell_type{1} = 'HO-AC';
cell_type{2} = 'BLA';
%psth_all contains the smoothed (50 ms window) and z-scored (wrt. a 1s spontaneous window) psth for each unit (includes RS and FS) in these two regions across all 30 FM sweeps, both directions
%167 HO-AC units and 63 BLA units and (1s of baseline +10 s post sound onset)
%%
%Plot heatmaps of the Neurograms where each row is a unit and the heat map
%corresponds to the psth z score values
%Heat map is sorted based on average value of the psth

figure();
for c_t = 1:length(psth_all)%c_t - cell type
    subplot(1,2,c_t);
    psth = psth_all{c_t};
    %my function to sort the neurogram based on mean psth
    sorted_psth = myfunc_sort_based_on_mean(psth);
    colormap hot
    imagesc(sorted_psth);
    colorbar();
    caxis([0 4]);%heat map limits (in z scores)
    xlim([1 7000]);
    hold on
    vline([1000,2000,3000,4000,5000],'k--');
    hold on
    xlabel('Time (ms)');
    ylabel('Units');
    title(cell_type{c_t});
    set(gca,'fontsize',12);
end

%%
%FFT to obtain stimulus synchrony
%The stimulus is a 1Hz train for 5s - therefore the 1Hz power of the psth
%will be a proxy for stimulus synchrony
poi = cell(1,2);%power of interest - power at 1Hz
for c_t = 1:2
    poi_temp = [];
    for c = 1:size(psth_all{c_t},1)
        X = psth_all{c_t}(c,1000+(1:5000));%Only including the 5 s sound period for FFT
        L = length(X);
        Fs = 1000;
        Y = fft(X);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        % figure();
        % plot(f,P1) 
        % xlim([0 10])
        % title("Single-Sided Amplitude Spectrum of X(t)")
        % xlabel("f (Hz)")
        % ylabel("|P1(f)|")
        ind1 = find(f == 1);%power at 1 Hz
        poi_temp(c) = P1(ind1);
    end
    poi{c_t} = poi_temp;
end
%%
%violin plot to show fft amps
for c_t = 1:2
    in_idx = find(poi{1,c_t}<1);%limit for power for visualization
    data_poi{c_t} = poi{1,c_t}(in_idx);
end
%Only RS ACtx and BLA - half violin
%Half violin plot with box plot superimposed
figure();
for c_t = 1:2
subplot(1,2,c_t)
violinplot([data_poi{c_t}]',{'region'}, 'QuartileStyle','boxplot','HalfViolin','right',...
    'ShowData',false);
xlim([0.75 1.5]);
ylim([0 1]);
box off
xlabel(cell_type{c_t});
ylabel('Stimulus synchrony (power at 1Hz)');
set(gca,'fontsize',12);
end
%%
%Stats
stats_data = poi;
%%lillietest and ks test to test for normality of the data
ltest = [];
for i = 1:length(stats_data)
    ltest(i) = lillietest(stats_data{i})%0 => data are normal; 1 - not normal
    ktest(i) = kstest((stats_data{i} - mean(stats_data{i}))./std(stats_data{i}))%0 => data are normal; 1 => not normal
end
%data - not normal based on both the tests
%wilcoxon ranksum test (non-param equiv of 2 sample t test)
values1 = stats_data{1};
values2 = stats_data{2};
[p,h,stats] = ranksum(values2,values1)
z = stats.zval;
nx = length(stats_data{1});
ny = length(stats_data{2});
u_stat = stats.ranksum - (nx*(nx+1)/2)
cliff_d = (2*u_stat/(nx*ny))-1