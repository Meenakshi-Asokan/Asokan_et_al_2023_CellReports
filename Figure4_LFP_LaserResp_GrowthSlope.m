%Meenakshi M Asokan's code to optogenetically evoked LFP response growth slopes
%Related to Figure 4 from Asokan et al., Cell Rep 2023
%September 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;

data_filename = fullfile(currentFolder,'/Data/Figure4/lfp_laservaryV_allmice.mat');
data = load(data_filename);
lfp_laservaryV_allmice = data.lfp_laservaryV_allmice;
laser_powers = 0:5:20;%mW
CT{1}=cbrewer('seq', 'Greens', length(laser_powers));
CT{2}=cbrewer('seq', 'RdPu', length(laser_powers));
%%
%Plot laser response growth and comparison of growth slopes on hab vs recall 
num_mice = size(lfp_laservaryV_allmice,1);
baseline = 20;
interburst_duration = 240;
win = 80;

for m_n = 1:num_mice
    lfp_laservaryV = lfp_laservaryV_allmice(m_n,:);
    num_chans = size(lfp_laservaryV{1,1}{1,1},2);
    lfp_laservaryV_amps = cell(1,length(days));
    for day = 1:2
    for ch = 1:num_chans
    for ii = 1:length(laser_powers)
        lfp = lfp_laservaryV{day};
        dat = lfp{1,ii}(:,ch);
        %sign of LFP varies across depth of the probe after CMR,
        %essentially interested in the absolute value of the laser evoked deflection
        sign_check_data = lfp{1,length(laser_powers)}(:,ch);
        sc_delta_dat = (sign_check_data - mean(sign_check_data(1:baseline)));
        sc_data_window = sc_delta_dat(baseline+(1:win));
        sc_minima = min(sc_data_window);
        sc_maxima = max(sc_data_window);
        if abs(sc_maxima)>abs(sc_minima)
            dat = -dat;
        end
        delta_dat = (dat - mean(dat(1:baseline)));
        lfp_laservaryV_amps{day}(ii,ch) = abs(min(delta_dat(baseline+(1:win))));
    end
    end
    end
    %Plot the growth functions (mean across all channels)
    figure();
    for day = 1:2
    x = 1:length(laser_powers);
    y = lfp_laservaryV_amps{day}'*1e6;
    colorr = CT{day}(length(laser_powers)-2,:);
    lw = 1.5;
    fa = 0.5;
    line_sem_plot(x, y, colorr,lw,fa)
    hold on
    box off
    end
    
    xticks(1:1:length(laser_powers));
    xticklabels(laser_powers);
    xlabel('Laser power (mW)');
    ylabel('ChR2-LFP amplitude (\muV)');
    set(gca,'fontsize',12);
    set(gcf, 'Color', 'w');

    %scatter and box plots for the slopes     
    slope_dat = cell(1,2);
    auc_dat = cell(1,2);
    for day = 1:2
        dat = lfp_laservaryV_amps{day}'*1e6;
        for ch = 1:size(dat,1)
            slope_dat{day}(ch) = (dat(ch,size(dat,2)) - dat(ch,1))/(laser_powers(size(dat,2)) - laser_powers(1));
        end
    end
    figure();
    measure = slope_dat;
    box_colors = 'k';
    box_width = 0.5;
    myfunc_boxplot_customized(measure,box_colors,box_width);
    hold on
    for i = 1:length(measure)
    a = i-0.25;
    b = i+0.25;
    rand_x = (b-a).*rand(1,length(measure{i})) + a;
    scatter(rand_x,measure{i},20,[0.5 0.5 0.5],'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
    hold on
    end
    box off
    xticklabels({'Pre','Post'});
    set(gca,'fontsize',12);
end
%%
%Example of optogenetically evoked LFP responses
m_n = 1;%example mouse
lfp_laservaryV = lfp_laservaryV_allmice(m_n,:);
num_chans = size(lfp_laservaryV{1,1}{1,1},2);

for day = 2 %Hab - 1 vs Recall - 2; Note - acute recordings; same channel number does not mean anything for hab vs recall comparisons 
for ch = 16 %Channel number varies from 1 to 64
    figure();
for ii = 1:length(laser_powers)
    lfp = lfp_laservaryV{day};
    dat = lfp{1,ii}(:,ch);
    sign_check_data = lfp{1,length(laser_powers)}(:,ch);
    sc_delta_dat = (sign_check_data - mean(sign_check_data(1:baseline)));
    sc_data_window = sc_delta_dat(baseline+(1:win));
    sc_minima = min(sc_data_window);
    sc_maxima = max(sc_data_window);
    if abs(sc_maxima)>abs(sc_minima)
        dat = -dat;
    end
    delta_dat = (dat - mean(dat(1:baseline)));
        
    plot(1:interburst_duration,delta_dat*1e6,'color',CT{day}(ii,:),'linewidth',1.5);
    hold on
end
box off
xticks(20:40:220);
xticklabels(0:40:200);
vline(20,'k--');
% ylim([-100 30]);
legend('0 mW','5 mW', '10 mW', '15 mW','20 mW');
xlabel('Time (ms)');
ylabel('\Delta LFP (\muV)');
set(gca,'fontsize',12);
end
end
%%
%Stats
data_filename = fullfile(currentFolder,'/Data/Figure4/lfp_laser_resp_growth_allchs_data_allmice_habrecall.mat');
data = load(data_filename);
stats_data = data.data_allmice;

for r_num = 1:length(stats_data)
    values1 = stats_data{r_num};
    ktest(r_num) = kstest((values1 - mean(values1))./std(values1));
end
% %ktest 1 => data are not normal

%wilcoxon ranksum test (non-param equiv of 2 sample t test)
values1 = stats_data{1};
values2 = stats_data{2};
[p,h,stats] = ranksum(values2,values1);
z = stats.zval;
nx = length(stats_data{1});
ny = length(stats_data{2});
u_stat = stats.ranksum - (nx*(nx+1)/2)
cliff_d = (2*u_stat/(nx*ny))-1