%Meenakshi M Asokan's code for Figure 1 from Asokan et al., Cell Rep 2023
%August 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;
data_filename = fullfile(currentFolder,'/Data/Figure1/behav_meas_allmice_3sessProtocol_22fm20pd.mat');
data = load(data_filename);
behav_meas_allmice = data.behav_meas_allmice_22fm20pd;
fs = 30;%sampling rate is 30 frames per second 
%Data - 2*2*3 cell
%Dimension 1 - Sound - CS+, CS-
%Dimension 2 - Measure - Facial motion, Pupil dilation
%Dimension 3 - Session - Habitutaion, Conditioning and Recall 
%Each trial is 22 seconds (2 s baseline + 20 s after CS onset), hence 660 frames
%Facial motion data from 22 mice and Pupil dilation data from 20 mice in total
%%
%Plot pupil dilations and facial motion for an example mouse
close all
ex_mouse = 1;
meas = 1;%Measure 1 - facial motion
figure();
t = 1:fs*22;
CT_sound = cell(1,2);
CT_sound{1}=cbrewer('seq', 'RdPu', 7);
CT_sound{2}=cbrewer('seq', 'Greens', 7);
for sess = 1:3
    for sound = 1:2
        c = CT_sound{sound}(5,:);
        pupil_dilations = behav_meas_allmice{sound,meas,sess};
        for mouse = ex_mouse
            subplot(3,1,sess);
            plot(t,pupil_dilations(mouse,:),'color',c,'linewidth',1.5);
            hold on  
            xlabel('Time (s)');
            ylabel('Facial Motion');
            xlim([1 12*fs]);
            xticks(0:2*fs:12*fs);
            xticklabels(-2:2:10);
            box off
            set(gca,'fontsize',12);
        end
    end
    if (sess == 1)
        legend({'CS+','CS-'});
    end
end


meas = 2;%Measure 2 - pupil dilation
figure();
t = 1:fs*22;
CT_sound = cell(1,2);
CT_sound{1}=cbrewer('seq', 'RdPu', 7);
CT_sound{2}=cbrewer('seq', 'Greens', 7);
for sess = 1:3
    for sound = 1:2
        c = CT_sound{sound}(5,:);
        facial_motions = behav_meas_allmice{sound,meas,sess};
        for mouse = ex_mouse
            subplot(3,1,sess);
            plot(t,facial_motions(mouse,:),'color',c,'linewidth',1.5);
            hold on  
            xlabel('Time (s)');
            ylabel('Pupil dilation');
            xlim([1 12*fs]);
            xticks(0:2*fs:12*fs);
            xticklabels(-2:2:10);
            box off
            set(gca,'fontsize',12);
        end
    end
    if (sess == 1)
        legend({'CS+','CS-'});
    end
end

%%
%Figure 1D - Scatter and box plots for the summary across mice
t_baseline = 2;
t_preshock = t_baseline+4;%Each sweep is one second long and shock occurs at the onset of the 5th sweep
figure();
num_sess = 3;
stats_data = cell(1,2);%2 measures - fm and pd
for meas = 1:2
    %For normalization, first compute Habituation AOI (Area Of Interest - area under the curve during the first 4 secs after CS onset b4 shock
    temp_aoi = cell(2,2);
    for hab_sess =1
    for sound = 1:2
        temp = behav_meas_allmice(sound,meas,hab_sess);
        data = temp{1};
        for m_n = 1:size(data,1)
            cum_area = cumtrapz(data(m_n,:));
            aoi = cum_area(round(t_preshock*fs)) - cum_area(round(t_baseline*fs));    
            temp_aoi{hab_sess,sound} = [temp_aoi{hab_sess,sound} aoi];
        end
    end
    end
    hab_aoi = temp_aoi(1,:);%only sess 1 as hab (even though multiple hab sessions for some mice)
    
    for sess = 1:num_sess
        subplot(2,num_sess,(meas-1)*num_sess+sess);
        %Area under the curve for each session and each CS
        cs_aoi = cell(1,2);%2sounds
        cs_aoi_diff = cell(1,2);
        for sound = 1:2
            temp = behav_meas_allmice(sound,meas,sess);
            data = temp{1};
            for m_n = 1:size(data,1)
                cum_area = cumtrapz(data(m_n,:));
                aoi = cum_area(round(t_preshock*fs)) - cum_area(round(t_baseline*fs));    
                cs_aoi{sound} = [cs_aoi{sound} aoi];
            end
            cs_aoi_diff{sound} = (cs_aoi{sound}-hab_aoi{sound});%Delta wrt. hab (there is only one hab)
        end
        bp_data = cs_aoi_diff;%box plot data
        stats_data{meas}(:,sess,1) = bp_data{1};
        stats_data{meas}(:,sess,2) = bp_data{2};
        box_colors = 'k';
        box_width = 0.5;
        %Box plot
        myfunc_boxplot_customized(bp_data,box_colors,box_width);%My function for a customized box plot (with mean display as well)
        hold on
        %Scatter plot on top
        scatter_width = 0.15;marker_size = 30;marker_color = [0.5 0.5 0.5];marker_alpha = 0.25;
        myfunc_scatterplot(bp_data,scatter_width,marker_size,marker_color,marker_alpha);
        %Line plot on top
        for m = 1:length(bp_data{1})
            plot(1:2,[bp_data{1}(m),bp_data{2}(m)],'color',[0.5 0.5 0.5]);
            hold on
        end
        hline(0,'k--');
        xticklabels({'CS+','CS-'});
        if (meas ==1)
            ylabel('Facial Motion');
            ylim([-110 75]);
        else
            ylabel('Pupil Dilation');
            ylim([-10 25]);
        end
        box off
        set(gca,'fontsize',12);
    end
end
%%
%Figure 1E-F
%discriminative and non-discriminative learning
for meas = 1:2%1-fm; 2-pd
    temp_data = stats_data{meas};
    non_discrim_data{meas} = mean(temp_data,3);
    discrim_data{meas} = temp_data(:,:,1) - temp_data(:,:,2);
end
%box-plot and scatter
close all
for meas = 1:2
figure();
for d = 1:2%1 - non-discriminatory; 2 - discriminatory changes
    bp_data = cell(1,2);
    subplot(1,2,d);
    if d ==1
    temp = non_discrim_data{meas};
    else
    temp = discrim_data{meas};
    end
    for sess = 2:3
        bp_data{sess-1} = temp(:,sess)';
    end
    box_colors = 'k';
    box_width = 0.5;
    %box plot with mean overlay
    myfunc_boxplot_customized(bp_data,box_colors,box_width);
    hold on
    %scatter plot on top
    scatter_width = 0.15;marker_size = 30;marker_color = [0.5 0.5 0.5];marker_alpha = 0.25;
    myfunc_scatterplot(bp_data,scatter_width,marker_size,marker_color,marker_alpha);
    hline(0,'k--');
    xticklabels({'CS+','CS-'});
    if (meas ==1)
        title('Facial Motion');
        ylim([-80 60]);
    else
        title('Pupil Dilation');
        ylim([-10 25]);
    end
    if d ==1
        ylabel('Non-discriminative changes');
    else
        ylabel('Discriminatory learning');
    end
    box off
    set(gca,'fontsize',12);
end
end

%%
%Stats
%Discriminatory learning
for meas = 1:2%1-fm; 2-pd
temp = discrim_data{meas};
y = [];
i = 1;
for sess = [2,3]
    y(:,i) = temp(:,sess)';
    i = i+1;
end
data_mat = y;
%one way rm anova (session - repeated measure)
[tbl,rm] = simple_mixed_anova(data_mat,[],{'Session'})
%posthoc one sample t test for each session
h = [];
p = [];
t = [];
d = [];
for sess = 1:2
    values = data_mat(:,sess);
    [h(sess),p(sess),ci,stats] = ttest(values,0);
    t(sess) = stats.tstat;
    d(sess) = computeCohen_d(values, 0);%effect size
end
p
h
[cor_p, h]=bonf_holm(p,.05)%correction for multiple comparisons
end
%%
%Non-discriminative changes
for meas = 1:2%meas1 - fm; meas2 - pd;
temp = non_discrim_data{meas};
y = [];
i = 1;
for sess = [2,3]
    y(:,i) = temp(:,sess)';
    i = i+1;
end
data_mat = y;
%one way rm anova (session - repeated measure)
[tbl,rm] = simple_mixed_anova(data_mat,[],{'Session'})
%posthoc one sample t test for each session
h = [];
p = [];
t = [];
d = [];
for sess = 1:2
    values = data_mat(:,sess);
    [h(sess),p(sess),ci,stats] = ttest(values,0);
    t(sess) = stats.tstat;
    d(sess) = computeCohen_d(values, 0);%effect size
end
p
h
[cor_p, h]=bonf_holm(p,.05)%correction for multiple comparisons
end