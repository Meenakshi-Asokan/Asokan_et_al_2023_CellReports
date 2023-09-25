%Meenakshi M Asokan's code to analyze ACh data
%Related to Figure 6 from Asokan et al., Cell Rep 2023
%September 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;

data_filename = fullfile(currentFolder,'/Data/Figure6/photometry_df_f_allmice.mat');
data = load(data_filename);
df_f_allmice = data.df_f_allmice;

data_filename = fullfile(currentFolder,'/Data/Figure6/behav_meas_allmice_final.mat');
data = load(data_filename);
behav_meas_allmice = data.behav_meas_allmice;
%%
%Initialize variables
t_ach = 4+12;%4s pre CS onset and 12s post CS onset
t_behav = 2+20;%2s pre CS and 20 s post CS
fs_behav = 30;
fs_ach = size(df_f_allmice{1,1,1},2)/t_ach;
CT_sound{1} = cbrewer('seq','RdPu',7);
CT_sound{2} = cbrewer('seq','Greens',7);
%%
%plot 6-day conditioning paradigm behav data (facial motion and pupil dilation) + ACh data from HO-AC and LA

t_baseline = 2;
t_post = 10;
T = t_baseline + t_post;
t_preshock = t_baseline+4;

la_animals = [1,2,3,6,7,8];%4,5 - animals where fiber wasn't precisely over LA - data excluded
num_sess = size(df_f_allmice,3);
for row = 1:4
for sess = 1:num_sess
    for sound = 1:2
        if (row ==1 || row ==2)
            behav_meas = row;
            y = behav_meas_allmice{sound,behav_meas,sess}(:,1:round(T*fs_behav));
        else
            reg = row -2;
            y = df_f_allmice{sound,reg,sess}(:,round(2*fs_ach):round((T+2)*fs_ach));
        end
    for animal = 1:size(y,1)
        if row ==4
            figure(la_animals(animal));
        else
            figure(animal);
        end
    subplot(4,num_sess,(row-1)*num_sess+sess);
    x = 1:size(y,2);
    plot(x, y(animal,:), 'color',CT_sound{sound}(5,:),'linewidth',1);
    % line_sem_plot(x, y, CT_sound{sound}(5,:), 0.5,0.2);
    hold on
    if (row ==1 || row == 2)
    xlim([1 round(T*fs_behav)]);
    xticks(0:round(2*fs_behav):round(T*fs_behav));
    xticklabels(-2:2:12);
    else
    xlim([1 round(T*fs_ach)]);
    xticks(0:round(2*fs_ach):round(T*fs_ach)); 
    xticklabels(-2:2:12);
    end
    hold on
    hline(0,'k--');
    hold on
    box off

    if (sess == 1)
        if (row == 1)
            ylabel('FMR (df/f0)');
            ylim_lowbound = -0.5; ylim_upbound = 3.5;
        elseif (row ==2)
            ylabel('PDR (dp/p0)');
            ylim_lowbound = -0.05; ylim_upbound = 0.3;
        elseif (row ==3)
            ylabel('TeA ACh (df/f0)');  
            ylim_lowbound = -0.5; ylim_upbound = 2;
        else
            ylabel('BLA ACh (df/f0)');
            ylim_lowbound = -0.5; ylim_upbound = 1;
        end
    end


    ylim([ylim_lowbound ylim_upbound]); 
    if row ==1 || row == 2
        fs = fs_behav;
        t_b = 2;
    else
        fs = fs_ach;
        t_b = 2;
    end
    for burst = 1:5
    x_on = round((t_b+burst-1)*fs);
    x_off = round((t_b+burst-1+0.5)*fs);
    x_fill = [x_on, x_on, x_off, x_off ];
    y_fill = [ylim_lowbound, ylim_upbound, ylim_upbound, ylim_lowbound];
    patch(x_fill, y_fill, 'c','FaceAlpha',.1,'EdgeAlpha',0);
    end
    hold on
    end
    end

end
end