%Meenakshi M Asokan's code for the neural trajectories in Figure 2 from Asokan et al., Cell Rep 2023
%August 2023

close all
clearvars
%Download and save the 'Data' folder in the current working directory;
%The 'Data' folder should have subfolders containing data corresponding to each Figure
currentFolder = pwd;
data_filename = fullfile(currentFolder,'/Data/Figure2/Habituation_HOAC_BLA_units_RS.mat');
data = load(data_filename);
hab_units = data.pre_units_type;
data_filename = fullfile(currentFolder,'/Data/Figure2/Recall_HOAC_BLA_units_RS.mat');
data = load(data_filename);
recall_units = data.post_units_type;
%Add to path the folder Utils including all subfolders
addpath(genpath(fullfile(currentFolder,'Utils')));
%%
mouse_numbers = [9,10,11,12,13,14,15,17];
%Which sweep is which CS gets counterbalanced across mice
cs_plus = [1,2,1,2,2,1,2,1];
cs_minus = [2,1,2,1,1,2,1,2];
regions = cell(1,2);
regions{1} = 'HO-AC';
regions{2} = 'BLA';
region = 1;%Change region number to plot the neural trajectories for HO-AC and BLA separately
Zproj = cell(2,2);%pre_post,2 sounds
num_proj = [];
var_exp = 0.8;
close all
% Pre-processing prior to dimensionality reduction to capture shared
% variance and not the stochastic fluctuations inherent to each neuron -
% Averaging across all 15 trials and then smoothing
sm = 100;
for pp = 1:2 %pre-post: 1 - hab; 2 - recall
    if pp ==1
        units = hab_units{1,region};
    else
        units = recall_units{1,region};
    end
    cs_resp = cell(1,2);
    for sound = 1:2
        temp = [];
        for clu = 1:length(units)
            m_num = find(mouse_numbers == units(clu).mouse_num);
            cs_p = cs_plus(m_num);
            cs_m = cs_minus(m_num);
            cs = [cs_p cs_m];
            temp(clu,:) = smoothdata(mean(units(clu).sound_raster{1,cs(sound)},1),'gaussian',sm);
        end
        cs_resp{sound} = temp(:,1:7000);
    end
%     %If you want to plot the neurograms for each sound for hab and recall    
%     figure();
%     for sound = 1:2
%         subplot(1,2,sound);
%         imagesc(cs_resp{sound});
%         colorbar();
%     end

    %Dimensionality red - PCA   
    %Combine the response to the two sounds and then reduce dimensions
    data = [cs_resp{1} cs_resp{2}];
    Zproj_temp = [];%Projections in the transformed space
    [Zproj_temp,num_proj(pp)] = myfunc_pca_aap(data,var_exp);
    Zproj{pp,1} = Zproj_temp(:,1:7000);
    Zproj{pp,2} = Zproj_temp(:,7000+(1:7000));

end

figure();
for pp = 1:2
    subplot(1,2,pp);
    for sound = 1:2  %sound 1-CS+; 2-CS- 
    plot3(Zproj{pp,sound}(1,:),Zproj{pp,sound}(2,:),Zproj{pp,sound}(3,:));
    hold on
    end
    grid on
    legend({'CS+','CS-'});
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    view([210 30]);
    xlim([-0.25 0.25]);
    ylim([-0.15 0.15]);
    zlim([-0.15 0.15]);
    if (pp ==1)
        title('Habituation');
    else
        title('Recall');
    end
end

%%
%fancy plotting with color indicating timepoint in a 3d trace
CT_sound{1} = cbrewer('seq', 'RdPu', 10000,'pchip');%For CS+
CT_sound{2} = cbrewer('seq', 'Greens', 10000,'pchip');%For CS-

figure();
for pp = 1:2
    subplot(1,2,pp);
    for sound = 1:2   
        x = Zproj{pp,sound}(1,:)';
        y = Zproj{pp,sound}(2,:)';
        z = Zproj{pp,sound}(3,:)';
        xseg = [x(1:end-1),x(2:end)]; 
        yseg = [y(1:end-1),y(2:end)]; 
        zseg = [z(1:end-1),z(2:end)]; 
        % Plot all line segments (invisible for now unless you remove 'visible','off')
        % h = plot(xseg',yseg','-','LineWidth',1,'Visible','Off'); 
        h = plot3(xseg',yseg',zseg','-','LineWidth',1.5,'Visible','Off'); 
        hold on
        % axis equal;
        % xlim([min(x) max(x)]);
        % ylim([min(y) max(y)]);
         
        % Set all line segment colors 
        segColors = CT_sound{sound}(3000+(1:size(xseg,1)),:);
        set(h, {'Color'}, mat2cell(segColors,ones(size(xseg,1),1),3))
        % Make segements visible in loop
        %     for n = 1:numel(t)-1
        %         pause(.05)
        %         set(h(n),'Visible', 'on')
        %         drawnow(); 
        %     end
        % Or make it all visible at once
        set(h, 'Visible', 'on')
    end
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    view([210 30]);
    xlim([-0.25 0.25]);
    ylim([-0.15 0.15]);
    zlim([-0.15 0.15]);
    grid on
    if (pp ==1)
        title('Habituation');
    else
        title('Recall');
    end
end

%%
%Bootstrapping to obtain eucledian distances between the cs+ and cs- neural trajectories

var_exp = 0.8;%Consider as many dimensions needed to explain 80% of variance - 5 and 9 for ACtx, 3 and 4 for BLA
sm = 100;

n_boot = 500;%Number of bootstraps
euc = zeros(n_boot, 2, 7000);
tic
for boot = 1:n_boot
    Zproj = cell(2,2);%pre_post,2sounds
for pp = 1:2
    if pp ==1
        units = hab_units{1,region};
    else
        units = recall_units{1,region};
    end
    num_units = length(units);
    %Random sampling with replacement
    clu_resamp = datasample(1:num_units,num_units,'Replace',true);
    units_resamp = units(clu_resamp);
    cs_resp = cell(1,2);
    for sound = 1:2
        temp = [];
        for clu = 1:length(units_resamp)
            m_num = find(mouse_numbers == units_resamp(clu).mouse_num);
            cs_p = cs_plus(m_num);
            cs_m = cs_minus(m_num);
            cs = [cs_p cs_m];
            temp(clu,:) = smoothdata(mean(units_resamp(clu).sound_raster{1,cs(sound)},1),'gaussian',sm);
        end
        cs_resp{sound} = temp(:,1:7000);
    end
    
    %Combine the response to the two sounds and then reduce dimensions
    data = [cs_resp{1} cs_resp{2}];
    [Zproj_temp,numproj] = myfunc_pca_aap(data,var_exp);
    Zproj{pp,1} = Zproj_temp(:,1:7000);
    Zproj{pp,2} = Zproj_temp(:,7000+(1:7000));

    %Compute the eucledian distance between the two trajectories in a
    %n-dimensional space defined by the number of projections that explain 80% of variance
    for i = 1:numproj
        if i ==1
            temp_dist = (Zproj{pp,1}(i,:) - Zproj{pp,2}(i,:)).^2;
        else
            temp_dist = temp_dist + (Zproj{pp,1}(i,:) - Zproj{pp,2}(i,:)).^2;
        end
    end
    euc(boot,pp,:) = sqrt(temp_dist);
end
boot
end
toc
%%
%Volin plots and stats
data_filename = fullfile(currentFolder,'/Data/Figure2/stats_data_pca_traj_maingroup.mat');
% data_filename =
% fullfile(currentFolder,'/Data/Figure2/stats_data_pca_traj_pseudocond.mat');%Uncomment this and Comment
% the line above to plot the viloin plots for pseudoconditioning data

stats_data = load(data_filename);
data = stats_data.stats_data;
%Only 500 bootstraps taken for consistency with other analyses

%normal as per ktest

%2 sample/ unpaired t-test
stats = cell(1,length(data));
h = [];
p = [];
t = [];
d = [];
for r_num = 1:length(data)
    values1 = data{r_num}(1,1:500);
    values2 = data{r_num}(2,1:500);
    [h(r_num),p(r_num),ci,stats{r_num}] = ttest2(values2,values1);
    t(r_num) = stats{r_num}.tstat;
    d(r_num) = computeCohen_d(values2, values1,'independent');
end
h
p
d

%%
%violin plot
figure();
for r_num = 1:2
    subplot(1,2,r_num);
    data_mat = [];
for pp = 1:2
    data_mat(:,pp) = data{r_num}(pp,1:500);%500 bootstraps, not 1000 (to keep it consistent)
end
violinplot(data_mat,{'hab','recall'}, 'QuartileStyle','boxplot','HalfViolin','right',...
    'ShowData',false);
xlim([0.5 2.5]);
box off
ylabel('Eucledian distance between trajectories');
title(regions{r_num});
set(gca,'fontsize',12);

end