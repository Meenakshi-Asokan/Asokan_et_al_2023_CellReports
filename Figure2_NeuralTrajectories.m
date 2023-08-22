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
mouse_numbers = [9,10,11,12,13,14,15,17];
cs_plus = [1,2,1,2,2,1,2,1];
cs_minus = [2,1,2,1,1,2,1,2];
%%
%region 1 - HOAC (RS); region 2 - BLA (RS)
region = 1;
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
    
    figure();
%     CT=cbrewer('seq', 'RdPu', 10);
    for sound = 1:2
        subplot(1,2,sound);
%         colormap(CT);
        imagesc(cs_resp{sound});
        colorbar();
    %     caxis([0 0.2]);
    end

    %Dimensionality red - PCA
    
%     for sound = 1:2
%     data = [cs_resp{sound}];
%     [Zproj{pp,sound},num_proj(pp,sound)] = myfunc_pca_aap(data,var_exp);
%     end
    
%Combine the response to the two sounds and then reduce dimensions
    data = [cs_resp{1} cs_resp{2}];
    Zproj_temp = [];
    [Zproj_temp,num_proj(pp)] = myfunc_pca_aap(data,var_exp);
    Zproj{pp,1} = Zproj_temp(:,1:7000);
    Zproj{pp,2} = Zproj_temp(:,7000+(1:7000));

end

figure();
for pp = 1:2
    subplot(1,2,pp);
    for sound = 1:2   
    plot3(Zproj{pp,sound}(1,:),Zproj{pp,sound}(2,:),Zproj{pp,sound}(3,:));
    hold on
    end
    grid on
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    view([210 30]);
    xlim([-0.25 0.25]);
    ylim([-0.15 0.15]);
    zlim([-0.15 0.15]);

end
% %%
% figure();
% for pp = 1:2
%     subplot(1,2,pp);
%     for sound = 1:2   
%     plot(Zproj{pp,sound}(1,:),Zproj{pp,sound}(2,:));
%     hold on
%     end
% end
%%
%fancy plotting with color indicating timepoint in a 3d trace
CT_sound{1} = cbrewer('seq', 'RdPu', 10000,'pchip');
CT_sound{2} = cbrewer('seq', 'Greens', 10000,'pchip');

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
 
% Set all line segment colors (FOR DATA PROVIDED BY OP)
%     segColors = [v./max(v),zeros(size(v)),1-v./max(v)]; 
%     set(h, {'Color'}, mat2cell(segColors,ones(size(v)),3))
% Set all line segment colors (GENERAL SOLUTION)
% segColors = jet(size(xseg,1)); % Choose a colormap
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
end

%%
%consider how many ever dimensions needed to explain 80% of the variance
%5 and 9 for ACtx, 3 and 4 for BLA
euc = [];
for pp = 1:2
    numproj = num_proj(pp);
    for i = 1:numproj
        if i ==1
            temp = (Zproj{pp,1}(i,:) - Zproj{pp,2}(i,:)).^2;
        else
            temp = temp + (Zproj{pp,1}(i,:) - Zproj{pp,2}(i,:)).^2;
        end
    end
    euc(pp,:) = sqrt(temp);
end
figure();
for pp = 1:2
    plot(1:size(euc,2),euc(pp,:));
    hold on
end
box off
legend('hab','recall');
xticks(0:1000:7000);
xticklabels(-1:1:6);
ylim([0 0.2]);
xlabel('Time wrt. sound onset (s)');
ylabel('Euclidean distance between CS+ and CS- trajectories');
set(gca,'fontsize',12);
% %%
% %Just for the colorbars
% figure();
% colormap(CT_sound{1}(3000+(1:7000),:));
% colorbar();
% caxis([3000 10000]);
% figure();
% colormap(CT_sound{2}(3000+(1:7000),:));
% colorbar();
% caxis([3000 10000]);
% %%
% %quantify trajectory distances 
% spont_win = 1000;
% sound_win = 5000;
% euc_sound = mean(euc(:,spont_win+(1:sound_win)),2);
% euc_spont = mean(euc(:,1:spont_win),2);
% euc_sound_spont_diff = euc_sound - euc_spont;
%%
%%
%Bootstrapping for traj distances
%Takes a bit to run
var_exp = 0.8;
sm = 100;

n_boot = 500;
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
    
%     figure();
% %     CT=cbrewer('seq', 'RdPu', 10);
%     for sound = 1:2
%         subplot(1,2,sound);
% %         colormap(CT);
%         imagesc(cs_resp{sound});
%         colorbar();
%     %     caxis([0 0.2]);
%     end

    %Dimensionality red - PCA
    
%     for sound = 1:2
%     data = [cs_resp{sound}];
%     [Zproj{pp,sound},num_proj(pp,sound)] = myfunc_pca_aap(data,var_exp);
%     end
    
%Combine the response to the two sounds and then reduce dimensions
    data = [cs_resp{1} cs_resp{2}];
    [Zproj_temp,numproj] = myfunc_pca_aap(data,var_exp);
    Zproj{pp,1} = Zproj_temp(:,1:7000);
    Zproj{pp,2} = Zproj_temp(:,7000+(1:7000));
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
