function [Zproj,num_proj] = myfunc_pca_aap(data,var_exp)
%INPUT - data - each row is a neuron and each column is the
%spikecount/rate/response in a time bin
[num_units, num_timepoints] = size(data);

%mean z_fr for each units across all time points
mu = mean(data,2);
sd = std(data,0,2);

%subtract mean
Z = data - mu;
%SVD
[U,S,V] = svd(Z/sqrt(num_timepoints),'econ');

lambdas = diag(S).^2;
% Variance explained per eigenvalue
var_explained = cumsum(lambdas)./sum(lambdas);
%A heuristic for number of projections to be chosen
num_proj = max(find(var_explained<var_exp));
% % Scree plots
% figure();
% yyaxis right
% plot(lambdas)
% ylabel('Variance per eigenvalue')
% hold on
% yyaxis left
% plot(var_explained,'-o')
% xlabel('Eigenvalue number')
% ylabel('Cumulative fraction of variance')
% hold on 
% vline(num_proj,'k--');
% box off
% set(gca, 'Fontsize',12);

%Projections in the transformed space
Zproj = U'*Z;

% figure();
% subplot(1,2,1);
% CT=cbrewer('div', 'RdBu', 100);
% colormap(flipud(CT));
% imagesc(Z);
% % imagesc(data);
% colorbar();
% % caxis([-5 5]);
% subplot(1,2,2);
% CT=cbrewer('div', 'RdBu', 100);
% colormap(flipud(CT));
% imagesc(Zproj(:,:));
% colorbar();
end