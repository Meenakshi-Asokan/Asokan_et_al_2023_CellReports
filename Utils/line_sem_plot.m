function h = line_sem_plot(x, y, color,lw,fa)
% INPUT y: a matrix, each row represents an observation, each collumne
%          represent a variable(e.g. time)
%       x: variable (e.g time)
%       lw:linewidth
%       fa:facealpha

y_mean = mean(y, 1, 'omitnan');
y_sem = std(y, 0, 1, 'omitnan')/sqrt(size(y, 1));

% y_sem(isnan(y_sem)) = 0; % add by ke to avoid nan
% lowBound = y_mean - y_sem;
% upBound  = y_mean;

lowBound = y_mean;
upBound  = y_mean + y_sem;
x2 = [x, fliplr(x)];
inBetween = [upBound, fliplr(lowBound)];
h1 = fill(x2, inBetween, color);
set(h1, 'facealpha', fa)
set(h1, 'EdgeColor', 'none')
hold on
h = plot(x, y_mean, 'color',color,'linewidth',lw);
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])