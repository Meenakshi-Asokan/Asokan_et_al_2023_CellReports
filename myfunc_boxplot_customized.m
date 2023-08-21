%%Here's some sample input
% data = cell(1,3);
% for i = 1:3
% data{i} = i*rand(1,i*10);
% end
% box_colors = 'k';
% box_width = 0.5;
% figure();
%%
function myfunc_boxplot_customized(data,box_colors,box_width)
%Function to create a boxplot with mean display as well
% - INPUT: data should be a cell array with each cell corresponding to the
% elements in one box, each cell could be of a diff size
% This function also adds the mean as a magenta circle
mean_data = [];
stderr_data = [];
for i = 1:length(data)
    mean_data(i) = mean(data{i});
    stderr_data(i) = std(data{i})/sqrt(length(data{i}));
end
plot(1:length(data),mean_data,'+k','MarkerFaceColor', [0.5 0.5 0.5]); %Mean addition
% errorbar(1:length(data),mean_data,stderr_data,'ok','MarkerFaceColor', [0.5 0.5 0.5]); %Mean-stderr addition
hold on
mat = [];
xaxiss = [];
for ii = 1:length(data)
    for i = 1:length(data{ii})
        xaxiss{1,length(mat)+i} = (sprintf('%d',ii));
    end
    mat = [mat data{ii}];
end
boxplot(mat,xaxiss,'colors',box_colors,'widths',box_width,'symbol','k+');
hold on
end

