%Function to create a scatter plot with randomly jittered x
% - INPUT: data should be a cell array with each cell corresponding to the
% elements in one scatter, each cell could be of a diff size

function myfunc_scatterplot(data)
    for i = 1:length(data)
        a = i-0.15;
        b = i+0.15;
        rand_x = (b-a).*rand(1,length(data{i})) + a;
        scatter(rand_x,data{i},20,[0.5 0.5 0.5],'filled',...
            'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25);
        hold on
    end
end