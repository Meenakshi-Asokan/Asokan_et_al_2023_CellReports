%Function to create a scatter plot with randomly jittered x
% - INPUT: data should be a cell array with each cell corresponding to the
% elements in one scatter, each cell could be of a diff size

function myfunc_scatterplot(data,scatter_width,marker_size,marker_color,marker_alpha)
    for i = 1:length(data)
        a = i-scatter_width;
        b = i+scatter_width;
        rand_x = (b-a).*rand(1,length(data{i})) + a;
        scatter(rand_x,data{i},marker_size,marker_color,'filled',...
            'MarkerEdgeColor',marker_color,'MarkerFaceAlpha',marker_alpha,'MarkerEdgeAlpha',marker_alpha);
        hold on
    end
end