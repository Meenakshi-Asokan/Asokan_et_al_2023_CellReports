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
        % % if you rather want each point colored by channel/depth for eg.
        % CT=cat(1,cbrewer('seq', 'Blues', 32),cbrewer('seq', 'Greens', 32));%Left shank or top half - blue; other one - green
        % scatter(rand_x,[units.cs_amp_assym],30,CT([units.channel],:),'filled','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.5);
    end
end