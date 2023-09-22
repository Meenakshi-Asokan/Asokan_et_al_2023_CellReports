%Function to plot a regular histogram while customizing the color and bandwidth
% - INPUT: data should be a 1D vector
function myfunc_histogram(data_histplot,colr,bw)
    h1 = histogram(data_histplot,length(data_histplot));
    % h1.DisplayStyle = 'stairs';
    h1.FaceColor = colr;
    h1.Normalization = 'probability';%can change to just counts
    h1.BinWidth = bw;  
end