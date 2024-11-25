function RH_Boxchart(data,color)
% Input: 
%   data = double or cell, cf. RH_Barplot
%   color = Nx3 array;

% Settings
    facealpha = 1; %non-transparent
    boxwidth = 0.5; %default
    linewidth = 1.5; %default = 1
    edgecolor = 'k'; %black
    markerstyle = '.'; %outliers are dots
    markersize = 8; %slim outliers
% Draw box chart
    if isa(data,'double')
        groupdata = repmat(1:size(data,2),size(data,1),1);
        groupdata = reshape(groupdata,numel(groupdata),1);
        plotdata = reshape(data,numel(data),1);
        b = boxchart(groupdata, plotdata);
    elseif isa(data,'cell')
        celldata = cell(length(data),1);
        for c = 1:length(data)
            cdata = data{c};
            groupdata = repmat(1:size(cdata,2),size(cdata,1),1);
            groupdata = reshape(groupdata,numel(groupdata),1);
            plotdata = reshape(cdata,numel(cdata),1);
            celldata{c} = [plotdata,groupdata,c*ones(size(plotdata))];
        end
        celldata = cell2mat(celldata);
        plotdata = celldata(:,1);
        groupdata = celldata(:,2);
        categorydata = celldata(:,3);
        b = boxchart(groupdata,plotdata,'GroupByColor',categorydata);
        if size(color,1) == 1
            color = repmat(color,length(data),1);
        end
    end
% Styling
    for i = 1:length(b)
        b(i).BoxFaceColor = color(i,:);
        b(i).BoxEdgeColor = edgecolor;
        b(i).BoxFaceAlpha = facealpha;
        b(i).BoxWidth = boxwidth;
        b(i).LineWidth = linewidth;
        b(i).MarkerStyle = markerstyle;
        b(i).MarkerSize = markersize;
        b(i).MarkerColor = color(i,:);
    end
    
% X ticks: special case if there are >1 colour groups but 1 primary group
    if length(unique(groupdata))==1 && length(unique(categorydata))>1 %#ok<ISCL>
        ncols = length(unique(categorydata)); %number of color groups
        x = 1 + ((1:ncols) - (ncols+1)/2)/ncols;
        xticks(x);
    end
end