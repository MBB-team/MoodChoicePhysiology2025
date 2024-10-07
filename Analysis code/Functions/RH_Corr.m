function [R_Pearson,P_Pearson,R_Spearman,P_Spearman] = RH_Corr(data,data2)
% Correlate variables, but remove NaNs.
% input: a matrix of variables, each column is a separate variable; 
%         OR: two vectors data and data2

if exist('data2','var')
    if size(data,1) < size(data,2)
        data = data';
    end
    if size(data2,1) < size(data2,2)
        data2 = data2';
    end
    X = data;
    Y = data2;
    XY = [X,Y];
    XY = XY(~any(isnan(XY),2),:);
    [R,P] = corr(XY(:,1),XY(:,2),'type','Pearson');
    R_Pearson = R;
    P_Pearson = P;
    [R,P] = corr(XY(:,1),XY(:,2),'type','Spearman');
    R_Spearman = R;
    P_Spearman = P;
    
else
    % Prepare output
        R_Pearson = NaN(size(data,2));
        P_Pearson = R_Pearson; R_Spearman = R_Pearson; P_Spearman = P_Pearson;

    % Set infinites as NaN
        data(isinf(data)) = NaN;

    % Loop through combinations of variables
        for i = 1:size(data,2)
            for j = 1:size(data,2)
                X = data(:,i);
                Y = data(:,j);
                XY = [X,Y];
                XY = XY(~any(isnan(XY),2),:);
                if size(XY,1)<2
                    R_Pearson(i,j) = NaN;
                    P_Pearson(i,j) = NaN;
                    R_Spearman(i,j) = NaN;
                    P_Spearman(i,j) = NaN;
                else
                    [R,P] = corr(XY(:,1),XY(:,2),'type','Pearson');
                    R_Pearson(i,j) = R;
                    P_Pearson(i,j) = P;
                    [R,P] = corr(XY(:,1),XY(:,2),'type','Spearman');
                    R_Spearman(i,j) = R;
                    P_Spearman(i,j) = P;
                end
            end
        end
        if size(data,2) == 2
            R_Pearson = R_Pearson(2);
            P_Pearson = P_Pearson(2);
            R_Spearman = R_Spearman(2);
            P_Spearman = P_Spearman(2);
        end
end

end