function [data] = RH_Normalize(data,dim)
% Normalize 1D or 2D data

if ~exist('dim','var')
    [~,dim] = max(size(data));
end

data = data-nanmin(data,[],dim);
data = data./nanmax(data,[],dim);
data(isinf(data)) = NaN; %Set Inf to NaN in case you're dividing by zero

end