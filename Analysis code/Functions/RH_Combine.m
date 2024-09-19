function [result] = RH_Combine(vectors)
% This function creates combinations of different sorts.
% At the moment (19-2-2019) it produces all possible unique combinations of vectors of different
% sizes.

% Input:
%   vectors: cell array with 1×N vectors to combine, e.g. vectors = {1:2, 3:5, 6:7, 8:10};

% First, determine whether the inserted vectors are matrices of type double, or cells
if isa(vectors{1},'double')
    
%% Combinations of vectors of different sizes
% From https://fr.mathworks.com/matlabcentral/answers/98191-how-can-i-obtain-all-possible-combinations-of-given-vectors-in-matlab
% Answer by stewpend0us on 30 Jan 2017
    combinations = cell(1, numel(vectors)); %set up the varargout result
    [combinations{:}] = ndgrid(vectors{:});
    combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
    result = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
    
elseif isa(vectors{1},'cell')
%% Combinations of cell inputs
    %First create a numeric equivalent (cfr. supra)
        num_vectors = cell(size(vectors));
        for i = 1:length(vectors)
            num_vectors{i} = 1:numel(vectors{i});
        end
        num_combinations = cell(1, numel(num_vectors)); %set up the varargout result
        [num_combinations{:}] = ndgrid(num_vectors{:});
        num_combinations = cellfun(@(x) x(:), num_combinations,'uniformoutput',false); %there may be a better way to do this
        num_result = [num_combinations{:}]; % NumberOfCombinations by N matrix. Each row is unique.
    %Now fill in the values from each cell
        result = cell(size(num_result));
        for i = 1:size(num_result,1)
            for j = 1:size(num_result,2)
                i_vectors = vectors{j};
                result{i,j} = i_vectors(num_result(i,j));
            end
        end
end
    
end