function [Tree] = MakeRandomTree(Max_Offspring, Max_Gens)

% This function uses the Galton-Watson branching process to generate a
% random tree.
%
% INPUTS:
%   Max_Offspring - maximum number of offspring any node can have
%   Max_Gens      - maximum number of generations allowed in the tree
%
% OUTPUT:
%   Tree - scalar value indicating the number of leaves (terminal nodes)
%          produced in the final tree
%
% The function samples a random number of generations (up to Max_Gens),
% then recursively generates offspring counts for each node using a
% randomly sampled probability distribution. Offspring are assigned
% according to a branching process that continues until extinction or
% generation depth is reached.
%
% The implementation is a modified version of Galton–Watson branching,
% adapted from Ingemar Kaj and Raimundas Gaigalas (2024) and available at:
% https://www.mathworks.com/matlabcentral/fileexchange/2516-random-trees
%
% NOTE: This function returns only the number of leaves in the tree,
% not the full tree structure — it is intended for integration with
% multifractal tree construction where scale values are tracked separately.
% To visualize such trees see MakeRandomTreeForVisual.m at 
% https://github.com/KevinAndrewHudnall/the-living-tree-of-life/tree/main/Functions

Done = 0;
while(Done == 0)
    
    % Generate offspring probabilies
    Probs = rand(1, Max_Offspring + 1); 
    
    Cum_Probs = [cumsum(Probs) 1]; %Cumulative probabilities

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the branching %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    Parent = 0;
    Gens = randi(Max_Gens); % Determine number of generations
    Child_Count = 1;
    Next_Gen = cell(1, Gens);
    Gens_To_Extinction = 0;
    
    for i = 1:Gens % Branching at each generation
        
        % Determine the linear indices of parents in previous generation.
        Parent_Index = length(Parent) - Child_Count + 1:length(Parent);
        
        % Update Parent_Count to continue getting linear indices of parents
        Parent_Count = Child_Count;
        
        % Generate the offspring distribution
        Child_Dist = rand(1, Parent_Count);

        Child_Count = 0;
        for j = 1:Max_Offspring
            
            % Set Index to contain all values of Parent_Index that are
            % greater than the jth value of Cum_Probs but less than the
            % (j+1)th value of Cum_Probs.
            Index = Parent_Index((Child_Dist > Cum_Probs(j)) &...
                (Child_Dist <= Cum_Probs(j + 1)));

            if (~isempty(Index)) % If at least one parent birthed j kids
                
                % Get the number of offspring for each parent by making a
                % vector of Parent concatenated with the vector Index
                % repeated j times.
                Parent = [Parent repmat(Index, 1, j)];
                
                % Update Child_Count with additional children
                Child_Count = Child_Count + length(Index) * j;
            end
        end
        Next_Gen{i} = Parent;

        % If the tree has gone extinct, record the generations to
        % extinction and break from the loop.
        if (Child_Count == 0)
            Gens_To_Extinction = i;
            break;
        end
    end

    % If the tree didn't go extinct in the first generation
    if (Gens_To_Extinction ~= 1) 
                        
        %Determine leaves and set it to Tree
        Parent = rem(Parent + length(Parent), length(Parent) + 1) + 1;
        Isaleaf = ones(1, length(Parent) + 1);
        Isaleaf(Parent) = zeros(length(Parent), 1);
        
        [~, Tree] = size(Isaleaf(Isaleaf == 1));
        
        Done = 1;
    end
end

