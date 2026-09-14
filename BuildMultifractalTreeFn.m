function [S,P, Leaves] = BuildMultifractalTreeFn(MaxOffspring, MaxGens, ITERATIONS)

% The function generates a multifractal tree through iterative
% branching and scale assignment, returning matrices of progeny structure
% (P), scale assignments (S), and the recursive lineage scales (Leaves).
%
% INPUTS:
%   MaxOffspring  - maximum number of offspring any node can have
%   MaxGens       - maximum number of generations in each finite tree
%   ITERATIONS    - number of recursive replacement iterations
%
% OUTPUTS:
%   S      - scale matrix reconstructed along complete lineage paths
%   P      - progeny-count matrix reconstructed along lineage paths
%   Leaves - cell array containing the scale of every recursive lineage
%            present at each iteration, with lineage multiplicity preserved
%
% The function recursively constructs a stochastic multifractal tree using
% a Galton–Watson process with randomly assigned scales. Each subtree is 
% expanded in parallel, and the full branching structure is recovered to 
% construct explicit ancestor-descendant matrices. Entropy is computed as 
% the log-change in scale across levels for each path.

%%%%%%%%%%%%%%%%%%%%% LEVEL 1: generate the first tree %%%%%%%%%%%%%%%%%%%%

Scale = rand(1, 1); % Randomly select a scale to generate tree
T = MakeRandomTree(MaxOffspring, MaxGens); % Generate random tree
Leaves = cell(ITERATIONS, 1); % Initialize cell array to hold leaf scales
Progenies = cell(ITERATIONS, 1); % Initialize cell array to hold progenies
Leaves{1} = ones(1, T) * Scale; % Store the scale of the leaves in tree 1
Progenies{1} = T; % Store the number of progeny in tree 1

% Parallelize tree generation for each leaf scale in further iterations
for n = 1:(ITERATIONS - 1)
    % Initialize temporary storage for this iteration
    Temp_Leaves = gpuArray([]);
    Temp_Progenies = gpuArray([]);

    Leaves_n_Chunk = Leaves{n};
    % Use parfor to parallelize tree generation for each leaf scale
    parfor i = 1:length(Leaves{n}) % loop over the leaves in the nth iterate
    
        % Randomly select a scale for each leaf in the nth iterate
        Scale = Leaves_n_Chunk(i) * rand(1, 1);
    
        % Generate a random tree for each one of these scales
        T = MakeRandomTree(MaxOffspring, MaxGens);

        % Use temporary arrays to avoid issues with parfor and also
        % avoid gpuArray operations inside the parfor loop
        Temp_Leaves = [Temp_Leaves, ones(1, T) * Scale];
        Temp_Progenies = [Temp_Progenies, T];
    end

    % Store the results back into the cell arrays
    Leaves{n + 1} = gather(Temp_Leaves);
    Progenies{n + 1} = gather(Temp_Progenies);
end

%%%%%%%%%%% Organize Progenies and Leaves into matrices P and S %%%%%%%%%%%

% The number of offspring of every leaf at every iterate has been 
% recorded in Progenies, and likewise the scale of every leaf at every 
% iterate has been recorded in Leaves. These data structures need to be 
% converted to matrices that have N rows and as many columns as there are
% leaves in the final iterate. Doing so will make ancestral relationships
% explicit. Columns of the matrices will then correspond to distinct paths
% through the tree of life. Constructing these matrices requires
% determining how many times each member of every cell except the last of
% Progenies (and Leaves) must be repeated and populated in the matrix.
% Basically, the ancestral relationships are being back constructed using
% the information in Progenies, and then P and S are constructed via these
% ancestral relationships.

Prog_Size = size(Progenies, 1);
% tally holds the number of times each member of each cell of Progenies
% must be repeated in order to form the matrix P
tally = cell(Prog_Size, 1);
% Ensure value is initialized properly
parfor i = 1:size(Progenies, 1) % Loop over all rows in Progenies

    tally_i = zeros(1, size(Progenies{i}, 2));
    value_i = cell(1, Prog_Size);

    for j = 1:size(Progenies{i}, 2) % Loop over all columns in Progenies

        % Count where we start column-wise in our repititions
        count_1 = sum(Progenies{i}(1:(j - 1))) + 1; 

        % Count where we stop column-wise in our repititions relative to
        % our start location
        count_2 = sum(Progenies{i}(1:j));

        % Count the total number of repitions we need; initialize to 1
        count_3 = 1;

        % Loop over all cells of Progenies in iterate (i + 1)
        for k = (i + 1):(Prog_Size)

            if (j == 1) % If we're on the first column

                % The total number of repeats equals the stop location
                count_3 = count_2;

                % Determine the next start-relative stop location
                count_2 = sum(Progenies{k}(count_1:count_2));

                % Start location is 1
                count_1 = 1;

            else
                % Number of repeats equals the start-relative stop location
                count_3 = count_2; 

                temp = count_1; % Set a temporary variable 

                if (k < Prog_Size) % If we're not at the last row

                    % Update start location
                    count_1 = value_i{1, k + 1} + 1; 

                    % Update start-relative stop location
                    count_2 = count_1 + sum(Progenies{k}(temp:count_2)) - 1;
                end
            end
            % Store the repeats required in a cell array
            value_i{1, k} = count_3; 
        end
        % Store the repeats required in a cell array
        tally_i(1, j) = count_3; 
    end

    % Assign the results back to the main cell array
    tally{i} = tally_i;
end

% Now use tally to build P
P = zeros(size(Progenies, 1), size(Progenies{end}, 2)); % Initialize P

% Generate the first row of P
P(1, 1:size(Progenies{end}, 2)) = ones(1, size(Progenies{end}, 2)) * Progenies{1}(1, 1);

% Preallocate temporary storage for Begin and Continue
Begin_All = cell(Prog_Size - 1, 1);
Continue_All = cell(Prog_Size - 1, 1);

% Parallelize the generation of rows 2 through (last - 1) of P
parfor i = 2:(Prog_Size - 1)
    tally_i = tally{i}; % Extract the current tally
    Progenies_i = Progenies{i}; % Extract the current progenies
    Continue = [];
    Begin = ones(1, tally_i(1)) * Progenies_i(1, 1);

    for j = 2:length(tally_i)
        
        Continue = [Continue, ones(1, (tally_i(j) - tally_i(j - 1))) *...
            Progenies_i(1, j)];
    end
    Begin_All{i - 1} = Begin;
    Continue_All{i - 1} = Continue;
end

% Combine the results to form P
for i = 2:(size(Progenies, 1) - 1)
    P(i, 1:size(Progenies{end}, 2)) = [Begin_All{i - 1}, Continue_All{i - 1}]; 
end

% Generate the last row of P
P(end, 1:size(Progenies{end}, 2)) = Progenies{end}(1, :);

% Use Progenies to build S
S = cell(size(Leaves, 1), 1); % Initialize S
Temp_1 = Leaves;
for p = 0:(size(Leaves, 1) - 2)
    Temp_2 = cell(size(Leaves, 1) - p, 1);
    for i = (size(Leaves, 1) - p):-1:2
        Progenies_Sizes = Progenies{i + p, 1};
        Values_To_Repeat = Temp_1{i - 1, 1};
        Repeated_Values = cellfun(@(sz, val) repmat(val, 1, sz), ...
            num2cell(Progenies_Sizes), num2cell(Values_To_Repeat), 'UniformOutput', false);
        Temp_2{i - 1} = [Repeated_Values{:}];
    end
    Temp_1 = Temp_2;
    S{(size(Leaves, 1) - (p + 1)), 1} = Temp_2{(size(Leaves, 1) - (p + 1)), 1};
end
S{size(Leaves, 1), 1} = Leaves{end, 1};
S = cell2mat(S);

end