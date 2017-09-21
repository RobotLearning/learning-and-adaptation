%% Sokoban environment

% define the grid
side_len = 9;
len_actions = 4;
num_grids = side_len*side_len;
num_blocks = 7;
num_walls = 10;
num_max_walks = 1e6;
num_walks = 0;

grid_vec = 1:grid_len;
% locations of the walls
walls = randperm(num_grids,num_walls);
% area not including the walls as a vector of indices
env = setdiff(grid_vec,walls);
% put in the blocks
blocks = env(randperm(length(env),num_blocks));
free_area = setdiff(env,blocks);
guy = randi(length(free_area),1);
guys_loc(1) = rem(guy,side_len);
guys_loc(2) = ceil(guy/side_len);
% put in the crosses
crosses = env(randperm(length(env),num_blocks));
solved = false;

while num_walks < num_max_walks && ~solved
    % sample an action
    action = randi(len_actions,1);
    % 1 is up, 2 is right, 3 is down, 4 is left
    % find index of next location
    switch action
        case 1
            next_cand_loc(2) = min(side_len,guys_loc(2) + 1);
        case 2
            next_cand_loc(1) = min(side_len,guys_loc(1) + 1);
        case 3
            next_cand_loc(2) = max(1,guys_loc(2) - 1);
        case 4
            next_cand_loc(1) = max(1,guys_loc(1) - 1);
        otherwise
            error('Shouldnt happen!');
    end
    % check if it is possible to move
    guys_next_cand_loc = next_cand_loc(1) * side_len + next_cand_loc(2);
    if ~isempty(find(free_area==guys_next_cand_loc,1))
        guys_loc = next_cand_loc;
    else if ~isempty(find(blocks==guys_next_cand_loc,1)) % there is a block
        % check if it can be moved along the action
        %next_block(1) = 
        end
    end
    
    solved = ~any(blocks - crosses);
end