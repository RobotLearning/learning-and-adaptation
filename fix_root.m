% Sets the root (Learning-and-Adaptation)

function root = fix_root(rootName)

    % set the root
    if ispc
        root = [cd, '\']; % in windows
    else
        if nargin == 0
            error('Please provide the root name');
        end
        root = [rootName, '/Learning-and-Adaptation/'];
    end