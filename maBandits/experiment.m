function experiment(game, n, N, policy, tsave, fname)
    % Run a number of experiments and save the results subsampling in time
    % Use: experiment(game, n, N, policy, tsave, fname),
    % n is the horizon and N the number of played games.
    %
    % authors: Olivier Cappé and Aurélien Garivier

    % $Id: experiment.m,v 1.17 2012-06-06 09:42:04 cappe Exp $
    
    if nargin<5, tsave = 1:n; end % Times at which the results are saved
    if nargin<6, fname = 'results/exp'; end % Defaut file name
    K = length(tsave); % tsave contains the time indices at which the
                       % results will be saved
    
    % Cumulated rewards: N repeats x K subsanpled times between 1 and n
    cumReward = zeros(N, K);
    % Number of times each arm has been played: N repeats x K subsanpled
    % times between 1 and n x number of arms
    cumNbPlayed = zeros(N, K, game.nbActions);
    
    fprintf('%s %d:', class(policy), N);
    for j = 1:N, 
        [reward, action] = game.play(policy, n);
        cr = cumsum(reward);
        cumReward(j, :) = cr(tsave);
        for a = 1:game.nbActions
            ca = cumsum(action == a);
            cumNbPlayed(j, :, a) = ca(tsave);
        end
        % Once every N/50 runs, display something and save current state
        % of variables
        if (rem(j, floor(N/50))==0) | (j == N)
            fprintf(' %d', j);
            % Expectations of the arms
            mu = game.mu;
            save([fname '_n_' num2str(n) '_N_' num2str(N) '_' class(policy)],...
              'mu', 'n', 'N', 'tsave', 'cumReward', 'cumNbPlayed');
        end
    end
    fprintf('\n');   
end
