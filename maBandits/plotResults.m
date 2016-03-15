function plotResults(game, n, N, policies, fname)
    % Plot the results of different algorithms in a given scenario.
    % When there are more than 4 algoprithms the figure will look
    % best when viewed in full screen.
    %
    % Use: plotResults(game, n, N, policies, fname)
    %
    % authors: Olivier Cappé and Aurélien Garivier

    % $Id: plotResults.m,v 1.5 2012-06-06 09:42:04 cappe Exp $
    
    l = length(policies);
    % Find best arm (assuming there is one only...)
    mu = game.mu;
    [~, iBest] = max(mu);
    if (length(iBest) > 1)
        warning('Found %d best arms, plots will be incorrect', length(iBest));
    end
    
    clf;
    
    % Common y-axis values
    ar = -Inf;
    ad = -Inf;
    for i = 1:l
        % Load saved results and compute regret
        load([fname '_n_' num2str(n) '_N_' num2str(N) '_' class(policies{i})]);
        policyName = class(policies{i}); policyName = policyName(7:end);
        regret = (mu(iBest(1))*ones(N,1)*tsave)-cumReward;

        % Regret plot
        risk1 = 25; % Plot quartiles (dark grey)
        risk2 = 5;  % and upper 5 percents quantile (light grey)
        subplot(2,l,i);
        h = area(tsave, [prctile(regret, risk1); prctile(regret, 100-risk1)- ...
                         prctile(regret, risk1); prctile(regret, 100-risk2)- ...
                         prctile(regret, 100-risk1)]');
        set(h(1),'Visible', 'off');
        set(h(2),'FaceColor', [0.86 0.86 0.86]);
        set(h(3),'FaceColor', [0.94 0.94 0.94]);
        hold on;  
        h = plot(tsave, mean(regret), 'k');
        set(h, 'LineWidth', 2);
        xlabel('time', 'FontSize', 8);
        if (i == 1), ylabel('regret'); end
        title(policyName, 'FontSize', 9);
        set(gca, 'FontSize', 8);
        tmp = axis; ar = max(ar, tmp(4));
      
        % Suboptimal draws plot
        subplot(2,l,l+i);
        subOpt = (ones(N,1)*tsave)-cumNbPlayed(:, :, iBest(1));
        clear cumNbPlayed; % Play it is easy on memory
        h = area(tsave, [prctile(subOpt, risk1); prctile(subOpt, 100-risk1)- ...
                         prctile(subOpt, risk1); prctile(subOpt, 100-risk2)- ...
                         prctile(subOpt, 100-risk1)]');
        set(h(1),'Visible', 'off');
        set(h(2),'FaceColor', [0.86 0.86 0.86]);
        set(h(3),'FaceColor', [0.94 0.94 0.94]);
        hold on;
        h = plot(tsave, mean(subOpt), 'k');
        set(h, 'LineWidth', 2);
        xlabel('time', 'FontSize', 8);
        if (i == 1), ylabel('suboptimal draws'); end
        title(policyName, 'FontSize', 9);
        set(gca, 'FontSize', 8);
        tmp = axis; ad = max(ad, tmp(4));
    end
    % Fix comparable y-axis
    for i = 1:l
      subplot(2,l,i);
      axis([0 tsave(end) 0 ar]);
      subplot(2,l,l+i);
      axis([0 tsave(end) 0 ad]);
    end
end
