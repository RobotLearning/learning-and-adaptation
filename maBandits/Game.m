classdef Game<handle
    % Generic bandit game interface
    %
    % authors: Olivier Cappé and Aurélien Garivier

    % $Id: Game.m,v 1.11 2012-06-05 13:26:38 cappe Exp $

    properties
        nbActions % number of actions available
        mu % expectations of arms
    end
    
    methods
        function [reward, action] = play(self, policy, n), end % starts a new game 
        function K = initRewards(self,n), end % to be called before a game
    end    
end
