classdef Policy<handle
    % Generic policy interface
    %
    % authors: Olivier Cappé, Aurélien Garivier

    % $Id: Policy.m,v 1.6 2012-06-05 13:26:38 cappe Exp $
    
    properties
    end
    
    methods
        function init(self, nbActions, horizon), end % to be called before a new game
        function a = decision(self), end % chooses the next action
        function getReward(self, reward), end % update after new observation
    end
    
end
