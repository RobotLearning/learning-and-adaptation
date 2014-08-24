% Controller superclass for holding performances, plotting, etc...

classdef (Abstract) Controller < handle
    
    properties (Abstract)
        
        % number of total episodes so far
        epi
        % color of particular controller
        color
        % name of the particular controller
        name
        % costs incurred (Q-SSE)
        sse
        
    end
    
    methods (Abstract)
        
        % apply control signal - feedforward or feedback
        control(i,trj,model,x)        
        
    end
    
    methods
        
        % get sse costs for the controller
        function record(obj,cost)
            obj.epi = obj.epi + 1;
            obj.sse(obj.epi) = cost;
            
        end
        
    end
    
    
    
end