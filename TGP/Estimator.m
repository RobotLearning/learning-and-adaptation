% Class for performing parameter estimation
% TODO: For now a wrapper for the minimize function of GPML toolbox

classdef Estimator < handle
    
    properties
        % covariance (kernel) structure [as function]
        cov
        % mean structure [as function]
        mean
        % likelihood function
        lik
        % inference function
        inf
        % hyperparameters (structure with mean/covar values as fields)
        hp
    end
    
    methods       
        % set covariance structure
        function set.cov(obj, handl)
            obj.cov = handl;
        end
        % set mean structure
        function set.mean(obj, handl)
            obj.mean = handl;
        end
        % set hyperparameters structure
        function set.hp(obj, STR)
            % initalize all fields to 0
            obj.hp.mean = []; obj.hp.cov = 0; obj.hp.lik = log(0.1);
            
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.hp), fieldnames(STR))));
            obj.hp = STR;
        end
    end
    
    methods (Access = public)
        
        % Constructor for estimator
        function obj = Estimator(model)
            
            dim_ctx = model.dim_x * 2;
            dim_u = model.dim_u;
            
            % Most general case
            %maskCtx = [ones(1,dim_ctx), zeros(1,dim_u)]; 
            %maskAct = [zeros(1,dim_ctx), ones(1,dim_u)];
            
            % Gravity mismatch
            maskCtx = [0,0,0,0,0,0,0,0,0,1,zeros(1,dim_u)];
            maskAct = [zeros(1,dim_ctx),1,0];
            covlen = sum(maskCtx + maskAct) + 1;
            
            c1 = {@covSEard};
            u1  = {@covLINard};
            covCtx = {'covMask',{maskCtx,c1{:}}}; % context kernel
            covAct  = {'covMask',{maskAct,u1{:}}}; % action kernel
            
            obj.cov = {@covProd,{covCtx,covAct}};
            obj.lik = @likGauss;
            obj.inf = @infExact;
            obj.mean = [];
            
            hyp.mean = [];
            hyp.cov = zeros(covlen,1);
            hyp.lik = log(0.1);
            obj.hp = hyp;
        end
        
        % process control experiment (trajectory runs)
        % trjs is an array of Trajectory classes
        function [x,y] = processCtrlExp(~,trjs)
            
            ctxs = arrayfun(@(trj) trj.PERF.context, [trjs{:}], ...
                            'UniformOutput', false);
            us = arrayfun(@(trj) trj.PERF.u, [trjs{:}], ...
                            'UniformOutput', false);
            yfulls = arrayfun(@(trj) trj.PERF.cost, [trjs{:}], ...
                            'UniformOutput', false);
            ypreds = arrayfun(@(trj) trj.PERF.nom_cost, [trjs{:}], ...
                            'UniformOutput', false);
            ctx = horzcat(ctxs{:});
            u = horzcat(us{:});                        
            y0 = vertcat(yfulls{:});
            y1 = vertcat(ypreds{:});
            
            x = [ctx;u];
            y = y0 - y1;
        end
        
        % cheats using the known function model
        % to see if the x-y's are correctly generated
        function cheat(obj,x,y)
            
            %TODO: process x and y's 
            
        end
        
        % estimate using maximum likelihood
        % fit hyperparameters
        % function handle gp is used in training mode
        % training: [nlZ dnlZ] = gp(hyp, inf, mean, cov, lik, x, y)
        % nlZ: negative log probability of the training data
        % dnlZ: derivative of above
        function mle(obj,x,y)
            Ncg = 100;
            addpath('../../gpml-matlab-v3.1-2010-09-27/util/');
            obj.hp = minimize(obj.hp, @gp, -Ncg, ...
                           obj.inf, obj.mean, obj.cov, obj.lik, x, y);
            %inferred noise standard deviation
            fprintf('Inferred noise standard deviation : %f \n', ...
                     exp(obj.hp.lik));
        end
        
        % estimate using restricted mle
        function rmle(obj,x,y)
            %TODO:
            obj.hp = [];
            obj.beta = 0;
        end
        
    end

    
    
end