function out = nlls(GPSTR, x, y, varargin)

% Function that accepts vectors of hyperparameter values
%
% Inputs : GPSTR - structure containing gp information
%          x - data points
%          y - noisy observations
%          varargin - hyperparameters [variable length to include mean
%          case]
%
% Output : array returning negative log likelihoods for each setting of
%          parameters

% grid of hyperparameters
inpsize = nargin - 3;
[varargout{1:inpsize}] = ndgrid(varargin{:});
fun = @(varargin) gps(GPSTR,x,y,varargin);
out = arrayfun(fun,varargout{:});

end