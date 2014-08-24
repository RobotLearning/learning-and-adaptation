function nll = gps(GPSTR, x, y, varargin)

covfunc = GPSTR.covfunc;
meanfunc = GPSTR.meanfunc;
likfunc = GPSTR.likfunc;
inf = GPSTR.inf;
hyp = GPSTR.hyp;
arg = varargin{1};

switch length(arg)
    case 3    
        % just sampling
        ell = arg(1);
        sf = arg(2);
        sn = arg(3);
        % each combination is considered
        hyp.cov = log([ell; sf]);
        hyp.lik = log(sn);
    case 5
        % called by direct
        ell = arg(1);
        sf = arg(2);
        sn = arg(3);
        beta0 = arg(4);
        beta1 = arg(5);
        hyp.cov = log([ell; sf]);
        hyp.lik = log(sn);
        hyp.mean = [beta0; beta1];
    case 4
        % called by direct
        ell1 = arg(1);
        ell2 = arg(2);
        sf = arg(3);
        sn = arg(4);
        hyp.cov = log([ell1; ell2; sf]);
        hyp.lik = log(sn);
    case 7
        ell1 = arg(1);
        ell2 = arg(2);
        sf = arg(3);
        sn = arg(4);
        beta0 = arg(5);
        beta1 = arg(6);
        beta2 = arg(7);
        hyp.cov = log([ell1; ell2; sf]);
        hyp.lik = log(sn);
        hyp.mean = [beta0; beta1; beta2];
end

nll = gp(hyp,inf,meanfunc,covfunc,likfunc,x,y);
end