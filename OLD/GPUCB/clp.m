function [x,lambda,status] = clp(Q,c,A,b,Aeq,beq,lb,ub,options)
% CLP Interface to LP/QP solver MEXCLP
%
% [x,z,status] = clp(Q,c,A,b,Aeq,beq,lb,ub,options)
%
% min c'*x
% s.t   A*x <  b, Aeq*x == beq, lb < x < ub
%
% Options structure (see CLP user manual for details)
%    options.solver           [1 (primal), 2 (dual), 3 (interior point) (default 1)].
%    options.maxnumiterations [int>=0 (default 99999999)]
%    options.maxnumseconds    [int>=0 (default 3600)]
%    options.primaltolerance  [double>=0 (default 1e-7)]
%    options.dualtolerance    [double>=0 (default 1e-7)]
%    options.primalpivot      [1 (steepest) | 2 (Dantzig) (default 1)]
%    options.dualpivot        [1 (steepest) | 2 (Dantzig) (default 1)]
%    options.verbose          [0|1|... (default 0)]
%    options.maxabs           [double>0 (default Inf)]
%    options.boundstol        [double>0 (default 1e-8)]
%    options.extractbounds    [0 (off) or 1 (on) (default = 0, use to extract bounds from A, b, Aeq, Beq)]
%    options.qfactor          [double>=0 (default 1, use 2 for optimizing x'Qx instead of 0.5x'Qx)]
%
% please make sure to extract bounds if A or Aeq contain rows with a single
% nonzero element
%
% output
%  x      : primal
%  z      : dual
%  status : 0 - optimal, 1 - infeasible, 2- unbounded

% Author of the original code: Johan Löfberg ETH Zürich.
% Extended by Jan Poland, ABB.CHCRC.C1

% **************************
% Check input
% **************************
%clear mexclp
if nargin<9
    options.solver = 1;              
    if nargin < 8
        ub = [];
        if nargin < 7
            lb = [];
            if nargin < 6
                beq = [];
                if nargin < 5
                    Aeq = [];
                    if nargin < 4
                        b = [];
                        if nargin < 3
                            A = [];
                            if nargin < 2
                                Q = [];
                                if nargin < 1
                                    help clp;return
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if isempty(c)
    c = ones(size(A,2),1);
end

A = [Aeq;A];
b = [beq;b];
neq = length(beq);

if isfield(options,'qfactor'), qfactor = options.qfactor; else qfactor = 1; end

if isfield(options,'extractbounds') && options.extractbounds
  nx = length(c);
  if isfield(options,'maxabs'), maxabs = options.maxabs; else maxabs=inf; end
  if isempty(lb), lb=ones(1,nx)*(-maxabs); else lb=max(lb,-maxabs); end
  if isempty(ub), ub=ones(1,nx)*(maxabs); else ub=min(ub,maxabs); end
  lb = lb(:);
  ub = ub(:);
  
  nzA = sum(A~=0,2);
  ii = find(nzA==1);
  
  [ii1,jj,aa]=find(A(ii,:));
  ii2 = ii(ii1);
  vv = b(ii2)./aa;
  jlb = find((ii2<=neq) | (aa<0));
  jub = find((ii2<=neq) | (aa>0));
  [vvlb,jjlb] = sort(vv(jlb));
  [vvub,jjub] = sort(vv(jub),1,'descend');
  j1lb = jj(jlb(jjlb));
  j1ub = jj(jub(jjub));
  lb(j1lb) = max(lb(j1lb), vvlb);
  ub(j1ub) = min(ub(j1ub), vvub);
  ii = find(nzA<2);
  A(ii,:)=[];
  b(ii)=[];
  neq = neq-sum(ii<=neq);

  if isfield(options,'boundstol'), boundstol = options.boundstol; else boundstol=1e-8; end
  if any(lb>ub+boundstol)
    error('conflicting bounds!');
  end
  j = find(lb>ub);
  if ~isempty(j), lb(j)=ub(j); end
end

% Bug in mexclp...
if isempty(b)
    if isempty(Q) || (nnz(Q)==0)
        x = zeros(length(c),1);
        lambda = [];
        status = 2;
        return
    end
    x = -qfactor*(Q\c);
    lambda = [];
    status = 0;
    return;
end

% **************************
% CLP sparse format
% **************************
cmatcnt = sum(A ~= 0,1);
cmatbeg = full(cumsum([0 cmatcnt]));
cmatbeg = cmatbeg(:)';
nzA = find(A);
cmatind = full(rem(nzA-1,size(A,1))');
cmatind = cmatind(:)';
cmatval = full(A(nzA));
cmatval = cmatval(:)';

if nnz(Q)==0
    cmatbegQ = [];
    cmatindQ = [];
    cmatvalQ = [];
else    
    Q = tril(Q);
    cmatcntQ = full(sum(Q ~= 0,1));
    cmatbegQ = full(cumsum([0 cmatcntQ]));
    cmatbegQ = cmatbegQ(:)';
    nzQ = find(Q);
    cmatindQ = full(rem(nzQ-1,size(Q,1))');
    cmatindQ = cmatindQ(:)';
    cmatvalQ = full(Q(nzQ));
    cmatvalQ = qfactor*cmatvalQ(:)';
end

c = full(c(:))';
b = full(b(:))';
lb = full(lb(:)');
ub = full(ub(:)');

% **************************
% CALL MEX FILE
% **************************
%try
if isfield(options,'x'), xwarmstart=options.x; else xwarmstart=[]; end
[x,lambda,status] = mexclp(cmatbeg,cmatind,cmatval,...
  c,b,neq,lb,ub,cmatbegQ,cmatindQ,cmatvalQ,options,xwarmstart);   
%catch
%end

