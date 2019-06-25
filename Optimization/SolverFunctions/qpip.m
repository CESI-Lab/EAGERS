function [x,feasible] = qpip(qp)
%function [vars,status,stats] = qpip(Q,c,C,d,A,b,options)
% QPIP: Implementation of Mehrotra-like Primal-Dual QP solver
% for general quadratic programs.  Based on the C++ solver OOQP.
%
% Usage:
%   [x,Feasible] = qpip(QP)   
%
%   This finds a solution to the problem:
%
%   min 1/2 x'Hx + f'x
%   s.t. Aeq*x = beq, A*x <= b.
%

% Original qpip Author : P. Goulart, Univ. Oxford.  Last modified 24/07/2017
%
% NB: this code is follows closely the the implemention of the C++
% solver OOQP, but in a pure matlab implementation

%% convert QP form into problem structure (p) used by qpip
n = length(qp.f);
if ~isfield(qp,'Organize')|| ~isfield(qp.Organize,'IC')
    ic = 0;
else
    ic = max(qp.Organize.IC); %number of initial conditions. Don't need upper/lower bounds because there is an equality constraint
end
if ~isempty(qp.lb)
    I = [zeros(n-ic,ic),eye(n-ic)];
    if isempty(qp.A)
        p.data.C = [I;-I];
        p.data.d = [qp.ub(ic+1:end);-qp.lb(ic+1:end);];
    else
        p.data.C = [qp.A; I;-I];
        p.data.d = [qp.b;qp.ub(ic+1:end);-qp.lb(ic+1:end);];
    end
else
    p.data.C = qp.A;
    p.data.d = qp.b;
end
            
p.data.Q = qp.H;
p.data.c = qp.f;
p.data.A = qp.Aeq;
p.data.b = qp.beq;
p.data.nx = size(p.data.A,2);
p.data.ny = size(p.data.A,1);
p.data.nz = size(p.data.C,1);
p.data.ns = p.data.nz;

p.options.max_iter = 100;
p.options.gamma = 0.99;
p.options.tolMu = 1e-7;
p.options.tolR = 1e-7;
p.options.minPhi = 1e10;
p.options.verbose = true;

%create vector variables of the right size
p.variables.x = zeros(p.data.nx,1);  
p.variables.y = zeros(p.data.ny,1);  
p.variables.z = ones(p.data.nz,1);  
p.variables.s = ones(p.data.ns,1);  

%create matrices of the right size to hold the residual components
p.residuals.rQ  = zeros(p.data.nx,1);
p.residuals.rA  = zeros(p.data.ny,1);
p.residuals.rC  = zeros(p.data.nz,1);
p.residuals.rS  = zeros(p.data.ns,1);

%% initialize problem structure and populate  with problem data
%Create the Jacoban matrix with a dummy Sigma
S  = -speye(p.data.ns);
Z1 = spalloc(p.data.ny,p.data.ny,0);
Z2 = spalloc(p.data.nz,p.data.ny,0);
Z3 = Z2';
ZA = spalloc(p.data.nx,p.data.ny,0);
ZC = spalloc(p.data.nx,p.data.ns,0);

%construct the jacobian.  Since ldl is used for factorization, only the lower part is needed
p.linsys.J = [p.data.Q  ZA  ZC; p.data.A  Z1  Z3; p.data.C  Z2   S];
        
%get the indices for the entries of S
idx = (p.data.nx + p.data.ny) + (1:p.data.ns);
p.linsys.idxSigma = sub2ind(size(p.linsys.J),idx,idx);

%configure structure for factors
p.linsys.factors.L = [];
p.linsys.factors.D = [];
p.linsys.factors.perm  = [];
p.linsys.factors.scale = [];

%construct the norm
normD = norm([p.data.A(:);p.data.b(:);p.data.C(:);p.data.d(:);p.data.Q(:);p.data.c(:)],inf);

%problem variable initialization, default first step, find some interior point (large z and s)
sdatanorm = sqrt(normD);
p.variables.z(:) = sdatanorm;
p.variables.s(:) = sdatanorm;

p.residuals.rQ  = p.data.Q*p.variables.x + p.data.A'*p.variables.y + p.data.C'*p.variables.z + p.data.c;
p.residuals.rA  = p.data.A*p.variables.x - p.data.b;
p.residuals.rC  = p.data.C*p.variables.x + p.variables.s - p.data.d;
p.residuals.rS(:) = p.variables.z.*p.variables.s + 0;  

%factor and solve
p.linsys.J(p.linsys.idxSigma) = -(p.variables.s./p.variables.z);%update the jacobian
[p.linsys.factors.L,p.linsys.factors.D,p.linsys.factors.perm,S] = ldl(p.linsys.J,'vector');%LDL' factorization
p.linsys.factors.scale = diag(S);
[step,~] = qpip_solve(p);   

%take the full affine scaling step
p.variables.x = p.variables.x + step.x;
p.variables.y = p.variables.y + step.y;
p.variables.z = p.variables.z + step.z;
p.variables.s = p.variables.s + step.s;

%shift the bound variables
shift = 1e3 + max(0,max(-[p.variables.z;p.variables.s]));
p.variables.z = p.variables.z + shift;
p.variables.s = p.variables.s + shift;

%get the complementarity measure mu
muval = 0;
if ~isempty(p.variables.z)
    muval = sum(abs(p.variables.z.*p.variables.s))./length(p.variables.z);
end

%% -------------------------------------------
for iter = 1:p.options.max_iter      
    %Update the right hand side residuals
    p.residuals.rQ  = p.data.Q*p.variables.x + p.data.A'*p.variables.y + p.data.C'*p.variables.z + p.data.c;
    p.residuals.rA  = p.data.A*p.variables.x - p.data.b;
    p.residuals.rC  = p.data.C*p.variables.x + p.variables.s - p.data.d;
    
    %test for convergence or infeasibility
    %dualitygap: Calculate duality gap defined as gap = x'Qx + c'x + b'y + d'z
    gap      = abs(p.variables.x'*(p.data.Q*p.variables.x) + p.data.c'*p.variables.x + p.data.b'*p.variables.y + p.data.d'*p.variables.z);
    normR    = norm([p.residuals.rQ(:);p.residuals.rA(:);p.residuals.rC(:)],inf);
    phi      = (normR + gap)./normD;
    minPhi   = min(p.options.minPhi,phi);
    if(muval <= p.options.tolMu && normR <= p.options.tolR*normD) 
        feasible = 1;
        break
    elseif(phi > 10-8 && phi > 10^4*minPhi)
        feasible = -1;
    else
        feasible = 0;
    end
    
    %PREDICTOR STEP   
    %find the RHS for this step
    p.residuals.rS(:) = p.variables.z.*p.variables.s + 0;   
    %factor and solve
    p.linsys.J(p.linsys.idxSigma) = -(p.variables.s./p.variables.z);%update the jacobian
    [p.linsys.factors.L,p.linsys.factors.D,p.linsys.factors.perm,S] = ldl(p.linsys.J,'vector');%LDL' factorization
    p.linsys.factors.scale = diag(S);
    [stepAff,dir] = qpip_solve(p);    

    % Calculate centering parameter
    %determine the largest step that preserves consistency of the multiplier constraint  
    if(any(dir<0))
        tmp   = dir./[p.variables.z;p.variables.s];
        alphaAff = min(max(0,1/max(-tmp)),1);
    else
        alphaAff = 1;
    end
    %mu_step: calculate the value of z's/m given a step in this input direction
    muAff = sum(abs(p.variables.z+alphaAff.*stepAff.z).*(p.variables.s+alphaAff.*stepAff.s))./length(p.variables.s);   
    sigma    = (muAff/muval)^3;
        
    %CENTERING-CORRECTOR STEP
    p.residuals.rS(:) = p.residuals.rS(:) + stepAff.z.*stepAff.s - sigma*muval;  
    %solve for the corrected system
    [stepCC,dir] = qpip_solve(p);          
        
    %determine the largest step that preserves consistency of the multiplier constraint  
    if(any(dir<0))
        tmp   = dir./[p.variables.z;p.variables.s];
        alphaMax = min(max(0,1/max(-tmp)),1);
    else
        alphaMax = 1;
    end 
    
    %a simple blocking mechanism
    stepSize = alphaMax*p.options.gamma;    
    
    %take the step and update mu
    p.variables.x = p.variables.x + stepSize.*stepCC.x;
    p.variables.y = p.variables.y + stepSize.*stepCC.y;
    p.variables.z = p.variables.z + stepSize.*stepCC.z;
    p.variables.s = p.variables.s + stepSize.*stepCC.s;
    %calculate complementarity measure 
    if ~isempty(p.variables.z)
        muval = sum(abs(p.variables.z.*p.variables.s))./length(p.variables.z);
    end
    
end
x = p.variables.x;%extract the full solution variable set

%calculate the residuals
p.residuals.rQ  = p.data.Q*p.variables.x + p.data.A'*p.variables.y + p.data.C'*p.variables.z + p.data.c;
p.residuals.rA  = p.data.A*p.variables.x - p.data.b;
p.residuals.rC  = p.data.C*p.variables.x + p.variables.s - p.data.d;

% %collect solve statistics
% stats.iterations = iter;
% stats.status     = status;
% stats.value      = 1/2*(x'*(p.data.Q*x)) + p.data.c'*x;
% stats.dualityGap = dualityGap(p);
% stats.solveTime  = toc(tstart);
% stats.residuals  = p.residuals;
end%ends function qpip

function [variables,dir] = qpip_solve(p)
%Solve the Newton system for a given set of residuals, using the current Jacobian factorization
%get the Jacobian factorization
L = p.linsys.factors.L;
D = p.linsys.factors.D;
perm  = p.linsys.factors.perm;
scale = p.linsys.factors.scale;

%construct the rhs to be solved (including any pre-elimination)
rC = p.residuals.rC - p.residuals.rS./p.variables.z;
rhs = [p.residuals.rQ;p.residuals.rA;rC];

%solve it
lhs(perm) = L'\(D\(L\(rhs(perm).*scale(perm)))).*scale(perm);
lhs = lhs(:); %force a column

%parse the solution (including any post-solving)back into a *copy* of the variables
variables = p.variables;
variables.x(:) =  -lhs(1:p.data.nx);
variables.y(:) =  -lhs((1:p.data.ny) + p.data.nx);                                
variables.z(:) =  -lhs((1:p.data.nz) + p.data.nx + p.data.ny);  
variables.s(:) = -(p.residuals.rS + p.variables.s.*variables.z)./p.variables.z;  
dir = [variables.z;variables.s];%the proposed z and s directions
end%ends function qpip_solve