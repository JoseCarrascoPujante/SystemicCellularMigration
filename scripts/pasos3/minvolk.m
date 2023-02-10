
 function [u,R,factor,improv,mxv,mnv,flagstep,lamhist,var,time,iter]...
     =  minvolk(X,k,tol,KKY,maxit,print,u)

% Finds the ellipsoidal cylinder with minimal k-dimensional cross-sectional
% area containing the columns of X:=[Z;Y] with Z in R^(r*m) 
% and Y in R^(k*m) using the Fedorov-Wynn-Frank-Wolfe method, 
% with Wolfe-Atwood away steps if KKY = 0.
%
% The algorithm returns an ellipsoidal cylinder providing a 
% (1+tol)n-rounding of the convex hull of the columns of X in n-space. 
%
%%%%%%%%%%%%%%%%%%%%  INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% X:=[Z;Y] where Z in R^(r*m) and Y in R^(k*m) is the input data set;
%
% k is the dimension of the space where the minimum-area cross-section
%   is desired;
%
% tol is the tolerance (measure of duality gap), set to 10^-6 by default;
%
% KKY is:
%     0 (default) for the Wolfe-Atwood method using Wolfe's away steps
%     (sometimes decreasing the weights) with the Kumar-Yildirim start;
%     1 if using the Fedorov-Wynn-Frank-Wolfe algorithm 
%     (just increasing the weights) with the Khachiyan initialization;
%     2 for the Wolfe-Atwood method with the Khachiyan initialization;
%
% maxit is the maximum number of iterations, default value 100000;
%
% print is the desired level of printing (default 1); and
%
% u is the initial value for the weights (default as above).
%
%%%%%%%%%%%%%%%%%%%%%  OUTPUT PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% u is the dual solution vector;
%
% R and factor are such that factor^-1/2 * R is the (upper
%     triangular) Cholesky factor of the optimal
%     M := XUX^T, and the trailing k x k submatrix Rbar
%     of R and factor are such that factor^-1/2 * Rbar is the (upper
%     triangular) Cholesky factor of the optimal
%     K(u) := YUY^T-YUZ^T(ZUZ^T)^{-1}ZUY^T;
%
% improv(i) holds the improvement in obj. function value at iteration i;
%
% mxv(i) holds the maximum variance over all points at iteration i;
%
% mnv(i) holds the minimum variance over all points with positive
%     weight at iteration i;
%
% flagstep(i) identifies the type of step taken at iteration i: 1(drop),
%     2(decrease), 3(add), and 4(increase);
%
% lamhist(i) holds the step length lam at iteration i;
%
% var gives the variances of all points w.r.t. the optimal u
%     (var_i(u)=(y_i+Ez_i)^T K(u)^{-1} (y_i+Ez_i));
%
% iterno is the total number of iterations taken; and
%
% time is the total cputime spent in order to obtain the optimal solution.
%
% Calls initwt, updatevar, and updateR.

%%%%%%%%%%%%%%%%%  INITIALIZE INPUT PARAMETERS IF NOT DEFINED  %%%%%%%%%%%%

 if (nargin < 2), error('Please input X and k'); end
 [n,m] = size(X);
 if (nargin < 3), tol = 1e-06; end;
 if (nargin < 4), KKY = 0; end;
 if (nargin < 5), maxit = 100000; end;
 if (nargin < 6), print = 1; end;
 if print,
    fprintf('\n Dimension = %5.0f, Number of points = %5.0f',n,m)
    fprintf(', Tolerance = %5.1e \n',tol);
    fprintf('\n Dimension of y-space = %5.0f \n',k);
 end;
 if (nargin < 7),
    if (KKY >= 1),
       u = (1/m) * ones(m,1);
       fprintf('\n Using Khachiyan initialization \n');
    else
        u = initwt(X,print);
    end;
 end;


%%%%%%%%%%%%%%%%%  INITIALIZE NECESSARY PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%

 st = cputime;
 r = n - k;
 iter = 1;
 n100 = max([n,100]);
 n50000 = max([n,50000]);
 ximult = 1;
 gamma = 1000;
 tol2 = 10^-8;
 mxv = zeros(1,maxit);    % pre-allocate memory for output vectors
 mnv = zeros(1,maxit);
 flagstep = zeros(1,maxit);
 lamhist = zeros(1,maxit);
 improv = zeros(1,maxit);

%%%%%%%%%%%%%%%%%  PARTITION X INTO BLOCKS Y AND Z  %%%%%%%%%%%%%%%%%%%%%%%

Z = X(1:r,:);

 %Y = X(r+1:n,:);     % Y is not used, but is helpful for explaining some
                      % computed quantities.
 
%%%%%%%%%%%%%%%%%  INITIALIZE FACTORIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 upos = find(u > 0);
 lupos = length(upos);

 % M = X(:,upos)*spdiags(u(upos),0,lupos,lupos)*X(:,upos)';
                                     % M   = XUX' (X * diag(u) * X')
 % MZZ = M(1:r,1:r);                 % MZZ = ZUZ' (Z * diag(u) * Z')
                                     % M and MZZ are not used, but are
                                     % helpful in explaining some
                                     % computed quantities.

 A = spdiags(sqrt(u(upos)),0,lupos,lupos)*X(:,upos)'; % A'A = M.
 factor = 1;
 [Q,R]  = qr(A,0);                   % M    = factor^-1 * R' * R
 RZ     = R(1:r,1:r);                % MZZ  = factor^-1 * RZ'* RZ
 
%%%%%%%%%%%%%%%%%%%%%%%%  INITIALIZE VARIANCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 RX   = R' \ X;                  % RX = R^-T * X;
 RZZ  = RZ' \ Z;                 % RZZ = RZ^-T * Z

 zeta = sum(RZZ .* RZZ,1);       % zeta(i) = z_i' * (ZUZ')^-1 * z_i
 xi   = sum(RX .* RX,1);         % xi(i) = x_i' * (XUX^T)^-1 * x_i
 var  = xi - zeta;               % var(i) = (y_i+Ez_i)'K(u)^{-1}(y_i+Ez_i)

%%%%%%%%%%%%%%%%%%%%%%%%  FIND "FURTHEST" AND "CLOSEST" POINTS %%%%%%%%%%%%

 upos = find(u > 0);
 [maxvar,maxj] = max(var);   % maxj has greatest variance among all points
 [minvar,ind] = min(var(upos)); minj = upos(ind); mnvup = minvar;
 if (maxvar > k*(m-1)), 
    fprintf('\n Initialization worse than Khachiyan''s');
    fprintf(', maxvar = %5.0f \n',maxvar);
 end;

% minj has smallest variance among points with positive weight

 mxv(iter) = maxvar; mnv(iter) = minvar;

 if (KKY == 1), fprintf('\n Using KKY \n'); mnvup = k; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%  START ITERATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

 while ((maxvar > (1+tol) * k) || (mnvup < (1-tol) * k)) && (iter < maxit),

%%%%%%%%%%%%%%%   SELECT THE COMPONENT TO INCREASE OR DECREASE  %%%%%%%%%%%

    if (maxvar + mnvup < 2*k), 
         j = minj;
         mvar = mnvup;
    else 
         j = maxj;
         mvar = maxvar;
    end
    uj = u(j);

    % Recompute var(j) and zeta(j).

    flag_recompute = 0;
    xj = X(:,j); zj = Z(:,j);
    RZzj = RZ' \ zj;                    % RZzj  = RZ^-T * z(j)
    MZZzj = factor * (RZ \ RZzj);       % MZZzj = (ZUZ')^{-1} * z(j)
    Rxj = R' \ xj;
    Mxj = factor * (R \ Rxj);		    % Mxj = (XUX')^{-1} * xj
    zetaj = factor * (RZzj' * RZzj);    % zeta(j) recomputed
    xij = factor * (Rxj' * Rxj);        % xi(j) recomputed
    mvarn = xij - zetaj;                % var(j) recomputed
    mvarerror = abs(mvarn - mvar)/max([1,mvar]);
    mvarerrhist(iter) = mvarerror;
    if (mvarerror > tol2),
       flag_recompute = 1;
    end;
    mvar = mvarn;
    
%%%%%%    COMPUTE STEPSIZE LAM (MAY BE NEGATIVE) %%%%%%%%%%%

    % The derivative of the improvement in the objective w.r.t lambda is a
    % negative quantity times a * lambda^2 - 2 b * lambda +c with
    % a, b, and c as below.

    % We need to find the roots of this quadratic; different
    % cases are investigated in the following if clauses such as no
    % real roots, single root, double real roots.

    b = - mvar / 2 + mvar / (2*k) - zetaj;
    a = zetaj * (mvar + zetaj);
    c = 1 - mvar / k;

    % Identify various cases that can arise when solving the quadratic
    % equation a*lambda^2 - 2b*lambda + c =0.

    if (abs(a) < 1e-15),
        if (abs(b) < 1e-15),
            if (c < 0),
                lam = 1e10;
            else
                lam = -uj;
            end
        else
            lam = c / (2*b);
        end
    else
        if b*b > a*c,
            lam = c / (b - sqrt(b*b - a*c));
        else
            if c < 0, lam = 1e10; else lam = -uj; end
        end
    end
    lam = max(lam,-uj);

%   If the step would make ZUZ' close to singular, take a slightly
%   shorter step.

    if (lam < -.9999*uj) & (1 + lam*zetaj < .0001),
        fprintf('\n Truncating step making ZUZ'' singular \n');
        lam = -.9999*uj;
    end

%   If the increase in xi is too great compared to lam, reject the decrease
%   or drop step, and perform an increase or add step

    if (lam < 0) && (lam > -uj),
       ximultold = ximult;
       ximult = (1+lam)*ximult / ((1+lam*xij)*(1-gamma*lam));
    end;
    if (lam == -uj),
       ximultold = ximult;
       ximult = (1+lam)*ximult / (1+lam*xij); 
    end;
    if (lam > 0), ximult = (1+lam)*ximult / (1+gamma*lam); end;

    if ximult > gamma,
       
       fprintf('\n Rejecting decrease or drop \n');
       
       % Reject drop or decrease step and do add or increase step

       ximult = ximultold;
       j = maxj;
       mvar = maxvar;
       uj = u(j);
   
       % Recompute var(j) and zeta(j).

       xj = X(:,j); zj = Z(:,j);
       RZzj = RZ' \ zj;                    % RZzj  = RZ^-T * z(j)
       MZZzj = factor * (RZ \ RZzj);       % MZZzj = (ZUZ')^{-1} * z(j)
       Rxj = R' \ xj;
       Mxj = factor * (R \ Rxj);
       zetaj = factor * (RZzj' * RZzj);    % zeta(j) recomputed
       xij = factor * (Rxj' * Rxj);        % xi(j) recomputed
       mvarn = xij - zetaj;                 % var(j) recomputed
       mvarerror = abs(mvarn - mvar)/max([1,mvar]);
       mvarerrhist(iter) = mvarerror;
       if (mvarerror > 1e-08),
          if print, mvarerror, end;
          flag_recompute = 1;
       end;
       mvar = mvarn;
   
       % The derivative of the improvement in the objective w.r.t lambda is a
       % negative quantity times a * lambda^2 - 2 b * lambda +c with
       % a, b, and c as below.
   
       % We need to find the roots of this quadratic; different
       % cases are investigated in the following if clauses such as no
       % real roots, single root, double real roots.
   
       b = - mvar / 2 + mvar / (2*k) - zetaj;
       a = zetaj * (mvar + zetaj);
       c = 1 - mvar / k;
   
       % Identify various cases that can arise when solving the quadratic
       % equation a*lambda^2 - 2b*lambda + c =0.
   
       if (abs(a) < 1e-15),
           if (abs(b) < 1e-15),
               if (c < 0),
                   lam = 1e10;
               else
                   lam = -uj;
               end
           else
               lam = c / (2*b);
           end
       else
           if b*b > a*c,
               lam = c / (b - sqrt(b*b - a*c));
           else
               if c < 0, lam = 1e10; else lam = -uj; end
           end
       end
       lam = max(lam,-uj);
       if (lam > 0), ximult = (1+lam)*ximult / (1+gamma*lam); end;
    end
    if (j == maxj) && (mvar < (1 + tol)*k),
       fprintf('\n Terminating with epsilon-primal feasible u \n');
       break;
    end;
    lamhist(iter)=lam;    % record the step size taken

    if lam < -uj+tol2, flagstep(iter) = 1;            % drop steps
    elseif lam < 0, flagstep(iter) = 2;               % decrease steps
    elseif uj < tol2, flagstep(iter) = 3;             % add steps
    else flagstep(iter) = 4;                          % increase steps
    end

    % Update u and make sure it stays nonnegative,

    imp = - k*log(1 + lam) + log(1 + lam*mvar/(1 + lam*zetaj));
    improv(iter) = imp;

    uold = u;
    u(j) = max(uj + lam,0); u = (1/(1 + lam)) * u;
    if print && (iter > 1) && (iter-1 == floor((iter-1)/n100) * n100),

%       Print statistics.

       fprintf('\n At iteration %6.0f, maxvar %9.5f',iter-1,maxvar)
       fprintf(', minvar %9.5f',minvar)
    end;

%%%%%%%%%    UPDATE (OR RECOMPUTE) CHOLESKY FACTOR AND VAR   %%%%%%%%%%%%%%

    if (iter > 1) && ((iter-1 == floor((iter-1)/n50000) * n50000) ...
                      || (flag_recompute && print)),
        upos = find(uold > 0);
        lupos = length(upos);
        if (k > 0.5*n),
            M = X(:,upos) * spdiags(uold(upos),0,lupos,lupos) * X(:,upos)';
            normdiff = norm(factor*M - R'*R) / (factor*norm(M));
        else
            MZZ = Z(:,upos) * spdiags(uold(upos),0,lupos,lupos) * Z(:,upos)';
            normdiff = norm(factor*MZZ - RZ'*RZ) / (factor*norm(MZZ));
        end;
        if (normdiff > tol2),
            flag_recompute = 1;
        end;
        if (flag_recompute && print)
            fprintf('\n Relative error in mvar = %8.1e', mvarerror);
            fprintf(' and in XUX'' = %8.1e; reinverting \n', normdiff);
       end;
    end;

    if flag_recompute == 1,
        upos = find(u > 0);
        lupos = length(upos);
        % M = X(:,upos) * spdiags(u(upos),0,lupos,lupos) * X(:,upos)';
                                            % M   = XUX' (X * diag(u) * X')
        A = spdiags(sqrt(u(upos)),0,lupos,lupos) * X(:,upos)';
        factor = 1;
        [Q,R]  = qr(A,0);               % M    = factor * R' * R
        RZ     = R(1:r,1:r);            % MZZ  = factor * RZ'* RZ
        RX     = R' \ X;
        RZZ    = RZ' \ Z;
        zeta   = sum(RZZ .* RZZ,1);
        xi     = sum(RX .* RX,1);
        var    = xi - zeta;
    else

        % Update factorizations.

        [R,factor,down_err]  = updateR(R,factor,xj,lam);
        if down_err, fprintf('\n Error in downdating Cholesky \n');break;end;
        RZ   = R(1:r,1:r);              % MZZ  = factor * RZ'* RZ
        mu = lam / (1 + lam*zetaj);
        zeta = updatevar(zeta,lam,mu,MZZzj,Z);
        nu = lam / (1 + lam*xij);
        xi = updatevar(xi,lam,nu,Mxj,X);
        var = xi - zeta;
    end

    %%%%% FIND "FURTHEST" AND "CLOSEST" POINTS USING UPDATED VAR %%%%%%%%%%%

    upos = find(u > tol2);
    [maxvar,maxj] = max(var);
    [minvar,ind] = min(var(upos)); minj = upos(ind); mnvup = minvar;
    iter = iter + 1;
    mxv(iter) = maxvar;
    mnv(iter) = minvar;
    if (KKY == 1), mnvup = k; end;

 end

%%%%%%%%%%%%%%%%% CALCULATE AND PRINT SOME OF THE OUTPUT VARIABLES %%%%%%%%%

 mxv = mxv(1:iter); mnv = mnv(1:iter);
 flagstep = flagstep(1:iter); improv=improv(1:iter);
 lamhist = lamhist(1:iter);
 iter = iter - 1;

 if print,
    fu = find(u > 1e-12);           % indices of points with positive weight
    for i=1:4, cases(i) = length(find(flagstep==i)); end
    fprintf('\n \n maxvar - k = %4.3e', max(var) - k)
    fprintf(', k - minvar = %4.3e \n', k - min(var(fu)));
    fprintf('\n Drop, decrease, add, increase cases: %6.0f', cases(1));
    fprintf('%6.0f %6.0f %6.0f \n',cases(2),cases(3),cases(4)),
    fprintf('\n Number of positive weights = %7.0f \n', length(fu));
    fprintf('\n Number of iterations       = %7.0f \n', iter);
    fprintf('\n Time taken                 = %7.2f \n \n', cputime - st);
 end;

 return
