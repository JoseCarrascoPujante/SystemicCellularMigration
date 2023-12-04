function [] = ellipse_switch(X,G,precision,colr,precisionType)
%# Modified from:
%# https://stackoverflow.com/a/3419973
%# PrecisionType==1 for Confidence Interval, PrecisionType==2 for STD
if precisionType == 2
    STD = precision;
end
for k=1:2
    %# indices of points in this group
    idx = ( G == k );

    %# substract mean
    Mu = mean( X(idx,:) );
    X0 = bsxfun(@minus, X(idx,:), Mu);

    %# eigen decomposition [sorted by eigen values]
%     [V, D] = eig( X0'*X0 ./ (sum(idx)-1) );     %#' cov(X0)
%     [D, order] = sort(diag(D), 'descend');
%     D = diag(D);
%     V = V(:, order);
    
    %# Represent a specific level of standard deviation or confidence
    %# interval
    if precisionType == 1
        scale = chi2inv(precision,2);     %# inverse chi-squared with dof=#dimensions
    elseif precisionType == 2
%         STD = 2;                  %# set n standard deviations
        precision = 2*normcdf(STD)-1;     %# covers 95.45% of population when STD=2, 68.27% when STD=1
        scale = chi2inv(precision,2);     %# inverse chi-squared with dof=#dimensions
    end
    
    Covar = cov(X0) * scale;
    [V, D] = eig(Covar);
    
    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space

    %# plot cov and major/minor axes
    plot(e(1,:), e(2,:), 'Color',colr,'LineWidth', .75,'LineStyle','--');
    %#quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
    %#quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
end