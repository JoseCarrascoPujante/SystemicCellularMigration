function [] = ellipse_gscatter(h,X,G,conf,colr)

%# Modified from:
%# https://stackoverflow.com/a/3419973
%# PrecisionType==1 for Confidence Interval, PrecisionType==2 for STD
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
    
    %# manual selection of standard deviation level to represent
%     STD = precision;                  %# set n standard deviations
%     conf = 2*normcdf(STD)-1;     %# covers 95.45% of population when STD=2, 68.27% when STD=1
    scale = chi2inv(conf/100,2);     %# inverse chi-squared with dof=#dimensions
    
    Covar = cov(X0) * scale;
    [V, D] = eig(Covar);
    
    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        %# unit circle
    VV = V*sqrt(D);               %# scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); %#' project circle back to orig space
    
    %# plot cov and major/minor axes
    plot(h,e(1,:), e(2,:),'Color',colr,'LineWidth', 1.25,'LineStyle','--');
    %#quiver(h,Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
    %#quiver(h,Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
end