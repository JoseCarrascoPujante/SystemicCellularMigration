function [] = ellipse(X,G,precision,colr,precisionType)

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
    if strcmp('CI',precisionType)            % If precisionType is CIs...
        scale = chi2inv(precision,2);        % inverse chi-squared with dof=#dimensions
    elseif strcmp('STD',precisionType)       % If precisionType is STDs...
        precision = 2*normcdf(precision)-1;  % covers 95.45% of population when STD=2, 68.27% when STD=1
        scale = chi2inv(precision,2);        % inverse chi-squared with dof=#dimensions
    end
    
    Cov = cov(X0) * scale;
    [V, D] = eig(Cov);
    
    t = linspace(0,2*pi,100);
    e = [cos(t) ; sin(t)];        % unit circle
    VV = V*sqrt(D);               % scale eigenvectors
    e = bsxfun(@plus, VV*e, Mu'); % project circle back to orig space

    % plot cov and major/minor axes
    plot(e(1,:), e(2,:), 'Color',colr,'LineWidth', 1,'LineStyle','--');
    % quiver(Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
    % quiver(Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
end