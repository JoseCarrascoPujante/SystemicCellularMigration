function [] = ellipse_scatter(h,X,conf,colr)
%# Modified from:
%# https://stackoverflow.com/a/3419973
%# substract mean
Mu = mean( X );
X0 = X - Mu; % substract to X its mean

%# eigen decomposition [sorted by eigen values]
%     [V, D] = eig( X0'*X0 ./ (sum(idx)-1) );     %#' cov(X0)
%     [D, order] = sort(diag(D), 'descend');
%     D = diag(D);
%     V = V(:, order);

%# select standard deviation level to represent
% STD = precision;            %# set # standard deviations
% conf = normcdf([-STD STD]); %# calculate confidence interval. Covers 95.45% of population when STD=2, 68.27% when STD=1
% conf = conf(2) - conf(1);
scale = chi2inv(conf/100,2);    %# inverse chi-squared with dof=#dimensions. If you generate random numbers from this chi-square distribution, you would observe numbers greater than scale only (1-conf)*100 of the time.

Covar = cov(X0) * scale;      %# Multiply variance vector by inverse chi-squared function output
[V, D] = eig(Covar);

t = linspace(0,2*pi,100);
e = [cos(t) ; sin(t)];        %# unit circle
VV = V*sqrt(D);               %# scale eigenvectors (ellipse axes)
e = (VV*e) + Mu'; %# project circle back to orig space

%# plot cov and major/minor axes
plot(h,e(1,:), e(2,:),'Color',colr,'LineWidth', 1.25,'LineStyle','--');
% quiver(h,Mu(1),Mu(2), VV(1,1),VV(2,1), 'Color','k')
% quiver(h,Mu(1),Mu(2), VV(1,2),VV(2,2), 'Color','k')
