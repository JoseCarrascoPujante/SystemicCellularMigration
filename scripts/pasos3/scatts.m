
function figures = scatts(field_names, results, figures, type)

% Plots systemic movement parameters in 2D and 3D scatterplots

% Merge all metrics in a single table
results.full = [];

for index = 1:length(field_names)

        results.full = cat(1,results.full, results.(field_names{index}));
end

if contains(type, 'normalized')

    % Normalized version of each metric
    
    rmsfAlpha = (results.full(:,1)-min(results.full(:,1)))/(max(results.full(:,1))-min(results.full(:,1)));
    srmsfAlpha = (results.full(:,2)-min(results.full(:,2)))/(max(results.full(:,2))-min(results.full(:,2)));
    rmsfR2 = (results.full(:,3)-min(results.full(:,3)))/(max(results.full(:,3))-min(results.full(:,3)));
    srmsfR2 = (results.full(:,4)-min(results.full(:,4)))/(max(results.full(:,4))-min(results.full(:,4)));
    rmsfTimeMax = (results.full(:,5)-min(results.full(:,5)))/(max(results.full(:,5))-min(results.full(:,5)));
    srmsfTimeMax = (results.full(:,6)-min(results.full(:,6)))/(max(results.full(:,6))-min(results.full(:,6)));
    dfaGamma = (results.full(:,7)-min(results.full(:,7)))/(max(results.full(:,7))-min(results.full(:,7)));
    sdfaGamma = (results.full(:,8)-min(results.full(:,8)))/(max(results.full(:,8))-min(results.full(:,8)));
    msdBeta = (results.full(:,9)-min(results.full(:,9)))/(max(results.full(:,9))-min(results.full(:,9)));
    smsdBeta = (results.full(:,10)-min(results.full(:,10)))/(max(results.full(:,10))-min(results.full(:,10)));
    ApEn = (results.full(:,11)-min(results.full(:,11)))/(max(results.full(:,11))-min(results.full(:,11)));
    sApEn = (results.full(:,12)-min(results.full(:,12)))/(max(results.full(:,12))-min(results.full(:,12)));

elseif contains(type, 'standardized')

    % Standardized version of each metric

    rmsfAlpha = (results.full(:,1)-mean(results.full(:,1)))/std(results.full(:,1));
    srmsfAlpha = (results.full(:,2)-mean(results.full(:,2)))/std(results.full(:,2));
    rmsfR2 = (results.full(:,3)-mean(results.full(:,3)))/std(results.full(:,3));
    srmsfR2 = (results.full(:,4)-mean(results.full(:,4)))/std(results.full(:,4));
    rmsfTimeMax = (results.full(:,5)-mean(results.full(:,5)))/std(results.full(:,5));
    srmsfTimeMax = (results.full(:,6)-mean(results.full(:,6)))/std(results.full(:,6));
    dfaGamma = (results.full(:,7)-mean(results.full(:,7)))/std(results.full(:,7));
    sdfaGamma = (results.full(:,8)-mean(results.full(:,8)))/std(results.full(:,8));
    msdBeta = (results.full(:,9)-mean(results.full(:,9)))/std(results.full(:,9));
    smsdBeta = (results.full(:,10)-mean(results.full(:,10)))/std(results.full(:,10));
    ApEn = (results.full(:,11)-mean(results.full(:,11)))/std(results.full(:,11));
    sApEn = (results.full(:,12)-mean(results.full(:,12)))/std(results.full(:,12));

elseif contains(type, 'original')
    
    % Original, untouched version of each metric
    
    rmsfAlpha = results.full(:,1);
    srmsfAlpha = results.full(:,2);
    rmsfR2 = results.full(:,3);
    srmsfR2 = results.full(:,4);
    rmsfTimeMax = results.full(:,5);
    srmsfTimeMax = results.full(:,6);
    dfaGamma = results.full(:,7)-1; % substract 1 to turn it into a Hurst exponent
    sdfaGamma = results.full(:,8);
    msdBeta = results.full(:,9);
    smsdBeta = results.full(:,10);
    ApEn = results.full(:,11);
    sApEn = results.full(:,12);

end

% Scatter plot 3D and surface plus sphericity measurement
figures.full3DScatter.(type) = figure('Name',strcat('Scatter3D_allConditions_', type),'NumberTitle','off') ;
ax3D = gca;

scatter3(rmsfAlpha, dfaGamma, msdBeta, 6, rmsfTimeMax,'filled','MarkerEdgeColor','k','LineWidth',0.001) % timeMax as colormap
hold on

colormap(gca,"jet")
xlabel('RMSF Alpha')
ylabel('DFA Gamma')
zlabel('MSD Beta')
cb = colorbar;                                     
cb.Label.String = 'RMSF Max Time';
cb.Label.FontSize = 12;

[b3d, volume] = convhull(rmsfAlpha, dfaGamma, msdBeta);

trisurf(b3d, rmsfAlpha, dfaGamma, msdBeta,... %  timeMax can be colormap 
    'FaceColor','red','FaceAlpha', 0.044, 'LineWidth', 0.001,'LineStyle',':');
hold on


% Compute MVEE (minimal volume enclosing ellipsoid) of the 3D scatterplot

b3d = unique(b3d(:));
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([rmsfAlpha(b3d) ...
    dfaGamma(b3d) msdBeta(b3d)]', 0.00001);
m = Ellipse_plot(ax3D, ellipsoid_equation, ellipsoid_center, 1000);

set( m, 'FaceColor', 'g', 'FaceAlpha', .03, 'EdgeColor', 'none' );
% view( -70, 40 );
% axis vis3d equal;
camlight
lighting gouraud;

% Compute MVES (minimal volume enclosing sphere) for the 3D scatterplot

hold on
[R,C,Xb] = ExactMinBoundSphere3D([rmsfAlpha(b3d) ...
    dfaGamma(b3d) msdBeta(b3d)]);
VisualizeBoundSphere([rmsfAlpha(b3d) ...
    dfaGamma(b3d) msdBeta(b3d)],R,C,Xb,ax3D);

% 2D scatters

figures.full2DScatters.(type) = figure('Name',strcat('Scatter2D_allConditions', type),'NumberTitle','off') ;
ax2D = gca;
figures.full2DScatters.(type).Position(1:4) = [10 100 2250 500];
tiledlayout(2,6,'TileSpacing','tight','Padding','tight')

nexttile

scatter(rmsfAlpha, dfaGamma, 1, 'w','filled')
hold on
xlabel('RMSF Alpha')
ylabel('DFA Gamma')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(rmsfAlpha, dfaGamma, 0); % s=0 gives the 2D convex hull
plot(rmsfAlpha(bds),dfaGamma(bds),'w');
hold on
fit_ellipse(rmsfAlpha, dfaGamma, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([rmsfAlpha(bds)'; ...
    dfaGamma(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(rmsfAlpha, msdBeta, 1, 'w','filled')
hold on
xlabel('RMSF Alpha')
ylabel('MSD Beta')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(rmsfAlpha, msdBeta, 0);
plot(rmsfAlpha(bds),msdBeta(bds),'w');
hold on
fit_ellipse(rmsfAlpha, msdBeta, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([rmsfAlpha(bds)'; ...
    msdBeta(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(dfaGamma, msdBeta, 1, 'w','filled')
hold on
xlabel('DFA Gamma')
ylabel('MSD Beta')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(dfaGamma, msdBeta, 0);
plot(dfaGamma(bds),msdBeta(bds),'w');
hold on
fit_ellipse(dfaGamma, msdBeta, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([dfaGamma(bds)'; ...
    msdBeta(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(msdBeta, ApEn, 1, 'w','filled')
hold on
xlabel('MSD Beta')
ylabel('ApEn')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(msdBeta, ApEn, 0);
plot(msdBeta(bds), ApEn(bds),'w');
hold on
fit_ellipse(msdBeta, ApEn, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([msdBeta(bds)'; ...
    ApEn(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(rmsfAlpha, ApEn, 1, 'w','filled')
hold on
xlabel('RMSF Alpha')
ylabel('ApEn')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(rmsfAlpha, ApEn, 0);
plot(rmsfAlpha(bds), ApEn(bds),'w');
hold on
fit_ellipse(rmsfAlpha, ApEn, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([rmsfAlpha(bds)'; ...
    ApEn(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(dfaGamma, ApEn, 1, 'w','filled')
hold on
xlabel('DFA Gamma')
ylabel('ApEn')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(dfaGamma, ApEn, 0);
plot(dfaGamma(bds),ApEn(bds),'w');
hold on
fit_ellipse(dfaGamma, ApEn, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([dfaGamma(bds)'; ...
    ApEn(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(srmsfAlpha, sdfaGamma, 1, 'w','filled')
hold on
xlabel('Shuffled RMSF Alpha')
ylabel('Shuffled DFA Gamma')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(srmsfAlpha, sdfaGamma, 0);
plot(srmsfAlpha(bds), sdfaGamma(bds),'w');
hold on
fit_ellipse(srmsfAlpha, sdfaGamma, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([srmsfAlpha(bds)'; ...
    sdfaGamma(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(srmsfAlpha, smsdBeta, 1, 'w','filled')
hold on
xlabel('Shuffled RMSF Alpha')
ylabel('Shuffled MSD Beta')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(srmsfAlpha, smsdBeta, 0);
plot(srmsfAlpha(bds), smsdBeta(bds),'w');
hold on
fit_ellipse(srmsfAlpha, smsdBeta, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([srmsfAlpha(bds)'; ...
    smsdBeta(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(sdfaGamma, smsdBeta, 1, 'w','filled')
hold on
xlabel('Shuffled DFA Gamma')
ylabel('Shuffled MSD Beta')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(sdfaGamma, smsdBeta, 0);
plot(sdfaGamma(bds),smsdBeta(bds),'w');
hold on
fit_ellipse(sdfaGamma, smsdBeta, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([sdfaGamma(bds)'; ...
    smsdBeta(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(smsdBeta, sApEn, 1, 'w','filled')
hold on
xlabel('Shuffled MSD Beta')
ylabel('Shuffled ApEn')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(smsdBeta, sApEn, 0);
plot(smsdBeta(bds), sApEn(bds),'w');
hold on
fit_ellipse(smsdBeta, sApEn, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([smsdBeta(bds)'; ...
    sApEn(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(srmsfAlpha, sApEn, 1, 'w','filled')
hold on
xlabel('Shuffled RMSF Alpha')
ylabel('Shuffled ApEn')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(srmsfAlpha, sApEn, 0);
plot(srmsfAlpha(bds), sApEn(bds),'w');
hold on
fit_ellipse(srmsfAlpha, sApEn, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([srmsfAlpha(bds)'; ...
    sApEn(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

nexttile

scatter(sdfaGamma, sApEn, 1, 'w','filled')
hold on
xlabel('Shuffled DFA Gamma')
ylabel('Shuffled ApEn')
set(gca, 'Color','k', 'XColor','w', 'YColor','w')
set(gcf, 'Color','k')
bds = boundary(sdfaGamma, sApEn, 0);
plot(sdfaGamma(bds),sApEn(bds),'w');
hold on
fit_ellipse(sdfaGamma, sApEn, gcf);
[ellipsoid_equation,ellipsoid_center] = MinVolEllipse([sdfaGamma(bds)'; ...
    sApEn(bds)'], 0.00001);
Ellipse_plot(ax2D, ellipsoid_equation, ellipsoid_center, 1000);

% % Create textbox
% annotation(figures.full2DScatters,'textbox',...
%     [0.01 0.4 0 0.1],...
%     'Color',[0.9 0.9 0.9],...
%     'String',['Original',newline,'data'],...
%     'FontSize',22,...
%     'FitBoxToText','off');
% 
% % Create textbox
% annotation(figures.full2DScatters,'textbox',...
%     [0.01 0.85 0 0.1],...
%     'Color',[0.9 0.9 0.9],...
%     'String',['Shuffled',newline,'data'],...
%     'FontSize',22,...
%     'FitBoxToText','off');



end