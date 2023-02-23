
function figures = scatts(field_names, results, figures, type)

% Plots systemic movement parameters in 2D and 3D scatterplots

% Group the 12 track-wise statistics in one table

results.full = [];

for index = 1:length(field_names)

        results.full = cat(1,results.full, results.(field_names{index}));
end

stat_names = ["rmsfAlpha" "srmsfAlpha" "rmsfR2" "srmsfR2" "rmsfTimeMax" ...
    "srmsfTimeMax" "dfaGamma" "sdfaGamma" "msdBeta" "smsdBeta" "ApEn" "sApEn"];
indexes = 1:length(stat_names);
pairs = nchoosek(indexes(1:2:end),2);

% 2D scatters

figures.full2DScatters.(type) = figure('Name',strcat('Scatter2D_allConditions', type),'NumberTitle','off') ;
figures.full2DScatters.(type).InvertHardcopy = 'off';
figures.full2DScatters.(type).Position(1:4) = [10 10 1250 850];
tiledlayout(3,5,'TileSpacing','tight','Padding','tight')

for ej=1:length(pairs)
    % if dfaGamma substract 1 else do not
    nexttile
    disp(ej)
    scatter(results.full(:,pairs(ej,1)), results.full(:,pairs(ej,2)), 1, 'g','filled') ;
    disp(strcat(stat_names(pairs(ej,1)),' & ',stat_names(pairs(ej,2))))
    set(gca, 'Color','k', 'XColor','w', 'YColor','w')
    set(gcf, 'Color','k') 
    hold on
    scatter(results.full(:,pairs(ej,1)+1), results.full(:,pairs(ej,2)+1), 1, 'y','filled') ;
    
    hold on
    EllipseDirectFit(cat(2, results.full(:,pairs(ej,1)), results.full(:,pairs(ej,2))));
    hold on
    EllipseDirectFit(cat(2, results.full(:,pairs(ej,1)+1), results.full(:,pairs(ej,2)+1)));
    xlabel(stat_names(pairs(ej,1)))
    ylabel(stat_names(pairs(ej,2)))
%     xlim([-inf inf]) 
%     ylim([-inf inf])
end

leg = legend('Original', 'Shuffled', 'Best-fit ellipse','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.TextColor = 'w';
leg.FontSize = 17 ;

versions = dir('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\Tracks 3600 frames, matlab files and tables\GraphicalAbstract') ;
count = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'GraphicalAbstract')
        count = count + 1 ;
    end
end
saveas(gcf,strcat('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\Tracks 3600 frames, matlab files and tables\GraphicalAbstract\GraphicalAbstract(',num2str(count),').png'))

end