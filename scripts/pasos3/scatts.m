
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
figures.full2DScatters.(type).Position(1:4) = [0 0 1300 900];
tiledlayout(3,5,'TileSpacing','tight','Padding','tight')

for ej=1:length(pairs)
    % if dfaGamma substract 1 else do not
    nexttile
    disp(ej)
    metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
    metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
    G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
    gscatter(metric1,metric2, G,'gy','..',1.75,'off') ;
    disp(strcat(stat_names(pairs(ej,1)),' & ',stat_names(pairs(ej,2))))
    set(gca, 'Color','k', 'XColor','w', 'YColor','w')
    set(gcf, 'Color','k') 
    conf_intervals = [.25, .5, 1, 1.5];
    colr = ['w','c','r','m'];
    for c = 1:length(conf_intervals)
        hold on
        simpler_ellipse(cat(2,metric1,metric2),G,conf_intervals(c),colr(c))
    end
%     fit_ellipse(results.full(:,pairs(ej,1)), results.full(:,pairs(ej,2)), gcf);
%     hold on
%     fit_ellipse(results.full(:,pairs(ej,1)+1), results.full(:,pairs(ej,2)+1), gcf);
    xlabel(stat_names(pairs(ej,1)))
    ylabel(stat_names(pairs(ej,2)))
end

h(1) = plot(nan, nan, 'go', 'MarkerFaceColor','g', 'MarkerSize', 10, 'DisplayName', 'Original');
h(2) = plot(nan, nan, 'yo', 'MarkerFaceColor','y', 'MarkerSize', 10, 'DisplayName', 'Shuffled');
h(3) = plot(nan, nan, 'm--', 'MarkerSize', 25, 'DisplayName', '0.25 STD');
h(4) = plot(nan, nan, 'c--', 'MarkerSize', 25, 'DisplayName', '0.5 STD');
h(5) = plot(nan, nan, 'r--', 'MarkerSize', 25, 'DisplayName', '1 STD');
h(6) = plot(nan, nan, 'w--', 'MarkerSize', 25, 'DisplayName', '1.5 STD');
leg=legend(h,Orientation='Horizontal',TextColor='w',FontSize=15);
leg.Layout.Tile = 'north';

versions = dir('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\Tracks 3600 frames, matlab files and tables\GraphicalAbstract') ;
count = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'GraphicalAbstract')
        count = count + 1 ;
    end
end
saveas(gcf,strcat('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\Tracks 3600 frames, matlab files and tables\GraphicalAbstract\GraphicalAbstract(',num2str(count),').png'))

end