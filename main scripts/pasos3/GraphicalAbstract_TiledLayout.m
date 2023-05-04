
function figures = GraphicalAbstract_TiledLayout(field_names, results, figures)

% Plots systemic movement parameters in 2D and 3D scatterplots

% Group the 12 track-wise statistics in one table

results.full = [];

for index = 1:length(field_names)

        results.full = cat(1,results.full, results.(field_names{index}));
end

stat_names = ["rmsfAlpha" "srmsfAlpha" "rmsfR2" "srmsfR2" "rmsfTimeMax" ...
    "srmsfTimeMax" "dfaGamma" "sdfaGamma" "msdBeta" "smsdBeta" "ApEn" "sApEn"];
indexes = 1:length(stat_names);
pairs = nchoosek(indexes([1,7:2:end]),2);

% 2D scatters

precision = [[.50, .75, .90, .95]; [.5, .75, 1, 1.5]];
ellipseFitType = [{'%Confidence'},{'xSTD'}];
colr = ['w','c','r','m']; % there must be one color per precision level
for precisionType=1:height(precision)

    figures.full2DScatters = figure('Name',strcat('GraphicalAbstract_',...
        ellipseFitType{precisionType},'fit'),'NumberTitle','off') ;
    figures.full2DScatters.InvertHardcopy = 'off';
    figures.full2DScatters.Position(1:4) = [0 0 1300 900];
    tiledlayout(2,3,'TileSpacing','tight','Padding','tight')
    
    for ej=1:length(pairs)
        nexttile;
        disp(ej)
        if pairs(ej,1) == 7
            metric1 = cat(1,results.full(:,pairs(ej,1))-1, results.full(:,pairs(ej,1)+1));
        else 
            metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
        end
        if pairs(ej,2) == 7
            metric2 = cat(1,results.full(:,pairs(ej,2))-1, results.full(:,pairs(ej,2)+1));
        else
            metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
        end
        G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
        gscatter(metric1,metric2, G,'gy','..',1.75,'off') ;
        disp(strcat(stat_names(pairs(ej,1)),' & ',stat_names(pairs(ej,2))))
        set(gca, 'Color','k', 'XColor','w', 'YColor','w')
        set(gcf, 'Color','k') 
        for c = 1:length(precision(precisionType,:))
            disp(strcat(ellipseFitType{precisionType},': ',num2str(precision(precisionType,c))))
            hold on
            ellipse_switch(cat(2,metric1,metric2),G,precision(precisionType,c),colr(c),precisionType)
        end
        xlabel(stat_names(pairs(ej,1)))
        ylabel(stat_names(pairs(ej,2)))
    end
    
    h(1) = plot(nan, nan, 'go', 'MarkerFaceColor','g', 'MarkerSize', 10, 'DisplayName', 'Original');
    h(2) = plot(nan, nan, 'yo', 'MarkerFaceColor','y', 'MarkerSize', 10, 'DisplayName', 'Shuffled');
    h(3) = plot(nan, nan, 'w--', 'MarkerSize', 25, 'DisplayName',...
        strcat(num2str(precision(precisionType,1)), ellipseFitType{precisionType}));
    h(4) = plot(nan, nan, 'c--', 'MarkerSize', 25, 'DisplayName',...
        strcat(num2str(precision(precisionType,2)), ellipseFitType{precisionType}));
    h(5) = plot(nan, nan, 'r--', 'MarkerSize', 25, 'DisplayName',...
        strcat(num2str(precision(precisionType,3)), ellipseFitType{precisionType}));
    h(6) = plot(nan, nan, 'm--', 'MarkerSize', 25, 'DisplayName',...
        strcat(num2str(precision(precisionType,4)), ellipseFitType{precisionType}));
    leg=legend(h,Orientation='Horizontal',TextColor='w',FontSize=15);
    leg.Layout.Tile = 'north';
    
    versions = dir('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\GraphicalAbstract') ;
    count = 0 ;
    for v = 1:length(versions)
        if  contains(versions(v).name, 'GraphicalAbstract')
            count = count + 1 ;
        end
    end
    saveas(gcf,strcat('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\GraphicalAbstract\GraphicalAbstract(',num2str(count),')',ellipseFitType{precisionType},'.png'))
    hold off
end
end