% Plots cell movement statistics in 2D scatterplots
function figures = GraphAbs_def(field_names, results, figures)

% Merge all 12 stats in one struct

results.full = [];

for index = 1:length(field_names)
        results.full = cat(1,results.full, results.(field_names{index}));
end

stat_names = ["rmsf\alpha" "srmsf\alpha" "rmsfR2" "srmsfR2" "rmsfTimeMax" ...
    "srmsfTimeMax" "dfa\gamma" "sdfa\gamma" "msd\beta" "smsd\beta" "ApEn" "sApEn"];

% Presets
precision = .65;
ellipseFitType = '%confidence';

figures.full2DScatters = figure('Name',strcat('GraphicalAbstract_',...
    ellipseFitType,'fit'),'NumberTitle','off') ;
figures.full2DScatters.InvertHardcopy = 'off';
figures.full2DScatters.Position(1:4) = [0 0 625 625];
T = tiledlayout(2,2,'TileSpacing','tight','Padding','tight'); % Outer layout

% pair stats to plot
indexes = 1:length(stat_names);
pairs = nchoosek(indexes([1,7:2:end]),2); % first column indexes go to X axis, 2nd column's to Y axis
pairs([5,6],:) = []; % remove two last combinations
hurstDFA = results.full(:,7)-1;

%# Plot
count = 0;
for ej=1:length(pairs)
    count = count + 1;
    disp(ej)
    if pairs(ej,1) == 7
        metric1 = cat(1,hurstDFA, results.full(:,pairs(ej,1)+1));
    else 
        metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
    end
    if pairs(ej,2) == 7
        metric2 = cat(1,hurstDFA, results.full(:,pairs(ej,2)+1));
    else
        metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
    end
    G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];

    switch true
        case count == 2
            t2=tiledlayout(T,'flow','TileSpacing','tight','Padding','tight'); % Inner layout
            t2.Layout.Tile = 2;
            nexttile(t2,[5 5])
            gscatter(metric1,metric2, G,'gy','..',1.75,'off') ;
            disp(strcat(stat_names(pairs(ej,1)),' & ',stat_names(pairs(ej,2))))
            hold on
            ellipse_gscatter(gca,cat(2,metric1,metric2),G,precision,'r')
            xlabel(stat_names(pairs(ej,1)))
            ylabel(stat_names(pairs(ej,2)))
            set(gca, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
            nexttile(t2,'east');
            scatter(results.full(:,pairs(ej,1)+1),results.full(:,pairs(ej,2)+1),1.75,'y','filled','o') ;
            hold on
            ellipse_scatter(gca,cat(2,results.full(:,pairs(ej,1)+1),results.full(:,pairs(ej,2)+1)),precision,'r')
            box on
            ax = gca;
            ax.LineWidth = .3;
            set(gca, 'Color','none', 'XColor','r', 'YColor','r','FontSize',8)
            ax.XAxis.TickLabelColor = 'w';
            ax.YAxis.TickLabelColor = 'w';
            alpha(gca,.5)
        case count == 3
            t3=tiledlayout(T,'flow','TileSpacing','tight','Padding','tight'); % Inner layout
            t3.Layout.Tile = 3;
            nexttile(t3,[5 5])
            gscatter(metric1,metric2, G,'gy','..',1.75,'off') ;
            disp(strcat(stat_names(pairs(ej,1)),' & ',stat_names(pairs(ej,2))))
            hold on
            ellipse_gscatter(gca,cat(2,metric1,metric2),G,precision,'r')
            xlabel(stat_names(pairs(ej,1)))
            ylabel(stat_names(pairs(ej,2)))
            set(gca, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
            nexttile(t3,'east')
            scatter(results.full(:,pairs(ej,1)),results.full(:,pairs(ej,2)),1.75,'g','filled','o') ;
            hold on
            ellipse_scatter(gca,cat(2,results.full(:,pairs(ej,1)),results.full(:,pairs(ej,2))),precision,'r')
            box on
            ax = gca;
            ax.LineWidth = .3;
            set(gca, 'Color','none', 'XColor','r', 'YColor','r','FontSize',8)
            ax.XAxis.TickLabelColor = 'w';
            ax.YAxis.TickLabelColor = 'w';
            alpha(gca,.5)
        case count == 4
            t4=tiledlayout(T,'flow','TileSpacing','tight','Padding','tight'); % Inner layout
            t4.Layout.Tile = 4;
            nexttile(t4,[5 5])
            gscatter(metric1,metric2, G,'gy','..',1.75,'off') ;
            disp(strcat(stat_names(pairs(ej,1)),' vs ',stat_names(pairs(ej,2))))
            hold on
            ellipse_gscatter(gca,cat(2,metric1,metric2),G,precision,'r')
            xlabel(stat_names(pairs(ej,1)))
            ylabel(stat_names(pairs(ej,2)))
            set(gca, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
            nexttile(t4,'east');
            scatter(results.full(:,pairs(ej,1)+1),results.full(:,pairs(ej,2)+1),1.75,'y','filled','o') ;
            hold on
            ellipse_scatter(gca,cat(2,results.full(:,pairs(ej,1)+1),results.full(:,pairs(ej,2)+1)),precision,'r')
            box on
            ax = gca;
            ax.LineWidth = .3;
            set(gca, 'Color','none', 'XColor','r', 'YColor','r','FontSize',8)
            ax.XAxis.TickLabelColor = 'w';
            ax.YAxis.TickLabelColor = 'w';
            alpha(gca,.5)
        otherwise
            nexttile(T)
            gscatter(metric1,metric2, G,'gy','..',1.75,'off') ;
            disp(strcat(stat_names(pairs(ej,1)),' & ',stat_names(pairs(ej,2))))
            hold on
            ellipse_gscatter(gca,cat(2,metric1,metric2),G,precision,'r')
            xlabel(stat_names(pairs(ej,1)))
            ylabel(stat_names(pairs(ej,2)))
            set(gca, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
    end
end

set(gcf, 'Color','k')
nexttile(T,1)

h(1) = plot(nan, nan, 'go', 'MarkerFaceColor','g', 'MarkerSize', 10, 'DisplayName', 'Original');
h(2) = plot(nan, nan, 'yo', 'MarkerFaceColor','y', 'MarkerSize', 10, 'DisplayName', 'Shuffled');
h(3) = plot(nan, nan, 'r--', 'MarkerSize', 25, 'DisplayName',...
    strcat(num2str(precision), ellipseFitType));
leg=legend(h,Orientation='Horizontal',TextColor='w',FontSize=15);
leg.Layout.Tile = 'north';

versions = dir('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\GraphicalAbstract') ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'GraphicalAbstract')
        gabs = gabs + 1 ;
    end
end
disp(strcat(num2str(gabs),' Graphical Abstract files found'))
exportgraphics(gcf,strcat(['E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\Graphical' ...
    'Abstract\GraphicalAbstract('],num2str(gabs),')',ellipseFitType,'.png'),"Resolution",600,'BackgroundColor','k')
saveas(gcf,strcat('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\GraphicalAbstract\GraphicalAbstract(',num2str(gabs),')',ellipseFitType,'.svg'))
end