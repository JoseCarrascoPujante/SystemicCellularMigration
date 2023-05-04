function figures = GraphAbs_subaxis_def(field_names, results, figures, stat_names)
%# create dataset
results.full = [];
for index = 1:length(field_names)
        results.full = cat(1,results.full, results.(field_names{index}));
end

%# Presets
conf = 68.27; %# set to either STD or confidence %
ellipseFitType = '% confidence interval';  %# set to either STD or confidence %

figures.GraphicalAbstract = figure('Name',strcat('GraphicalAbstract_',...
    ellipseFitType,'fit'),'NumberTitle','off','Visible','on') ;
figures.GraphicalAbstract.InvertHardcopy = 'off';
figures.GraphicalAbstract.Position(1:4) = [0 0 735 735];

% pair stats
indexes = 1:length(stat_names);
pairs = nchoosek(indexes([1,7:2:end]),2); % choose metrics to compare
pairs([5,6],:) = []; % leave best 4 metrics to display on the GraphicalAbstract
results.full(:,7) = results.full(:,7)-1; %convert DFAgamma values to Hurst exponent values by substracting 1 to them

%# build axes positions
props = {'sh', 0.02, 'sv', 0.03, 'padding', 0.03 'margin', 0.03};
hBig = [subaxis(2,2,1, props{:}) subaxis(2,2,2, props{:}) subaxis(2,2,3, props{:}) subaxis(2,2,4, props{:})]; %# create subplots
posBig = get(hBig, 'Position');             %# record their positions
delete(hBig)                                %# delete them
posSmall{1} = [0.83 0.61 0.13 0.13];
posSmall{2} = [0.339 0.18 0.13 0.13];
posSmall{3} = [0.83 0.16 0.13 0.13];

%# Create axes (big/small)
hAxB(1) = axes('Position',posBig{1});
hAxB(2) = axes('Position',posBig{2});
hAxB(3) = axes('Position',posBig{3});
hAxB(4) = axes('Position',posBig{4});
hAxS(1) = axes('Position',posSmall{1});
hAxS(2) = axes('Position',posSmall{2});
hAxS(3) = axes('Position',posSmall{3});

%# Plot on big axes
for ej=1:length(pairs)
    disp(strcat('Plot nº',num2str(ej),': ',stat_names{pairs(ej,1)},'_vs_',stat_names{pairs(ej,2)}))
    metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
    metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
    G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
    if ej == 1 || ej == 3
        gscatter(hAxB(ej),metric1,metric2, G,'gy','..',2.3,'off')
    else
        gscatter(hAxB(ej),metric1,metric2, G,'gy','..',1.75,'off') ;
    end
    hold(hAxB(ej),'on')
    ellipse_gscatter(hAxB(ej),cat(2,metric1,metric2),G,conf,'r')
    xlabel(hAxB(ej),stat_names{pairs(ej,1)})
    ylabel(hAxB(ej),stat_names{pairs(ej,2)})
end

%# plot on small axes
for ek=1:length(pairs)-1 % for number of small axes do...
    if ek == 2
        scatter(hAxS(ek),results.full(:,pairs(ek+1,1)),results.full(:,pairs(ek+1,2)),.7,'g','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek+1,1)),results.full(:,pairs(ek+1,2))),conf,'r')
    else
        scatter(hAxS(ek),results.full(:,pairs(ek+1,1)+1),results.full(:,pairs(ek+1,2)+1),.7,'y','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek+1,1)+1),results.full(:,pairs(ek+1,2)+1)),conf,'r')
    end
    box(hAxS(ek),"on")
    xl= xlim(hAxS(ek));
    yl= ylim(hAxS(ek));
    rPos = [xl(1)-(((xl(2)-xl(1))*1.1)-((xl(2)-xl(1))))/2 yl(1)-(((yl(2)-yl(1))*7)-((yl(2)-yl(1))))/2 (xl(2)-xl(1))*1.05 (yl(2)-yl(1))*7];
    hold(hAxB(ek+1),'on')
    rectangle(hAxB(ek+1),'Position', rPos,'EdgeColor','r','FaceColor','none','LineWidth',.5);
    [Xor,Yor] = ds2nfu(hAxB(ek+1),rPos(1),rPos(2));
    [Xfr,Yfr] = ds2nfu(hAxB(ek+1),rPos(1)+rPos(3),rPos(2)+rPos(4));
%     disp([Xor,Yor; ... %# rectangle bottom left
%         Xfr,Yfr; ... %# rectangle top right
%         posSmall{ek}(1),posSmall{ek}(2); ... %# smallAxis{ek} bottom left
%         posSmall{ek}(1)+posSmall{ek}(3),posSmall{ek}(2)+posSmall{ek}(4)]); %# smallAxis{ek} top right
    annotation(figures.GraphicalAbstract,'line',[Xfr posSmall{ek}(1)+posSmall{ek}(4)], [Yor posSmall{ek}(2)], ...
        'Color','w','LineStyle','--','LineWidth',.5);
    annotation(figures.GraphicalAbstract,'line',[Xor posSmall{ek}(1)], [Yfr posSmall{ek}(2)+posSmall{ek}(4)], ...
        'Color','w','LineStyle','--','LineWidth',.5);
    alpha(hAxS(ek),1)
end

%# set axes properties
set(hAxB, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
set(hAxS, 'Color','none', 'XColor','r', 'YColor','r','LineWidth',.5,'FontSize',6, ...
    'XAxisLocation','bottom', 'YAxisLocation','left');
hAxS(1).XAxis.TickLabelColor = 'w';
hAxS(1).YAxis.TickLabelColor = 'w';
hAxS(2).XAxis.TickLabelColor = 'w';
hAxS(2).YAxis.TickLabelColor = 'w';
hAxS(3).XAxis.TickLabelColor = 'w';
hAxS(3).YAxis.TickLabelColor = 'w';
set(gcf, 'Color','k')

[h,objh] = legend(hAxB(1),'Systemic Cell Migrations','Non-Systemic Cell Migrations',...
    '', Orientation='Horizontal',TextColor='w',FontSize=12);
objhl = findobj(objh, 'type', 'line'); %// objects of legend #1 of type line
set(objhl, 'Markersize', 30); %// set marker size as desired

[h2,objh2] = legend(hAxB(4),'','',strcat(num2str(conf),ellipseFitType), Orientation='Vertical', ...
    TextColor='w',FontSize=8.5);
objhl2 = findobj(objh2, 'type', 'line'); %// objects of legend #2 of type line
set(objhl2, 'Markersize', 15); %// set marker size as desired

h.Position(1:2)=[0.15,0.953];
h2.Position(1:2)=[0.35,0];

%# save as .png and .svg

versions = dir('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\GraphicalAbstract') ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'GraphicalAbstract')
        gabs = gabs + 1 ;
    end
end
disp(strcat(num2str(gabs),' Graphical Abstract files found'))
exportgraphics(gcf,strcat(['E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\Graphical' ...
    'Abstract\GraphicalAbstract('],num2str(gabs),')','.jpg'),"Resolution",600,'BackgroundColor','k')
saveas(gcf,strcat('E:\Doctorado\Amebas\Pavlov 2 y 3\Resultados movimiento sistémico\GraphicalAbstract\GraphicalAbstract(',num2str(gabs+1),')','.svg'))

end