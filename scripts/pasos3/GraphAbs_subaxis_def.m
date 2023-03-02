function figures = GraphAbs_subaxis_def(field_names, results, figures)
%# create dataset
results.full = [];
for index = 1:length(field_names)
        results.full = cat(1,results.full, results.(field_names{index}));
end
stat_names = ["RMSF\alpha" "sRMSF\alpha" "RMSFR2" "sRMSFR2" "RMSFtimeMax" ...
    "sRMSFTimeMax" "DFA\gamma" "sDFA\gamma" "MSD\beta" "sMSD\beta" "Approximate Entropy" "sApproximate Entropy"];
%# Presets
precision = .6827; %Equivalent to 1xSTD
ellipseFitType = '% confidence interval';

figures.full2DScatters = figure('Name',strcat('GraphicalAbstract_',...
    ellipseFitType,'fit'),'NumberTitle','off','Visible','on') ;
figures.full2DScatters.InvertHardcopy = 'off';
figures.full2DScatters.Position(1:4) = [0 0 700 700];

% pair stats to plot
indexes = 1:length(stat_names);
pairs = nchoosek(indexes([1,7:2:end]),2); % first column indexes go to X axis, 2nd column's to Y axis
pairs([5,6],:) = []; % remove two last combinations
results.full(:,7) = results.full(:,7)-1;

%# build axes positions
props = {'sh', 0.02, 'sv', 0.03, 'padding', 0.03 'margin', 0.03};
hBig = [subaxis(2,2,1, props{:}) subaxis(2,2,2, props{:}) subaxis(2,2,3, props{:}) subaxis(2,2,4, props{:})]; %# create subplots
posBig = get(hBig, 'Position');             %# record their positions
delete(hBig)                                %# delete them
posSmall{1} = [0.83 0.61 0.13 0.13];
posSmall{2} = [0.3 0.18 0.13 0.13];
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
    disp(strcat('Plot nº',num2str(ej)))
    metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
    metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
    G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
    gscatter(hAxB(ej),metric1,metric2, G,'gy','..',1.75,'off') ;
    disp(strcat(stat_names(pairs(ej,1)),' vs ',stat_names(pairs(ej,2))))
    hold(hAxB(ej),'on')
    ellipse_gscatter(hAxB(ej),cat(2,metric1,metric2),G,precision,'r')
    xlabel(hAxB(ej),stat_names(pairs(ej,1)))
    ylabel(hAxB(ej),stat_names(pairs(ej,2)))
end

%# plot on small axes
for ek=1:length(pairs)-1
    if ek == 2
        scatter(hAxS(ek),results.full(:,pairs(ek+1,1)),results.full(:,pairs(ek+1,2)),.75,'g','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek+1,1)),results.full(:,pairs(ek+1,2))),precision,'r')
        box(hAxS(ek),"on")
        xl= xlim(hAxS(ek));
        yl= ylim(hAxS(ek));
        rPos = [xl(1)-(((xl(2)-xl(1))*1.1)-((xl(2)-xl(1))))/2 yl(1)-(((yl(2)-yl(1))*7)-((yl(2)-yl(1))))/2 (xl(2)-xl(1))*1.1 (yl(2)-yl(1))*7];
        hold(hAxB(ek+1),'on')
        rectangle(hAxB(ek+1), 'Position',rPos, 'EdgeColor','r','FaceColor','none', 'LineWidth',.5);
        [Xor,Yor] = ds2nfu(hAxB(ek+1),rPos(1),rPos(2));
        [Xfr,Yfr] = ds2nfu(hAxB(ek+1),rPos(1)+rPos(3),rPos(2)+rPos(4));
        disp([Xor,Yor; ... %# rectangle bottom left
            Xfr,Yfr; ... %# rectangle top right
            posSmall{ek}(1),posSmall{ek}(2); ... %# small axis ek bottom left
            posSmall{ek}(1)+posSmall{ek}(3),posSmall{ek}(2)+posSmall{ek}(4)]) %# small axis ek top right
        annotation(figures.full2DScatters,'line',[Xfr posSmall{ek}(1)+posSmall{ek}(4)], [Yor posSmall{ek}(2)], ...
            'Color','w','LineStyle','--','LineWidth',.5);
        annotation(figures.full2DScatters,'line',[Xor posSmall{ek}(1)], [Yfr posSmall{ek}(2)+posSmall{ek}(4)], ...
            'Color','w','LineStyle','--','LineWidth',.5);
        alpha(hAxS(ek),1)
    else
        scatter(hAxS(ek),results.full(:,pairs(ek+1,1)+1),results.full(:,pairs(ek+1,2)+1),.75,'y','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek+1,1)+1),results.full(:,pairs(ek+1,2)+1)),precision,'r')
        box(hAxS(ek),"on")
        xl= xlim(hAxS(ek));
        yl= ylim(hAxS(ek));
        rPos = [xl(1)-(((xl(2)-xl(1))*1.1)-((xl(2)-xl(1))))/2 yl(1)-(((yl(2)-yl(1))*7)-((yl(2)-yl(1))))/2 (xl(2)-xl(1))*1.1 (yl(2)-yl(1))*7];
        hold(hAxB(ek+1),'on')
        rectangle(hAxB(ek+1),'Position', rPos,'EdgeColor','r','FaceColor','none','LineWidth',.5);
        [Xor,Yor] = ds2nfu(hAxB(ek+1),rPos(1),rPos(2));
        [Xfr,Yfr] = ds2nfu(hAxB(ek+1),rPos(1)+rPos(3),rPos(2)+rPos(4));
        disp([Xor,Yor; ... %# rectangle bottom left
            Xfr,Yfr; ... %# rectangle top right
            posSmall{ek}(1),posSmall{ek}(2); ... %# small axis ek bottom left
            posSmall{ek}(1)+posSmall{ek}(3),posSmall{ek}(2)+posSmall{ek}(4)]) %# small axis ek top right
        annotation(figures.full2DScatters,'line',[Xfr posSmall{ek}(1)+posSmall{ek}(4)], [Yor posSmall{ek}(2)], ...
            'Color','w','LineStyle','--','LineWidth',.5);
        annotation(figures.full2DScatters,'line',[Xor posSmall{ek}(1)], [Yfr posSmall{ek}(2)+posSmall{ek}(4)], ...
            'Color','w','LineStyle','--','LineWidth',.5);
        alpha(hAxS(ek),1)
    end
end

%# set axes properties
set(hAxB, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
set(hAxS, 'Color','none', 'XColor','r', 'YColor','r','LineWidth',.3,'FontSize',6, ...
    'XAxisLocation','bottom', 'YAxisLocation','left');
hAxS(1).XAxis.TickLabelColor = 'w';
hAxS(1).YAxis.TickLabelColor = 'w';
hAxS(2).XAxis.TickLabelColor = 'w';
hAxS(2).YAxis.TickLabelColor = 'w';
hAxS(3).XAxis.TickLabelColor = 'w';
hAxS(3).YAxis.TickLabelColor = 'w';
set(gcf, 'Color','k')

[h,objh] = legend(hAxB(1),'Original','Shuffled',strcat(num2str(precision),ellipseFitType), ...
    Orientation='Horizontal',TextColor='w',FontSize=12);
objhl = findobj(objh, 'type', 'line'); %// objects of legend of type line
set(objhl, 'Markersize', 30); %// set marker size as desired

h.Position(1:2)=[0.15,0.95];

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