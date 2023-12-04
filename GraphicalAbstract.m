
clear
close all
load 'coordinates.mat' destination_folder
load numerical_results.mat results
fig = figure('Position',[0 0 650 650],'NumberTitle','off','Visible','off') ;
fig.InvertHardcopy = 'off';  % To saveas white axes

%% Presets

% create dataset
field_names = ...
    {'SinEstimuloProteus11_63'
    'GalvanotaxisProteus11_63'
    'QuimiotaxisProteus11_63'
    'InduccionProteus11_63'
    'SinEstimuloLeningradensis11_63'
    'GalvanotaxisLeningradensis11_63'
    'QuimiotaxisLeningradensisVariosPpmm'
    'InduccionLeningradensis11_63'
    'SinEstimuloBorokensis23_44'
    'GalvanotaxisBorokensis11_63'
    'QuimiotaxisBorokensis23_44'
    'InduccionBorokensis11_63'
    };

full = [];
for index = 1:length(field_names)
        full = cat(1,full, results.(field_names{index}));
end

full(:,[5,6]) = full(:,[5,6])/120; % convert frames to minutes

stat_names = {'RMSF\alpha', 'sRMSF\alpha', 'RMSF_R2', 'sRMSF_R2', 'RMSFCorrelationTime', ...
    'sRMSFCorrelationTime', 'DFA\gamma', 'sDFA\gamma', 'MSD\beta', 'sMSD\beta',...
    'Approximate Entropy', 'sApproximate Entropy'} ;
conf = 68.27; %# set to either a confidence or STD value. 68.27% CI equates to 1xSTD
ellipseFitType = '% confidence interval';  %# set to either STD or confidence %

% pair stats
indexes = 1:length(stat_names);
pairs = nchoosek(indexes([1,7:2:end]),2); % choose metrics to compare
pairs([5,6],:) = [];                      % remove the two least worst metrics from the GraphicalAbstract

%# build axes positions
props = {'sh', 0.02, 'sv', 0.03, 'padding', 0.03 'margin', 0.03};
hBig = [subaxis(2,2,1, props{:})...         %# create big axes
    subaxis(2,2,2, props{:})...
    subaxis(2,2,3, props{:})...
    subaxis(2,2,4, props{:})];
posBig = get(hBig, 'Position');             %# record their positions
delete(hBig)                                %# delete them
% establish small axes' position
posSmall{1} = [0.85 0.61 0.13 0.13];
posSmall{2} = [0.335 0.27 0.13 0.13];
posSmall{3} = [0.61 0.27 0.13 0.13];

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
    metric1 = cat(1,full(:,pairs(ej,1)), full(:,pairs(ej,1)+1));
    metric2 = cat(1,full(:,pairs(ej,2)), full(:,pairs(ej,2)+1));
    G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
    if ej == 1
        gscatter(hAxB(ej),metric1,metric2, G,'gy','..',2.4,'off')
    elseif ej == 3
        gscatter(hAxB(ej),metric1,metric2, G,'gy','..',2.2,'off')
    else
        gscatter(hAxB(ej),metric1,metric2, G,'gy','..',2.1,'off') ;
    end
    hold(hAxB(ej),'on')
    ellipse_gscatter(hAxB(ej),cat(2,metric1,metric2),G,conf,'r')
    xlabel(hAxB(ej),stat_names{pairs(ej,1)})
    ylabel(hAxB(ej),stat_names{pairs(ej,2)})
end

%# plot on small axes
for ek=1:length(pairs)-1 % for number of small axes do...
    if ek == 2
        scatter(hAxS(ek),full(:,pairs(ek+1,1)),full(:,pairs(ek+1,2)),.7,'g','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,full(:,pairs(ek+1,1)),full(:,pairs(ek+1,2))),conf, 'r')
    else
        scatter(hAxS(ek),full(:,pairs(ek+1,1)+1),full(:,pairs(ek+1,2)+1),.7,'y','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,full(:,pairs(ek+1,1)+1),full(:,pairs(ek+1,2)+1)),conf, 'r')
    end
    box(hAxS(ek),'on')
    xl = xlim(hAxS(ek));
    yl = ylim(hAxS(ek));
    rPos = [xl(1)-(((xl(2)-xl(1))*1.05)-((xl(2)-xl(1))))/2 yl(1)-(((yl(2)-yl(1))+.07)-((yl(2)-yl(1))))/2 (xl(2)-xl(1))*1.05 (yl(2)-yl(1))+.07];
    hold(hAxB(ek+1),'on')
    rectangle(hAxB(ek+1),'Position', rPos,'EdgeColor', 'r','FaceColor','none','LineWidth',.5);
    [Xor,Yor] = ds2nfu(hAxB(ek+1),rPos(1),rPos(2));
    [Xfr,Yfr] = ds2nfu(hAxB(ek+1),rPos(1)+rPos(3),rPos(2)+rPos(4));
%     disp([Xor,Yor; ... %# rectangle bottom left
%         Xfr,Yfr; ... %# rectangle top right
%         posSmall{ek}(1),posSmall{ek}(2); ... %# smallAxis{ek} bottom left
%         posSmall{ek}(1)+posSmall{ek}(3),posSmall{ek}(2)+posSmall{ek}(4)]); %# smallAxis{ek} top right
    annotation(fig,'line',[Xfr posSmall{ek}(1)+posSmall{ek}(4)], [Yor posSmall{ek}(2)], ...
        'Color','w','LineStyle','--','LineWidth',.5);
    annotation(fig,'line',[Xor posSmall{ek}(1)], [Yfr posSmall{ek}(2)+posSmall{ek}(4)], ...
        'Color','w','LineStyle','--','LineWidth',.5);
    alpha(hAxS(ek),1) % make small axes' background invisible
end

%# set axes properties
set(hAxB, 'Color','k', 'XColor','w', 'YColor','w','FontSize',9.5,'FontWeight','bold')
set(hAxS, 'Color','none', 'XColor', 'r', 'YColor','r','LineWidth',.5,'FontSize',6, ...
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

[h2,objh2] = legend(hAxB(4),'','','1x STD', Orientation='Vertical', ...
    TextColor='w',FontSize=8.5);
objhl2 = findobj(objh2, 'type', 'line'); %// objects of legend #2 of type line
set(objhl2, 'Markersize', 15); %// set marker size as desired

h.Position(1:2)=[0.15,0.953];
h2.Position(1:2)=[0.35,0];

%# save as .svg
versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'GraphicalAbstract'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs-1),' Graphical Abstract files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig, strcat(destination_folder,'\Figures\GraphicalAbstract(',num2str(gabs+1),').svg'))
