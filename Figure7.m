%% Figure 7

clear 
close all
set(groot,'defaultFigurePaperPositionMode','manual')
load 'numerical_results.mat' results

fig = figure('Visible','off','Position', [0 0 1070 650]);

% create dataset
field_names = fieldnames(results);
results.full = [];
for index = 1:length(field_names)
    results.full = cat(1,results.full, results.(field_names{index})); % field name order does not matter here
end

results.full(:,[5,6]) = results.full(:,[5,6])/120; % convert frames to minutes

% Presets, use non-kynetic metrics only
stat_names = {'RMSF\alpha', 'sRMSF\alpha', 'RMSF_R2', 'sRMSF_R2', 'RMSFCorrelationTime', ...
    'sRMSFCorrelationTime', 'DFA\gamma', 'sDFA\gamma', 'MSD\beta', 'sMSD\beta',...
    'Approximate Entropy', 'sApproximate Entropy'} ;
conf = 68.27;                              % set to either a confidence or STD value. 68.27% CI equates to 1xSTD
ellipseFitType = 'x STD';  % set to either STD or confidence %

% pair stats
indexes = 1:length(stat_names);
pairs = nchoosek(indexes([1,7:2:end]),2); % choose metrics to compare
pairs([3,5],:) = pairs([5,3],:);
pairs([4,5],:) = pairs([5,4],:);

% build axes positions
props = {'sh', 0.02, 'sv', 0.03, 'padding', 0.03 'margin', 0.03};
hBig = [subaxis(2,3,1, props{:}) subaxis(2,3,2, props{:})...
    subaxis(2,3,3, props{:}) subaxis(2,3,4, props{:}) ...
    subaxis(2,3,5, props{:}) subaxis(2,3,6, props{:})]; % create subplots
posBig = get(hBig, 'Position');                         % record their positions
delete(hBig)                                            % delete them
ratio = fig.Position(4)/fig.Position(3);
% Assign small axes to corresponding big axes by following the natural tile
% ordinal index progression: when [row(1),column(end)] is reached -> [row(2),column(1)]
posSmall{1} = [0.55 0.615 0.13*ratio 0.13];
posSmall{2} = [0.89 0.68 0.13*ratio 0.13];
posSmall{3} = [0.24 0.18 0.13*ratio 0.13];
posSmall{4} = [0.438 0.163 0.13*ratio 0.13];
posSmall{5} = [0.76 0.34 0.13*ratio 0.13];
posSmall{6} = [0.89 0.18 0.13*ratio 0.13];

% Create axes (big/small)
hAxB(1) = axes('Position',posBig{1});
hAxB(2) = axes('Position',posBig{2});
hAxB(3) = axes('Position',posBig{3});
hAxB(4) = axes('Position',posBig{4});
hAxB(5) = axes('Position',posBig{5});
hAxB(6) = axes('Position',posBig{6});
hAxS(1) = axes('Position',posSmall{1});
hAxS(2) = axes('Position',posSmall{2});
hAxS(3) = axes('Position',posSmall{3});
hAxS(4) = axes('Position',posSmall{4});
hAxS(5) = axes('Position',posSmall{5});
hAxS(6) = axes('Position',posSmall{6});

% Plot on big axes
for ej=1:length(pairs)
    disp(strcat('Plot nº',num2str(ej),': ',stat_names{pairs(ej,1)},'_vs_',stat_names{pairs(ej,2)}))
    metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
    metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
    G = [ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
    gscatter(hAxB(ej),metric1,metric2, G,[.1 .5 0;0 0 1],'..',2.9,'off',...
        stat_names{pairs(ej,1)},stat_names{pairs(ej,2)})
    hold(hAxB(ej),'on')
    ellipse_gscatter(hAxB(ej),cat(2,metric1,metric2),G,conf,'r')
    axis(hAxB(ej),'padded')
    % xlabel(hAxB(ej),stat_names{pairs(ej,1)})
    % ylabel(hAxB(ej),stat_names{pairs(ej,2)})
end

% Plot on small axes, use ek+n when small panels start at big panel #n, use 
% (ek[+n])+1 when plotting shuffled data
for ek=1:length(pairs) % for # small axes do...
    if ek == 2 || ek == 3
        scatter(hAxS(ek),results.full(:,pairs(ek+1,1)),results.full(:,pairs(ek+1,2)),.7,[.1 .5 0],'filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek+1,1)),results.full(:,pairs(ek+1,2))),conf,'r')
        axis(hAxS(ek),'padded')
    elseif ek == 6 % last small axis corresponds to original (green) data
        scatter(hAxS(ek),results.full(:,pairs(ek,1)),results.full(:,pairs(ek,2)),.7,[.1 .5 0],'filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek,1)),results.full(:,pairs(ek,2))),conf,'r')
        axis(hAxS(ek),'padded')
    else
        scatter(hAxS(ek),results.full(:,pairs(ek+1,1)+1),results.full(:,pairs(ek+1,2)+1),.7,'b','filled','o') ;
        hold(hAxS(ek),'on')
        ellipse_scatter(hAxS(ek),cat(2,results.full(:,pairs(ek+1,1)+1),results.full(:,pairs(ek+1,2)+1)),conf,'r')
        axis(hAxS(ek),'padded')
    end
    box(hAxS(ek),'on')
    xl = xlim(hAxS(ek));
    yl = ylim(hAxS(ek));
    rPos = [xl(1)-(((xl(2)-xl(1))+.02)-((xl(2)-xl(1))))/2 yl(1)-(((yl(2)-yl(1))+.04)-((yl(2)-yl(1))))/2 (xl(2)-xl(1))+.02 (yl(2)-yl(1))+.04];
    if ek ~= 6
        hold(hAxB(ek+1),'on')
        rectangle(hAxB(ek+1),'Position', rPos,'EdgeColor','k','FaceColor','none','LineWidth',.5);
        [Xor,Yor] = ds2nfu(hAxB(ek+1),rPos(1),rPos(2));
        [Xfr,Yfr] = ds2nfu(hAxB(ek+1),rPos(1)+rPos(3),rPos(2)+rPos(4));
    else
        hold(hAxB(ek),'on')
        rectangle(hAxB(ek),'Position', rPos,'EdgeColor','k','FaceColor','none','LineWidth',.5);
        [Xor,Yor] = ds2nfu(hAxB(ek),rPos(1),rPos(2));
        [Xfr,Yfr] = ds2nfu(hAxB(ek),rPos(1)+rPos(3),rPos(2)+rPos(4));
    end
%     disp([Xor,Yor; ...                                                     % rectangle's bottom left
%         Xfr,Yfr; ...                                                       % rectangle's top right
%         posSmall{ek}(1),posSmall{ek}(2); ...                               % smallAxis{ek}'s bottom left
%         posSmall{ek}(1)+posSmall{ek}(3),posSmall{ek}(2)+posSmall{ek}(4)]); % smallAxis{ek}'s top right
    annotation(fig,'line',[Xfr posSmall{ek}(1)+posSmall{ek}(4)*ratio], [Yor posSmall{ek}(2)], ...
        'Color','k','LineStyle','--','LineWidth',.5);
    annotation(fig,'line',[Xor posSmall{ek}(1)], [Yfr posSmall{ek}(2)+posSmall{ek}(4)], ...
        'Color','k','LineStyle','--','LineWidth',.5);
    % alpha(hAxS(ek),1) % make small axes' background invisible
end

%# set axes properties
set(hAxB, 'Color','w', 'XColor','k', 'YColor','k','FontSize',9.5)
set(hAxS, 'Color','w', 'XColor','k', 'YColor','k','LineWidth',.5,'FontSize',6, ...
    'XAxisLocation','bottom', 'YAxisLocation','left');

[h,objh] = legend(hAxB(1),'Systemic Cell Migrations','Non-Systemic Cell Migrations',...
    '', Orientation='Horizontal',TextColor='k',FontSize=12);
objhl = findobj(objh, 'type', 'line'); %// objects of legend #1 of type line
set(objhl, 'Markersize', 30); %// set marker size as desired

[h2,objh2] = legend(hAxB(4),'','',strcat(num2str(1),' ',ellipseFitType), Orientation='Vertical', ...
    TextColor='k',FontSize=8.5);
objhl2 = findobj(objh2, 'type', 'line'); %// objects of legend #2 of type line
set(objhl2, 'Markersize', 15); %// set marker size as desired

h.Position(1:2)=[0.25,0.96];
h2.Position(1:2)=[0.4,0];

%% Export as vector graphics svg

load 'coordinates.mat' destination_folder 

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 1 ;
for v = 1:length(versions)
    if contains(versions(v).name, 'Fig7'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs-1),' Fig7 files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig7(',num2str(gabs),')'),'svg')
