
clear
close all
load 'coordinates.mat' coordinates destination_folder

set(groot,'defaultFigurePaperPositionMode','manual')

fig = figure('NumberTitle','off','Visible','off','Position',[0 0 1000 1500]);

%% Layouts

layout0 = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','tight') ;
layout1 = tiledlayout(layout0,2,2,'TileSpacing','compact','Padding','tight') ;
layout1.Layout.Tile = 1;
layout1.Layout.TileSpan = [2 1];
layout2 = tiledlayout(layout0,3,3,'TileSpacing','none','Padding','tight');
layout2.Layout.Tile = 3;
layouts = struct;

%% Panels 1-4

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

panels = {'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'};
indexes = {
{%No stimulus
{1,3,7,8,9,10,11,12,13,14,15,16,18,19,20,24,25,28,33,39}%Proteus
{1,4,6,7,9,10,14,16,17,18,20,21,22,23,24,26,27,28,29,31}%Leningradensis
{1,2,4,6,12,13,14,15,17,18,20,22,24,27,28,29,30,32,34,35}%Borokensis
};
{%Galvanotaxis
{1,2,3,4,5,7,8,11,14,15,16,17,19,21,22,23,24,25,26,27}%Proteus
{1:20}%Leningradensis
{1:20}%Borokensis
};
{%Chemotaxis
{2,3,5,6,7,13,14,16,17,18,19,20,21,22,25,26,27,28,29,38}%Proteus
{1,2,3,4,6,12,14,18,19,20,21,22,23,26,27,28,29,42,45,59}%Leningradensis
{8,9,11,13,15,19,22,25,29,31,32,33,34,44,45,47,50,51,53,55}%Borokensis
};
{%Induction
{1,4,5,7,8,9,11,16,20,22,24,25,26,28,30,32,33,35,36,70}%Proteus
{1,8,18,24,30,38,39,40,41,42,43,44,48,49,50,51,52,53,54,57}%Leningradensis
{2,32,37,38,40,42,44,47,50,52,57,59,52,68,69,72,75,76,77,78}%Borokensis
}
};
dim = [[.1 .3];[.1 .2];[.1 .1]]; % text label locations on left side of chart
dim2 = [[.75 .3];[.75 .2];[.75 .1]]; % text label locations on right side of chart
for i = 1:length(panels)
    layouts.(strcat('x',num2str(i))) = tiledlayout(layout1,10,3,'TileSpacing','tight','Padding','tight');
    layouts.(strcat('x',num2str(i))).Layout.Tile = i;
    speciesScenarioIdx = find(contains(field_names,panels{i}));
    ax = nexttile(layouts.(strcat('x',num2str(i))),[8,3]);
    for species = [length(speciesScenarioIdx),1,2] % specify plotting zorder
        listx=[];
        if species == 1
            colr = [0, 0, 0];
        elseif species == 2
            colr = [1, 0, 0];
        elseif species == 3
            colr = [0, 0, 1];
        end
        for track = cell2mat(indexes{i}{species})
            %# Plot trajectory and a 'ko' marker at its tip      
            plot(coordinates.(field_names{speciesScenarioIdx(species)}).scaled_x(:,track), ...
                coordinates.(field_names{speciesScenarioIdx(species)}).scaled_y(:,track), 'Color', colr) ;
            hold on;
            plot(coordinates.(field_names{speciesScenarioIdx(species)}).scaled_x(end,track), ...
                coordinates.(field_names{speciesScenarioIdx(species)}).scaled_y(end,track), ...
                'o', 'MarkerFaceColor', colr, 'MarkerEdgeColor', colr, 'MarkerSize', 1.75) ;

            listx = [listx,coordinates.(field_names{speciesScenarioIdx(species)}).original_x(end,track) ...
                - coordinates.(field_names{speciesScenarioIdx(species)}).original_x(1,track)];
        end
        if i ~= 1
            text(ax,dim(species,1),dim(species,2),...
                [num2str(round((sum(listx<0)/length(listx))*100)),'%'],...
                'Units','normalized','Color',colr,...
                'HorizontalAlignment','left','FontSize',12)
            text(ax,dim2(species,1),dim2(species,2),...
                [num2str(round((sum(listx>0)/length(listx))*100)),'%'],...
                'Units','normalized','Color',colr,...
                'HorizontalAlignment','left','FontSize',12)
        end
    end
    axis tight
    axis square
    text(ax,.1,.9,["N=60(20-\color{red}20\color{black}-\color{blue}20)","\color{black}t=30'"],'Units','normalized',...
        'HorizontalAlignment','left','FontSize',12)
    ax.FontSize = 7;
    MaxX = max(abs(ax.XLim));    MaxY = max(abs(ax.YLim));
    % Add x-line
    x = 0; 
    xl = plot([x,x],ylim(ax), 'k-', 'LineWidth', .5);
    % Add y-line
    y = 0; 
    yl = plot(xlim(ax), [y, y], 'k-', 'LineWidth', .5);
    % Send x and y lines to the bottom of the stack
    uistack([xl,yl],'bottom')
    % Update the x and y line bounds any time the axes limits change
    ax.XAxis.LimitsChangedFcn = @(ruler,~)set(xl, 'YData', ylim(ancestor(ruler,'axes')));
    ax.YAxis.LimitsChangedFcn = @(ruler,~)set(yl, 'XData', xlim(ancestor(ruler,'axes')));
    if MaxX > MaxY
        axis([-MaxX MaxX -MaxY*(MaxX/MaxY) MaxY*(MaxX/MaxY)]);
    elseif MaxY > MaxX
        axis([-MaxX*(MaxY/MaxX) MaxX*(MaxY/MaxX) -MaxY MaxY]);
    elseif MaxY == MaxX
        axis([-MaxX MaxX -MaxY MaxY]);
    end
    cn=0;
    for sp = 1:3
        if sp == 1
            colr = [0, 0, 0];
        elseif sp == 2
            colr = [1, 0, 0];
        elseif sp == 3
            colr = [0, 0, 1];
        end
        
        if i == 1
            pax = polaraxes(layouts.(strcat('x',num2str(i))));
            pax.Layout.Tile = 25+cn; %Initialize tile
            pax.Layout.TileSpan = [2 1];
            thetaIn360 = mod(coordinates.(field_names{speciesScenarioIdx(sp)}).theta(...
                end,:) + 2*pi, 2*pi);
            pol = polarhistogram(pax,thetaIn360,'Normalization','probability','LineWidth',1,...
                'FaceColor','None','EdgeColor',colr,'DisplayStyle','bar','BinEdges',linspace(-pi, pi, 9));
            rlim([0 .39])
            x = pol.BinEdges ;
            y = pol.Values ;
            text(pax,x(1:end-1)+pi/9,zeros(length(y),1) + .3,...
                strcat(num2str(round(y'*100,1)),'%'),'vert','bottom','horiz',...
                'center','FontSize',4); %Add labels as percentages
            rticks([])
            thetaticks([])
            thetaticklabels([])
            thetaticks(0:45:360)
            pax.GridColor = [0,0,0];
            pax.GridAlpha = 1;
            % pax.GridLineWidth = 2;
        else
            axh = nexttile(layouts.(strcat('x',num2str(i))),25+cn,[2,1]); %Initialize tile
            b = histogram(axh,cos(coordinates.(field_names{speciesScenarioIdx(sp)}).theta(...
                end,:)),10,"BinEdges",-1:.2:1,'FaceColor',colr);
            xticks([-1 0 1])
            xlim([-1,1])
            yticks([min(b.Values) max(b.Values)])
            ylim([0 max(b.Values)])
            axh.FontSize = 5;
            axis square
        end
        cn=cn+1;
    end
end

%% Panel 5

ax1 = nexttile(layout2,[3 3]);

% import and plot trajectory
original_x = readmatrix('C:\Users\pc\Desktop\Doctorado\Papers publicados\mov_sist\Tracking panel E fig1.xlsx', "Range", "A:A");
original_y = readmatrix('C:\Users\pc\Desktop\Doctorado\Papers publicados\mov_sist\Tracking panel E fig1.xlsx', "Range", "B:B");

% Do not center this trajectory
scaled_x = original_x/11.63;
scaled_y = original_y/11.63;

plot(scaled_x, scaled_y, 'Color', 'b') ;

ax1.FontSize = 13;
xlabel('x(mm)','FontSize',17)
ylabel('y(mm)','FontSize',17)
ax1.XAxis.TickLength = [0 0];
ax1.YAxis.TickLength = [0 0];
axis image
axis padded

% Define bounds of the first rectangle
left = 20.4;
bottom = 15.4;
width = 0.2;
height = 0.23;

% Display the first rectangle
hold(ax1,'on');
rectangle('Position',[left bottom width height], ...
    'EdgeColor','red','LineWidth',0.75);

% Create first axes for zoomed-in view
ax2 = axes(layout2);
ax2.Layout.Tile = 1;
plot(scaled_x, scaled_y, 'Color', 'b')

% Adjust first axis limits and remove ticks
axis padded
axis equal
ax2.XLim = [left left+width];
ax2.YLim = [bottom bottom+height];
ax2.XAxis.TickLength = [0 0];
ax2.YAxis.TickLength = [0 0];

% Set other properties on the first axes
ax2.Box = 'on';
ax2.XAxis.Color = 'k';
ax2.YAxis.Color = 'k';
ax2.FontSize = 6;
xlabel('x(mm)','FontSize',7)
ylabel('y(mm)','FontSize',7)

% Define bounds of the second rectangle
left = 20.422;
bottom = 15.54;
width = .07;
height = .08;

% Display the second rectangle
hold(ax2,'on');
rectangle('Position',[left bottom width height], ...
    'EdgeColor','red','LineWidth',0.5);

% Create second axes for zoomed-in view
ax3 = axes(layout2);
ax3.Layout.Tile = 7;
plot(scaled_x, scaled_y, 'Color', 'b')

% Adjust second axis limits and remove ticks
axis image
axis padded
ax3.XLim = [left left+width];
ax3.YLim = [bottom bottom+height];
ax3.XAxis.TickLength = [0 0];
ax3.YAxis.TickLength = [0 0];

% Set other properties on the second axes
ax3.Box = 'on';
ax3.XAxis.Color = 'k';
ax3.YAxis.Color = 'k';
ax3.FontSize = 6;
xlabel('x(mm)','FontSize',7)
ylabel('y(mm)','FontSize',7)




%% Export as jpg, tiff and vector graphics pdf

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig1'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs-1),' Fig1 files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),')'),'svg')
% exportgraphics(fig,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),').jpg') ...
%     ,"Resolution",600)
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),').tiff') ...
%     ,"Resolution",600)
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),').pdf') ...
%     ,'BackgroundColor','white', 'ContentType','vector')
