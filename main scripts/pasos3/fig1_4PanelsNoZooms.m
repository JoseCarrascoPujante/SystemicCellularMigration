function [fig1,cosines] = fig1_4PanelsNoZooms(coordinates, destination_folder)

field_names = fieldnames(coordinates) ;
fig1 = figure('Position',[10 10 650 1000]);
layout0 = tiledlayout(2,2,'TileSpacing','tight','Padding','tight') ;
layouts = struct;
cosines = struct;
%% Panels 1-4

panels = {'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'};
indexes = {
{%No stimulus
{1,2,4,6,12,13,14,15,17,18,20,22,24,27,28,29,30,32,34,35}%Borokensis
{1,3,7,8,9,10,11,12,13,14,15,16,18,19,20,24,28,33,38,39}%Proteus
{1,4,6,7,9,10,14,16,17,18,20,21,22,23,24,26,27,28,29,31}%Leningradensis
};
{%Galvanotaxis
{1:20}%Borokensis
{1,2,3,4,5,6,8,9,12,15,16,17,18,20,22,23,24,25,26,27}%Proteus
{1:20}%Leningradensis
};
{%Chemotaxis
{8,9,11,13,15,19,22,25,29,31,32,33,34,44,45,47,50,51,53,55}%Borokensis
{2,3,5,6,7,13,14,16,17,18,19,20,21,22,25,26,27,28,29,38}%Proteus
{1,2,3,4,6,12,14,18,19,20,21,22,23,26,27,28,29,42,45,59}%Leningradensis
};
{%Induction
{21:40}%Borokensis
{1,2,3,4,5,7,8,9,11,20,22,24,25,26,28,32,33,35,36,43}%Proteus
{1,8,18,24,30,38,39,40,41,42,43,44,48,49,50,51,52,53,54,57}%Leningradensis
}
};
dim = [[.1 .1];[.1 .3];[.1 .2]];
dim2 = [[.75 .1];[.75 .3];[.75 .2]];
for i = 1:length(panels)
    layouts.(strcat('x',num2str(i))) = tiledlayout(layout0,6,3,'TileSpacing','tight','Padding','tight');
    layouts.(strcat('x',num2str(i))).Layout.Tile = i;
    cond = find(contains(field_names,panels{i}));
    ax = nexttile(layouts.(strcat('x',num2str(i))),[4,3]);
    for species = 1:length(cond)
        listx=[];
        if species == 1
            colr = [0, 0, 1];
        elseif species == 2
            colr = [0, 0, 0];
        elseif species == 3
            colr = [1, 0, 0];
        end
        for track = cell2mat(indexes{i}{species})
            %# Plot trajectory and a 'ko' marker at its tip      
            plot(coordinates.(field_names{cond(species)}).scaled_x(:,track), ...
                coordinates.(field_names{cond(species)}).scaled_y(:,track), 'Color', colr) ;
            hold on;
            plot(coordinates.(field_names{cond(species)}).scaled_x(end,track), ...
                coordinates.(field_names{cond(species)}).scaled_y(end,track), ...
                'o', 'MarkerFaceColor', colr, 'MarkerEdgeColor', colr, 'MarkerSize', 1.75) ;
%             rng('default')
            listx = [listx,coordinates.(field_names{cond(species)}).original_x(end,track) ...
                - coordinates.(field_names{cond(species)}).original_x(1,track)];
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
    axis equal
    text(ax,.1,.9,['N=60(20-\color{red}20\color{black}-\color{blue}20)',"\color{black}t=30'"],'Units','normalized',...
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
            colr = [0, 0, 1];
        elseif sp == 2
            colr = [0, 0, 0];
        elseif sp == 3
            colr = [1, 0, 0];
        end
        
        if i == 1
            pax = polaraxes(layouts.(strcat('x',num2str(i))));
            pax.Layout.Tile = 13+cn; %Initialize tile
            pax.Layout.TileSpan = [2 1];
            thetaIn360 = mod(coordinates.(field_names{cond(sp)}).theta(...
                end,cell2mat(indexes{i}{sp})) + 2*pi, 2*pi);
            pol = polarhistogram(pax,thetaIn360,'Normalization','probability','LineWidth',1,...
                'FaceColor','None','DisplayStyle','bar','BinEdges',linspace(-pi, pi, 9));
            rlim([0 .35])
            x = pol.BinEdges ;
            y = pol.Values ;
            text(pax,x(1:end-1)+pi/9,zeros(length(y),1) + .3,...
                strcat(num2str(round(y'*100,0)),'%'),'vert','bottom','horiz',...
                'center','FontSize',7.5); %Add labels as percentages
            rticks([])
            thetaticks([])
            thetaticklabels([])
            thetaticks(0:45:360)
            pax.GridColor = [0,0,0];
            pax.GridAlpha = 1;
            % pax.GridLineWidth = 2;
        else
            ax2 = nexttile(layouts.(strcat('x',num2str(i))),13+cn,[2,1]); %Initialize tile
            b = histogram(ax2,cos(coordinates.(field_names{cond(sp)}).theta(...
                end,cell2mat(indexes{i}{sp}))),10,"BinEdges",[-1:.2:1],'FaceColor',colr);
            xticks([-1 0 1])
            xlim([-1,1])
            yticks([min(b.Values) max(b.Values)])
            ylim([min(b.Values) max(b.Values)])              
        end
        cn=cn+1;
    end
end

%% Export as jpg, tiff and vector graphics pdf

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures')) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig1'+wildcardPattern+'.jpg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig1 files found'))

exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),').jpg') ...
    ,"Resolution",600)
exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),').tiff') ...
    ,"Resolution",400)
exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig1(',num2str(gabs),').pdf') ...
    ,'BackgroundColor','white', 'ContentType','vector')

end