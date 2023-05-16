function figure1 = fig1_version1(coordinates, destination_folder)
field_names = fieldnames(coordinates) ;
figure
figure1 = tiledlayout(3,1,'TileSpacing','tight','Padding','tight') ;

%# Panels 1-4
fig1_1 = tiledlayout(figure1,2,2,'TileSpacing','tight','Padding','tight');
fig1_1.Layout.Tile = 1;
fig1_1.Layout.TileSpan = [2 1];
panels = {'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'};
for i = 1:length(panels)
    idx = find(contains(field_names,panels{i}));
    ax = nexttile(fig1_1);
    for species = 1:length(idx)
        if species == 1
            colr = [0, 0, 0];
        elseif species == 2
            colr = [0, 0, 1];
        elseif species == 3
            colr = [1, 0, 0];
        end
        rng('default')
        for track = randperm(length(coordinates.(field_names{idx(species)}).scaled_x(1,:)), 30)
            %# Plot trajectory and a 'ko' marker at its tip      
            plot(coordinates.(field_names{idx(species)}).scaled_x(:,track),...
                coordinates.(field_names{idx(species)}).scaled_y(:,track), 'Color', colr) ;
            hold on;
            plot(coordinates.(field_names{idx(species)}).scaled_x(end,track),...
                coordinates.(field_names{idx(species)}).scaled_y(end,track),...
                'ko',  'MarkerFaceColor',  colr, 'MarkerSize', 2) ;
            rng('default')
        end
    end
%     set(ax,'YTickLabel',[],'YTick',[]);
    MaxX = max(abs(ax.XLim));    MaxY = max(abs(ax.YLim));
    axis([-MaxX MaxX -MaxY MaxY]);
    xline(0,'-');  yline(0,'-');
end

% Panels 5-6
fig1_2 = tiledlayout(figure1,1,2,'TileSpacing','none','Padding','tight');
fig1_2.Layout.Tile = 3;

% Panel 5
fig1_3_1 = tiledlayout(fig1_2,3,3,'TileSpacing','none','Padding','tight');
fig1_3_1.Layout.Tile = 1;
ax1 = nexttile(fig1_3_1,[3 3]);
plot(coordinates.('SinEstimuloProteus11_63').scaled_x(:,1), ...
    coordinates.('SinEstimuloProteus11_63').scaled_y(:,1), 'Color', 'b') ;

% Define bounds of the rectangle
left = -0.5;
bottom = -0.7;
width = 0.4;
height = 0.4;

% Display the rectangle
hold(ax1,'on');
r = rectangle('Position',[left bottom width height], ...
    'EdgeColor','red','LineWidth',1.5);

% Set properties on the axes
ax1.FontSize = 10;
ax1.XLim = [-4.5 4.5];
ax1.YLim = [-4.5 4.5];
% grid(ax1,'on')

% Create axes for zoomed-in view
ax2 = axes(fig1_3_1);
ax2.Layout.Tile = 3;
scatter(ax2,x,y,10,'filled');

% Adjust axis limits and remove ticks
ax2.XLim = [left left+width];
ax2.YLim = [bottom bottom+height];
ax2.XTick = [];
ax2.YTick = [];

% Set other properties on the axes
ax2.Box = 'on';
ax2.XAxis.Color = 'red';
ax2.YAxis.Color = 'red';
title(ax2,'100x Magnification','Color','red');

% Export as .jpg and .svg
versions = dir(destination_folder) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig1')
        gabs = gabs + 1 ;
    end
end
disp(strcat(num2str(gabs),' Fig1 files found'))
exportgraphics(gcf,strcat(destination_folder, 'Fig1(',num2str(gabs),')','.jpg'),"Resolution",600,'BackgroundColor','k')
saveas(gcf,strcat(destination_folder, 'Fig1(',num2str(gabs+1),')','.svg'))

end