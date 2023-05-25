function figure1 = fig1_4PanelsPlusZooms(coordinates, destination_folder)

field_names = fieldnames(coordinates) ;
figure('Position', [770 10 500 1500])
figure1 = tiledlayout(3,1,'TileSpacing','tight','Padding','tight') ;

%% Panels 1-4
fig1_1 = tiledlayout(figure1,2,2,'TileSpacing','tight','Padding','tight');
fig1_1.Layout.Tile = 1;
fig1_1.Layout.TileSpan = [2 1];
panels = {'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'};
indexes = {
{
{1 2}
{1}
{1 2}
};
{
{1:20}
{1}
{1:20}
};
{
{1:20}
{1}
{1}
};
{
{1:20}
{1}
{37 38 39 40 41 42 43}
}
};
for i = 1:length(panels)
    idx = find(contains(field_names,panels{i}));
    ax = nexttile(fig1_1);
    for species = 1:length(idx)
        if species == 1
            colr = [0, 0, 1];
        elseif species == 2
            colr = [0, 0, 0];
        elseif species == 3
            colr = [1, 0, 0];
        end
%         rng()
%         for track = randperm(length(coordinates.(field_names{idx(species)}).scaled_x(1,:)), 20)
        for track = cell2mat(indexes{i}{species})
            %# Plot trajectory and a 'ko' marker at its tip      
            plot(coordinates.(field_names{idx(species)}).scaled_x(:,track),...
                coordinates.(field_names{idx(species)}).scaled_y(:,track), 'Color', colr) ;
            hold on;
            plot(coordinates.(field_names{idx(species)}).scaled_x(end,track),...
                coordinates.(field_names{idx(species)}).scaled_y(end,track),...
                'ko', 'MarkerFaceColor', colr, 'MarkerEdgeColor', colr, 'MarkerSize', 1.75) ;
%             rng('default')
        end
    end
%     set(ax,'YTickLabel',[],'YTick',[]);
    ax.FontSize = 7;
    MaxX = max(abs(ax.XLim));    MaxY = max(abs(ax.YLim)+.75);
    axis([-MaxX MaxX -MaxY MaxY]);
    xline(0,'-','Alpha',1,'Color',[0 0 0]); % xline and yline cannot be sent to back
    yline(0,'-','Alpha',1,'Color',[0 0 0]);
%     daspect([1 1 1])
end

%% Panels 5-6
fig1_2 = tiledlayout(figure1,1,2,'TileSpacing','tight','Padding','tight');
fig1_2.Layout.Tile = 3;

%% Panel 5
fig1_3 = tiledlayout(fig1_2,3,3,'TileSpacing','none','Padding','tight');
fig1_3.Layout.Tile = 1;
ax1 = nexttile(fig1_3,[3 3]);
scenario1 = 'SinEstimuloProteus11_63';
Ameba_wo_stim = 9;
plot(coordinates.(scenario1).scaled_x(:,Ameba_wo_stim), ...
    coordinates.(scenario1).scaled_y(:,Ameba_wo_stim), 'Color', 'b') ;

% Define bounds of the rectangle
left = 2.97;
bottom = 0.2;
width = 0.12;
height = 0.12;

% Display the rectangle
hold(ax1,'on');
rectangle('Position',[left bottom width height], ...
    'EdgeColor','red','LineWidth',0.75);

% Create axes for zoomed-in view
ax2 = axes(fig1_3);
ax2.Layout.Tile = 7;
plot(coordinates.(scenario1).scaled_x(:,Ameba_wo_stim), ...
    coordinates.(scenario1).scaled_y(:,Ameba_wo_stim), 'Color', 'b') 

% Adjust axis limits and remove ticks
ax2.XLim = [left left+width];
ax2.YLim = [bottom bottom+height];
% ax2.XTick = [];
% ax2.YTick = [];

% Set other properties on the axes
ax2.Box = 'on';
ax2.XAxis.Color = 'k';
ax2.YAxis.Color = 'k';
ax1.FontSize = 7;
ax2.FontSize = 6;
% title(ax1,'Trajectory in the absence of stimuli','Color','red','FontSize',9);

%% Panel 6
fig1_3_1 = tiledlayout(fig1_2,3,3,'TileSpacing','none','Padding','tight');
fig1_3_1.Layout.Tile = 2;
ax1 = nexttile(fig1_3_1,[3 3]);
scenario2 = 'SinEstimuloLeningradensis11_63';
Ameba_w_stim = 14;
plot(coordinates.(scenario2).scaled_x(:,Ameba_w_stim), ...
    coordinates.(scenario2).scaled_y(:,Ameba_w_stim), 'Color', 'b') ;

% Define bounds of the rectangle
left = -1.35;
bottom = -2.9;
width = 0.3;
height = 0.15;

% Display the rectangle
hold(ax1,'on');
rectangle('Position',[left bottom width height], ...
    'EdgeColor','red','LineWidth',0.75);

% Create axes for zoomed-in view
ax2 = axes(fig1_3_1);
ax2.Layout.Tile = 9;
plot(coordinates.(scenario2).scaled_x(:,Ameba_w_stim), ...
    coordinates.(scenario2).scaled_y(:,Ameba_w_stim), 'Color', 'b') 
% Adjust axis limits
ax2.XLim = [left left+width];
ax2.YLim = [bottom bottom+height];
% ax2.XTick = [];
% ax2.YTick = [];

% Set other properties on the axes
ax2.Box = 'on';
ax2.XAxis.Color = 'k';
ax2.YAxis.Color = 'k';
ax1.FontSize = 7;
ax2.FontSize = 6;

% Export as .jpg and .svg
versions = dir(destination_folder) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig1')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig1 files found'))

exportgraphics(gcf,strcat(destination_folder, '\Fig1(',num2str(gabs),').jpg') ...
    ,"Resolution",400)
exportgraphics(gcf,strcat(destination_folder, '\Fig1(',num2str(gabs),').tiff') ...
    ,"Resolution",400)
exportgraphics(gcf,strcat(destination_folder, '\Fig1(',num2str(gabs),').pdf') ...
    ,'BackgroundColor','white', 'ContentType','vector')

end