function figure1 = fig1_version2(coordinates)
N = 40 ;
field_names = fieldnames(coordinates) ;
figure
figure1 = tiledlayout(3,1,'TileSpacing','tight','Padding','tight') ; % outer layout

%# Panels 1-4
panels = {'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'};
fig1_1 = tiledlayout(figure1,2,2,'TileSpacing','tight','Padding','tight'); % inner layout
fig1_1.Layout.Tile = 1;
fig1_1.Layout.TileSpan = [2 1];
pans = {'a','b','c','d'};
for i = 1:length(panels)
    fig1_2.(pans{i}) = tiledlayout(fig1_1,3,1,'TileSpacing','none','Padding','none'); % innermost layouts
    fig1_2.(pans{i}).Layout.Tile = i;
    idx = find(contains(field_names,panels{i}));
    for species=1:length(idx)
        nexttile(fig1_2.(pans{i}))
        hold on
        colr = [0, 0, 0];
        colr(species) = colr(species)+1;
        for track = 1:N
            %# Plot trajectory and 'ko' marker at its tip
            plot(coordinates.(field_names{idx(species)}).scaled_x(:,track),...
                coordinates.(field_names{idx(species)}).scaled_y(:,track), 'Color', colr) ;
            plot(coordinates.(field_names{idx(species)}).scaled_x(end,track),...
                coordinates.(field_names{idx(species)}).scaled_y(end,track),...
                'ko',  'MarkerFaceColor',  colr, 'MarkerSize', 2) ;
        end
        ax = gca;
        set(ax,'YTickLabel',[],'YTick',[]);
        ax.YAxisLocation = 'origin';
        axis([-15 15 -15 15])
        if species ~= 3
            ax.XAxis.Visible = 'off';
        end
    end
    if pans{i} == 'a' || pans{i} == 'b'
        ax.XAxis.Visible = 'off';
    end
end

%# Panels 5-6
fig1_3=tiledlayout(figure1,1,3,'TileSpacing','none','Padding','tight'); % Inner layout
fig1_3.Layout.Tile = 3;
nexttile(fig1_3)
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
ax2 = axes(t);
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

end