% Initialize condition track figure
figures.('InduccionBorokensis11_63').tracks = figure('Name',strcat('Tracks_Induccion Borokensis 11.63'),...
    'Visible','on','NumberTitle','off') ;
hTracks = gca;

% Tracks loop
A = 78 ;
for i = 1:A

% Plot trajectory and place black dot marker at its tip      
    plot(hTracks,coordinates.('InduccionBorokensis11_63').scaled_x(:,i),...
        coordinates.('InduccionBorokensis11_63').scaled_y(:,i), 'Color', [0, 0, 0]) ;
    hold on;
    plot(hTracks,coordinates.('InduccionBorokensis11_63').scaled_x(end,i),...
        coordinates.('InduccionBorokensis11_63').scaled_y(end,i),...
        'ko',  'MarkerFaceColor',  [0, 0, 0], 'MarkerSize', 2) ;
    
end
% Adjust track plot axes' proportions
hold on
box on
MaxX = max(abs(hTracks.XLim))+1;   MaxY = max(abs(hTracks.YLim))+1;
% Add x-line
x = 0; 
xl = plot([x,x],ylim(hTracks), 'k-', 'LineWidth', .5);
% Add y-line
y = 0; 
yl = plot(xlim(hTracks), [y, y], 'k-', 'LineWidth', .5);
% Send x and y lines to the bottom of the stack
uistack([xl,yl],'bottom')
% Update the x and y line bounds any time the axes limits change
hTracks.XAxis.LimitsChangedFcn = @(ruler,~)set(xl, 'YData', ylim(ancestor(ruler,'axes')));
hTracks.YAxis.LimitsChangedFcn = @(ruler,~)set(yl, 'XData', xlim(ancestor(ruler,'axes')));
axis equal
if MaxX > MaxY
    axis([-MaxX MaxX -MaxY*(MaxX/MaxY) MaxY*(MaxX/MaxY)]);
elseif MaxY > MaxX
    axis([-MaxX*(MaxY/MaxX) MaxX*(MaxY/MaxX) -MaxY MaxY]);
elseif MaxY == MaxX
    axis([-MaxX MaxX -MaxY MaxY]);
end
hold off