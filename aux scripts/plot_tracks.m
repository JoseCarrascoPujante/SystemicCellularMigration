fig = figure('Position', [130 10 800 800]);
f_n = fieldnames(coordinates);
uicontrol('String','Next','Callback','uiresume(fig)',...
    'Position',[30 10 90 30]);
for f=11:length(f_n)
    for i = 1:size(coordinates.(f_n{f}).original_x,2)
        plot(coordinates.(f_n{f}).scaled_x(:,i), ...
            coordinates.(f_n{f}).scaled_y(:,i),'Color','b')
        hold on
        plot(coordinates.(f_n{f}).scaled_x(end,i), ...
            coordinates.(f_n{f}).scaled_y(end,i),...
            'ko','MarkerFaceColor','k','MarkerSize', 1.75)
        text(.05,.95,[f_n{f},'\_nº',num2str(i)],'Units','normalized')
        ax = gca;
        MaxX = max(abs(ax.XLim)+1);    MaxY = max(abs(ax.YLim)+1);
        axis equal
        if MaxX > MaxY
            axis([-MaxX MaxX -MaxY*(MaxX/MaxY) MaxY*(MaxX/MaxY)]);
        elseif MaxY > MaxX
            axis([-MaxX*(MaxY/MaxX) MaxX*(MaxY/MaxX) -MaxY MaxY]);
        elseif MaxY == MaxX
            axis([-MaxX MaxX -MaxY MaxY]);
        end
        uiwait(fig)
        hold off
    end
end