function figure1 = fig1_12Panels4x3NoZooms(coordinates, destination_folder)

field_names = {
    'SinEstimuloProteus11_63';
    'SinEstimuloLeningradensis11_63';
    'SinEstimuloBorokensis23_44';
    'GalvanotaxisProteus11_63';
    'GalvanotaxisLeningradensis11_63';
    'GalvanotaxisBorokensis11_63';
    'QuimiotaxisProteus11_63';
    'QuimiotaxisLeningradensisVariosPpmm';
    'QuimiotaxisBorokensis23_44';
    'InduccionProteus11_63';    
    'InduccionLeningradensis11_63';    
    'InduccionBorokensis11_63'               
    };
figure('Position', [10 10 600 800])
figure1 = tiledlayout(4,3,'TileSpacing','tight','Padding','tight') ; % outer layout

%# Panels 1-4
for i = 1:length(field_names)
    ax = nexttile;
    box on
    ax.LineWidth = .4;
    hold on
    if ismember(i,[1,4,7,10])
        colr = 'k';
    elseif ismember(i,[2,5,8,11])
        colr = 'r';
    elseif ismember(i,[3,6,9,12])
        colr = 'b';
    end
    for track = 1:length(coordinates.(field_names{i}).scaled_x(1,:))
        %# Plot trajectory and 'ko' marker at its tip
        plot(coordinates.(field_names{i}).scaled_x(:,track),...
            coordinates.(field_names{i}).scaled_y(:,track), 'Color',colr);
        plot(coordinates.(field_names{i}).scaled_x(end,track),...
            coordinates.(field_names{i}).scaled_y(end,track),'o','MarkerFaceColor', ...
            colr,'MarkerEdgeColor',colr,'MarkerSize',1.5) ;
    end
    ax.FontSize = 7;
    MaxX = max(abs(ax.XLim))+1;   MaxY = max(abs(ax.YLim))+1;
    xline(0,'-','Alpha',1,'Color',[0 0 0]); % xline and yline cannot be sent to plot back
    yline(0,'-','Alpha',1,'Color',[0 0 0]);
    axis equal
    if MaxX > MaxY
        axis([-MaxX MaxX -MaxY*(MaxX/MaxY) MaxY*(MaxX/MaxY)]);
    elseif MaxY > MaxX
        axis([-MaxX*(MaxY/MaxX) MaxX*(MaxY/MaxX) -MaxY MaxY]);
    elseif MaxY == MaxX
    axis([-MaxX MaxX -MaxY MaxY]);
    end
end

% Export as .jpg and .svg
versions = dir(destination_folder) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig1_4x3')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig1_4x3 files found'))

exportgraphics(gcf,strcat(destination_folder, '\Fig1_4x3(',num2str(gabs),').jpg') ...
    ,"Resolution",400)
exportgraphics(gcf,strcat(destination_folder, '\Fig1_4x3(',num2str(gabs),').tiff') ...
    ,"Resolution",400)
exportgraphics(gcf,strcat(destination_folder, '\Fig1_4x3(',num2str(gabs),').pdf') ...
    ,'BackgroundColor','white', 'ContentType','vector')

end