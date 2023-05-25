topLevelFolder = uigetdir('C:\') ; %choose folder containing the excel files
AllDirs = dir(fullfile(topLevelFolder, '**\*\*\*.xlsx')) ;
AllsubFolderNames = {AllDirs.folder} ;
UsefulSubFolderNames = unique(AllsubFolderNames, 'sorted') ;
for f = 1:length(UsefulSubFolderNames)
    cd(UsefulSubFolderNames{f})
    files = dir('*.xlsx') ;
    [~,name2,name3] = fileparts(pwd) ;
    condition = strcat(name2, name3) ;
    conditionValidName = matlab.lang.makeValidName(condition) ;
    files = natsortfiles(files) ;
    figures.(conditionValidName).tracks = figure('Name',strcat('Tracks_',condition),...
        'Visible','off','NumberTitle','off') ;
    hTracks = gca;
    A = coordinates.(conditionValidName).scaled_x(1,:) ; 
    for i = 1:length(A)  
        plot(hTracks,coordinates.(conditionValidName).scaled_x(:,i),...
            coordinates.(conditionValidName).scaled_y(:,i), 'Color', [0, 0, 0]) ;
        hold on;
        plot(hTracks,coordinates.(conditionValidName).scaled_x(end,i),...
            coordinates.(conditionValidName).scaled_y(end,i),...
            'ko',  'MarkerFaceColor',  [0, 0, 0], 'MarkerSize', 2) ;
    end
    hold on
    box on
    MaxX = max(abs(hTracks.XLim))+1;   MaxY = max(abs(hTracks.YLim))+1;
    xline(0,'-','Alpha',1,'Color',[0 0 0]); % xline and yline cannot be sent to plot's back
    yline(0,'-','Alpha',1,'Color',[0 0 0]);
    axis equal
    if MaxX > MaxY
        axis([-MaxX MaxX -MaxY*(MaxX/MaxY) MaxY*(MaxX/MaxY)]);
    elseif MaxY > MaxX
        axis([-MaxX*(MaxY/MaxX) MaxX*(MaxY/MaxX) -MaxY MaxY]);
    elseif MaxY == MaxX
        axis([-MaxX MaxX -MaxY MaxY]);
    end
    hold off
end