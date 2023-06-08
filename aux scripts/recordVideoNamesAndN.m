% Run date & time
run_date = char(datetime('now','Format','yyyy-MM-dd_HH.mm''''ss''''''''')) ;

% Select dir containing .xlsx track files
topLevelFolder = uigetdir('C:\') ;

% Create destination folder
destination_folder = strcat(fileparts(topLevelFolder), '\', run_date) ;
mkdir(destination_folder) ;

% Initialize log
diary off
diary_filename = strcat(destination_folder,'\VideoTrackCount_', run_date,'.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on

% Get a list of all .xlsx files and their containing subfolders in this
% directory
AllDirs = dir(fullfile(topLevelFolder, '**\*\*\*.xlsx')) ;

% Dump folder names from .folder field of "AllDirs" struct into a cell
% array
AllsubFolderNames = {AllDirs.folder} ;

% Filter unique instances of names
UsefulSubFolderNames = unique(AllsubFolderNames, 'sorted') ;

% Display subfolders' name on the console
for k = 1 : length(UsefulSubFolderNames)
	fprintf('Sub folder #%d = %s\n', k, UsefulSubFolderNames{k}) ;
end
for f = 12:length(UsefulSubFolderNames)

    % Find all subfolders containing xlsx files
    thisfoldertic = tic;
    cd(UsefulSubFolderNames{f})
    files = dir('*.xlsx') ;
    
    % Store valid condition names as a variable
    [~,name2,name3] = fileparts(pwd) ;
    disp(name2(1:end-3))
    condition = strcat(name2, name3) ;
    conditionValidName = matlab.lang.makeValidName(condition) ;
    
    % Sort file names in natural order 
    files = natsortfiles(files) ;
     
    % Retrieve tracks from condition's xlsx files
    condition_track_n = 0 ;
    for l = 1:length(files)
        thisxlsx = files(l).name ;
        temp_xlsx = readmatrix(thisxlsx, "Range", "E:F");
        startRow = 1;  %starting row number
        nRows = 3600;  %number of rows in segment
        totalRows = size(temp_xlsx,1) ;  %number of total rows
        nSegments = floor((totalRows - startRow + 1) / nRows);  %number of complete segments
        nLeftover = totalRows - startRow + 1 - (nSegments * nRows);  %number of rows at end that will be ignored
        [name2(1:end-3) ' ' thisxlsx ' ' 'contains' ' ' num2str(nSegments) ' ' 'tracks' ' ' 'with' ' ' num2str(nLeftover) ' ' 'leftover']
        if nLeftover ~= 0
            disp(nLeftover)
            error('File contains one or more tracks of length not equal to 3600')
        end
        % Compute start and end indices of each segment
        segStart = startRow : nRows : totalRows ;
        segStop = segStart(2)-1 : nRows : totalRows ;
        segment = cell(nSegments, 1) ;
        % Plot each track separately
        for s = 1:nSegments
            ['This is the' ' ' num2str(s) 'th' ' ' 'track' ' ' 'out of' ' ' ...
                num2str(nSegments) ' ' 'in' ' ' thisxlsx ' ' 'of' ' ' condition]
            condition_track_n = condition_track_n + 1 ;
            figure('Position',[5 440 315 200],'Name',strcat(thisxlsx(1:end-5),'_nº', ...
                num2str(s),'/',num2str(nSegments)),'NumberTitle','off')
            plot(temp_xlsx(segStart(s) : segStop(s), 1),temp_xlsx(segStart(s) : segStop(s), 2),'Color','k')
            hold on
            plot(temp_xlsx(segStop(s), 1),temp_xlsx(segStop(s),2),'ko','MarkerFaceColor','k','MarkerSize', 1.75)
            hold off
            axis equal
        end
    end
end