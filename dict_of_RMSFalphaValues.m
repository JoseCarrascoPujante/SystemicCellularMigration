close all
clear

load 'coordinates.mat' destination_folder
% Select dir containing .xlsx track files
topLevelFolder = uigetdir('C:\') ;

% Get a list of all .xlsx files and their containing subfolders in this
% directory
AllDirs = dir(fullfile(topLevelFolder, '**\*\*\*.xlsx')) ;

% Dump folder names from .folder field of "AllDirs" struct into a cell
% array
AllsubFolderNames = {AllDirs.folder} ;

% Filter unique instances of names
UsefulSubFolderNames = unique(AllsubFolderNames, 'sorted') ;

% Initialize bulk data structures
RMSFalphaDictionary = dictionary ;

% Custom pixel/mm ratio list for Metamoeba leningradensis chemotaxis
ratio_list = [23.44,23.44,23.44,23.44,23.44,23.44,11.63,11.63,11.63] ;

bar1 = waitbar(0,'In progress...','Name','Condition') ;
bar2 = waitbar(0,'In progress...','Name','File') ;
bar3 = waitbar(0,'In progress...','Name','Segment') ;

fg=figure(Visible="off") ;
h = gca;

%% Track extraction and plotting

tImportSec = tic;

for f = 1:length(UsefulSubFolderNames)

    % Find all subfolders containing xlsx files
    thisfoldertic = tic;
    cd(UsefulSubFolderNames{f})
    files = dir('*.xlsx') ;
    
    % Store valid condition names as a variable
    [~,name2,name3] = fileparts(pwd) ;
    condition = strcat(name2, name3) ;
    conditionValidName = matlab.lang.makeValidName(condition) ;
    
    % update condition/subfolder progress bar
    bar1 = waitbar(f/length(UsefulSubFolderNames), bar1, condition) ;

    % Sort file names in natural order 
    files = natsortfiles(files) ;
      
    % Retrieve tracks from condition's xlsx files
    condition_track_n = 0 ;
    for l = 1:length(files)
        thisxlsx = files(l).name ;
        bar2 = waitbar(l/length(files), bar2, thisxlsx) ;
        temp_xlsx = readmatrix(thisxlsx, "Range", "E:F");
        startRow = 1;  %starting row number
        nRows = 3600;  %number of rows in segment
        totalRows = size(temp_xlsx,1) ;  %number of total rows
        nSegments = floor((totalRows - startRow + 1) / nRows);  %number of complete segments
        nLeftover = totalRows - startRow + 1 - (nSegments * nRows);  %number of rows at end that will be ignored
        % [name2(1:end-3) ' ' thisxlsx ' ' 'contains' ' ' num2str(nSegments) ' ' 'tracks' ' ' 'with' ' ' num2str(nLeftover) ' ' 'leftover']
        if nLeftover ~= 0
            disp(nLeftover)
            error('File contains one or more tracks of length not equal to 3600')
        end
        % Compute start and end indices of each segment
        segStart = startRow : nRows : totalRows ;
        segStop = segStart(2)-1 : nRows : totalRows ;
        segment = cell(nSegments, 1) ;
        % Save each track separately
        for s = 1:nSegments
            % This is the s_th section:
            bar3 = waitbar(s/nSegments, bar3, strcat('Track_nº: ',num2str(s))) ;
            condition_track_n = condition_track_n + 1 ;
            raw_x = temp_xlsx(segStart(s) : segStop(s), 1) ;
            raw_y = temp_xlsx(segStart(s) : segStop(s), 2) ;
            
            centered_x = raw_x(:) - raw_x(1);
            centered_y = raw_y(:) - raw_y(1);
            
            [theta,rho] = cart2pol(centered_x, centered_y) ;

            if contains(condition, '11.63')
                scaled_rho = rho/11.63;
            elseif contains(condition, '23.44')
                scaled_rho = rho/23.44;
            else
                scaled_rho = rho/ratio_list(l);
            end
            
            [RMSFalphaDictionary(strcat(conditionValidName,'_',thisxlsx,'_global_',num2str(condition_track_n),...
                '_segment_',num2str(s))),~,~,~] = ...           
                amebas5(scaled_rho, h,'orig') ;            
        end
    end
end

zip = [RMSFalphaDictionary.values,RMSFalphaDictionary.keys];
sort(zip)
[sortedzip,I] = sort(zip(:,1));
sortedKeys = zip(I,2);
RMSFalphaRanking = [sortedKeys,sortedzip];

toc(tImportSec)

close(fg)

save("C:\Users\pc\Desktop\mov_sist\2023-06-07_14.16'19''_200000_shuffles\RMSFalphaDictionary700Cells.mat",'RMSFalphaDictionary','RMSFalphaRanking')