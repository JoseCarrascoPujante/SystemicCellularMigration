tic
% Select folder containing .xlsx track files
topLevelFolder = uigetdir('C:\') ;

% Get a list of all .xlsx files and their containing subfolders in this
% directory
AllDirs = dir(fullfile(topLevelFolder, '**\*\*\*.xlsx')) ;

% Dump folder names from .folder field of "AllDirs" struct into a cell
% array
AllsubFolderNames = {AllDirs.folder} ;

% Filter unique instances of names
UsefulSubFolderNames = unique(AllsubFolderNames, 'sorted') ;

% Display to-be-used folder names on the console
for k = 1 : length(UsefulSubFolderNames)
	fprintf('Sub folder #%d = %s\n', k, UsefulSubFolderNames{k}) ;
end

bar1 = waitbar(0,'In progress...','Name','Processing condition...') ;
bar2 = waitbar(0,'In progress...','Name','Processing file...') ;

% Iterate over UsefulSubfolders to extract X and Y
% coordinate vectors from the excel files
for f=1:length(UsefulSubFolderNames)
    folder=UsefulSubFolderNames{f};
    cd(folder)
    files = dir('*.xlsx') ;
    
    % Save condition name as a variable
    [~,name2,name3]=fileparts(pwd) ;
    condition = strcat(name2, name3) ;
    bar1 = waitbar(f/length(UsefulSubFolderNames), bar1, condition) ;
    condition_ValidName = matlab.lang.makeValidName(condition) ;
    
    % Sort file names alphanumerically to make their future index match
    % their name (e.g., Conditioned31.xlsx indexed as i=31) 
    files = natsortfiles(files) ;
    
    N = length(files) ;
    for i = 1:N
        thisfile = files(i).name ;
        bar2 = waitbar(i/N, bar2, thisfile) ;
        temp_x = readmatrix(thisfile,'Range','E1:E3600') ;
        temp_y = readmatrix(thisfile,'Range','F1:F3600') ;
        
        % Check for decimal separator problems in the data imported from the excel files
        % (run in two separate loops in case X and Y have different lengths)
        if mod(temp_x,1) == 0
            for k1 = 1:length(temp_x)
                if temp_x(k1) >= 1000000
                    msgbox(strcat(temp_x(k1),' X value in ', ...
                        thisfile,' is >1000'),"Error","error")
                    error(strcat('Program terminated because ', ...
                        temp_x(k1),' X value in ',thisfile,' is >1000')')
                end
                temp_x(k1) = str2double(separatethousands(temp_x(k1),'.',0)) ;
            end
        end
        
        if mod(temp_y,1) == 0
            for k1 = 1:length(temp_y)
                if temp_y(k1) >= 1000000
                    msgbox(strcat(temp_y(k1),' Y value in ', ...
                        thisfile,' is >1000'),"Error","error")
                    error(strcat('Program terminated because ', ...
                        temp_y(k1),' Y value in ',thisfile,' is >1000')')
                end
                temp_y(k1) = str2double(separatethousands(temp_y(k1),'.',0)) ;
            end
        end
    end
end

elapsed_time = datevec(toc./(60*60*24)) ;


