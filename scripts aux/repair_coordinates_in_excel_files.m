tic
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

% Display to-be-used folders' name on the console
for k = 1 : length(UsefulSubFolderNames)
	fprintf('Sub folder #%d = %s\n', k, UsefulSubFolderNames{k}) ;
end

% Iterate over UsefulSubfolders to extract X and Y
% coordinate vectors from the excel files
for f=1:length(UsefulSubFolderNames)
    folder=UsefulSubFolderNames{f};
    cd(folder)
    files = dir('*.xlsx') ;
    
    % Sort file names alphanumerically so their index/position
    % matches their name (e.g. Conditioned31.xlsx shall be indexed as i=31) 
    files = natsortfiles(files) ;
    N = length(files); % length(files) ;
    for i = 1:N
        thisfile = files(i).name ;
        temp_x = readmatrix(thisfile,'Range','E1:E3600') ;
        temp_y = readmatrix(thisfile,'Range','F1:F3600') ;
        
        % Check for decimal separator problems in the data imported from the excel files
        % (run in two separate loops in case X and Y have different lengths)
        if mod(temp_x,1) == 0 % If the modulus of every number/1 is 0, all numbers are integers and vector needs to be repaired
            for k1 = 1:length(temp_x)
                if temp_x(k1) >= 1000000
                    msgbox(strcat(temp_x(k1),' X coordinate in ', ...
                        thisfile,' is >1.000.000'),"Error","error")
                    error(strcat('Program terminated because ', ...
                        temp_x(k1),' X coordinate in ',thisfile,' is >1.000.000')')
                end
                temp_x(k1) = str2double(separatethousands(temp_x(k1),'.',0)) ;
            end
            writematrix(temp_x,thisfile,'Sheet',1,'Range','E1:E3600')
        end
        
        if mod(temp_y,1) == 0 % If the modulus of every number/1 is 0, all numbers are integers and vector needs to be repaired
            for k1 = 1:length(temp_y)
                if temp_y(k1) >= 1000000
                    msgbox(strcat(temp_y(k1),' Y coordinate in ', ...
                        thisfile,' is >1.000.000'),"Error","error")
                    error(strcat('Program terminated because ', ...
                        temp_y(k1),' Y coordinate in ',thisfile,' is >1.000.000')')
                end
                temp_y(k1) = str2double(separatethousands(temp_y(k1),'.',0)) ;
            end
            writematrix(temp_y,thisfile,'Sheet',1,'Range','F1:F3600')
        end
    end
end
toc