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

% List the parameters to be calculated by the script
parameters = {'alpha','r2','time_max','dfa_gamma','msd_beta','ApEn'} ;

% Initialize bulk data structures
coordinates = struct ;
results = struct ;

% Custom pixel/mm ratio list for Chemotaxis of Metamoeba leningradensis
ratio_list = [23.44,23.44,23.44,23.44,23.44,23.44,23.44,11.63,11.63,11.63,...
    11.63,11.63,11.63,11.63,11.63,23.44,23.44,23.44,23.44,23.44,23.44,...
    23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44,...
    23.44,23.44,23.44,23.44,11.63,11.63,11.63,11.63,23.44,23.44,23.44...
    ,11.63,11.63,11.63,11.63,11.63,11.63,11.63,23.44,23.44] ;

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
        
        % Save original X and Y coordinates as x(i) and y(i)
        coordinates.(condition_ValidName).original_x(:,i) = temp_x ;
        coordinates.(condition_ValidName).original_y(:,i) = temp_y ;
        
        % Center X and Y coordinates
        coordinates.(condition_ValidName).centered_x(:,i) = temp_x-temp_x(1) ;
        coordinates.(condition_ValidName).centered_y(:,i) = temp_y-temp_y(1) ;
        
        % Convert centered X and Y coordinates to polar coordinates
        [coordinates.(condition_ValidName).teta(:,i),...
            coordinates.(condition_ValidName).rho(:,i)] = ...
            cart2pol(coordinates.(condition_ValidName).centered_x(:,i),...
            coordinates.(condition_ValidName).centered_y(:,i)) ;
        
        % Scale Rho, X and Y
        if contains(condition, '11.63')
            coordinates.(condition_ValidName).scaled_rho(:,i) = coordinates.(condition_ValidName).rho(:,i)/11.63 ;
            coordinates.(condition_ValidName).scaled_x(:,i) = coordinates.(condition_ValidName).centered_x(:,i)/11.63 ;
            coordinates.(condition_ValidName).scaled_y(:,i) = coordinates.(condition_ValidName).centered_y(:,i)/11.63 ;
        elseif contains(condition, '23.44')
            coordinates.(condition_ValidName).scaled_rho(:,i) = coordinates.(condition_ValidName).rho(:,i)/23.44 ;
            coordinates.(condition_ValidName).scaled_x(:,i) = coordinates.(condition_ValidName).centered_x(:,i)/23.44 ;
            coordinates.(condition_ValidName).scaled_y(:,i) = coordinates.(condition_ValidName).centered_y(:,i)/23.44 ;
        elseif contains(condition, 'en nombre de cada')
            coordinates.(condition_ValidName).scaled_rho(:,i) = coordinates.(condition_ValidName).rho(:,i)/ratio_list(i) ;
            coordinates.(condition_ValidName).scaled_x(:,i) = coordinates.(condition_ValidName).centered_x(:,i)/ratio_list(i) ;
            coordinates.(condition_ValidName).scaled_y(:,i) = coordinates.(condition_ValidName).centered_y(:,i)/ratio_list(i) ;
        end
        
        % RMSF
        [results.(condition_ValidName)(i,strcmp(parameters(:), 'alpha'))...
            ,results.(condition_ValidName)(i,strcmp(parameters(:), 'r2')),...
            results.(condition_ValidName)(i,strcmp(parameters(:), 'time_max'))] = ...
            amebas5(coordinates.(condition_ValidName).scaled_rho(:,i)) ;
        
        % DFA
        results.(condition_ValidName)(i,strcmp(parameters(:), 'dfa_gamma')) = ...
            DFA_main2(coordinates.(condition_ValidName).scaled_rho(:,i)) ;

        % MSD
        results.(condition_ValidName)(i,strcmp(parameters(:), 'msd_beta')) =...
            msd(coordinates.(condition_ValidName).scaled_x(:,i),coordinates.(condition_ValidName).scaled_y(:,i)) ;
        
        % Approximate entropy (Kolmogorov-Sinai entropy)
        results.(condition_ValidName)(i,strcmp(parameters(:), 'ApEn')) = ...
            ApEn(2, 0.2*std(coordinates.(condition_ValidName).scaled_rho(:,i)),...
            coordinates.(condition_ValidName).scaled_rho(:,i)) ;
    end
end

elapsed_time = datevec(toc./(60*60*24)) ;


