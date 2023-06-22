%% Preprocessing

% Run date & time
run_date = char(datetime('now','Format','yyyy-MM-dd_HH.mm''''ss''''''''')) ;

% Select dir containing .xlsx track files
topLevelFolder = uigetdir('C:\') ;

% Number of times Rho(s) will be randomly permuted to test against the
% null hypothesis that results are random
shuffles = 200000;

% Create destination folder
destination_folder = strcat(fileparts(topLevelFolder), '\', run_date, '_', ...
    num2str(shuffles), '_shuffles') ;
mkdir(destination_folder) ;

% Initialize log
diary off
diary_filename = strcat(destination_folder,'\MatlabCommandWindowlog_', run_date,'.txt') ;
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

% List the parameters to be calculated by the script
stat_names = {'RMSF_alpha', 'sRMSF_alpha', 'RMSF_R2', 'sRMSF_R2', 'RMSFCorrelationTime', ...
    'sRMSFCorrelationTime', 'DFA_gamma', 'sDFA_gamma', 'MSD_beta', 'sMSD_beta', 'AppEn', 'sAppEn'} ;

% Initialize bulk data structures
tracks = struct ;
coordinates = struct ;
results = struct ;

% Custom pixel/mm ratio list for Metamoeba leningradensis chemotaxis
ratio_list = ...
[
    [23.44,23.44,23.44,23.44,23.44],...
    [23.44,23.44,23.44,23.44,23.44,23.44,23.44, 23.44],...
    [23.44,23.44,23.44,23.44,23.44],...
    [23.44,23.44,23.44,23.44,23.44],...
    [23.44,23.44,23.44,23.44,23.44,23.44,23.44,23.44],...
    [23.44,23.44,23.44,23.44,23.44,23.44,23.44],...
    [11.63,11.63,11.63,11.63,11.63,11.63,11.63,11.63],...
    [11.63,11.63,11.63,11.63,11.63,11.63,11.63,11.63],...
    [11.63,11.63,11.63,11.63,11.63,11.63],...
] ;

bar1 = waitbar(0,'In progress...','Name','Reading condition files...') ;
bar2 = waitbar(0,'In progress...','Name','Reading file...') ;

%% Track extraction and plotting

tImportSec = tic;

for f = 1:length(UsefulSubFolderNames)

    % Find all subfolders containing xlsx files
    thisfoldertic = tic;
    cd(UsefulSubFolderNames{f})
    files = dir('*.xlsx') ;
    
    % Store valid condition names as a variable
    [~,name2,name3] = fileparts(pwd) ;
    disp(name2(1:end-3))
    condition = strcat(name2, name3) ;
    conditionValidName = matlab.lang.makeValidName(condition) ;
    
    % update condition/subfolder progress bar
    bar1 = waitbar(f/length(UsefulSubFolderNames), bar1, condition) ;

    % Sort file names in natural order 
    files = natsortfiles(files) ;
    
    % Initialize condition track figure
    figures.(conditionValidName).tracks = figure('Name',strcat('Tracks_',condition),...
        'Visible','off','NumberTitle','off') ;
    hTracks = gca;
    
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
        % Save each track separately
        for s = 1:nSegments
            % This is the s_th section:
            condition_track_n = condition_track_n + 1 ;
            tracks.(conditionValidName).(genvarname(num2str(condition_track_n))) = ...
                temp_xlsx(segStart(s) : segStop(s), :) ; 
        end
    end

    % Tracks loop
    A = fieldnames(tracks.(conditionValidName)) ;
    for i = 1:length(A)
        
        thisfiletic = tic;
      
        % update track progress bar
        bar2 = waitbar(i/length(A), bar2, A{i}) ;
             
        % Save original X and Y coordinates as x(i) and y(i)
        coordinates.(conditionValidName).original_x(:,i) = tracks.(conditionValidName).(A{i})(:,1) ;
        coordinates.(conditionValidName).original_y(:,i) = tracks.(conditionValidName).(A{i})(:,2) ;
        
        % Center X and Y coordinates (substract first value to all series)
        coordinates.(conditionValidName).centered_x(:,i) = ...
            coordinates.(conditionValidName).original_x(:,i) - coordinates.(conditionValidName).original_x(1,i) ;
        coordinates.(conditionValidName).centered_y(:,i) = ...
            coordinates.(conditionValidName).original_y(:,i) - coordinates.(conditionValidName).original_y(1,i) ;
        
        % Polar coordinate conversion
        [coordinates.(conditionValidName).theta(:,i),...
            coordinates.(conditionValidName).rho(:,i)] = ...
            cart2pol(coordinates.(conditionValidName).centered_x(:,i),...
            coordinates.(conditionValidName).centered_y(:,i)) ;
        
        % Scale Rho, X, and Y
        if contains(condition, '11.63')
            coordinates.(conditionValidName).scaled_rho(:,i) = ...
                coordinates.(conditionValidName).rho(:,i)/11.63 ;
            coordinates.(conditionValidName).scaled_x(:,i) = ...
                coordinates.(conditionValidName).centered_x(:,i)/11.63 ;
            coordinates.(conditionValidName).scaled_y(:,i) = ...
                coordinates.(conditionValidName).centered_y(:,i)/11.63 ;
        elseif contains(condition, '23.44')
            coordinates.(conditionValidName).scaled_rho(:,i) = ...
                coordinates.(conditionValidName).rho(:,i)/23.44 ;
            coordinates.(conditionValidName).scaled_x(:,i) = ...
                coordinates.(conditionValidName).centered_x(:,i)/23.44 ;
            coordinates.(conditionValidName).scaled_y(:,i) = ...
                coordinates.(conditionValidName).centered_y(:,i)/23.44 ;
        else
            coordinates.(conditionValidName).scaled_rho(:,i) = ...
                coordinates.(conditionValidName).rho(:,i)/ratio_list(i) ;
            coordinates.(conditionValidName).scaled_x(:,i) = ...
                coordinates.(conditionValidName).centered_x(:,i)/ratio_list(i) ;
            coordinates.(conditionValidName).scaled_y(:,i) = ...
                coordinates.(conditionValidName).centered_y(:,i)/ratio_list(i) ;
        end
        
        % Shuffle X and Y trajectories
       coordinates.(conditionValidName).shuffled_x(:,i) = ...
            coordinates.(conditionValidName).scaled_x(:,i) ; % initialize shuffled X

        coordinates.(conditionValidName).shuffled_y(:,i) = ...
            coordinates.(conditionValidName).scaled_y(:,i) ; % initialize shuffled Y
            
        for k=1:shuffles
          coordinates.(conditionValidName).shuffled_x(:,i) = ...
              shuff(coordinates.(conditionValidName).shuffled_x(:,i)) ;

          coordinates.(conditionValidName).shuffled_y(:,i) = ...
              shuff(coordinates.(conditionValidName).shuffled_y(:,i)) ;
        end
        
        % Shuffle Scaled_Rho
        coordinates.(conditionValidName).shuffled_rho(:,i) = ...
            coordinates.(conditionValidName).scaled_rho(:,i) ;
        
        for k=1:shuffles
          coordinates.(conditionValidName).shuffled_rho(:,i) = ...
              shuff(coordinates.(conditionValidName).shuffled_rho(:,i)) ;
        end
        
        % Plot trajectory and place black dot marker at its tip      
        plot(hTracks,coordinates.(conditionValidName).scaled_x(:,i),...
            coordinates.(conditionValidName).scaled_y(:,i), 'Color', [0, 0, 0]) ;
        hold on;
        plot(hTracks,coordinates.(conditionValidName).scaled_x(end,i),...
            coordinates.(conditionValidName).scaled_y(end,i),...
            'ko',  'MarkerFaceColor',  [0, 0, 0], 'MarkerSize', 2) ;
        
        [A{i} ' runtime was ' num2str(toc(thisfiletic)) ' seconds']
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
    
    [condition ' runtime was ' num2str(toc(thisfoldertic)) ' seconds']
end

tImportSec = num2str(toc(tImportSec)) ;

% Save data
save(strcat(destination_folder, '\', run_date, '_coordinates.mat'),...
    'coordinates','figures','stat_names','shuffles','tImportSec','destination_folder') ;

['Coordinate section FINISHED in ', tImportSec, ' seconds']

%% Parameter calculation

tCalcSec=tic;

field_names = fieldnames(coordinates) ;

bar3 = waitbar(0,'In progress...','Name','Condition...') ;
bar4 = waitbar(0,'In progress...','Name','Track number...') ;

for i = 1:length(field_names)
    
    bar3 = waitbar(i/length(field_names), bar3, field_names{i}) ;
    
    foldertime = tic;

    figures.(field_names{i}).msd = figure('Name',strcat('MSD_',field_names{i}),...
        'Visible','off','NumberTitle','off') ;
    xlabel('Log(MSD(\tau))');
    ylabel('Log(\tau(s))');
    msdhandle = gca;

    figures.(field_names{i}).msd_shuff = figure('Name',strcat('MSD_Shuffled_',...
        field_names{i}),'Visible','off','NumberTitle','off') ;
    xlabel('Log(MSD(\tau))');
    ylabel('Log(\tau(s))');
    Shuffmsdhandle = gca;
    
    N = length(coordinates.(field_names{i}).original_x(1,:)); % Trajectories in condition
    for j = 1:N 
        
        tic

        bar4 = waitbar(j/N, bar4, strcat('Track number', ' ', num2str(j))) ;

        % RMSFalpha
		figures.(field_names{i}).rmsf.(strcat('number', num2str(j))) = ...
		figure('Name',strcat(field_names{i}, '_amoeba_number_',...
		num2str(j)),'NumberTitle','off','Visible','off');
		
		rmsfhandle = gca;
        set(rmsfhandle,'xscale','log')
        set(rmsfhandle,'yscale','log')
		
        [results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_alpha')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_R2')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'RMSFCorrelationTime')), tc2] = ...
            amebas5(coordinates.(field_names{i}).scaled_rho(:,j), rmsfhandle) ;
        
        % Shuffled RMSFalpha
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSF_alpha')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSF_R2')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSFCorrelationTime')), ~] = ...
            amebas5(coordinates.(field_names{i}).shuffled_rho(:,j), rmsfhandle) ;

        legend('Original RMSF','Shuffled RMSF','Location','best')

        % DFAgamma
        figures.(field_names{i}).dfa_original.(strcat('number', num2str(j))) = ...
            figure('Name',strcat('DFA_', field_names{i}, '_cell_no_', num2str(j)),...
            'NumberTitle','off','Visible','off') ;

        dfahandle = gca;

        [results.(field_names{i})(j,strcmp(stat_names(:), 'DFA_gamma'))] = ...
            DFA_main2(coordinates.(field_names{i}).scaled_rho(:,j),...
            'Original_DFA_', dfahandle) ;

        % Shuffled DFAgamma
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sDFA_gamma'))] = ...
            DFA_main2(coordinates.(field_names{i}).shuffled_rho(:,j),...
            'Shuffled_DFA_', dfahandle) ;

        legend('Original data','Shuffled data (Gaussian noise)',...
            'Original DFA \gamma','Shuffled DFA \gamma',...
            'Location','northwest')

        % MSDbeta
        [results.(field_names{i})(j,strcmp(stat_names(:), 'MSD_beta'))] = ...
            msd(coordinates.(field_names{i}).scaled_x(:,j),...
            coordinates.(field_names{i}).scaled_y(:,j),msdhandle) ;
        
        % Shuffled MSDbeta
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sMSD_beta'))] = ...
            msd(coordinates.(field_names{i}).shuffled_x(:,j),...
            coordinates.(field_names{i}).shuffled_y(:,j),Shuffmsdhandle) ;

        % Approximate entropy (Kolmogorov-Sinai entropy)
        results.(field_names{i})(j,strcmp(stat_names(:), 'AppEn')) = ...
            ApEn(2, 0.2*std(coordinates.(field_names{i}).scaled_rho(:,j)),...
            coordinates.(field_names{i}).scaled_rho(:,j)) ;

        % Shuffled Approximate entropy (Kolmogorov-Sinai entropy)
        results.(field_names{i})(j,strcmp(stat_names(:), 'sAppEn')) = ...
            ApEn(2, 0.2*std(coordinates.(field_names{i}).shuffled_rho(:,j)),...
            coordinates.(field_names{i}).shuffled_rho(:,j)) ;

        [field_names{i} ' amoeba ' num2str(j) ' runtime was ' num2str(toc) ' seconds']

    end

    [field_names{i} ' runtime was ' num2str(toc(foldertime)) ' seconds']

end

%# Save data

tCalcSec = num2str(toc(tCalcSec)) ;

save(strcat(destination_folder, '\', run_date ,'_numeric_results.mat'),...
    'tCalcSec', 'results') ;

['Calculations section runtime FINISHED in ' tCalcSec ' seconds']

%% Statistical analyses

tStatsSec=tic;

% Kolmogorov-Smirnov test against the standard normal distribution of parameters

for field = 1:length(field_names)
    n_hypotheses = size(results.(field_names{field}),2) ; % N hypotheses being tested, needed by Bonferroni correction)
    for column_index = 1:n_hypotheses

        x = (results.(field_names{field})(:,column_index)...
            -mean(results.(field_names{field})(:,column_index)))...
            /std(results.(field_names{field})(:,column_index)) ;
        
        [hypothesis, p_value, ksstat, cv] = kstest(x,'Alpha',(0.05/n_hypotheses)) ; % with Bonferroni correction

        kolmogorov_smirnov.(field_names{field})(2,column_index) = p_value ;
        kolmogorov_smirnov.(field_names{field})(1,column_index) = hypothesis ;        
        kolmogorov_smirnov.(field_names{field})(3,column_index) = ksstat ;
        kolmogorov_smirnov.(field_names{field})(4,column_index) = cv ;
        
    end
end

['K-S test: done']

% Create groups to test data by condition and species under
% Kruskal-Wallis test

kw_conds = {{'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'},...
    {'Proteus','Leningradensis','Borokensis'}} ; % Conditions for the Kruskal-Wallis test

for group = 1:length(kw_conds)
    count = 0 ;
    for condit = 1:length(kw_conds{group})
        count = 0 ;
        for field = 1:length(field_names)
            if contains(field_names{field},kw_conds{group}{condit})
                count = count + 1 ;
                for ind = 1:length(stat_names) % Shuffled parameters can be excluded by tweaking this line
                    kruskal_groups.(kw_conds{group}{condit}).(stat_names{ind}){:,count} = ...
                        results.(field_names{field})(:,ind) ;
                end
            end
        end
    end
end

['Created K-W groups']

% Kruskal-Wallis test to determine if the samples come from the same population/different
% populations with the same distribution).

fields_krusk = fieldnames(kruskal_groups) ;
for kwfield = 1:length(fields_krusk)
    if kwfield <= 4
        subgroup = 2;
    elseif kwfield > 4
        subgroup = 1;
    end
    for params=1:length(stat_names)
        
        [kruskal_results.(fields_krusk{kwfield}){1, params}, ... % p_value
            kruskal_results.(fields_krusk{kwfield}){2, params}, ... % anova_table
            kruskal_results.(fields_krusk{kwfield}){3, params}] = ... % test_stats
        kruskalwallis(padcat(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){:}),kw_conds{subgroup}',"off") ;
        
        % Run multicomparison post-hoc test with bonferroni correction
        % on the 'stats' structure output by the Kruskal-Wallis function
        [multcomp_results.(fields_krusk{kwfield}){1, params},... % matrix of multiple comparison results
            multcomp_results.(fields_krusk{kwfield}){2, params},... % matrix of estimates
            multcomp_results.(fields_krusk{kwfield}){3, params},... % handle to the figure
            multcomp_results.(fields_krusk{kwfield}){4, params}] = ... % group names
        multcompare(kruskal_results.(fields_krusk{kwfield}){3, params},...
            "Alpha",0.05,"Display","off","CriticalValueType","bonferroni");
        
        % Run Dunn's post-hoc test to benchmark "multcompare"
        if kwfield <= 4
        [dunn_results.(fields_krusk{kwfield}){1, params},...
            dunn_results.(fields_krusk{kwfield}){2, params},...
            dunn_results.(fields_krusk{kwfield}){3, params}] = ...
        dunn(vertcat(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){:}).',...
            [ones(1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){1}))...
            repmat(2,1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){2}))...
            repmat(3,1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){3}))],0);

        elseif kwfield > 4
            [dunn_results.(fields_krusk{kwfield}){1, params},...
                dunn_results.(fields_krusk{kwfield}){2, params},...
                dunn_results.(fields_krusk{kwfield}){3, params}] = ...
            dunn(vertcat(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){:}).',...
                [ones(1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){1}))...
                repmat(2,1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){2}))...
                repmat(3,1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){3}))...
                repmat(4,1,length(kruskal_groups.(fields_krusk{kwfield}).(stat_names{params}){4}))],0);
        end
    end
end

['Done performing statistical tests']

tStatsSec = num2str(toc(tStatsSec)) ;

% Save data
save(strcat(destination_folder, '\', run_date, '_statistics.mat'), ...
    'kolmogorov_smirnov','kruskal_groups','kruskal_results','multcomp_results','dunn_results','tStatsSec') ;


%% File backup and plot exporting

% Back up 'scripts' folder

copyfile('C:\Users\pc\Desktop\mov_sist\code', strcat(destination_folder, '\Scripts_backup'));

% Export figures as vector graphic files (.svg)

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig) ;
%   FigName = get(FigHandle, 'Name') ;
%   set(0, 'CurrentFigure', FigHandle) ;
%   export_fig(strcat(destination_folder,'\',FigName),'-svg')
% %   exportgraphics(FigHandle, fullfile(destination_folder, [FigName '.pdf']), 'Resolution', 600) ;
% end

%% Generate publishable figures

% Navigate through coordinates struct sequentially plotting each track
plot_tracks_iterator(coordinates)

% Graphical Abstract
figures.GraphAbs = GraphAbs_subaxis_def(field_names, results, figures, destination_folder) ;

% Figure 1
figures.figure1_4Panels = Figure1(coordinates, destination_folder) ;

% Figure 2
Figure2

% Figure 3
Figure3


% Save data
save(strcat(destination_folder, '\', run_date, '_figures.mat'), ...
    'figures') ;

diary off
