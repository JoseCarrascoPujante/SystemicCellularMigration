%% Preprocessing section

% Select dir containing .xlsx track files
topLevelFolder = uigetdir('C:\') ;

% Number of times Rho(s) will be randomly permuted to test against the
% null hypothesis that DFA results are random
shuffles = 200000;

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

% List the parameters to be calculated by the script

parameters = {'alpha','shuff_alpha','r2','shuff_r2','correlation_time','shuff_correlation_time',...
    'dfa_gamma','shuff_dfa_gamma','msd_beta','shuff_msd_beta','ApEn','shuff_ApEn'} ;

% Initialize bulk data structures
results = struct ;
coordinates = struct ;

% Custom pixel/mm ratio list for Metamoeba leningradensis chemotaxis
ratio_list = [[23.44,23.44,23.44,23.44,23.44,23.44,23.44],...
    [11.63,11.63,11.63,11.63,11.63,11.63,11.63,11.63],...
    [23.44,23.44,23.44,23.44,23.44,23.44,23.44],...
    [23.44,23.44,23.44,23.44,23.44],...
    [23.44,23.44],...
    [23.44,23.44,23.44,23.44,23.44,23.44,23.44],...
    [11.63,11.63,11.63,11.63],...
    [23.44,23.44,23.44],...
    [11.63,11.63,11.63,11.63,11.63,11.63,11.63],...
    [23.44,23.44],...
    [38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651],...
    [38.7651,38.7651,38.7651,38.7651,38.7651],...
    [38.7651,38.7651,38.7651,38.7651,38.7651,38.7651,38.7651]] ;

bar1 = waitbar(0,'In progress...','Name','Reading condition files...') ;
bar2 = waitbar(0,'In progress...','Name','Reading file...') ;

% Create destination folder

run_date = datestr(now,'YYYY-mm-dd hhºMM''ss''''') ;

destination_folder = strcat(topLevelFolder, '\', run_date, '_', ...
    num2str(shuffles), '_shuffles') ;

mkdir(destination_folder) ;


%% Coordinate extraction section

ImportTime = tic;

for f=1:length(UsefulSubFolderNames)

    % Find all subfolders containing xlsx files
    thisfoldertic = tic;
    folder=UsefulSubFolderNames{f};
    cd(folder)
    files = dir('*.xlsx') ;
    
    % Store valid condition names as a variable

    [~,name2,name3] = fileparts(pwd) ;
    condition = strcat(name2, name3) ;
    conditionValidName = matlab.lang.makeValidName(condition) ;
    
    % update condition/subfolder progress bar

    bar1 = waitbar(f/length(UsefulSubFolderNames), bar1, condition) ;

    % Sort file names in natural alphanumerical order so index/position
    % matches naming (e.g. Conditioned31.xlsx shall be indexed at i=31) 

    files = natsortfiles(files) ;
    
    % Initialize condition track figure

    figures.(conditionValidName).tracks = figure('Name',strcat('Tracks_',condition),...
        'Visible','off','NumberTitle','off') ;
    hTracks = gca;

    % Tracks loop
    N = 1:40 ;
    for i = N
        
        thisfile = files(i).name ;
        thisfiletic = tic;

        % update track progress bar
        bar2 = waitbar(i/length(N), bar2, thisfile) ;
        
        % Read individual amoeba trajectory from excel file
        temp_x = readmatrix(thisfile,'Range','E1:E3600') ;
        temp_y = readmatrix(thisfile,'Range','F1:F3600') ;
     
        % Save original X and Y coordinates as x(i) and y(i)
        coordinates.(conditionValidName).original_x(:,i) = temp_x ;
        coordinates.(conditionValidName).original_y(:,i) = temp_y ;
        
        % Center X and Y coordinates
        coordinates.(conditionValidName).centered_x(:,i) = temp_x-temp_x(1) ;
        coordinates.(conditionValidName).centered_y(:,i) = temp_y-temp_y(1) ;
        
        % Polar coordinate conversion
        [coordinates.(conditionValidName).teta(:,i),...
            coordinates.(conditionValidName).rho(:,i)] = ...
            cart2pol(coordinates.(conditionValidName).centered_x(:,i),...
            coordinates.(conditionValidName).centered_y(:,i)) ;
        
        % Scale Rho, X and Y
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
        elseif contains(condition, 'en nombre de cada')
            coordinates.(conditionValidName).scaled_rho(:,i) = ...
                coordinates.(conditionValidName).rho(:,i)/ratio_list(i) ;
            coordinates.(conditionValidName).scaled_x(:,i) = ...
                coordinates.(conditionValidName).centered_x(:,i)/ratio_list(i) ;
            coordinates.(conditionValidName).scaled_y(:,i) = ...
                coordinates.(conditionValidName).centered_y(:,i)/ratio_list(i) ;
        end
        
        
        % Shuffle X and Y
       coordinates.(conditionValidName).shuffled_x(:,i) = ...
            coordinates.(conditionValidName).scaled_x(:,i) ; % make copy of scaled x

        coordinates.(conditionValidName).shuffled_y(:,i) = ...
            coordinates.(conditionValidName).scaled_y(:,i) ; % make copy of scaled y
        
        tic
        
        for k=1:shuffles

          coordinates.(conditionValidName).shuffled_x(:,i) = ...
              shuff(coordinates.(conditionValidName).shuffled_x(:,i)) ;

          coordinates.(conditionValidName).shuffled_y(:,i) = ...
              shuff(coordinates.(conditionValidName).shuffled_y(:,i)) ;

        end
        
        % Shuffle Rho
        coordinates.(conditionValidName).shuffled_rho(:,i) = ...
            coordinates.(conditionValidName).scaled_rho(:,i) ; % make a copy of scaled rho
        
        tic
        for k=1:shuffles
          coordinates.(conditionValidName).shuffled_rho(:,i) = ...
              shuff(coordinates.(conditionValidName).shuffled_rho(:,i)) ;
        end
        
        % Plot trajectory and 'ko' marker at its tip      
        plot(hTracks,coordinates.(conditionValidName).scaled_x(:,i),...
            coordinates.(conditionValidName).scaled_y(:,i), 'Color', [0, 0, 0]) ;
        hold on;
        plot(hTracks,coordinates.(conditionValidName).scaled_x(end,i),...
            coordinates.(conditionValidName).scaled_y(end,i),...
            'ko',  'MarkerFaceColor',  [0, 0, 0], 'MarkerSize', 2) ;
        
        [thisfile ' runtime was ' num2str(toc(thisfiletic)) ' seconds']
    end
    
    
    %Adjust 'tracks' figure's proportions
    divx=[-18 18];
    divy=[0 0];
    plot(hTracks,divx,divy,'k');
    hold on;
    plot(hTracks,divy,divx,'k');
    hold on;
    axis([-18 18 -18 18]);
    daspect([1 1 1]);
    box on;
    hold off
    
    [condition ' runtime was ' num2str(toc(thisfoldertic)) ' seconds']
end

CoordImportAndShuffleTime = datevec(toc(ImportTime)./(60*60*24)) ;

save(strcat(destination_folder, '\', run_date, '_coordinates'), 'N', 'parameters', 'shuffles', 'coordinates',...
    'CoordImportAndShuffleTime') ;

['Track section runtime was ' num2str(ImportTime) ' seconds']


%% Calculations section

tCalc=tic;

field_names = fieldnames(coordinates) ;

bar3 = waitbar(0,'In progress...','Name','Condition...') ;
bar4 = waitbar(0,'In progress...','Name','Track number...') ;

for i=1:length(field_names)
    
    bar3 = waitbar(i/length(field_names), bar3, field_names{i}) ;
    
    foldern = tic;

    figures.(field_names{i}).msd = figure('Name',strcat('MSD_',field_names{i}),...
        'Visible','off','NumberTitle','off') ;
    xlabel('Log(MSD(\tau))');
    ylabel('Log(\tau(s))');
    msdhandle = gca;

    figures.(field_names{i}).msd_shuff = figure('Name',strcat('MSD_Shuffled_',field_names{i}),...
        'Visible','off','NumberTitle','off') ;
    xlabel('Log(MSD(\tau))');
    ylabel('Log(\tau(s))');
    Shuffmsdhandle = gca;

    for j=N %number of trajectories per condition
        
        tic

        bar4 = waitbar(j/length(N), bar4, strcat('Track number', ' ', num2str(j))) ;

        % RMSF
		
		figures.(field_names{i}).rmsf.(strcat('number', num2str(j))) = ...
		figure('Name',strcat(field_names{i}, '_amoeba_number_',...
		num2str(j),'_','_steps_'),'NumberTitle','off','Visible','off');
		
		rmsfhandle = gca;
		
        [results.(field_names{i})(j,strcmp(parameters(:), 'alpha')),...
            results.(field_names{i})(j,strcmp(parameters(:), 'r2')),...
            results.(field_names{i})(j,strcmp(parameters(:), 'correlation_time')), tc2] = ...
            amebas5(coordinates.(field_names{i}).scaled_rho(:,j), rmsfhandle) ;
        
        % Shuffled RMSF

        [results.(field_names{i})(j,strcmp(parameters(:), 'shuff_alpha')),...
            results.(field_names{i})(j,strcmp(parameters(:), 'shuff_r2')),...
            results.(field_names{i})(j,strcmp(parameters(:), 'shuff_correlation_time')), ~] = ...
            amebas5(coordinates.(field_names{i}).shuffled_rho(:,j), rmsfhandle) ;

        legend('Original RMSF','Shuffled RMSF','Location','best')

        % DFA
		
        figures.(field_names{i}).dfa_original.(strcat('number', num2str(j))) = ...
            figure('Name',strcat('DFA_', field_names{i}, '_cell_no_', num2str(j)),...
            'NumberTitle','off','Visible','off') ;

        dfahandle = gca;

        [results.(field_names{i})(j,strcmp(parameters(:), 'dfa_gamma'))] = ...
            DFA_main2(coordinates.(field_names{i}).scaled_rho(:,j), 'Original_DFA_', dfahandle) ;

        % Shuffled DFA (using same shuffled data as for shuffled RMSF)

        [results.(field_names{i})(j,strcmp(parameters(:), 'shuff_dfa_gamma'))] = ...
            DFA_main2(coordinates.(field_names{i}).shuffled_rho(:,j), 'Shuffled_DFA_', dfahandle) ;

        legend('Original data','Shuffled data (Gaussian noise)',...
            'Original DFA \gamma','Shuffled DFA \gamma',...
            'Location','northwest')

        % MSD
        
        [results.(field_names{i})(j,strcmp(parameters(:), 'msd_beta'))] = ...
            msd(coordinates.(field_names{i}).scaled_x(:,j), coordinates.(field_names{i}).scaled_y(:,j), ...
            msdhandle) ;
        
        % Shuffled MSD
        
        [results.(field_names{i})(j,strcmp(parameters(:), 'shuff_msd_beta'))] = ...
            msd(coordinates.(field_names{i}).shuffled_x(:,j), coordinates.(field_names{i}).shuffled_y(:,j), ...
            Shuffmsdhandle) ;

        % Approximate entropy (Kolmogorov-Sinai entropy)
        results.(field_names{i})(j,strcmp(parameters(:), 'ApEn')) = ...
            ApEn(2, 0.2*std(coordinates.(field_names{i}).scaled_rho(:,j)),...
            coordinates.(field_names{i}).scaled_rho(:,j)) ;

        % Shuffled Approximate entropy (Kolmogorov-Sinai entropy)
        results.(field_names{i})(j,strcmp(parameters(:), 'shuff_ApEn')) = ...
            ApEn(2, 0.2*std(coordinates.(field_names{i}).shuffled_rho(:,j)),...
            coordinates.(field_names{i}).shuffled_rho(:,j)) ;

        [field_names{i} ' amoeba ' num2str(j) ' runtime was ' num2str(toc) ' seconds']

    end

    [field_names{i} ' runtime was ' num2str(toc(foldern)) ' seconds']

end

% 3D and 2D scatter plots

figures = scatts(field_names,results,figures, 'original') ;
figures = scatts(field_names,results,figures, 'standardized') ;
figures = scatts(field_names,results,figures, 'normalized') ;

calcTime = datevec(toc(tCalc)./(60*60*24)) ;

save(strcat(destination_folder, '\', run_date, '_rmsf_', num2str(tc2), 'tmax_calculations&figures'),...
    'calcTime', 'results', 'figures') ;

old_dest_folder = destination_folder;
destination_folder = strcat(old_dest_folder, '_', num2str(tc2), 'tmax');

movefile(old_dest_folder, destination_folder)
clear old_dest_folder

['Calculations section runtime was ' num2str(toc(tCalc)) ' seconds']

%% Statistical analyses

% Kolmogorov-Smirnov test against the standard normal distribution of parameters

for field = 1:length(field_names)
    n_hypotheses = size(results.(field_names{field}),2) ; % N hypotheses being tested, needed for Bonferroni correction)
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

['Finished K-S test']

% Create groups to test data by condition and species under the
% Kruskal-Wallis test

kw_conds = {{'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion','Comprobacion'},...
    {'Proteus','Leningradensis','Borokensis'}} ; % Conditions for the Kruskal-Wallis test

for group = 1:length(kw_conds)
    count = 0 ;
    for condit = 1:length(kw_conds{group})
        count = 0 ;
        for field = 1:length(field_names)
            if contains(field_names{field},kw_conds{group}{condit})
                count = count + 1 ;
                for ind = 1:length(parameters) % Shuffled parameters can be excluded by tweaking this line
                    kruskal_groups.(kw_conds{group}{condit}).(parameters{ind})(:,count) = ...
                        results.(field_names{field})(:,ind) ;
                end
            end
        end
    end
end

['Created K-W groups']

% Kruskal-Wallis test to determine if the samples come from the same population
% (or, equivalently, from different populations with the same distribution).

fields_krusk = fieldnames(kruskal_groups) ;
for kwfield = 1:length(fields_krusk)
    if kwfield <= 5
        subgroup = 2;
    elseif kwfield > 5
        subgroup = 1;
    end
    for params =  1:length(parameters)
        
        [kruskal_results.(fields_krusk{kwfield}){1, params}, ... % p_value
            kruskal_results.(fields_krusk{kwfield}){2, params}, ... % anova_table
            kruskal_results.(fields_krusk{kwfield}){3, params}] = ... % test_stats
        kruskalwallis(kruskal_groups.(fields_krusk{kwfield}).(parameters{params}),kw_conds{subgroup}',"off") ;
        
        % Run multicomparison post-hoc test with bonferroni correction
        % on the 'stats' structure created by Kruskal-Wallis
        [multcomp_results.(fields_krusk{kwfield}){1, params},... % matrix of multiple comparison results
            multcomp_results.(fields_krusk{kwfield}){2, params},... % matrix of estimates
            multcomp_results.(fields_krusk{kwfield}){3, params},... % handle to the figure
            multcomp_results.(fields_krusk{kwfield}){4, params}] = ... % group names
        multcompare(kruskal_results.(fields_krusk{kwfield}){3, params},...
            "Alpha",0.05,"Display","off","CriticalValueType","bonferroni");
        
        if kwfield <= 5

        % Run Dunn's post-hoc test to benchmark "multcompare"
        [dunn_results.(fields_krusk{kwfield}){1, params},...
            dunn_results.(fields_krusk{kwfield}){2, params},...
            dunn_results.(fields_krusk{kwfield}){3, params}] = ...
        dunn(reshape(kruskal_groups.(fields_krusk{kwfield}).(parameters{params}),1,[]),...
            [ones(1,40) repmat(2,1,40) repmat(3,1,40)],0);

        elseif kwfield > 5
            [dunn_results.(fields_krusk{kwfield}){1, params},...
                dunn_results.(fields_krusk{kwfield}){2, params},...
                dunn_results.(fields_krusk{kwfield}){3, params}] = ...
            dunn(reshape(kruskal_groups.(fields_krusk{kwfield}).(parameters{params}),1,[]),...
                [ones(1,40) repmat(2,1,40) repmat(3,1,40) repmat(4,1,40) repmat(5,1,40)],0);
        end
    end
end

['Performed statistical testing']

save(strcat(destination_folder, '\', run_date, '_statistics'), ...
    'kolmogorov_smirnov', 'kruskal_groups', 'kruskal_results','multcomp_results','dunn_results') ;


%% Back up files and export figures as vector files

% Backing up 'scripts' folder

copyfile('C:\Users\pc\Desktop\Paper amebas 3 (movimiento sistémico)\_Tracks 3600 frames, matlab files and tables\scripts', strcat(destination_folder, '\Scripts_backup'));

% Export figures as vector graphic files (.svg)

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig) ;
%   FigName = get(FigHandle, 'Name') ;
%   set(0, 'CurrentFigure', FigHandle) ;
%   export_fig(strcat(destination_folder,'\',FigName),'-svg')
% %   exportgraphics(FigHandle, fullfile(destination_folder, [FigName '.pdf']), 'Resolution', 600) ;
% end


