clear
close all

load coordinates.mat coordinates destination_folder stat_names
load numerical_results.mat results

field_names = ...
    {'SinEstimuloProteus11_63'
    'GalvanotaxisProteus11_63'
    'QuimiotaxisProteus11_63'
    'InduccionProteus11_63'
    'SinEstimuloLeningradensis11_63'
    'GalvanotaxisLeningradensis11_63'
    'QuimiotaxisLeningradensisVariosPpmm'
    'InduccionLeningradensis11_63'
    'SinEstimuloBorokensis23_44'
    'GalvanotaxisBorokensis11_63'
    'QuimiotaxisBorokensis23_44'
    'InduccionBorokensis11_63'
    };

diary off
diary_filename = strcat(destination_folder,'\RMSFvalues.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on
tic

bar1 = waitbar(0,'In progress...','Name','Condition...') ;
bar2 = waitbar(0,'In progress...','Name','Track number...') ;

for i = 1:length(field_names)
    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
    for j = 1:N
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        
        figure('Name',strcat('RMSF_Original_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','off');   
        hrmsfOr = gca;
        [results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_alpha')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_R2')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'RMSFCorrelationTime')), ~] =...
            amebas5(coordinates.(field_names{i}).scaled_rho(:,j), hrmsfOr, 'orig') ;
        
        
        figure('Name',strcat('RMSF_Shuffled_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','off');   
        hrmsfshuff = gca;
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSF_alpha')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSF_R2')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSFCorrelationTime')), tc2] =...
            amebas5(coordinates.(field_names{i}).shuffled_rho(:,j), hrmsfshuff,'shuff');

        % Print RMSF values
        [field_names{i} 'nº' num2str(j) ':' ...
            newline 'alpha:' num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_alpha')))...
            newline 'R2:' num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_R2')))...
            newline 'MaxCorrTime:' num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'RMSFCorrelationTime')))...
            newline 'Shuff alpha:' num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSF_alpha'))) ...
            newline 'Shuff R2:' num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSF_R2'))) ...
            newline 'Shuff MaxCorrTime:' num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'sRMSFCorrelationTime')))]
    end
end

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig) ;
%   FigName = get(FigHandle, 'Name') ;
%   set(0, 'CurrentFigure', FigHandle) ;
%   exportgraphics(FigHandle,fullfile(destination_folder,'\Figures', [FigName '.jpg']), ...
%     "Resolution",300)
% end

toc

diary off

