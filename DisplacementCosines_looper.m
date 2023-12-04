clear
close all

load numerical_results.mat
load coordinates.mat destination_folder coordinates

diary off
diary_filename = strcat(destination_folder,'\DispCosvalues.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on

tic
field_names = fieldnames(coordinates) ;
% List the parameters to be calculated by the script
stat_names = {'RMSF_alpha', 'sRMSF_alpha', 'RMSF_R2', 'sRMSF_R2', 'RMSFCorrelationTime', ...
    'sRMSFCorrelationTime', 'DFA_gamma', 'sDFA_gamma', 'MSD_beta', 'sMSD_beta', 'ApEn', ...
    'sApEn','Intensity','sIntensity','DR','sDR','AvgSpeed','sAvgSpeed','DispCos'} ;

bar1 = waitbar(0,'In progress...','Name','Condition...') ;
bar2 = waitbar(0,'In progress...','Name','Track number...') ;

for i = 1:length(field_names)
    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
    for j = 1:N
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
      
        results.(field_names{i})(j,strcmp(stat_names(:), 'DispCos')) = ...
            cos(coordinates.(field_names{i}).theta(end,j)) ;
        
        [field_names{i} 'nº' num2str(j) 'displacement cosine:' ...
            num2str(results.(field_names{i})(j,strcmp(stat_names(:), 'DispCos')))]
    end
end

toc

diary off

% save(strcat(destination_folder, '\', 'numerical_results.mat'),...
% 'tCalcSec', 'results', 'field_names') ;
