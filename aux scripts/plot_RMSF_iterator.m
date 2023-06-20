diary off
diary_filename = strcat(destination_folder,'\rmsfGoodness.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on
tic
field_names = fieldnames(coordinates) ;
% List the parameters to be calculated by the script
stat_names = {'RMSF_alpha', 'sRMSF_alpha', 'RMSF_R2', 'sRMSF_R2', 'RMSFCorrelationTime', ...
    'sRMSFCorrelationTime', 'DFA_gamma', 'sDFA_gamma', 'MSD_beta', 'sMSD_beta', 'AppEn', 'sAppEn'} ;

bar1 = waitbar(0,'In progress...','Name','Condition...') ;
bar2 = waitbar(0,'In progress...','Name','Track number...') ;

for i = 1:length(field_names)
    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
    for j = 1:N
		disp(strcat(field_names{i},'nº',num2str(j)))
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        figure('Name',strcat('RMSF_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','off');
        
        rmsfhandle = gca;
        set(rmsfhandle,'xscale','log')
        set(rmsfhandle,'yscale','log')
        
        [results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_alpha')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'RMSF_R2')),...
            results.(field_names{i})(j,strcmp(stat_names(:), 'RMSFCorrelationTime')), tc2] = ...
        amebas5(coordinates.(field_names{i}).scaled_rho(:,j), rmsfhandle) ;
    end
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig) ;
  FigName = get(FigHandle, 'Name') ;
  set(0, 'CurrentFigure', FigHandle) ;
  exportgraphics(FigHandle,fullfile(destination_folder, [FigName '.jpg']), ...
    "Resolution",300)
end
diary off

toc