field_names = fieldnames(coordinates) ;

% List the parameters to be calculated by the script
stat_names = {'RMSF_alpha', 'sRMSF_alpha', 'RMSF_R2', 'sRMSF_R2', 'RMSFCorrelationTime', ...
    'sRMSFCorrelationTime', 'DFA_gamma', 'sDFA_gamma', 'MSD_beta', 'sMSD_beta', 'AppEn', 'sAppEn'} ;

for i = 1:length(field_names)
     N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
     for j = 1:N     
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