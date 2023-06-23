diary off
diary_filename = strcat(destination_folder,'\RMSFvalues.txt') ;
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
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        
        figure('Name',strcat('RMSF_Original_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','off');   
        hrmsfOr = gca;
        [alpha,R2,CorrelationTime, ~] = amebas5(coordinates.(field_names{i}).scaled_rho(:,j), hrmsfOr) ;
        
        
        figure('Name',strcat('RMSF_Shuffled_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','off');   
        hrmsfshuff = gca;
        [Salpha,SR2,SCorrelationTime, ~] = amebas5(coordinates.(field_names{i}).shuffled_rho(:,j), hrmsfOr) ;

        [field_names{i} 'nº' num2str(j) ':' newline 'alpha:' num2str(alpha)...
            newline 'R2:' num2str(R2) newline 'MaxCorrTime:' num2str(CorrelationTime)...
            newline 'Shuff alpha:' num2str(Salpha) newline 'Shuff R2:' num2str(SR2) newline...
            'Shuff MaxCorrTime:' num2str(SCorrelationTime)]
    end
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig) ;
  FigName = get(FigHandle, 'Name') ;
  set(0, 'CurrentFigure', FigHandle) ;
  exportgraphics(FigHandle,fullfile(destination_folder,'\Figures', [FigName '.jpg']), ...
    "Resolution",300)
end
toc

diary off

