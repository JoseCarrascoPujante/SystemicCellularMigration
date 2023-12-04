clear
close all

diary off
diary_filename = strcat(destination_folder,'\DFAvalues.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on
tic
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
bar1 = waitbar(0,'In progress...','Name','Condition...') ;
bar2 = waitbar(0,'In progress...','Name','Track number...') ;

for i = 1:length(field_names)
    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
    for j = 1:N
		disp(strcat(field_names{i},'nº',num2str(j)))
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        figure('Name',strcat('DFA_Original_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','on');
        dfahandle = gca;       
        gammaO = DFA_main2(coordinates.(field_names{i}).scaled_rho(:,j),'Original_DFA_', dfahandle) ;
        disp(strcat('Original Gamma=',num2str(gammaO)))
        figure('Name',strcat('DFA_Shuffled_',field_names{i}, '_amoeba_number_',...
            num2str(j)),'NumberTitle','off','Visible','on');
        dfahandle = gca;       
        gammaSh = DFA_main2(coordinates.(field_names{i}).shuffled_rho(:,j),'Shuffled_DFA_', dfahandle) ;
        disp(strcat('Shuffled Gamma=',num2str(gammaSh)))
    end
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig) ;
  FigName = get(FigHandle, 'Name') ;
  set(0, 'CurrentFigure', FigHandle) ;
  exportgraphics(FigHandle,fullfile(destination_folder, '\Figures', [FigName '.jpg']), ...
    "Resolution",200)
end
toc

diary off
