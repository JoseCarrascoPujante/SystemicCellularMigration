clear
close all

load coordinates.mat destination_folder coordinates

diary off
diary_filename = strcat(destination_folder,'\MSDvalues.txt') ;
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

for i = [1,5,9]%1:length(field_names)
    figure('Name',strcat('MSD_Original_',field_names{i}),'Visible','on','NumberTitle','off') ;
    xlabel('Log(MSD(\tau))');
    ylabel('Log(\tau(s))');
    msdOr = gca;

    % figure('Name',strcat('MSD_Shuffled_',field_names{i}),'Visible','on','NumberTitle','off') ;
    % xlabel('Log(MSD(\tau))');
    % ylabel('Log(\tau(s))');
    % msdSh = gca;

    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition

    for j = 1:N
		disp(strcat(field_names{i},'nº',num2str(j)))
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        
        betaO = msd(coordinates.(field_names{i}).scaled_x(:,j),coordinates.(field_names{i}).scaled_y(:,j),msdOr) ;
        disp(strcat('Original Beta=',num2str(betaO)))

        % betaSh = msd(coordinates.(field_names{i}).shuffled_x(:,j),coordinates.(field_names{i}).shuffled_y(:,j),msdSh,'shuff') ;
        % disp(strcat('Shuffled Beta=',num2str(betaSh)))
    end
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig) ;
  FigName = get(FigHandle, 'Name') ;
  set(0, 'CurrentFigure', FigHandle) ;
  exportgraphics(FigHandle, fullfile(destination_folder, '\Figures', [FigName '.jpg']), ...
    "Resolution", 200)
end
toc

diary off