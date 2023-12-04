
function fig = all2DscatterplotCombinations(field_names, results)

% Plots systemic movement parameters in 2D and 3D scatterplots

% Group the 12 track-wise statistics in one table

results.full = [];

for index = 1:length(field_names)

        results.full = cat(1,results.full, results.(field_names{index})); % field name order does not matter
end

results.full(:,[5,6]) = results.full(:,[5,6])/120; % convert frames to minutes

stat_names = ["RMSF\alpha" "sRMSF\alpha" "rmsfR2" "srmsfR2" "rmsfTimeMax" ...
    "srmsfTimeMax" "DFA\gamma" "sdfaGamma" "MSD\beta" "smsdBeta"...
    "Approximate Entropy" "sApproximate Entropy"];
indexes = 1:length(stat_names);
pairs = nchoosek(indexes(1:2:end),2);

% 2D scatters

fig = figure('Visible','off','Position',[0 0 1300 750],'NumberTitle','off') ;
fig.InvertHardcopy = 'off';
tiledlayout(3,5,'TileSpacing','tight','Padding','tight')

for ej=1:length(pairs)
    nexttile;
    disp(ej)
    
    % if pairs(ej,1) == 7   % if X-axis metric is dfa\Gamma substract 1, else do not
    %     metric1 = cat(1,results.full(:,pairs(ej,1))-1, results.full(:,pairs(ej,1)+1));
    % else 
    metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
    % end
    % if pairs(ej,2) == 7   % if Y-axis metric is dfa\Gamma substract 1, else do not
    %     metric2 = cat(1,results.full(:,pairs(ej,2))-1, results.full(:,pairs(ej,2)+1));
    % else
    metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
    % end
    G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
    gscatter(metric1,metric2, G,'gy','..',2,'off') ;
    disp(strcat(stat_names(pairs(ej,1)),' ','vs',' ',stat_names(pairs(ej,2))))
    set(gca, 'Color','k', 'XColor','w', 'YColor','w')
    set(gcf, 'Color','k') 
    precision_vals = 1;
    precision_colrs = 'r';
    precision_type = 'STD';
    for c = 1:length(precision_vals)
        hold on
        ellipse(cat(2,metric1,metric2),G,precision_vals(c),precision_colrs(c),precision_type)
    end
    xlabel(stat_names(pairs(ej,1)))
    ylabel(stat_names(pairs(ej,2)))
end

% h(1) = plot(nan, nan, 'go', 'MarkerFaceColor','g', 'MarkerSize', 10, 'DisplayName', 'Original');
% h(2) = plot(nan, nan, 'yo', 'MarkerFaceColor','y', 'MarkerSize', 10, 'DisplayName', 'Shuffled');
% h(3) = plot(nan, nan, 'w--', 'MarkerSize', 25, 'DisplayName', '0.5xSTD');
% h(4) = plot(nan, nan, 'c--', 'MarkerSize', 25, 'DisplayName', '0.75xSTD');
h(1) = plot(nan, nan, 'r--', 'MarkerSize', 25, 'DisplayName', strcat(num2str(precision_vals),'xSTD'));
% h(6) = plot(nan, nan, 'm--', 'MarkerSize', 25, 'DisplayName', '1.5xSTD');
leg=legend(h,Orientation='Horizontal',TextColor='w',FontSize=15);
leg.Layout.Tile = 'north';

%% Export as jpg and vector graphics svg

load '2023-06-07_14.16''19''''_coordinates.mat' destination_folder 

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures')) ;
gabs = 1 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'All2Dscats'+wildcardPattern+'.jpg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs-1),' All2Dscats files found'))

fig.Units = 'centimeters';          % set figure units to cm
fig.PaperUnits = 'centimeters';     % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\All2Dscats(',num2str(gabs),')'),'jpg')

end