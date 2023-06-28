%% Figure 6
clear 
close all
set(groot,'defaultFigurePaperPositionMode','manual')

fig = figure('Visible','off','Position', [0 0 990 700]);
layout0 = tiledlayout(3,2,'TileSpacing','compact','Padding','none') ;


%% Boxcharts with scatterplots

% Intensity of response (mm)
ax = nexttile(layout0,1);
kyneticBoxplots(ax,13)

% Shuffled Intensity of response (mm)
ax = nexttile(layout0,2);
kyneticBoxplots(ax,14)

% Directionality ratio (straightness)
ax = nexttile(layout0,3);
kyneticBoxplots(ax,15)

% Shuffled Directionality ratio (straightness)
ax = nexttile(layout0,4);
kyneticBoxplots(ax,16)

% Average speed (mm/s)
ax = nexttile(layout0,5);
kyneticBoxplots(ax,17)

% Shuffled Average speed (mm/s)
ax = nexttile(layout0,6);
kyneticBoxplots(ax,18)


%% Export as jpg and vector graphics pdf

load '2023-06-07_14.16''19''''_coordinates.mat' destination_folder 

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig6'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig6 files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig6(',num2str(gabs),')'),'svg')


%% Define functions
function kyneticBoxplots(ax,parameter_index)
    
    load '2023-06-07_14.16''19''''_numerical_results.mat' results
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
    
    species_keyword = {'Proteus','Leningradensis','Borokensis'};
    values = {[],[],[]};
    for i=1:length(species_keyword) % species    
        for f = find(contains(field_names(:),species_keyword(i)))' % conditions
            disp(field_names{f})
            values{i} = [values{i}; results.(field_names{f})(:,parameter_index)];
        end
    end
    values_conds = {{[],[],[],[]},{[],[],[],[]},{[],[],[],[]}};
    for i=1:length(species_keyword) % main boxes (species)
        f = find(contains(field_names(:),species_keyword(i)))'; % condition indexes
        for j = 1:length(f) % secondary boxes (conditions)
            disp(field_names{f(j)})
            values_conds{i}{j} = results.(field_names{f(j)})(:,parameter_index);
        end
    end
    
    % Plot boxchart
    boxplot_values=cat(1,values{1},values{2},values{3});
    
    groups{1} = [repmat({'Amoeba proteus'},length(values{1}),1); ...
        repmat({'Metamoeba leningradensis'},length(values{2}),1); ...
        repmat({'Amoeba borokensis'},length(values{3}),1)];
    
    groups{2} = [repmat({'Sc1'},length(values_conds{1}{1}),1);
        repmat({'Sc2'},length(values_conds{1}{2}),1);
        repmat({'Sc3'},length(values_conds{1}{3}),1);
        repmat({'Sc4'},length(values_conds{1}{4}),1);
        repmat({'Sc1'},length(values_conds{2}{1}),1);
        repmat({'Sc2'},length(values_conds{2}{2}),1);
        repmat({'Sc3'},length(values_conds{2}{3}),1);
        repmat({'Sc4'},length(values_conds{2}{4}),1);
        repmat({'Sc1'},length(values_conds{3}{1}),1);
        repmat({'Sc2'},length(values_conds{3}{2}),1);
        repmat({'Sc3'},length(values_conds{3}{3}),1);
        repmat({'Sc4'},length(values_conds{3}{4}),1)];
    
    idxX = [repmat(1.6,length(values_conds{1}{1}),1);
        repmat(2.8,length(values_conds{1}{2}),1);
        repmat(4,length(values_conds{1}{3}),1);
        repmat(5.2,length(values_conds{1}{4}),1);
        repmat(7.4,length(values_conds{2}{1}),1);
        repmat(8.6,length(values_conds{2}{2}),1);
        repmat(9.8,length(values_conds{2}{3}),1);
        repmat(11,length(values_conds{2}{4}),1);
        repmat(13.15,length(values_conds{3}{1}),1);
        repmat(14.4,length(values_conds{3}{2}),1);
        repmat(15.6,length(values_conds{3}{3}),1);
        repmat(16.8,length(values_conds{3}{4}),1)];
    
    c = [repmat([.1 .1 .1],length(values_conds{1}{1}),1);
        repmat([.25 .25 .25],length(values_conds{1}{2}),1);
        repmat([.5 .5 .5],length(values_conds{1}{3}),1);
        repmat([.6 .6 .6],length(values_conds{1}{4}),1);
        repmat([1 0 0],length(values_conds{2}{1}),1);
        repmat([1 .25 .25],length(values_conds{2}{2}),1);
        repmat([1 .5 .5],length(values_conds{2}{3}),1);
        repmat([1 .6 .6],length(values_conds{2}{4}),1);
        repmat([0 0 1],length(values_conds{3}{1}),1);
        repmat([.25 .25 1],length(values_conds{3}{2}),1);
        repmat([.5 .5 1],length(values_conds{3}{3}),1);
        repmat([.6 .6 1],length(values_conds{3}{4}),1)];
    
    %%% Box plot
    hold(ax,"on")
    b = boxplot(boxplot_values,groups,'Notch','off','LabelVerbosity','majorminor',...
        'FactorGap',[10,2],'colors',[.1 .1 .1;.25 .25 .25;.5 .5 .5;.6 .6 .6;1 0 0;1 .25 .25;1 .5 .5;...
        1 .6 .6;0 0 1;.25 .25 1;.5 .5 1;.6 .6 1],'Symbol','');
    set(b(5:7:end),'LineWidth',.8);
    set(b(6:7:end),'LineWidth',.8);
    % h=findobj('LineStyle','--'); set(h, 'LineStyle','-'); % Make whiskers a solid line
    
    scatter(idxX',boxplot_values,5,c,'jitter','off','jitterAmount',0);
    
    meanVal = groupsummary(boxplot_values,idxX,'mean');
    
    %%% Plot mean values
    plot(unique(idxX)-.6,meanVal,'--og',...
        'LineWidth',.8,...
        'MarkerSize',6,...
        'MarkerEdgeColor','g')
    axis auto
    yl = ylim;
    xlim(ax,[0 17])
    ylim(ax,[-yl(2)*.05 yl(2)*1.01])
    box off
end