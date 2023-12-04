%% Figure 6
clear 
close all
set(groot,'defaultFigurePaperPositionMode','manual')

fig = figure('Visible','off','Position', [0 0 550 750]);
layout0 = tiledlayout(3,1,'TileSpacing','loose','Padding','none') ;


%% Boxcharts with scatterplots

% Intensity of response (mm)
ax = nexttile(layout0);
kyneticBoxplots(ax,13,'Intensity of response (mm)')

% % Shuffled Intensity of response (mm)
% ax = nexttile(layout0,2);
% kyneticBoxplots(ax,14,'Shuffled Intensity of response (mm)')

% Directionality ratio (straightness)
ax = nexttile(layout0);
kyneticBoxplots(ax,15,'Directionality ratio')

% % Shuffled Directionality ratio (straightness)
% ax = nexttile(layout0,4);
% kyneticBoxplots(ax,16,'Shuffled Directionality ratio')

% Average speed (mm/s)
ax = nexttile(layout0);
kyneticBoxplots(ax,17,'Average speed (mm/s)')

% % Shuffled Average speed (mm/s)
% ax = nexttile(layout0,6);
% kyneticBoxplots(ax,18,'Shuffled average speed (mm/s)')


%% Export as jpg and vector graphics svg

load 'coordinates.mat' destination_folder 

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
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
saveas(fig,strcat(destination_folder, '\Figures\Fig6(',num2str(gabs+1),')'),'svg')


%% Define functions
function kyneticBoxplots(ax,parameter_index,string)
% Grouped boxplot chart of any of the 18 metrics
    
    load 'numerical_results.mat' results
    
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
    
    Xidx = [repmat(1.6,length(values_conds{1}{1}),1);
        repmat(2.8,length(values_conds{1}{2}),1);
        repmat(4.04,length(values_conds{1}{3}),1);
        repmat(5.25,length(values_conds{1}{4}),1);
        repmat(7.35,length(values_conds{2}{1}),1);
        repmat(8.6,length(values_conds{2}{2}),1);
        repmat(9.8,length(values_conds{2}{3}),1);
        repmat(11.03,length(values_conds{2}{4}),1);
        repmat(13.12,length(values_conds{3}{1}),1);
        repmat(14.35,length(values_conds{3}{2}),1);
        repmat(15.58,length(values_conds{3}{3}),1);
        repmat(16.8,length(values_conds{3}{4}),1)];
    
    colours = [repmat([.1 .1 .1],length(values_conds{1}{1}),1);
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
    set(findobj(ax,'Type','text'),'FontSize',9);
    set(b(5:7:end),'LineWidth',.8);
    set(b(6:7:end),'LineWidth',.8);
    % h=findobj('LineStyle','--'); set(h, 'LineStyle','-'); % Make whiskers a solid line
    
    scatter(Xidx',boxplot_values,5,colours,'jitter','off','jitterAmount',0);
    
    meanVal = groupsummary(boxplot_values,Xidx,'mean');
    
    %%% Plot mean values
    plot(unique(Xidx)-.6,meanVal,'--og',...
        'LineWidth',.8,...
        'MarkerSize',6,...
        'MarkerEdgeColor','g')
    axis auto
    yl = ylim;
    xlim(ax,[0 17])
    ylim(ax,[-yl(2)*.05 yl(2)])
    box off
    ylabel(string,"FontSize",9)
end