%% Figure 2

clear
close all
load('2023-06-07_14.16''19''''_coordinates.mat')
load('2023-06-07_14.16''19''''_numerical_results.mat')

%% Layouts

set(groot,'defaultFigurePaperPositionMode','manual')
fig = figure('Visible','off','Position',[0 0 900 750]);
layout0 = tiledlayout(3,1,'TileSpacing','loose','Padding','none') ;
layout1 = tiledlayout(layout0,1,3,'TileSpacing','loose','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,3,'TileSpacing','compact','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - RMSF max_correlations

fields = {"InduccionProteus11_63","InduccionLeningradensis11_63","QuimiotaxisBorokensis23_44"};
amoebas = {8,55,44};
for i=1:3 % subpanels (species)
    nexttile(layout1)   
    rmsfhandle = gca;
    set(rmsfhandle,'xscale','log')
    set(rmsfhandle,'yscale','log')
    amebas5(coordinates.(fields{i}).scaled_rho(:,amoebas{i}),rmsfhandle,'orig') ;
end

%% Panel 2 - RMSF \alpha

species = {'Proteus','Leningradensis','Borokensis'};
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

dataSpecies = {[] [] []};
dataSpeciesShuff = {[] [] []};
tiles = {
[1,7,13,19;4,10,16,22]
[2,8,14,20;5,11,17,23]
[3,9,15,21;6,12,18,24]
};

for i = 1:length(species)
    idx = find(contains(field_names(:),species{i}))';
    for f = 1:length(idx)
        t = nexttile(layout2,tiles{i}(1,f));
        hold on
               
        exes = zeros(size(results.(field_names{idx(f)}),1));
        plot(results.(field_names{idx(f)})(:,1),exes,'ro','MarkerSize',7)
        plot(results.(field_names{idx(f)})(:,2),exes,'bo','MarkerSize',7)
        ylim([0 eps]) % minimize y-axis height
        xlim([.35 .87])
        t.YAxis.Visible = 'off'; % hide y-axis
        t.Color = 'None';
        hold off

        t2 = nexttile(layout2,tiles{i}(2,f));
        hold on
        datamean = mean(results.(field_names{idx(f)})(:,1));
        datastd = std(results.(field_names{idx(f)})(:,1));
        datameanshuff = mean(results.(field_names{idx(f)})(:,2));
        datastdshuff = std(results.(field_names{idx(f)})(:,2));
        line([datamean-datastd datamean+datastd],[0 0],'Color','red',...
            'LineWidth',.5)
        line([datamean-datastd+.001 datamean-datastd],[0 0],'Color','red',...
            'LineWidth',5)
        line([datamean+datastd datamean+datastd+.001],[0 0],'Color','red',...
            'LineWidth',5)
        text(t2,datamean,-1.3,[num2str(round(datamean,2)) ' ' char(177) ' '...
            num2str(round(datastd,2))],'HorizontalAlignment','center','FontSize',8)
        line([datameanshuff-datastdshuff datameanshuff+datastdshuff],[0 0],'Color','blue',...
            'LineWidth',.5)
        line([datameanshuff-datastdshuff+.001 datameanshuff-datastdshuff],[0 0],'Color','blue',...
            'LineWidth',5)
        line([datameanshuff+datastdshuff datameanshuff+datastdshuff+.001],[0 0],'Color','blue',...
            'LineWidth',5)
        text(t2,datameanshuff,-1.3,[num2str(round(datameanshuff,2)) ' ' char(177)...
            ' ' num2str(datastdshuff,'%.e')],'HorizontalAlignment','center','FontSize',8)
        ylim([-0.08 0]) % minimize y-axis height
        xlim([.35 .87])
        t2.YAxis.Visible = 'off'; % hide y-axis
        t2.XAxis.Visible = 'off'; % hide y-axis
        t2.Color = 'None';
        hold off
    end
end

%% Panel 3 - RMSF Violin plots
ax=nexttile(layout0,3);
hold on

h=gca;
xlim([.5 3.5])
ylim([-8.35 30])
h.XTick = [1 2 3];
xticklabels([{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])
h.XAxis.TickLength = [0 0];
ylabel('Memory persistence (min)')

%%% RainCloud plots
data = cell(3,4);
for i=1:length(species) % species    
    f = find(contains(field_names(:),species(i))); % conditions
    for j = 1:length(f)
        data{i,j} = [data{i}; results.(field_names{f(j)})(:,5)/120];
    end
end
cb = [.2,.2,.2;.4,.4,.4;.6,.6,.6;.8,.8,.8;
    1,0,0;1,.25,.25;1,.5,.5; 1,.75,.75;
    0,0,1;.25,.25,1;.5,.5,1;.75,.75,1];

% plot
count = 0;
for p = 1:(size(data,1)) % species
    for q = 1:size(data,2) % conditions
        count = count+1;
        plot_rainclouds(data{p,q},cb,count); % impossible to plot more than one of these plots in the same axes
    end
end

%% Export figures as jpg and svg

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig2'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig2 files found'))

% Save main figure (the tiled layout)
fig.Units = 'centimeters';                % set figure units to cm
fig.PaperUnits = 'centimeters';           % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig2(',num2str(gabs+1),')'),'svg')
close(fig)

% Save the remaining figures
FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
for iFig = length(FigList):-1:1
    FigHandle = FigList(iFig) ;
    FigName = get(FigHandle, 'Name') ;
    set(0, 'CurrentFigure', FigHandle) ;
    FigHandle.Units = 'centimeters';                % set figure units to cm
    FigHandle.PaperUnits = 'centimeters';           % set pdf printing paper units to cm
    FigHandle.PaperSize = FigHandle.Position(3:4);  % assign to the printing paper the size of the figure
    FigHandle.PaperPosition = [0 0 FigHandle.Position(3:4)];
    set(FigHandle, 'Renderer', 'painters');
    saveas(FigHandle,strcat(destination_folder, '\Figures\Fig2(',num2str(iFig+gabs),')'),'svg')
end

%% Define functions

function plot_rainclouds(data,cb,count)
    figure('Visible','off','Position', [0, 0, 110 250])
    raincloud_plot(data,'box_on',1,'box_dodge',1,...
        'box_dodge_amount',-.08,'dot_dodge_amount',.2,'alpha',0.2,...
        'bxcl',[0 0 0],'color',cb(count,:),'line_width',1);
    view([90 -90])
    xlim([-8.35 30])
    ylim([-.03 .095])
    box off
    AxesH = gca;
    set(AxesH, 'Units', 'pixels', 'Position', [0, 0, 110 250]);
    AxesH.XTick = [];
    AxesH.YTick = [];
    % axis padded
end