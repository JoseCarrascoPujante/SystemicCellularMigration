%% Figure 2

clear
close all
load('coordinates.mat')
load('numerical_results.mat')

%% Layouts

set(groot,'defaultFigurePaperPositionMode','manual')
fig = figure('Visible','off','Position',[0 0 1000 833]);
layout0 = tiledlayout(3,1,'TileSpacing','tight','Padding','none') ;
layout1 = tiledlayout(layout0,1,3,'TileSpacing','loose','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,3,'TileSpacing','loose','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - RMSF max_correlations

fields = {"InduccionProteus11_63","InduccionLeningradensis11_63","QuimiotaxisBorokensis23_44"};
amoebas = {8,55,44};
for i=1:3 % subpanels (species)
    nexttile(layout1)   
    rmsfhandle = gca;
    set(rmsfhandle,'xscale','log')
    set(rmsfhandle,'yscale','log')
    conu.amebas5(coordinates.(fields{i}).scaled_rho(:,amoebas{i}),rmsfhandle,'orig') ;
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
        ylim([0 eps]) % minimize y-axis height
        xlim([.5 1])
        t.XRuler.TickLabelGapOffset = 2;
        t.YAxis.Visible = 'off'; % hide y-axis
        t.Color = 'None';
        hold off

        t2 = nexttile(layout2,tiles{i}(2,f));
        hold on
        datamean = mean(results.(field_names{idx(f)})(:,1));
        datastd = std(results.(field_names{idx(f)})(:,1));

        line([datamean-datastd datamean+datastd],[0 0],'Color','red',...
            'LineWidth',.5)
        line([datamean-datastd+.001 datamean-datastd],[0 0],'Color','red',...
            'LineWidth',5)
        line([datamean+datastd datamean+datastd+.001],[0 0],'Color','red',...
            'LineWidth',5)
        text(t2,datamean,-.75,[num2str(round(datamean,2)) ' ' char(177) ' '...
            num2str(round(datastd,2))],'HorizontalAlignment','center','FontSize',10)
        ylim([-0.08 0]) % minimize y-axis height
        xlim([.5 1])
        t2.XRuler.TickLabelGapOffset = 2;
        t2.YAxis.Visible = 'off'; % hide y-axis
        t2.XAxis.Visible = 'off'; % hide y-axis
        t2.Color = 'None';
        hold off
    end
end

%% Panel 3 - RMSF Violin plots
ax=nexttile(layout0,3);
rmsf_conds = {{[],[],[],[]},{[],[],[],[]},{[],[],[],[]}};

for i=1:length(species) % main boxes (species)
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        rmsf_conds{i}{j} = results.(field_names{f(j)})(:,5)/120;
    end
end

for i=1:length(species) % main boxes (species)
    conu.superviolin(rmsf_conds{i},'Parent',ax,'Xposition',i,'FaceAlpha',0.15,...
        'Errorbars','ci','Centrals','mean','LineWidth',0.1)
end
colorgroups = [ones(length(rmsf_conds{1}{1}),1);
    repmat(2,length(rmsf_conds{1}{2}),1);
    repmat(3,length(rmsf_conds{1}{3}),1);
    repmat(4,length(rmsf_conds{1}{4}),1);
    ones(length(rmsf_conds{2}{1}),1);
    repmat(2,length(rmsf_conds{2}{2}),1);
    repmat(3,length(rmsf_conds{2}{3}),1);
    repmat(4,length(rmsf_conds{2}{4}),1);
    ones(length(rmsf_conds{3}{1}),1);
    repmat(2,length(rmsf_conds{3}{2}),1);
    repmat(3,length(rmsf_conds{3}{3}),1);
    repmat(4,length(rmsf_conds{3}{4}),1)];

rmsfs = {[],[],[]};
for i=1:length(species) % species    
    for f = find(contains(field_names(:), species(i)))' % conditions
        rmsfs{i} = [rmsfs{i}; results.(field_names{f})(:,5)/120];
    end
end

boxChart_rmsf = cat(1, rmsfs{1}, rmsfs{2}, rmsfs{3});
boxchart([ones(length(rmsfs{1}), 1); repmat(2, length(rmsfs{2}), 1); ...
    repmat(3, length(rmsfs{3}), 1)], boxChart_rmsf, 'Notch', 'off',...
    'GroupByColor', colorgroups, 'BoxFaceAlpha',0) %Box charts whose notches do not overlap have different medians at the 5% significance level.
h=gca;
xlim([.5 3.5])
ylim([0 20])
% h.XAxisLocation = 'top';
box on
h.XTick = [1 2 3];
xticklabels(h,[{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])
h.XAxis.FontSize = 16;
h.XAxis.TickLength = [0 0];
ylabel('Correlations duration (min)')



%% Export figures as jpg and svg

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig2'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig2 files found'))

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
