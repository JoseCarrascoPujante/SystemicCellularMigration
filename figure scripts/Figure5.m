%% Figure 5
clear
close all
load('coordinates.mat')
load('numerical_results.mat')

%% Layouts
set(groot,'defaultFigurePaperPositionMode','manual')
fig = figure('Visible','off','Position',[0 0 900 1200]);
layout0 = tiledlayout(4,1,'TileSpacing','compact','Padding','none') ;
layout1 = tiledlayout(layout0,2,3,'TileSpacing','compact','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,12,'TileSpacing','loose','Padding','none') ;
layout2.Layout.Tile = 2;
% layout3 = tiledlayout(layout0,2,1,'TileSpacing','compact','Padding','none') ;
% layout3.Layout.Tile = 3;

%% Panel 1 - DFA correlations
scenarios = {"InduccionProteus11_63","QuimiotaxisLeningradensisVariosPpmm","GalvanotaxisBorokensis11_63"};
amoebas = {39,10,48};
for i=1:3 % subpanels (species)
    nexttile(layout1,i)
    box on
    dfahandle = gca;
    gamma = DFA_main2(coordinates.(scenarios{i}).scaled_rho(:,amoebas{i}),'Original_DFA_', dfahandle) ;
    yl = ylim();
    xl = xlim();
    text(xl(1)+1,yl(1)+.25,strcat('\gamma=',num2str(round(gamma,2))))
    nexttile(layout1,i+3)
    box on
    dfahandle = gca;
    gamma = DFA_main2(coordinates.(scenarios{i}).shuffled_rho(:,amoebas{i}),'Shuffled_DFA_', dfahandle) ;
    yl = ylim();
    xl = xlim();
    text(xl(1)+1,yl(1)+.25,strcat('\gamma=',num2str(round(gamma,2))),"FontSize",10)
end

%% Panel 2 - DFA \gamma original vs shuffled

field_names = fieldnames(results) ;
species = {'Proteus','Leningradensis','Borokensis'};
tiles = {
[25,73,49,1;37,85,61,13]
[29,77,53,5;41,89,65,17]
[33,81,57,9;45,93,69,21]};
for i = 1:length(species)
    idx = find(contains(field_names(:),species{i}))';
    for f = 1:length(idx)
        disp(field_names{idx(f)})
        t = nexttile(layout2,tiles{i}(1,f),[1,4]);
        hold on
        exes = zeros(size(results.(field_names{idx(f)}),1));
        plot(results.(field_names{idx(f)})(:,7),exes,'ro','MarkerSize',7)
        plot(results.(field_names{idx(f)})(:,8),exes,'bo','MarkerSize',7)
        box off
        ylim([0 eps]) % minimize y-axis height
        xlim([0 2])
        t.XRuler.TickLabelGapOffset = 2;
        t.YAxis.Visible = 'off'; % hide y-axis
        t.Color = 'None';
        hold off
        
        t2 = nexttile(layout2,tiles{i}(2,f),[1,4]);
        hold on
        datamean = mean(results.(field_names{idx(f)})(:,7));
        datastd = std(results.(field_names{idx(f)})(:,7));
        datameanshuff = mean(results.(field_names{idx(f)})(:,8));
        datastdshuff = std(results.(field_names{idx(f)})(:,8));
        
        % original
        line([datamean-datastd datamean+datastd],[0 0],'Color','red',...
            'LineWidth',.5)
        line([datamean-datastd+.005 datamean-datastd],[0 0],'Color','red',...
            'LineWidth',5)
        line([datamean+datastd datamean+datastd+.005],[0 0],'Color','red',...
            'LineWidth',5)
        text(t2,datamean,-.22,[num2str(round(datamean,2)) ' ' char(177) ' '...
            num2str(round(datastd,2))],'HorizontalAlignment', 'center','FontSize',9)
        
        % shuffled
        line([datameanshuff-datastdshuff datameanshuff+datastdshuff],[0 0],'Color','blue',...
            'LineWidth',.5)
        line([datameanshuff-datastdshuff+.005 datameanshuff-datastdshuff],[0 0],'Color','blue',...
            'LineWidth',5)
        line([datameanshuff+datastdshuff datameanshuff+datastdshuff+.005],[0 0],'Color','blue',...
            'LineWidth',5)
        text(t2,datameanshuff,-.22,[num2str(round(datameanshuff,2)) ' ' char(177)...
            ' ' num2str(round(datastdshuff,2))],'HorizontalAlignment',...
            'center','FontSize',9)

        ylim([-0.08 0]) % minimize y-axis height
        xlim([0 2])
        t2.YAxis.Visible = 'off'; % hide y-axis
        t2.XAxis.Visible = 'off'; % hide y-axis
        t2.Color = 'None';
        hold off
    end
end


%% Panel 3A - DFA \gamma Violin plots original

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

h = nexttile(layout0,3);

species = {'Proteus','Leningradensis','Borokensis'};
col = [.1,.1,.1;.3,.3,.3;.5,.5,.5;.7,.7,.7;1,0,0;1,.25,.25;1,.5,.5; 1,.7,.7;0,0,1;.25,.25,1;.5,.5,1;.7,.7,1];

count = 0;
c = 0;
for i=1:length(species) % main boxes (species)
    count = count + 1;
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        c = c+1;
        count = count+2;
        disp(field_names{f(j)})
        al_goodplot(results.(field_names{f(j)})(:,7), count, [], col(c,:),'bilateral', [], [], 0); %Shuffled
    end
end
xlim([1.5 28])
xticklabels([])
h.XAxis.TickLength = [0 0];
h.YAxis.FontSize = 8;
ylabel('DFA\gamma','FontSize',10)

%% Panel 3B - DFA \gamma Violin plots shuffled

h = nexttile(layout0,4);
count = 0;
c = 0;
for i=1:length(species) % main boxes (species)
    disp(species(i))
    count = count + 1;
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        c = c+1;
        count = count+2;
        disp(field_names{f(j)})
        al_goodplot(results.(field_names{f(j)})(:,8), count, [], col(c,:),'bilateral', [], [], 0); % original
    end
end
xlim([1.5 28])
xticks([5.9 15 24])
h.XAxis.TickLength = [0 0];
h.YAxis.FontSize = 8;
ylabel('DFA\gamma (Shuffled)','FontSize',10)
set(h,'XTickLabel',[{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])


%% Export as jpg, tiff and vector graphics pdf

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig5'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig5 files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig5(',num2str(gabs),')'),'svg')
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig5(',num2str(gabs),').jpg') ...
%   ,"Resolution",600)
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig5(',num2str(gabs),').tiff') ...
%   ,"Resolution",600)
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig5(',num2str(gabs),').pdf'), ...
%   'BackgroundColor','white', 'ContentType','vector')

