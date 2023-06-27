% Figure 3
%% Layouts
set(groot,'defaultFigurePaperPositionMode','manual')
fig = figure('Visible','off','Position',[0 0 900 1200]);
layout0 = tiledlayout(3,1,'TileSpacing','tight','Padding','none') ;
layout1 = tiledlayout(layout0,2,3,'TileSpacing','none','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,12,'TileSpacing','tight','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - DFA correlations
scenarios = {"InduccionProteus11_63","QuimiotaxisLeningradensisVariosPpmm","GalvanotaxisBorokensis11_63"};
amoebas = {39,10,51};
for i=1:3 % subpanels (species)
    nexttile(layout1,i)
    box on
    dfahandle = gca;
    gamma = DFA_main2(coordinates.(scenarios{i}).scaled_rho(:,amoebas{i}),'Original_DFA_', dfahandle) ;
    yl = ylim();
    xl = xlim();
    text(xl(1)+1,yl(1)+.5,strcat('\gamma=',num2str(round(gamma,2))))
    nexttile(layout1,i+3)
    box on
    dfahandle = gca;
    gamma = DFA_main2(coordinates.(scenarios{i}).shuffled_rho(:,amoebas{i}),'Shuffled_DFA_', dfahandle) ;
    yl = ylim();
    xl = xlim();
    text(xl(1)+1,yl(1)+.5,strcat('\gamma=',num2str(round(gamma,2))))
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
        t.YAxis.Visible = 'off'; % hide y-axis
        t.Color = 'None';
        hold off
        
        t2 = nexttile(layout2,tiles{i}(2,f),[1,4]);
        hold on
        datamean = mean(results.(field_names{idx(f)})(:,7));
        datastd = std(results.(field_names{idx(f)})(:,7));
        datameanshuff = mean(results.(field_names{idx(f)})(:,8));
        datastdshuff = std(results.(field_names{idx(f)})(:,8));
        line([datamean-datastd datamean+datastd],[12 12],'Color','red',...
            'LineWidth',.5)
        text(t2,datamean,-10,[num2str(round(datamean,2)) ' ' char(177) ' '...
            num2str(round(datastd,2))],'HorizontalAlignment', 'center','FontSize',9)
        line([datameanshuff-datastdshuff datameanshuff+datastdshuff],[12 12],'Color','blue',...
            'LineWidth',.5)
        text(t2,datameanshuff,-10,[num2str(round(datameanshuff,2)) ' ' char(177)...
            ' ' num2str(round(datastdshuff,2))],'HorizontalAlignment',...
            'center','FontSize',9)
        ylim([6 12]) % minimize y-axis height
        xlim([0 2])
        t2.YAxis.Visible = 'off'; % hide y-axis
        t2.XAxis.Visible = 'off'; % hide y-axis
        t2.Color = 'None';
        hold off
    end
end

%% Panel 3 - DFA \gamma Violin plots

h = nexttile(layout0,3);

col = [.2,.2,.2;.4,.4,.4;.6,.6,.6;.8,.8,.8;1,0,0;1,.25,.25;1,.5,.5; 1,.75,.75;0,0,1;.25,.25,1;.5,.5,1;.75,.75,1];
count = 0;
c = 0;
for i=1:length(species) % main boxes (species)ç
    count = count + 1;
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        c = c+1;
        count = count+1;
        disp(field_names{f(j)})
        al_goodplot(results.(field_names{f(j)})(:,7), count, [], col(c,:),'left', [], [], 0);
    end
end
xlim([.6 15.5])
xticks([3.5 8.5 13.5])
xticklabels([{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},{'\itAmoeba borokensis'}])
h.XAxis.TickLength = [0 0];
ylabel('DFA\gamma')

%% Export as jpg, tiff and vector graphics pdf

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig3'+wildcardPattern+'.jpg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig3 files found'))

FigList = findobj(allchild(0), 'flat', 'Type', 'figure') ;
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig) ;
  FigName = get(FigHandle, 'Name') ;
  set(0, 'CurrentFigure', FigHandle) ;
  set(FigHandle,'PaperSize',[16.5 22],'PaperPosition',[0 0 16.5 22]);
  set(FigHandle, 'Renderer', 'painters');
  saveas(FigHandle,strcat(destination_folder, '\Figures\Fig3(',num2str(iFig+gabs),')'),'svg')
  exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig3(',num2str(iFig+gabs),').jpg') ...
    ,"Resolution",600)
  % exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig3(',num2str(gabs),').tiff') ...
  %   ,"Resolution",600)
  % exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig3(',num2str(iFig)+gabs,').pdf'), ...
  %   'BackgroundColor','white', 'ContentType','vector')
end