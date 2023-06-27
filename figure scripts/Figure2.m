% Figure 2
%% Layouts
set(groot,'defaultFigurePaperPositionMode','manual')
fig = figure('Visible','off','Position',[0 0 900 1200]);
layout0 = tiledlayout(3,1,'TileSpacing','tight','Padding','none') ;
layout1 = tiledlayout(layout0,1,3,'TileSpacing','tight','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,8,'TileSpacing','none','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - RMSF max_correlations
fields = {"InduccionProteus11_63","InduccionLeningradensis11_63","QuimiotaxisBorokensis23_44"};
amoebas = {8,55,44};
for i=1:3 % subpanels (species)
    nexttile(layout1)   
    rmsfhandle = gca;
    set(rmsfhandle,'xscale','log')
    set(rmsfhandle,'yscale','log')
    amebas5(coordinates.(fields{i}).scaled_rho(:,amoebas{i}),rmsfhandle) ;
end

%% Panel 2 - RMSF \alpha
field_names = fieldnames(results) ;
t = gca;
tiles = [17,49,33,1;25,57,41,9];
idx = find(contains(field_names(:),'Leningradensis'))';
for f = 1:length(idx)
    t = nexttile(layout2,tiles(1,f),[1,5]);
    exes = zeros(size(results.(field_names{idx(f)}),1));
    plot(results.(field_names{idx(f)})(:,1),exes,'ro','MarkerSize',7)
    box off
    ylim([0 eps]) % minimize y-axis height
    xlim([0.5 0.9])
    t.YAxis.Visible = 'off'; % hide y-axis
    t.Color = 'None';
    
    t2 = nexttile(layout2,tiles(2,f),[1,5]);
    datamean = mean(results.(field_names{idx(f)})(:,1));
    datastd = std(results.(field_names{idx(f)})(:,1));
    line([datamean-datastd datamean+datastd],[0.05 0.05],'Color','red',...
        'LineWidth',.5)
    text(t2,datamean,-.35,[num2str(round(mean(results.(field_names{idx(f)})(:,1)),2))...
        ' ' char(177) ' ' num2str(round(datastd,2))],'HorizontalAlignment',...
        'center','FontSize',9)
    ylim([0 1]) % minimize y-axis height
    xlim([0.5 0.9])
    t2.YAxis.Visible = 'off'; % hide y-axis
    t2.XAxis.Visible = 'off'; % hide y-axis
    t2.Color = 'None';
end

%% Panel 3 - RMSF Violin plots
ax=nexttile(layout0,3);
hold on
species = {'Proteus','Leningradensis','Borokensis'};
rmsfs = {[],[],[]};
for i=1:length(species) % species    
    for f = find(contains(field_names(:),species(i)))' % conditions
        rmsfs{i} = [rmsfs{i}; results.(field_names{f})(:,5)/120];
    end
end

%%%Simple violin plots
% violin(rmsfs,'xlabel',...
% {'\itAmoeba proteus','\itMetamoeba leningradensis','\itAmoeba borokensis'}, ...
%     'ylabel',{'RMSF\alpha'},'FontSize',4,'facecolor',...
%     [[0,0,0];[1,0,0];[0,0,1]],'facealpha',0.15,'mc','y','medc','r')
% rmsfs_pad = padcat(rmsfs{1},rmsfs{2},rmsfs{3});
% boxplot(rmsfs_pad)

%%%More elaborate violin plots
% rmsfs_pad = padcat(rmsfs{1},rmsfs{2},rmsfs{3});
% violinplot(rmsfs_pad,...
%     {'\itAmoeba proteus','\itMetamoeba leningradensis','\itAmoeba borokensis'}, ...
%     'ViolinColor',[1,0,0;0,0,0;0,0,1],'ViolinAlpha',0.15,'ShowData',...
%     false,'ShowNotches',false,'ShowMean',true,'ShowMedian',true,'MedianColor',...
%     [1 1 0],'HalfViolin','left','BoxColor',[0 0 0],'BoxWidth',0.02)
% boxplot(rmsfs_pad)

%%%"Superviolin" plots
% rmsf_conds = {{[],[],[],[]},{[],[],[],[]},{[],[],[],[]}};
% for i=1:length(species) % main boxes (species)
%     f = find(contains(field_names(:),species(i)))'; % condition indexes
%     for j = 1:length(f) % secondary boxes (conditions)
%         rmsf_conds{i}{j} = results.(field_names{f(j)})(:,5)/120;
%     end
% end
% for i=1:length(species) % main boxes (species)
%     superviolin(rmsf_conds{i},'Parent',ax,'Xposition',i,'FaceAlpha',0.15,...
%         'Errorbars','ci','Centrals','mean','LineWidth',0.1)
% end
% colorgroups = [repmat({'Galvanotaxis'},length(rmsf_conds{1}{1}),1);
%     repmat({'Inducción'},length(rmsf_conds{1}{2}),1);
%     repmat({'Quimiotaxis'},length(rmsf_conds{1}{3}),1);
%     repmat({'Sin estímulo'},length(rmsf_conds{1}{4}),1);
%     repmat({'Galvanotaxis'},length(rmsf_conds{2}{1}),1);
%     repmat({'Inducción'},length(rmsf_conds{2}{2}),1);
%     repmat({'Quimiotaxis'},length(rmsf_conds{2}{3}),1);
%     repmat({'Sin estímulo'},length(rmsf_conds{2}{4}),1);
%     repmat({'Galvanotaxis'},length(rmsf_conds{3}{1}),1);
%     repmat({'Inducción'},length(rmsf_conds{3}{2}),1);
%     repmat({'Quimiotaxis'},length(rmsf_conds{3}{3}),1);
%     repmat({'Sin estímulo'},length(rmsf_conds{3}{4}),1)];
% boxChart_rmsf=cat(1,rmsfs{1},rmsfs{2},rmsfs{3});
% boxchart([ones(length(rmsfs{1}),1); repmat(2,length(rmsfs{2}),1); ...
%     repmat(3,length(rmsfs{3}),1)],boxChart_rmsf,'Notch','off',...
%     'GroupByColor',colorgroups,'BoxFaceAlpha',0) %Box charts whose notches do not overlap have different medians at the 5% significance level.
% h=gca;
% xlim([.5 3.5])
% h.XTick = [1 2 3];
% xticklabels([{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
%     {'\itAmoeba borokensis'}])
% h.XAxis.TickLength = [0 0];

%%%RainCloud plots
h=gca;
xlim([.5 3.5])
ylim([-8.35 30])
h.XTick = [1 2 3];
xticklabels([{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])
h.XAxis.TickLength = [0 0];

data = cell(3,4);
for i=1:length(species) % species    
    f = find(contains(field_names(:),species(i))); % conditions
    for j = 1:length(f)
        data{i,j} = [data{i}; results.(field_names{f(j)})(:,5)/120];
        disp(field_names{f(j)})
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
        plot_rainclouds(data{p,q},cb,count);
    end
end

%% Export as jpg, tiff and vector graphics pdf

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

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
saveas(fig,strcat(destination_folder, '\Figures\Fig2(',num2str(gabs),')'),'svg')
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig2(',num2str(gabs),').jpg') ...
%   ,"Resolution",600)
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig2(',num2str(gabs),').tiff') ...
%   ,"Resolution",600)
% exportgraphics(gcf,strcat(destination_folder, '\Figures\Fig2(',num2str(gabs),').pdf'), ...
%   'BackgroundColor','white', 'ContentType','vector')

function plot_rainclouds(data,cb,count)
    figure    
    raincloud_plot(data,'box_on',1,'box_dodge',1,...
        'box_dodge_amount',-.08,'dot_dodge_amount',.15,'alpha',0.2,...
        'bxcl',[0 0 0],'color',cb(count,:),'LineWidth',1);
    view([90 -90])
    xlim([-8.35 30])
    ylim([-.03 .095])
    box off
end