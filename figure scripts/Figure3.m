%% Figure 3

close all
clear
load('coordinates.mat')
load('numerical_results.mat')

%% Layouts
set(groot,'defaultFigurePaperPositionMode','manual')
fig = figure('Visible','off','Position', [0 0 1000 900]);

layout0 = tiledlayout(3,1,'TileSpacing','loose','Padding','none') ;
layout1 = tiledlayout(layout0,1,3,'TileSpacing','compact','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,3,'TileSpacing','compact','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - MSD \Beta plots
fields = {"SinEstimuloProteus11_63","SinEstimuloLeningradensis11_63","SinEstimuloBorokensis23_44"};
chosen = [1,2,4,5,8,13,27,29;3,8,19,28,37,44,46,50;2,3,4,5,32,35,36,45];

for i = 1:3

    nexttile(layout1,i)
    h = gca;
    for j = chosen(i,:)
        [~,deltat] = msd(coordinates.(fields{i}).scaled_x(:,j),...
          coordinates.(fields{i}).scaled_y(:,j), h) ;
    end
    plot(h, log(deltat), log(deltat)-10, 'k--')
    plot(h, log(deltat), log(deltat.^2)-11, 'k--')
    text(h, log(deltat(5)),0,['\beta=2, ballistic' newline 'diffusion']...
        ,'HorizontalAlignment', 'center','FontSize',10.5)
    text(h, log(deltat(5)),-4.5,'Super-diffusion'...
        ,'HorizontalAlignment', 'center','FontSize',10.5)
    text(h, log(deltat(5)),-9,['\beta=1, normal' newline 'diffusion']...
        ,'HorizontalAlignment', 'center','FontSize',10.5)
    xlabel('Log(MSD(\tau))');
    ylabel('Log(\tau(s))');
    axis padded
    xlim([-.8    6.11])
    ylim([-14    2.2185])

    % nexttile(layout1,i+3)
    % h = gca;
    % for j=1:8
    %     [~,deltat] = msd(coordinates.(fields{i}).shuffled_x(:,j),...
    %         coordinates.(fields{i}).shuffled_y(:,j), h) ;
    % end
    % plot(h, log(deltat), log(deltat)-1.5, 'k--')
    % plot(h, log(deltat), log(deltat.^2)-1.5, 'k--')
    % text(h, log(deltat(5)),8,['\beta=2, ballistic' newline 'diffusion']...
    %     ,'HorizontalAlignment', 'center','FontSize',8)
    % text(h, log(deltat(5)),5,'Superdiffusion'...
    %     ,'HorizontalAlignment', 'center','FontSize',8)
    % text(h, log(deltat(5)),2,['\beta=1, normal' newline 'diffusion']...
    %     ,'HorizontalAlignment', 'center','FontSize',8)
    % xlabel('Log(MSD(\tau))');
    % ylabel('Log(\tau(s))');
    % xlim([-1 6.2])
    % ylim([-2.3863    11.2185])
end

%% Panel 2 - MSD \Beta circles

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

species = {'Proteus','Leningradensis','Borokensis'};
tiles = {
[1,7,13,19;4,10,16,22]
[2,8,14,20;5,11,17,23]
[3,9,15,21;6,12,18,24]
};

for i = 1:length(species)
    idx = find(contains(field_names(:),species{i}))';
    for f = 1:length(idx)
        disp(field_names{idx(f)})
        t = nexttile(layout2,tiles{i}(1,f));
        hold on
        
        exes = zeros(size(results.(field_names{idx(f)}),1));
        plot(results.(field_names{idx(f)})(:,9),exes,'ro','MarkerSize',7)
        plot(results.(field_names{idx(f)})(:,10),exes,'bo','MarkerSize',7)
        ylim([0 eps]) % minimize y-axis height
        xlim([-0.1 2.1])
        t.XRuler.TickLabelGapOffset = 2;
        t.YAxis.Visible = 'off'; % hide y-axis
        t.Color = 'None';
        hold off

        t2 = nexttile(layout2,tiles{i}(2,f));
        hold on
        datamean = mean(results.(field_names{idx(f)})(:,9));
        datastd = std(results.(field_names{idx(f)})(:,9));
        datameanshuff = mean(results.(field_names{idx(f)})(:,10));
        datastdshuff = std(results.(field_names{idx(f)})(:,10));
        line([datamean-datastd datamean+datastd],[0 0],'Color','red',...
            'LineWidth',.5)
        line([datamean-datastd+.005 datamean-datastd],[0 0],'Color','red',...
            'LineWidth',5)
        line([datamean+datastd datamean+datastd+.005],[0 0],'Color','red',...
            'LineWidth',5)
        text(t2,datamean,-.075,[num2str(round(datamean,2)) ' ' char(177) ' '...
            num2str(round(datastd,2))],'HorizontalAlignment','center','FontSize',8)
        line([datameanshuff-datastdshuff datameanshuff+datastdshuff],[0 0],'Color','blue',...
            'LineWidth',.5)
        line([datameanshuff-datastdshuff+.005 datameanshuff-datastdshuff],[0 0],'Color','blue',...
            'LineWidth',5)
        line([datameanshuff+datastdshuff datameanshuff+datastdshuff+.005],[0 0],'Color','blue',...
            'LineWidth',5)
        text(t2,datameanshuff,-.075,[num2str(round(datameanshuff,2)) ' ' char(177)...
            ' ' num2str(datastdshuff,'%.e')],'HorizontalAlignment','center','FontSize',8)
        ylim([-0.08 0]) % minimize y-axis height
        xlim([-0.1 2.1])
        t2.YAxis.Visible = 'off'; % hide y-axis
        t2.XAxis.Visible = 'off'; % hide y-axis
        t2.Color = 'None';
        hold off
    end
end


%% "Superviolin" plots of MSD \beta
ax=nexttile(layout0,3);

rmsf_conds = {{[],[],[],[]},{[],[],[],[]},{[],[],[],[]}};

for i=1:length(species) % main boxes (species)
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        rmsf_conds{i}{j} = results.(field_names{f(j)})(:,9);
    end
end

for i=1:length(species) % main boxes (species)
    superviolin(rmsf_conds{i},'Parent',ax,'Xposition',i,'FaceAlpha',0.15,...
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
    for f = find(contains(field_names(:),species(i)))' % conditions
        rmsfs{i} = [rmsfs{i}; results.(field_names{f})(:,9)];
    end
end

boxChart_rmsf=cat(1,rmsfs{1},rmsfs{2},rmsfs{3});
boxchart([ones(length(rmsfs{1}),1); repmat(2,length(rmsfs{2}),1); ...
    repmat(3,length(rmsfs{3}),1)],boxChart_rmsf,'Notch','off',...
    'GroupByColor',colorgroups,'BoxFaceAlpha',0)
h=gca;
xlim([.5 3.5])
ylim([1 2.1])
% h.XAxisLocation = 'top';
box on
h.XTick = [1 2 3];
xticklabels(h,[{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])
h.XAxis.FontSize = 8;
h.XAxis.TickLength = [0 0];
ylabel('MSD\beta','FontSize',10)

%% Export as jpg, tiff and vector graphics pdf

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig3'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig3 files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig3(',num2str(gabs),')'),'svg')
% exportgraphics(FigHandle,strcat(destination_folder, '\Figures\Fig3(',num2str(gabs),').jpg') ...
    % ,"Resolution",600)
% exportgraphics(FigHandle,strcat(destination_folder, '\Figures\Fig3(',num2str(gabs),').tiff') ...
%     ,"Resolution",600)
% exportgraphics(FigHandle,strcat(destination_folder, '\Figures\Fig3(',num2str(gabs),').pdf'), ...
%   'Resolution',600,'ContentType','vector')

