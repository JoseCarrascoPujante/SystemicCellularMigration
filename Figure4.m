%% Figure 4
close all
clear
load('coordinates.mat')
load('numerical_results.mat')
load('ApEN_72heatmapColumns.mat')

%% Layouts
set(groot,'defaultFigurePaperPositionMode','manual')

fig = figure('Visible','off','Position', [0 0 850 850]);
layout0 = tiledlayout(2,1,'TileSpacing','tight','Padding','none') ;
layout1 = tiledlayout(layout0,2,3,'TileSpacing','tight','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,2,1,'TileSpacing','compact','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - ApEn heatmaps

% create "full" dataset
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

AE.AmoebaProteus = [];
AESh.AmoebaProteus = [];
AE.MetamoebaLeningradensis = [];
AESh.MetamoebaLeningradensis = [];
AE.AmoebaBorokensis = [];
AESh.AmoebaBorokensis = [];

for index = 1:4
        AE.AmoebaProteus = cat(1,AE.AmoebaProteus, AE.(field_names{index}));
        AESh.AmoebaProteus = cat(1,AESh.AmoebaProteus, AESh.(field_names{index}));
end

for index = 5:8
        AE.MetamoebaLeningradensis = cat(1,AE.MetamoebaLeningradensis, AE.(field_names{index}));
        AESh.MetamoebaLeningradensis = cat(1,AESh.MetamoebaLeningradensis, AESh.(field_names{index}));
end

for index = 9:12
        AE.AmoebaBorokensis = cat(1,AE.AmoebaBorokensis, AE.(field_names{index}));
        AESh.AmoebaBorokensis = cat(1,AESh.AmoebaBorokensis, AESh.(field_names{index}));
end

fields = {"AmoebaProteus","MetamoebaLeningradensis","AmoebaBorokensis"};

% Optionally use AE min and max as colormap range limits for non-shuffled
% heatmaps
Zmin = min([min(min(AE.AmoebaProteus)),min(min(AE.MetamoebaLeningradensis)),min(min(AE.AmoebaBorokensis))]);

Zmax = max([max(max(AE.AmoebaProteus)),max(max(AE.MetamoebaLeningradensis)),max(max(AE.AmoebaBorokensis))]);

% Optionally use AESh min and max as colormap range limits for shuffled
% heatmaps
sZmin = min([min(min(AESh.AmoebaProteus)),min(min(AESh.MetamoebaLeningradensis)),min(min(AESh.AmoebaBorokensis))]);

sZmax = max([max(max(AESh.AmoebaProteus)),max(max(AESh.MetamoebaLeningradensis)),max(max(AESh.AmoebaBorokensis))]);

columns = size(AESh.SinEstimuloProteus11_63,2)-6;

for i = 1:length(fields)
    nexttile(layout1,i)
    h = gca;
    dataExp = AE.(fields{i});
    imagesc(dataExp(:,7:end))
    set(h,'YDir','normal')
    colormap(flipud(jet))
    xticklabels(h,{});
    h.XAxis.TickLength = [0 0];
    if i == 1
        h.YAxis.FontSize = 8;
        ylabel(h,'Series','FontSize',10);
        clim([Zmin Zmax])
    elseif i == 3
        % yticklabels(h,{});
        % h.YAxis.TickLength = [0 0];
        a=colorbar;
        clim([Zmin Zmax])
        ylabel(a,'Approximate Entropy','FontSize',10,'Rotation',270);
    else
        % yticklabels(h,{});
        % h.YAxis.TickLength = [0 0];
        clim([Zmin Zmax])
    end

    nexttile(layout1,i+3)
    h = gca;
    dataShuf = AESh.(fields{i});
    imagesc(dataShuf(:,7:end))
    set(h,'YDir','normal')
    colormap(flipud(jet))
    xticks([1:columns/6:columns,columns])
    xticklabels(h,300:550:3600);
    h.XAxis.FontSize = 8;
    if i == 1
        h.YAxis.FontSize = 8;
        ylabel(h,'Series (shuffled)','FontSize',10);
        clim([sZmin sZmax]);
    elseif i == 2
        % yticklabels(h,{});
        % h.YAxis.TickLength = [0 0];
        xlabel(h,'Time(s)','FontSize',10);
        clim([sZmin sZmax]);
    elseif i == 3
        % yticklabels(h,{});
        % h.YAxis.TickLength = [0 0];
        a=colorbar;
        a.Label.Position(1) = 3.2;
        clim([sZmin sZmax]);
        ylabel(a,'Approximate Entropy','FontSize',10,'Rotation',270);
    end
end

%% Panel 2 - Violin plots

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
col = [.1,.1,.1;.3,.3,.3;.5,.5,.5;.7,.7,.7;1,0,0;1,.25,.25;1,.5,.5; 1,.7,.7;0,0,1;.25,.25,1;.5,.5,1;.7,.7,1];

h = nexttile(layout2,1);

count = 0;
c = 0;
for i=1:length(species) % main boxes (species)
    count = count + 1;
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        c = c+1;
        count = count+2;
        disp(field_names{f(j)})
        conu.al_goodplot(results.(field_names{f(j)})(:,12), count, [], col(c,:),'bilateral', [], [], 0); %Shuffled
    end
end
xlim([1.5 28])
xticklabels([])
h.XAxis.TickLength = [0 0];
h.YAxis.FontSize = 8;
ylabel('Approximate Entropy (Shuffled)','FontSize',10)


h = nexttile(layout2,2);
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
        conu.al_goodplot(results.(field_names{f(j)})(:,11), count, [], col(c,:),'bilateral', [], [], 0); % original
    end
end
xlim([1.5 28])
xticks([5.9 15 24])
h.XAxis.TickLength = [0 0];
set(gca,'XTickLabel',[{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])
h.YAxis.FontSize = 8;
ylabel('Approximate Entropy','FontSize',10)

%% Export as vector graphics svg

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures\')) ;
gabs = 0 ;
for v = 1:length(versions)
    if  contains(versions(v).name, 'Fig4'+wildcardPattern+'.svg')
        gabs = gabs + 1 ;
    end
end

disp(strcat(num2str(gabs),' Fig4 files found'))

fig.Units = 'centimeters';        % set figure units to cm
fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
fig.PaperPosition = [0 0 fig.Position(3:4)];
set(fig, 'Renderer', 'painters');
saveas(fig,strcat(destination_folder, '\Figures\Fig4(',num2str(gabs),').svg'),'svg')

