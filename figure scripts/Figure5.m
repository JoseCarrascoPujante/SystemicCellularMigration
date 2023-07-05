%% Figure5
close all
clear
load('2023-06-07_14.16''19''''_coordinates.mat')
load('2023-06-07_14.16''19''''_numerical_results.mat')
load('ApEn.mat')

%% Layouts
set(groot,'defaultFigurePaperPositionMode','manual')

fig = figure('Visible','off','Position', [0 0 850 800]);
layout0 = tiledlayout(2,1,'TileSpacing','tight','Padding','none') ;
layout1 = tiledlayout(layout0,2,3,'TileSpacing','tight','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,2,1,'TileSpacing','compact','Padding','none') ;
layout2.Layout.Tile = 2;

%% Panel 1 - ApEn heatmaps

fields = {"SinEstimuloProteus11_63","SinEstimuloLeningradensis11_63","SinEstimuloBorokensis23_44"};

% Calculate Approximate Entropy
% AE = struct;
% AESh = struct;
% for i = 1:length(fields)
% 
%     for j = 1:length(coordinates.(fields{i}).original_x(1,:))
%         for k = 72:72:3600
%             AE.(fields{i})(j,k/72) = ApEn(2, 0.2*std(coordinates.(fields{i}).scaled_rho(1:k,j)),...
%                 coordinates.(fields{i}).scaled_rho(1:k,j)) ;
%             AESh.(fields{i})(j,k/72) = ApEn(2, 0.2*std(coordinates.(fields{i}).shuffled_rho(1:k,j)),...
%                 coordinates.(fields{i}).shuffled_rho(1:k,j)) ;
%         end
%     end
% end

% Zmin = min([min(min(AE.SinEstimuloBorokensis23_44)),min(min(AE.SinEstimuloLeningradensis11_63)),min(min(AE.SinEstimuloProteus11_63)), ...
%     min(min(AESh.SinEstimuloBorokensis23_44)),min(min(AESh.SinEstimuloLeningradensis11_63)),min(min(AESh.SinEstimuloProteus11_63))]);
% 
% Zmax = max([max(max(AE.SinEstimuloBorokensis23_44)),max(max(AE.SinEstimuloLeningradensis11_63)),max(max(AE.SinEstimuloProteus11_63)), ...
%     max(max(AESh.SinEstimuloBorokensis23_44)),max(max(AESh.SinEstimuloLeningradensis11_63)),max(max(AESh.SinEstimuloProteus11_63))]);

for i = 1:length(fields)
    nexttile(layout1,i)
    h = gca;
    imagesc(AE.(fields{i}))
    colormap(jet)
    a=colorbar;
    ylabel(a,'Approximate Entropy','FontSize',7.5,'Rotation',270);
    xticklabels(h,{});
    if i == 1
        ylabel(h,'Series');
    else
        yticklabels(h,{});
    end


    nexttile(layout1,i+3)
    h = gca;
    imagesc(AESh.(fields{i}))
    colormap(jet)
    a=colorbar;
    ylabel(a,'Approximate Entropy','FontSize',7.5,'Rotation',270);
    % a.Label.Position(1) = 3.2;
    % clim([Zmin Zmax]);
    xticks(10:10:50)
    xticklabels(h,compose('%d',720:720:3600));
    if i == 1
        ylabel(h,'Series (shuffled)','FontSize',10);
    elseif i == 2
        yticklabels(h,{});
        xlabel(h,'time(s)');
    elseif i == 3
        yticklabels(h,{});
    end
end

%% Panel 2 - Violin plots

h = nexttile(layout2,1);

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

count = 0;
c = 0;
for i=1:length(species) % main boxes (species)
    count = count + 1;
    f = find(contains(field_names(:),species(i)))'; % condition indexes
    for j = 1:length(f) % secondary boxes (conditions)
        c = c+1;
        count = count+2;
        disp(field_names{f(j)})
        al_goodplot(results.(field_names{f(j)})(:,12), count, [], col(c,:),'bilateral', [], [], 0); %Shuffled
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
        al_goodplot(results.(field_names{f(j)})(:,11), count, [], col(c,:),'bilateral', [], [], 0); % original
    end
end
xlim([1.5 28])
xticks([5.9 15 24])
h.XAxis.TickLength = [0 0];
set(gca,'XTickLabel',[{'\itAmoeba proteus'},{'\itMetamoeba leningradensis'},...
    {'\itAmoeba borokensis'}])
h.YAxis.FontSize = 8;
ylabel('Approximate Entropy','FontSize',10)

%% Export as jpg and vector graphics pdf

if ~exist(strcat(destination_folder,'\Figures'), 'dir')
   mkdir(strcat(destination_folder,'\Figures'))
end

versions = dir(strcat(destination_folder,'\Figures')) ;
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

