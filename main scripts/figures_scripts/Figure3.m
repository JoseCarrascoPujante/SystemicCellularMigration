% Figure 2
%% Layouts
fig = figure('Position',[350 10 700 1500]);
layout0 = tiledlayout(4,1,'TileSpacing','tight','Padding','none') ;
layout1 = tiledlayout(layout0,2,3,'TileSpacing','tight','Padding','none') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,8,12,'TileSpacing','tight','Padding','none') ;
layout2.Layout.Tile = 2;
layout3 = tiledlayout(layout0,1,3,'TileSpacing','tight','Padding','none') ;
layout3.Layout.Tile = 3;

%% Panel 1 - DFA correlations
fields = {"InduccionProteus11_63","QuimiotaxisLeningradensisVariosPpmm","GalvanotaxisBorokensis11_63"};
amoebas = {39,10,51};
for i=1:3 % subpanels (species)
    nexttile(layout1,i)
    dfahandle = gca;
    gamma = DFA_main2(coordinates.(fields{i}).scaled_rho(:,amoebas{i}),'Original_DFA_', dfahandle) ;
    yl = ylim();
    xl = xlim();
    text(xl(1)+1,yl(1)+.5,strcat('\gamma=',num2str(gamma)))
    nexttile(layout1,i+3)
    dfahandle = gca;
    gamma = DFA_main2(coordinates.(fields{i}).shuffled_rho(:,amoebas{i}),'Shuffled_DFA_', dfahandle) ;
    yl = ylim();
    xl = xlim();
    text(xl(1)+1,yl(1)+.5,strcat('\gamma=',num2str(gamma)))
end

%% Panel 2 - DFA \gamma
field_names = fieldnames(results) ;
species = {'Proteus','Leningradensis','Borokensis'};
t = gca;
tiles = {
[25,73,49,1;37,85,61,13]
[29,77,53,5;41,89,65,17]
[33,81,57,9;45,93,69,21]};
for i = 1:length(species)
    idx = find(contains(field_names(:),species{i}))';
    for f = 1:length(idx)
        t = nexttile(layout2,tiles{i}(1,f),[1,4]);
        exes = zeros(size(results.(field_names{idx(f)}),1));
        plot(results.(field_names{idx(f)})(:,7),exes,'ro','MarkerSize',7)
        box off
        ylim([0 eps]) % minimize y-axis height
        xlim([1 2])
        t.YAxis.Visible = 'off'; % hide y-axis
        t.Color = 'None';
        
        t2 = nexttile(layout2,tiles{i}(2,f),[1,4]);
        datamean = mean(results.(field_names{idx(f)})(:,7));
        datastd = std(results.(field_names{idx(f)})(:,7));
        line([datamean-datastd datamean+datastd],[1 1],'Color','red',...
            'LineWidth',.5)
        text(t2,datamean,-20,[num2str(round(mean(results.(field_names{idx(f)})(:,7)),2))...
            ' ' char(177) ' ' num2str(round(datastd,2))],'HorizontalAlignment',...
            'center','FontSize',9)
        ylim([0 2]) % minimize y-axis height
        xlim([1 2])
        t2.YAxis.Visible = 'off'; % hide y-axis
        t2.XAxis.Visible = 'off'; % hide y-axis
        t2.Color = 'None';
    end
end