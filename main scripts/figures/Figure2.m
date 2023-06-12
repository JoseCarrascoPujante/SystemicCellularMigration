% Figure 2
%% Layouts
fig = figure('Position',[350 10 800 1000]);
layout0 = tiledlayout(4,1,'TileSpacing','tight','Padding','tight') ;
layout1 = tiledlayout(layout0,1,3,'TileSpacing','tight','Padding','tight') ;
layout1.Layout.Tile = 1;
layout2 = tiledlayout(layout0,3,1,'TileSpacing','tight','Padding','tight') ;
layout2.Layout.Tile = 2;
layout3 = tiledlayout(layout0,1,3,'TileSpacing','tight','Padding','tight') ;
layout3.Layout.Tile = 3;
layout4 = tiledlayout(layout0,1,3,'TileSpacing','tight','Padding','tight') ;
layout4.Layout.Tile = 4;

%% Panel 1 - RMSF max_correlation
fields = {"InduccionProteus11_63","InduccionLeningradensis11_63","QuimiotaxisBorokensis23_44"};
amoebas = {8,55,44};
for i=1:3 % subpanels (species)
    nexttile(layout1)   
    rmsfhandle = gca;
    set(rmsfhandle,'xscale','log')
    set(rmsfhandle,'yscale','log')
    amebas5(coordinates.(fields{i}).scaled_rho(:,amoebas{i}), rmsfhandle) ;
end

%% Panel 2 - RMSF \alpha
field_names = fieldnames(results) ;
species = {'Proteus','Leningradensis','Borokensis'};

for i=1:length(species) % subpanels (species)
    nexttile(layout2);
    t = gca;
    hold on
    for f = find(contains(field_names(:),species(i)))' % condition indexes
        exes = zeros(size(results.(field_names{f}),1));
        plot(results.(field_names{f})(:,1),exes,'ro','MarkerSize',5)
    end
    ylim([0 eps])
    t.YAxis.Visible = 'off'; % remove y-axis
    t.Color = 'None';
    hold off
end
