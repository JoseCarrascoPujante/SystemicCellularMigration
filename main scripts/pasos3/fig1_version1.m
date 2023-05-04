function figure1 = fig1_version1(coordinates)
N = 40 ;
field_names = fieldnames(coordinates) ;
figure
figure1 = tiledlayout(3,1,'TileSpacing','tight','Padding','tight') ; % outer layout

%# Panels 1-4
fig1_1=tiledlayout(figure1,2,2,'TileSpacing','tight','Padding','tight'); % inner layout
fig1_1.Layout.Tile = 1;
fig1_1.Layout.TileSpan = [2 1];
panels = {'SinEstimulo','Galvanotaxis','Quimiotaxis','Induccion'};
for i = 1:length(panels)
    idx = find(contains(field_names,panels{i}));
    nexttile(fig1_1)
    for species=1:length(idx)
        colr = [0, 0, 0];
        colr(species) = colr(species)+1;
        for track = 1:N
            %# Plot trajectory and 'ko' marker at its tip      
            plot(coordinates.(field_names{idx(species)}).scaled_x(:,track),...
                coordinates.(field_names{idx(species)}).scaled_y(:,track), 'Color', colr) ;
            hold on;
            plot(coordinates.(field_names{idx(species)}).scaled_x(end,track),...
                coordinates.(field_names{idx(species)}).scaled_y(end,track),...
                'ko',  'MarkerFaceColor',  colr, 'MarkerSize', 2) ;
        end
    end
end

%# Panels 5-6
fig1_3=tiledlayout(figure1,1,3,'TileSpacing','tight','Padding','tight'); % Inner layout
fig1_3.Layout.Tile = 3;
nexttile(fig1_3,[1 3])

end