
function figures = GraphicalAbstract_IndividualLayout(field_names, results, figures)

% Plots systemic movement parameters in 2D and 3D scatterplots

% Group the 12 track-wise statistics in one table

results.full = [];

for index = 1:length(field_names)

        results.full = cat(1,results.full, results.(field_names{index}));
end

stat_names = ["rmsfAlpha" "srmsfAlpha" "rmsfR2" "srmsfR2" "rmsfTimeMax" ...
    "srmsfTimeMax" "dfaGamma" "sdfaGamma" "msdBeta" "smsdBeta" "ApEn" "sApEn"];
indexes = 1:length(stat_names);
pairs = nchoosek(indexes(1:2:end),2);
disp(pairs)

% 2D scatters

precision = [[.50, .75, .90, .95]; [.5, .75, 1, 1.5]];
ellipseFitType = [{'%Confidence'},{'xSTD'}];
colr = ['w','c','r','m'];
for precisionType=1:height(precision)
    for ej=1:length(pairs)
        disp(strcat('Combination: ',stat_names(pairs(ej,1)), '_and_',stat_names(pairs(ej,2))))
        for c = 1:length(precision(precisionType,:))
            disp(strcat(ellipseFitType{precisionType},': ',num2str(precision(precisionType,c))))
            f = figure;
            f.InvertHardcopy = 'off';
            if pairs(ej,1) == 7
                metric1 = cat(1,results.full(:,pairs(ej,1))-1, results.full(:,pairs(ej,1)+1));
            else 
                metric1 = cat(1,results.full(:,pairs(ej,1)), results.full(:,pairs(ej,1)+1));
            end
            if pairs(ej,2) == 7
                metric2 = cat(1,results.full(:,pairs(ej,2))-1, results.full(:,pairs(ej,2)+1));
            else
                metric2 = cat(1,results.full(:,pairs(ej,2)), results.full(:,pairs(ej,2)+1));
            end
            G = [1*ones(length(metric1)/2,1) ; 2*ones(length(metric2)/2,1)];
            gscatter(metric1,metric2, G,'gy','..',2.25,'off') ;
            set(gca, 'Color','k', 'XColor','w', 'YColor','w')
            set(gcf, 'Color','k')
            hold on
            ellipse(cat(2,metric1,metric2),G,precision(precisionType,c),colr(c),precisionType)
            xlabel(stat_names(pairs(ej,1)))
            ylabel(stat_names(pairs(ej,2)))
            h(1) = plot(nan, nan, 'go', 'MarkerFaceColor','g', 'MarkerSize', 10, 'DisplayName', 'Original');
            h(2) = plot(nan, nan, 'yo', 'MarkerFaceColor','y', 'MarkerSize', 10, 'DisplayName', 'Shuffled');
            h(3) = plot(nan, nan, strcat(colr(c),'--'), 'MarkerSize', 25, ...
                'DisplayName', strcat(num2str(precision(precisionType,c)), ellipseFitType{precisionType}));
            leg=legend(h,Orientation='Horizontal',TextColor='w',FontSize=15);
            leg.Location = 'northoutside';
            saveas(gcf,strcat(stat_names(pairs(ej,1)),'_vs_',stat_names(pairs(ej,2)),...
                '_',num2str(precision(precisionType,c)),ellipseFitType{precisionType},'.png'))
            hold off
        end
    end
end



end