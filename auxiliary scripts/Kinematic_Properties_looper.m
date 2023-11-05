clear
close all
load numerical_results.mat results field_names
load coordinates.mat coordinates destination_folder stat_names

diary off
diary_filename = strcat(destination_folder,'\cytokineticValues.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on
elapsedCytoki = tic;


bar1 = waitbar(0,'In progress...','Name','Condition...') ;
bar2 = waitbar(0,'In progress...','Name','Track number...') ;

for i = 1:length(field_names)
    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
    for j = 1:N
		disp(strcat(field_names{i},'nº',num2str(j)))
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        
        % Intensity of response (mm)
        [results.(field_names{i})(j,strcmp(stat_names(:), 'Intensity'))] = ...
            norm([coordinates.(field_names{i}).scaled_x(end,j) coordinates.(field_names{i}).scaled_y(end,j)]...
            - [coordinates.(field_names{i}).scaled_x(1,j) coordinates.(field_names{i}).scaled_y(1,j)]);

        % Shuffled intensity of response (mm)
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sIntensity'))] = ...
            norm([coordinates.(field_names{i}).shuffled_x(1,j) coordinates.(field_names{i}).shuffled_y(1,j)]...
            - [coordinates.(field_names{i}).shuffled_x(end,j) coordinates.(field_names{i}).shuffled_y(end,j)]);

        % Directionality ratio (straightness)
        d = hypot(diff(coordinates.(field_names{i}).scaled_x(:,j)), diff(coordinates.(field_names{i}).scaled_y(:,j))) ;
        distTrav = sum(d);
        disp(strcat('distTrav:',' ',num2str(distTrav)))
        [results.(field_names{i})(j,strcmp(stat_names(:), 'DR'))] = ...
            results.(field_names{i})(j,strcmp(stat_names(:), 'Intensity'))/distTrav;
        
        % Shuffled Directionality ratio (straightness)
        d = hypot(diff(coordinates.(field_names{i}).shuffled_x(:,j)), diff(coordinates.(field_names{i}).shuffled_y(:,j))) ;
        SdistTrav = sum(d);
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sDR'))] = ...
            results.(field_names{i})(j,strcmp(stat_names(:), 'sIntensity'))/SdistTrav;

        % Average speed (mm/s)
        [results.(field_names{i})(j,strcmp(stat_names(:), 'AvgSpeed'))] = ...
            distTrav/1800;

        % Shuffled average speed (mm/s)
        [results.(field_names{i})(j,strcmp(stat_names(:), 'sAvgSpeed'))] = ...
            SdistTrav/1800;
    end
end

toc(elapsedCytoki)

diary off
