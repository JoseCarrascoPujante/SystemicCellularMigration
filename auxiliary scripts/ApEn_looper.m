clear
close all

load coordinates.mat coordinates destination_folder

elapsedApEn = tic;
step = 60;

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

bar1 = waitbar(0,'In progress...','Name','Condition...') ;
bar2 = waitbar(0,'In progress...','Name','Track number...') ;

AE = struct;
AESh = struct;

for i = 1:length(field_names)
    bar1 = waitbar(i/length(field_names), bar1, field_names{i}) ;
    N = length(coordinates.(field_names{i}).original_x(1,:)) ; % N trajectories in condition
    for j = 1:N
		disp(strcat(field_names{i},'nº',num2str(j)))
        bar2 = waitbar(j/N, bar2, strcat('Track number', ' ', num2str(j))) ;
        nframes = length(coordinates.(field_names{i}).original_x(:,1));
        for k = step:step:nframes % calculation step
            AE.(field_names{i})(j,k/step) = ApEn(2, 0.2*std(coordinates.(field_names{i}).scaled_rho(1:k,j)),...
                coordinates.(field_names{i}).scaled_rho(1:k,j)) ;
            % disp(strcat('Original ApEn=',num2str(AE.(field_names{i})(j,k/step))))
               
            AESh.(field_names{i})(j,k/step) = ApEn(2, 0.2*std(coordinates.(field_names{i}).shuffled_rho(1:k,j)),...
                coordinates.(field_names{i}).shuffled_rho(1:k,j)) ;
            % disp(strcat('Shuffled ApEn=',num2str(AESh.(field_names{i})(j,k/step))))
        end
    end
end

elapsedApEn = toc(elapsedApEn);

save(strcat(destination_folder,'\','ApEN_',num2str(nframes/step),'heatmapColumns.mat'), ...
    'AESh', 'AE', 'elapsedApEn') ;

diary off
