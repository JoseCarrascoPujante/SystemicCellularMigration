diary off
diary_filename = strcat(destination_folder,'\ApEnValues.txt') ;
set(0,'DiaryFile',diary_filename)
clear diary_filename
diary on
elapsedApEn = tic;
field_names = fieldnames(results) ;

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
        for k = 72:72:3600 % calculation step
            AE.(field_names{i})(j,k/72) = ApEn(2, 0.2*std(coordinates.(field_names{i}).scaled_rho(1:k,j)),...
                coordinates.(field_names{i}).scaled_rho(1:k,j)) ;
            % disp(strcat('Original ApEn=',num2str(AE.(field_names{i})(j,k/72))))
               
            AESh.(field_names{i})(j,k/72) = ApEn(2, 0.2*std(coordinates.(field_names{i}).shuffled_rho(1:k,j)),...
                coordinates.(field_names{i}).shuffled_rho(1:k,j)) ;
            % disp(strcat('Shuffled ApEn=',num2str(AESh.(field_names{i})(j,k/72))))
        end
    end
end

elapsedApEn = toc(elapsedApEn);

save(strcat(destination_folder,'\',run_date,'ApEN.mat'), ...
    'AESh', 'AE', 'elapsedApEn') ;

diary off
