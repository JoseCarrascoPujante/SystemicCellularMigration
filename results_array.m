function full = results_array(field_names,results)    
    full = [];
    for index = 1:length(field_names)
            full = cat(1,full, results.(field_names{index}));
    end
end