function [yt]=compute_y (u,t)

% u is the vector time series, t a specific time point

yt=0;

for i=1:t
    
    yt=yt+ u(i);
    
end

end

