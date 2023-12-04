function [Ft,cont]=compute_F(y,t)

% y is the vector time series, t a specific time point


dim=size(y,1);

y1=0;
y2=0;

cont=0;


 
t0=1;


while t+t0<dim
    
    cont=cont+1;
    
    dy=y(t+t0)-y(t0);  
    
    y1=y1+dy;
    y2=y2+dy*dy;
    
    
    t0=t0+1;

    
end

y1=y1/cont;

y2=y2/cont;

% y1Squared=y1^2;
% precise_substraction = vpa(y2-y1Squared);
% Ft=sqrt(precise_substraction);

Ft=sqrt(y2-y1*y1);



end

