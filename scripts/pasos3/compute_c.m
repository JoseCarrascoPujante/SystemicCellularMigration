function [Ct,cont]=compute_c(y,t)

% u is the vector time series, t a specific time point

dim=length(y);

cont=0;


suma=0; 
resta=0;

for i=1:dim-t
    suma=suma+y(i)*y(i+t);
    cont=cont+1;
    resta=resta+y(i);
end


termino1=suma/cont;
termino2=(resta/cont)^2;



Ct= termino1 -termino2;

end

