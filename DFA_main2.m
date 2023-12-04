
function [Alpha1,D,F_n,z1,z2]=DFA_main2(DATA, input_type, axeshandle)
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
%A is the alpha in the paper
%D is the dimension of the time series
%n can be changed to your interest
n=[ 158 251 398 631 1000 1585 2512 ];% 3981 6310 10000 15849 25119 39811 63096 100000]; % % (tamaño óptimo, desde 100 hasta tamaño de la muestra entre 10 aproximadamente)
%n=[100 126 158 200 251 316 398 501 631 794 1000]; % de 0.1 en 0.1 en vez de 0.2 como arriba
%%elegimos n para que haga ventanas desde el 2 hasta el 4 en escala
%%logarítmica, separados de 0.2 en 0.2 (redondeados para ser un tamaño de
%%ventana entero)
N1=length(n);
F_n=zeros(N1,1);
for i=1:N1
 F_n(i)=DFA2(DATA,n(i),1);
end

n=n';

hold(axeshandle, 'on')

if contains(input_type, 'Original')
    scatter(log10(n),log10(F_n), 'filled', 'ro');
elseif contains(input_type, 'Shuffled')
    scatter(log10(n),log10(F_n), 'filled', 'bo');
end

hold on
xlabel('log10(n)','FontSize',10)
ylabel('log10(F(n))','FontSize',10)
A=polyfit(log10(n(1:end)),log10(F_n(1:end)),1);
B=polyval(A,log10(n(1:end)));

[z1,z2]=rsquare(log10(F_n(1:end)),B);

if contains(input_type, 'Original')
    plot(log10(n(1:end)),A(1)*log10(n(1:end))+A(2), 'Color', [1 0 0], 'LineStyle','-')
elseif contains(input_type, 'Shuffled')
    plot(log10(n(1:end)),A(1)*log10(n(1:end))+A(2), 'Color', [0 0 1], 'LineStyle','-')
end
Alpha1=A(1);
D=3-A(1);


return