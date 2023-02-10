
function [dev1,dev2,dev3,tc2]=amebas5(u_all, fighandle)

%clear all;     PARA MIRAR UNA SERIE EN VARIOS TMAX
%close all;
% clc;

%data=load('cortasrho.mat');

% microns per pixel at scale(i), cell i

scale_time=0.5; %tiempo en segundos entre cada paso

%%%%%%%%%%%%%%%%%
% calculating F %
%%%%%%%%%%%%%%%%%
tc=30;
numseries=1;

exponents=zeros(tc,numseries);
goodness=zeros(tc,numseries);
ks=zeros(tc,numseries);

% close all;

 
clear title xlabel ylabel;
tc2=tc*100;
zz=1;
%time_max=500;
for time_max=100:100:tc2    %ojo, este bucle no es compatible con el de varias series, no te vayas a liar!
    
numseries=1;
umin=1;
for jjj=1:numseries       
i=1;  
    
%for i=1:21
    
% %    filex=strcat('data.x',num2str(i));
% %    filey=strcat('data.c',num2str(i));
%     filey=strcat('data.c',num2str(i));
% %    x=eval(filex);
%     u_all=eval(filey);
% u_all=rho2;    
 
    nu=3500; %20000
    du=1;
    u=u_all(umin:du:umin+nu);
   
    dim=size(u,1);
    l=zeros(dim,1);
    
  
    for j=2:dim
        
        l(j)=sqrt(  ( u(j)-u(j-1) ) * ( u(j)-u(j-1)) );
        
    end
    
  
    
    % calculates y, net displacement for each t
    
    y=zeros(dim,1);
    
   % l=shuff(l);
    
    for j=1:dim
        
        y(j)=compute_y(l,j);
        
    end
    
    % root mean square fluctuation of the displacement
    
    F=zeros(dim,1);
    
%     tic
    for j=1:dim
        
        F(j)=compute_F(y,j);
        
    end
%     toc
    
    
    %%%%%%%%%%%%%
    % F fitting %
    %%%%%%%%%%%%%
%     
time=(1:dim)'.*scale_time*du;
idx=find( time(1:time_max)>0  & F(1:time_max)>0 );
[Rsq,slope,intercept]=rsquare (log10(time(idx)),log10(F(idx)));

alpha=slope;
v=time.^(alpha);
k=power(10,intercept/alpha);


    %%%%%%%%%%%%
    % Plotting %
    %%%%%%%%%%%%


% if time_max is in its last iteration...

if time_max == tc2
% ...Plot the NOT SHUFFLED (original) loglog plot of rmsf F versus l step...
        
    ylabelvar='rms fluctuation (mA)';
    xlabelvar='time (sec)';
    fontsizevar=20;
    
    hold(fighandle, 'on')
    loglog(fighandle, time(1:time_max), F(1:time_max), 'o', 'MarkerSize', 7);
    set(gca,'FontSize',fontsizevar);
    xlabel(xlabelvar);
    ylabel(ylabelvar);
    
% Plot regression fit for not-shuffled data
%     loglog(time(1:time_max),k*v(1:time_max),'k','LineWidth',3);
%     hold on       
    
    
% ...and plot the SHUFFLED loglog plot of rmsfF versus l step
       % l=shuff(l);
         
        % for j=1:dim
            
            % y(j)=compute_y(l,j);
            
        % end
        
        % F=zeros(dim,1);
        
        % for j=1:dim
            
            % F(j)=compute_F(y,j);
            
        % end

        % time=(1:dim)'.*scale_time*du;
        
        
        % loglog(time(1:time_max),F(1:time_max),'bo','MarkerSize',7);
        % set(gca,'FontSize',fontsizevar);
        % legend('Original RMSF','Shuffled RMSF')
 		% Plot regression fit for shuffled data
        % idx2=find( time(1:time_max)>0  & F(1:time_max)>0 );
        % [Rsq2,slope2,intercept2]=rsquare (log10(time(idx2)),log10(F(idx2)));
        % alpha=slope2;
        % v=time.^(alpha);
        % hold on;
        % loglog(time(1:time_max),k*v(1:time_max),'k--','LineWidth',3);
     
     
end
    exponents(zz,jjj)=slope;
    goodness(zz,jjj)=Rsq;
    ks(zz,jjj)=power(10,intercept);
    junto(zz,:)=[slope Rsq time_max];  
%     
 %%%% DEVOLVER SOLO UNA FILA EN FUNCIÓN DEL RCUADRADO (bondad del ajuste)
 B=goodness>0.975;
 res=junto(B,:,:);
 dev=res(end,:);
 dev1=dev(:,1);
 dev2=dev(:,2);
 dev3=dev(:,3);
    
%end 
%     
   umin=umin+nu;     
end
    medias(zz)=mean(exponents(zz,:));
    standard(zz)=std(exponents(zz,:));
    zz=zz+1;
end
end
