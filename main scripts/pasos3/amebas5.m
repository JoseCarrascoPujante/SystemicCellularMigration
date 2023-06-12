
function [dev1,dev2,dev3,tc2]=amebas5(u_all, fighandle)

%clear all;     
%close all;
% clc;
% microns per pixel at scale(i), cell i

scale_time=0.5; %time between frames (seconds)

%%%%%%%%%%%%%%%%%
% calculating F %
%%%%%%%%%%%%%%%%%
tc=30;
numseries=1;

exponents=zeros(tc,numseries);
goodness=zeros(tc,numseries);
ks=zeros(tc,numseries);
 
clear title xlabel ylabel;
tc2=tc*100;
zz=1;

for time_max=100:100:tc2    % loop for one series
    
numseries=1;
umin=1;
for jjj=1:numseries       
% i=1;  
    
%for i=1:21
    
% %    filex=strcat('data.x',num2str(i));
% %    filey=strcat('data.c',num2str(i));
%     filey=strcat('data.c',num2str(i));
% %    x=eval(filex);
%     u_all=eval(filey);
% u_all=rho2;    
 
    nu=3500;
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
    
    for j=1:dim
        
        F(j)=compute_F(y,j);
        
    end
   
    
    %%%%%%%%%%%%%
    % F fitting %
    %%%%%%%%%%%%%
%     
time=(1:dim)'.*scale_time*du;
idx=find( time(1:time_max)>0  & F(1:time_max)>0 );
[Rsq,slope,intercept]=rsquare (log10(time(idx)),log10(F(idx)));

alpha=slope;
v=time.^(alpha);
k=F(1);

    %%%%%%%%%%%%
    % Plotting %
    %%%%%%%%%%%%

% if time_max for loop is in its last iteration
if time_max == tc2

    % loglog plot of rmsf F versus l step non-shuffled coordinates     
    hold(fighandle, 'on')
    loglog(fighandle, time(1:time_max), F(1:time_max), 'or', 'MarkerSize', 2);
    % set(gca,'FontSize',20);
    xlabel('Log({\itl}(s))','FontName','times new roman');
    ylabel('Log(F({\itl}))','FontName','times new roman');
    
    % Plot regression fit line for the original data
    loglog(time(1:time_max),k*v(1:time_max),'k--','LineWidth',1.5);
    text(time(1),k*v(1),strcat('\alpha=',num2str(alpha)),...
        'FontName','times new roman')

    % Plot alpha=0.5 line for the original data
    v2 = time.^(.5);
    k2 = .05;
    loglog(time(1:time_max),k2*v2(1:time_max),'k--','LineWidth',1.5);
    text(time(1),k2*v2(1),'\alpha=0.5, uncorrelated','FontName','times new roman')
    
end
    exponents(zz,jjj)=slope;
    goodness(zz,jjj)=Rsq;
    ks(zz,jjj)=power(10,intercept);
    junto(zz,:)=[slope Rsq time_max];  
%     
 %%%% Return only one row depending on R2 value (goodness of fit)
 B=goodness>0.975;
 res=junto(B,:,:);
 dev=res(end,:);
 dev1=dev(:,1);
 dev2=dev(:,2);
 dev3=dev(:,3); 
   umin=umin+nu;     
end
    medias(zz)=mean(exponents(zz,:));
    standard(zz)=std(exponents(zz,:));
    zz=zz+1;
end
end
