
function [dev1,dev2,dev3,tc2]=amebas5(u_all, rmsfhandle, type)

scale_time=0.5; %time between frames (seconds)

%%%%%%%%%%%%%%%%%
% calculating F %
%%%%%%%%%%%%%%%%%
tc=24;
numseries=1;
step=125;

exponents=zeros(tc,numseries);
goodness=zeros(tc,numseries);
ks=zeros(tc,numseries);
 
clear title xlabel ylabel;
tc2=tc*step;
zz=1;

for time_max=step:step:tc2   % plot one graph every "step" data points
        
    umin=1;     
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
            
    for j=1:dim
        
        y(j)=conu.compute_y(l,j);
        
    end
    
    % root mean square fluctuation of the displacement
    
    F=zeros(dim,1);
    
    for j=1:dim
        
        F(j)=conu.compute_F(y,j);
        
    end
   
    %%%%%%%%%%%%%
    % F fitting %
    %%%%%%%%%%%%%
     
    time = (1:dim)'.*scale_time*du;
    idx = find( time(1:time_max)>0 & F(1:time_max)>0 );
    [Rsq,slope,intercept] = conu.rsquare (log10(time(idx)),log10(F(idx)));
    
    % alpha=slope;
    % v=time.^(alpha);
    % k=power(10,intercept/slope);
    k=F(1);

    exponents(zz,1)=slope;
    goodness(zz,1)=Rsq;
    ks(zz,1)=power(10,intercept);
    junto(zz,:)=[slope Rsq time_max];  

    % Return only one row depending on R2 value (goodness of fit)
    B=goodness>0.98;
    res=junto(B,:);
    dev=res(end,:);
    dev1=dev(:,1);
    dev2=dev(:,2);
    dev3=dev(:,3); 
    % umin=umin+nu;

    %%%%%%%%%%%%
    % Plotting %
    %%%%%%%%%%%%

    % loglog plot of rmsf F versus l step non-shuffled coordinates
    if time_max == tc2
        max_corr = res(end,3); % Use max tc2 when R2>0.99   
        hold(rmsfhandle, 'on')
        if strcmp(type,'orig')
            loglog(time(1:max_corr), F(1:max_corr), 'or', 'MarkerSize', 2);
        elseif strcmp(type,'shuff')
            loglog(time(1:max_corr), F(1:max_corr), 'ob', 'MarkerSize', 2);
        end
        xlabel('Log({\itl}(s))');
        ylabel('Log(F({\itl}))');
        
        % Plot regression fit line for the original data
        vmax=time.^(res(end,1));
        loglog(time(1:max_corr),(k/vmax(1))*vmax(1:max_corr),'k--',...
            'LineWidth',1);
        text(time(7),(k/vmax(1))*vmax(200),strcat('\alpha=',...
            num2str(round(res(end,1),3))))
    
        % Plot alpha=0.5 line with the original data
        v2 = time.^(.5);
        loglog(time(1:max_corr),(k/v2(1))*v2(1:max_corr),'k--','LineWidth',1);
        text(time(4),(k/v2(1))*v2(2),'\alpha=0.5, uncorrelated')
        set(rmsfhandle,'xscale','log')
        set(rmsfhandle,'yscale','log')
        % disp(goodness(1))
    end
    % medias(zz)=mean(exponents(zz,:));
    % standard(zz)=std(exponents(zz,:));
    zz=zz+1;
end
end
