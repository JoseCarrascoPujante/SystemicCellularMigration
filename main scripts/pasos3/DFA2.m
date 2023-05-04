function output1=DFA2(DATA,win_length,order)
      
       N=length(DATA);   
       n=floor(N/win_length); %redondeo al menor de puntos/tamaño, esto es, nº de intervalos
       N1=n*win_length;  %tamaño real de muestra, porque puede que queden puntos sueltos a no considerar
       y=zeros(N1,1);   %variable para almacenar y(k), primer término
       Yn=zeros(N1,1);  %%variable para almacenar yn(k), línea de tendencia de intervalo k de tamaño n
       
       fitcoef=zeros(n,order+1); %%variable para almacenar los coeficientes de la regresión lineal (de la linea de tendencia)
       mean1=mean(DATA(1:N1));  %%media del tamaño utilizado de la muestra
       
       for i=1:N1   
           y(i)=sum(DATA(1:i)-mean1);  %%y(k)=suma de 0 a k de los datos hasta k, menos la media total (esto es, cada dato menos la media, acumulados, para utilizar agrupados en "intervalos") INTEGRAL DISCRETA DE TODA LA SERIE
       end
       y=y';
       
       for j=1:n
           fitcoef(j,:)=polyfit(1:win_length,y(((j-1)*win_length+1):j*win_length),order);
       end
       
       for j=1:n
           Yn(((j-1)*win_length+1):j*win_length)=polyval(fitcoef(j,:),1:win_length);
       end
       
%        sum1=sum((y'-Yn).^2)/N1;
%        sum1=sqrt(sum1);
%        output1=sum1;
       
            
       %%notas para los intervalos solapados, con el ejemplo de 1500 puntos, e intervalos de 100      
            
       N2=length(DATA);   
       n2=floor(N2/win_length);
       N12=n2*win_length;
       y2=zeros(N12,1);
       Yn2=zeros(N12,1);  
       fitcoef2=zeros(n2,order+1); %%variable para almacenar los coeficientes de la regresión lineal (de la linea de tendencia)
       mean2=mean(DATA(floor((win_length)/2:N12-floor((win_length)/2))));  %%media del tamaño utilizado de la muestra (hay que quitarle lo que no usamos para calcularla)
       
        for i=1:N2    
           y2(i)=sum(DATA(1:i)-mean2);
        end
        for i=1:floor((win_length)/2)   %%%los primeros valores no van a ser usados porque el solapamiento empieza después
            y2(i)=0;
        end
        for i=N12-floor((win_length)/2):N12   %%%los ultimos valores no van a ser usados porque el solapamiento empieza después
            y2(i)=0;
        end  
             y2=y2';   
  
       for j2=1:n2-1
           fitcoef2(j2,:)=polyfit(1:win_length,y2(((j2)*win_length-floor((win_length)/2)+1):(j2+1)*win_length-floor(win_length/2)),order);  %%coeficientes desde 1 hasta 149 de los coef de la regresión por intervalos.
       end
          
       for j2=1:n2-1
           Yn2(((j2-1)*win_length+1):j2*win_length)=polyval(fitcoef2(j2,:),1:win_length);  %sus últimos 100 valores son 0, porque están desfasados los primeros 50 hacia delante
       end
       
      resta1=y'-Yn;     %%calculamos por separado las restas para concatenarlas después
       
      for i=1:N12-(win_length)      %%la resta de los solapados
          resta2(i)=y2(floor((win_length)/2)+i)-Yn2(i);
      end
      if N12-(win_length)>=1
      restatotal=[resta1; resta2'];
      else
      restatotal=[resta1];
      end
       %%Calculamos la fluctuación del conjunto
       sum2=sum((restatotal).^2/length(restatotal));
       sum2=sqrt(sum2);
       output1=sum2;
       
       