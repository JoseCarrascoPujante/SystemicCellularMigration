function [Rsq,slope,intercept]=rsquare(x,y)

p = polyfit(x,y,1);  % linear coefficients, p(1) is the slope, p(2) the intercept 

slope=p(1);

intercept=p(2);

yfit = polyval(p,x);  % fully equivalent to yfit =  p(1) .* x + p(2);

yresid = y - yfit;  % Compute the residual

SSresid = sum(yresid.^2);    % Residual sum of squares:

SStotal = (length(y)-1) * var(y);  % Compute the total sum of squares of y 

Rsq= 1 - SSresid/SStotal;  % compute R2

end



