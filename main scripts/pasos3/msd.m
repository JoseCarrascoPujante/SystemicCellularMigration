function [slope,deltat] = msd(x, y, axis)
% msd - calculates the msd curve for 2-D dissusion 
%   x: vector of x positions
%   y: vector of y positions

% modified from
% http://stackoverflow.com/questions/7489048/calculating-mean-squared-displacement-msd-with-matlab
tau = 0.5; % s time between trajectory points
data = sqrt(x.^2 + y.^2);
nData = length(data); %# number of data points
numberOfDeltaT = floor(nData/4); %# for MSD, dt should be up to 1/4 of number of data points

sem = zeros(numberOfDeltaT, 1); 
deltat = zeros(numberOfDeltaT, 1); 
msdpts = zeros(numberOfDeltaT, 1); 
sqdisp = zeros(numberOfDeltaT, 1); 

%# calculate msd for all deltaT's

for dt = 1:numberOfDeltaT
   sqdisp = (data(1+dt:end) - data(1:end-dt)).^2;
   deltat(dt) = dt * tau;
   sem(dt) = std(sqdisp) / sqrt(length(msdpts) / dt);
   msdpts(dt) = mean(sqdisp); 
end
 [yy1,~]=fit(log10(deltat),log10(msdpts),'poly1');
 p=coeffvalues(yy1);
 slope=p(1);

%PLOTS

hold(axis, 'on')

plot(axis, log(deltat), log(msdpts),'color',rand(1,3))
% hold(axis, "off")
