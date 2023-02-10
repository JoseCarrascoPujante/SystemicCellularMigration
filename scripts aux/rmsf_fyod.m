scale = 11.63;
dt = 0.5;
x = original_x / scale;
y = original_y / scale;
r = [x - x(1) y - y(1)];
r2 = r(:,1).^2+r(:,2).^2;
% u = r(2:end,:) - r(1:end-1,:);
% u2 = u(:,1).^2+u(:,2).^2;

taus = 1:(size(r2)-1);
data = zeros(size(taus));

for tau = taus
    data(tau) = rmsf(sqrt(r2),tau) / dt;
end

% Rmsf vs time plot
figure();
plot(dt * taus, data);
xlabel('tau');
ylabel('rmsf');
coeffs = polyfit(dt * taus(1:end), data, 1);
yfit = polyval(coeffs,dt * taus(1:end));
res_err = data - yfit;

% Log-log Rmsf vs time plot and linear regression for alpha calculation
figure();
logtaus = log(dt * taus);
logdata = log(data);
xlabel('Log(tau)');
ylabel('Log(RMSF(tau))');
p1 = polyfit(logtaus,logdata,1);
fit = polyval(p1,logtaus);
plot(logtaus,fit,'o','Color', [0.85, 0.1, 0.1, 0.1])
hold on
plot(logtaus, logdata,'k')

% Plot derivative of Rmsf vs time
figure();
alphas = diff(logdata) ./ diff(logtaus);
plot(dt * taus(2:end), alphas);

function res = rmsf(r, tau)
    res = sqrt(mean((r(1+tau:end)-r(1:end-tau)).^2));
end
