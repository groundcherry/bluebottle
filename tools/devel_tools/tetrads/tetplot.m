clear all; close all; clc;
load tetrad_stats.mat
style = {'k', 'b', 'r', 'g', 'm', 'c'};

% Plot volume
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(tnum(time), avgVol(rr,:)./(4/3*pi*dom.r^3), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog(tnum(time), 0.01*tnum(time).^(2), 'k--')
xlabel('Time')
ylabel('<V>/(4/3 \pi r^3)')
title('Tetrad Volume')
leg = [leg {'t^{2}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg
print('vol', '-dpdf', '-r300')
error('\n')

% Plot radius of gyration
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  loglog(tnum(time), avgRsq(rr,:).^(1/2)./dom.r, style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
loglog(tnum(time), 1.25*tnum(time).^(3/4), 'k--')
ylim([2*10^0, 10^(1.5)]);
xlabel('Time')
ylabel('<R^2>^{1/2}/r')
title('Tetrad Radius of Gyration')
leg = [leg {'t^{3/4}'}];
legend(leg, 'Location', 'SouthEast')
clearvars leg
print('rsq', '-dpdf', '-r300')
 
% Plot lambda
figure
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(tnum(time), avgLambda(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
xlabel('Time [ms]')
ylabel('\Lambda')
title('\Lambda = V^{2/3}/R^2')
legend(leg)
clearvars leg
print('lambda', '-dpdf', '-r300')

% Plot I1
figure
subplot(3,1,1)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(tnum(time), avgI1(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('I1')
% Plot I2
subplot(3,1,2)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(tnum(time), avgI2(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
ylabel('I2')
% Plot I3
subplot(3,1,3)
for rr = 1:length(r0)
  if r0(rr) == -1
    continue;
  end
  semilogx(tnum(time), avgI3(rr,:), style{rr})
  hold on
  leg{rr} = ['r0 = ' num2str(r0(rr))];
end
xlabel('Time [ms]')
ylabel('I3')
legend(leg, 'Location', 'NorthEast')
print('ifactor', '-dpdf', '-r300')
