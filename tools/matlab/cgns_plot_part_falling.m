%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLUEBOTTLE-1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   Copyright 2012 - 2014 Adam Sierakowski, The Johns Hopkins University
% 
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%       http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
% 
%   Please contact the Johns Hopkins University to use Bluebottle for
%   commercial and/or for-profit applications.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cgns_plot_part_falling(casename, ts, te)
% CGNS_PLOT_PART_POSITION  Plot the particle positions from a BLUEBOTTLE-
%   generated CGNS file.
%
%   CGNS_PLOT_PART_POSITION(CASENAME, TS, TE, DT) plots the particle positions
%   from the simulation CASENAME at all times TS <= T <= TE, with DT matching
%   the output time interval set in record.config.
%
%   Example:
%     cgns_plot_part_position('simulation', 0, 1, 0.01) will plot the position
%     history for the CASENAME = 'simulation' for 0 <= T <= 1 outputted at an
%     interval of 0.01.

% read part time
[tstr tnum] = cgns_read_part_time(casename);
n = 1:1:length(tnum);

% 'cut' the time array between the ts and te values, inclusive
i = find(tnum<ts | tnum>te);
i(1) = []; % include the te value
n(i) = [];
tnum(i) = [];
ts = tstr{n(1)};

% calculate the number of timesteps
nt = length(n);

% find the number of particles using the initial timestep
[x, y, z] = cgns_read_part_position(casename, ts);
[u, v, w] = cgns_read_part_vel(casename, ts);
[fx, fy, fz] = cgns_read_part_force_total(casename, ts);
[ox, oy, oz] = cgns_read_part_omega(casename, ts);
[lx, ly, lz] = cgns_read_part_moment(casename, ts);

[np d] = size(x); % np = number of particles

% create data storage locations
T = zeros(nt, 1);
X = zeros(nt, np);
Y = zeros(nt, np);
Z = zeros(nt, np);
U = zeros(nt, np);
V = zeros(nt, np);
W = zeros(nt, np);
Fx = zeros(nt, np);
Fy = zeros(nt, np);
Fz = zeros(nt, np);
Ox = zeros(nt, np);
Oy = zeros(nt, np);
Oz = zeros(nt, np);
Lx = zeros(nt, np);
Ly = zeros(nt, np);
Lz = zeros(nt, np);

% read time and particle position data
T = tnum;
for i = 1:nt
  [X(i,:), Y(i,:), Z(i,:)] = cgns_read_part_position(casename, tstr{n(i)});
  [U(i,:), V(i,:), W(i,:)] = cgns_read_part_vel(casename, tstr{n(i)});
  [Fx(i,:), Fy(i,:), Fz(i,:)] = cgns_read_part_force_total(casename, tstr{n(i)});
  [Ox(i,:), Oy(i,:), Oz(i,:)] = cgns_read_part_omega(casename, tstr{n(i)});
  [Lx(i,:), Ly(i,:), Lz(i,:)] = cgns_read_part_moment(casename, tstr{n(i)});
end

% plot
subplot(5,3,1)
plot(T, X)
axis([0 6 -5 5])
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')
subplot(5,3,2)
plot(T, Y)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
subplot(5,3,3)
plot(T, Z)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
subplot(5,3,4)
plot(T, U)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$u$', 'interpreter', 'latex')
subplot(5,3,5)
plot(T, V)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$v$', 'interpreter', 'latex')
subplot(5,3,6)
plot(T, W)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$w$', 'interpreter', 'latex')
subplot(5,3,7)
plot(T, Fx)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_x$', 'interpreter', 'latex')
subplot(5,3,8)
plot(T, Fy)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_y$', 'interpreter', 'latex')
subplot(5,3,9)
plot(T, Fz)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$F_z$', 'interpreter', 'latex')
subplot(5,3,10)
plot(T, Ox)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Omega_x$', 'interpreter', 'latex')
subplot(5,3,11)
plot(T, Oy)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Omega_y$', 'interpreter', 'latex')
subplot(5,3,12)
plot(T, Oz)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$\Omega_z$', 'interpreter', 'latex')
subplot(5,3,13)
plot(T, Lx)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$L_x$', 'interpreter', 'latex')
subplot(5,3,14)
plot(T, Ly)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$L_y$', 'interpreter', 'latex')
subplot(5,3,15)
plot(T, Lz)
%hold on
xlabel('$t$', 'interpreter', 'latex')
ylabel('$L_z$', 'interpreter', 'latex')

dt = (T(end) - T(1)) / (nt-1);
[XX f] = pwelch( X, [], [], [], 1/dt);
[YY f] = pwelch( Y, [], [], [], 1/dt);
figure
semilogx(f,XX,f,YY)

%%%% interpolate force signal to constant time interval dt
%% (turns out it makes very little difference)

%Fx_int = zeros(nt,1);
%for i = 1:nt
%  if i == 1
%    Tdt = 0;
%    T0 = 0;
%    T1 = T(i);
%    T2 = T(i+1);
%    dT0 = T2-T1;
%    dT1 = T2-T1;
%    F0 = 0;
%    F1 = Fx(i);
%    F2 = Fx(i+1);
%    dF0 = F2-F1;
%    dF1 = F2-F1;
%    dFdT0 = dF0/dT0;
%    dFdT1 = dF1/dT1+1;
%  elseif i == nt
%    Tdt = (i-1) * dt;
%    T0 = T(i-1);
%    T1 = T(i);
%    T2 = 0;
%    dT0 = T1-T0;
%    dT1 = T1-T0;
%    F0 = Fx(i-1);
%    F1 = Fx(i);
%    F2 = 0;
%    dF0 = F1-F0;
%    dF1 = F1-F0;
%    dFdT0 = dF0/dT0+1;
%    dFdT1 = dF1/dT1;
%  else
%    Tdt = (i-1) * dt;
%    T0 = T(i-1);
%    T1 = T(i);
%    T2 = T(i+1);
%    dT0 = T1-T0;
%    dT1 = T2-T1;
%    F0 = Fx(i-1);
%    F1 = Fx(i);
%    F2 = Fx(i+1);
%    dF0 = F1-F0;
%    dF1 = F2-F1;
%    dFdT0 = dF0/dT0;
%    dFdT1 = dF1/dT1;
%  end
%  if abs(dFdT0) < abs(dFdT1)
%    Fx_int(i) = F1 + dFdT0*(T1-Tdt);
%  else
%    Fx_int(i) = F2 - dFdT1*(T2-Tdt);
%  end
%end
%TTdt = 0:dt:T(end);
%figure
%plot(T,Fx,TTdt,Fx_int)

[UX f] = pwelch( U, [], [], [], 1/dt);
[VY f] = pwelch( V, [], [], [], 1/dt);
figure
loglog(f,0.5*(UX+VY))
xlabel('$f$','interpreter','latex')
ylabel('$u_\perp$ SPD','interpreter','latex')

[FX f] = pwelch( Fx, [], [], [], 1/dt);
[FY f] = pwelch( Fy, [], [], [], 1/dt);
figure
loglog(f,0.5*(FX+FY))
xlabel('$f$','interpreter','latex')
ylabel('$F_\perp$ SPD','interpreter','latex')

figure
plot(f,0.5*(FX+FY) ./ (0.5*(UX+VY)))
xlabel('$f$','interpreter','latex')
ylabel('$F_\perp / u_\perp$ SPD','interpreter','latex')

%figure
%for i = 1:nt-1
%TT(i) = (T(i+1) - T(i))/dt;
%TT_int(i) = (TTdt(i+1) - TTdt(i))/dt;
%end
%plot(1:nt-1,TT,1:nt-1,TT_int)
