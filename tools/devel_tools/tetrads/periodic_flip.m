% X -- [np x nt] array of particle x position over time
% Y -- [np x nt] array of particle y position over time
% Z -- [np x nt] array of particle z position over time
% tstep -- number of time steps (could just get from size(X,2)
% xl, yl, zl -- domain extents
function [X Y Z] = periodic_flip(X, Y, Z, N, tstep, xl, yl, zl);

for nn = 1:N
  for tt = 1:tstep-1
    if abs(X(nn,tt) - X(nn,tt+1)) >= 0.5*xl
      X(nn,tt+1:end) = X(nn,tt+1:end) + xl*sign(X(nn,tt) - X(nn,tt+1));
    end
    if abs(Y(nn,tt) - Y(nn,tt+1)) >= 0.5*yl
      Y(nn,tt+1:end) = Y(nn,tt+1:end) + yl*sign(Y(nn,tt) - Y(nn,tt+1));
    end
    if abs(Z(nn,tt) - Z(nn,tt+1)) >= 0.5*zl
      Z(nn,tt+1:end) = Z(nn,tt+1:end) + zl*sign(Z(nn,tt) - Z(nn,tt+1));
    end
  end
end

