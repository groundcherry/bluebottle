function [T] = form_tetrads(r0, X, Y, Z, dom, tol);

T = -ones(4,1);

% Superpose periodicity to make search easier
% TODO: explicitly check periodicity
N = dom.N;
X = [X; X + dom.xl; X - dom.xl];
Y = [Y; Y + dom.yl; Y - dom.yl];
Z = [Z; Z + dom.zl; Z - dom.xl];
%Z = [Z; Z         ; Z         ];
Nper = 3*N;

% Loop over ~all particles and find EW,NS,TB coordinates
for n1 = 1:dom.N
  % Set target particle
  p1.n = n1;
  p1.x = X(n1);
  p1.y = Y(n1);
  p1.z = Z(n1);

  % Desired position for each direction -- correct for periodicity
  eDes = p1.x + r0;
  wDes = p1.x - r0;
  nDes = p1.y + r0;
  sDes = p1.y - r0;
  tDes = p1.z + r0;
  bDes = p1.z - r0;

  % Initialize comparision -- want to find closest to desired position
  eMin = 10000*dom.r*dom.xl;
  ePoss = -1;
  wMin = 10000*dom.r*dom.xl;
  wPoss = -1;
  nMin = 10000*dom.r*dom.yl;
  nPoss = -1;
  sMin = 10000*dom.r*dom.yl;
  sPoss = -1;
  tMin = 10000*dom.r*dom.zl;
  tPoss = -1;
  bMin = 10000*dom.r*dom.zl;
  bPoss = -1;


  for n2 = 1:Nper
    nmod = mod(n2,N);
    if nmod == 0
      nmod = N;
    end
    if nmod == n1
      continue;
    end

    % Pull target part info
    pTarg.n = nmod;
    pTarg.x = X(n2);
    pTarg.y = Y(n2);
    pTarg.z = Z(n2);

    % Calculate distance from desired position
    eDist = sqrt((eDes - pTarg.x)^2 + (p1.y - pTarg.y)^2 + (p1.z - pTarg.z)^2);
    wDist = sqrt((wDes - pTarg.x)^2 + (p1.y - pTarg.y)^2 + (p1.z - pTarg.z)^2);
    nDist = sqrt((p1.x - pTarg.x)^2 + (nDes - pTarg.y)^2 + (p1.z - pTarg.z)^2);
    sDist = sqrt((p1.x - pTarg.x)^2 + (sDes - pTarg.y)^2 + (p1.z - pTarg.z)^2);
    tDist = sqrt((p1.x - pTarg.x)^2 + (p1.y - pTarg.y)^2 + (tDes - pTarg.z)^2);
    bDist = sqrt((p1.x - pTarg.x)^2 + (p1.y - pTarg.y)^2 + (bDes - pTarg.z)^2);

    % If the distance is acceptable AND it's closer than another stored particle
    % save it
    if (eDist < tol & eDist < eMin)
      ePoss = pTarg.n;
      eMin = eDist;
    end
    if (wDist < tol & wDist < wMin)
      wPoss = pTarg.n;
      wMin = wDist;
    end
    if (nDist < tol & nDist < nMin)
      nPoss = pTarg.n;
      nMin = nDist;
    end
    if (sDist < tol & sDist < sMin)
      sPoss = pTarg.n;
      sMin = sDist;
    end
    if (tDist < tol & tDist < tMin)
      tPoss = pTarg.n;
      tMin = tDist;
    end
    if (bDist < tol & bDist < bMin)
      bPoss = pTarg.n;
      bMin = bDist;
    end
  end

  % Store tetrads for current particle n1
  % -- ENT, ENB, ESB, EST
  if (ePoss ~= -1)
    if (nPoss ~= -1)
      if (tPoss ~= -1)
        temp = [n1; ePoss; nPoss; tPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
      if (bPoss ~= -1)
        temp = [n1; ePoss; nPoss; bPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
    end
    if (sPoss ~= -1)
      if (tPoss ~= -1)
        temp = [n1; ePoss; sPoss; tPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
      if (bPoss ~= -1)
        temp = [n1; ePoss; sPoss; bPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
    end
  % -- WNT, WNB, WSB, WST
  elseif (wPoss ~= -1)
    if (nPoss ~= -1)
      if (tPoss ~= -1)
        temp = [n1; wPoss; nPoss; tPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
      if (bPoss ~= -1)
        temp = [n1; wPoss; nPoss; bPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
    end
    if (sPoss ~= -1)
      if (tPoss ~= -1)
        temp = [n1; wPoss; sPoss; tPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
      if (bPoss ~= -1)
        temp = [n1; wPoss; sPoss; bPoss];
        anyDuplicates = ~all(diff(sort(temp)));
        if anyDuplicates ~= 1
          T = [T, temp];
        end
      end
    end
  end
end

% Remove initialization
if size(T,2) > 1
  T = T(:,2:end);
end
% sort columns of T from small to large
T = sort(T, 1);
% find only unique entries
% array gymnastics required due to no 'columns' option
T = unique(T', 'rows')';
