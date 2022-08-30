1;

function res = rob (t, y, yp)
  res =[-(yp(1) + 0.04*y(1) - 1e4*y(2)*y(3));
        -(yp(2) - 0.04*y(1) + 1e4*y(2)*y(3) + 3e7*y(2)^2);
        y(1) + y(2) + y(3) - 1];
endfunction

function [DFDY, DFDYP] = jacfunsparse (t, y, yp)
  DFDY = sparse ([-0.04,           1e4*y(3),  1e4*y(2);
                   0.04, -1e4*y(3)-6e7*y(2), -1e4*y(2);
                      1,                  1,         1]);
  DFDYP = sparse ([-1,  0, 0;
                    0, -1, 0;
                    0,  0, 0]);
endfunction

opt = odeset ("Jacobian", @jacfunsparse, "AbsTol", 1e-7, "RelTol", 1e-7);
[t, y] = ode15i (@rob, [0, 100], [1; 0; 0], [-1e-4; 1e-4; 0], opt);