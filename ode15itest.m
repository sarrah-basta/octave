1;

function res = fun (t, y, yp, lambda1, lambda2)
  res =[(yp(1) - lambda1* y(1));
        (yp(2) - lambda2* y(2))];
endfunction

function [DFDY, DFDYP] = jacfunsparse (t, y, yp, lambda1, lambda2)
  DFDY = - sparse ([lambda1, 0; 0, lambda2]);
  DFDYP = sparse ([1, 0;
                  0, 1]);
endfunction

opt = odeset ("Jacobian", @(t,y,yp) jacfunsparse(t,y,yp,-1,-2), "AbsTol", 1e-7, "RelTol", 1e-7);
[t, y] = ode15i (@(t,y,yp) fun(t,y,yp,-1,-2), [0, 10], [1; 1], -[1;2], opt)