%% PROBLEM SIZE (NDOF) WILL BE N^3
N = floor (10000^(1/3));
printf ("ndofs = %d\n", N^3)

pkg load bim msh

m = bim3c_mesh_properties (msh3m_structured_mesh (1:N, 1:N, 1:N, 1, 1:6));
A = bim3a_laplacian (m, 1, 1);
M = bim3a_reaction (m, 1, 1);
b = bim3a_rhs (m, 1, 1);

function r = res (t, y, yp, A, M, b)
  t
  r = full (M * yp + (A+M) * y - ((1 - exp(-t+1)) .* (t>1) - (1 - exp(-t+5)) .* (t>5)) * b);
endfunction

function [dfdy, dfdyp] = jac (t, y, yp, A, M, b)
  dfdy = (A+M);
  dfdyp = M;
endfunction

options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'Jacobian', @(t, y, yp) jac(t, y, yp, A, M, b));

tic
[t, y] = ode15i (@(t, y, yp) res(t, y, yp, A, M, b), [0, 10],
                 zeros (N^3,1), zeros (N^3,1), options);
toc

plot (t, y(:,1), 'x-')
