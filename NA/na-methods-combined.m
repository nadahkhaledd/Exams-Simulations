1;
%%% T_newton_sys — NONLINEAR SYSTEM with provided newtonsys.m

function [x,res,niter,difv,x_vect] = newtonsys(Ffun, Jfun, x0, tol, kmax)
  k = 1; x_vect = x0; nd = tol + 1; difv=[]; x = x0;
  while nd >= tol && k < kmax
    J = Jfun(x);
    F = Ffun(x);
    delta = - J\F;
    x = x + delta;
    nd = norm(delta);
    difv = [difv; nd];
    k = k+1;
    x_vect(:,k) = x;
  end
  res = norm(Ffun(x));
  niter = k;
end

%%% write f as column vector returning function
function r = residual(x)
  r(1,1) = % TODO: f1(x1,x2);
  r(2,1) = % TODO: f2(x1,x2);
end
f = @residual;

%%% Jacobian
function J = jacobian(x)
  J(1,1) = % TODO: df1/dx1;
  J(1,2) = % TODO: df1/dx2;
  J(2,1) = % TODO: df2/dx1;
  J(2,2) = % TODO: df2/dx2;
end
J = @jacobian;

%%% stopping criterion
%% newtonsys uses the increment norm(delta) as an estimate for the error.
%% This is reliable when the method converges at least quadratically:
%% under quadratic convergence the increment and the true error are of the same order.

%%% run with given x0, tol, maxit
x0 = [% TODO];
tol = 1e-6;
maxit = % TODO: 100 or 1000;
[x, res, niter, difv, x_vect] = newtonsys(f, J, x0, tol, maxit);
x
f(x)
det(J(x))
difv
x_vect

%%% interpretation
%% If det(J) at the solution is non-zero → expect (at least) quadratic convergence.
%% Estimate convergence order a-posteriori with:
%%   diff(log(abs(difv(2:end)))) ./ diff(log(abs(difv(1:end-1))))
%% Different x0 may converge to a different root (Newton is locally convergent).


%%% T_newton_scalar — SCALAR NONLINEAR EQUATION f(x)=0
function [x, niter, difv] = newton(f, df, x0, tol, kmax)
  x = x0; k = 1; nd = tol + 1; difv = [];
  while nd >= tol && k < kmax
    delta = -f(x)/df(x);
    x = x + delta;
    nd = abs(delta);
    difv = [difv; nd];
    k = k + 1;
  end
  niter = k;
end

function [x, niter] = bisection(f, a, b, tol, kmax)
  k = 0;
  while (b - a)/2 > tol && k < kmax
    x = (a + b)/2;
    if f(a)*f(x) < 0
      b = x;
    else
      a = x;
    end
    k = k + 1;
  end
  x = (a + b)/2; niter = k;
end


%%% T_ode — CAUCHY PROBLEM y' = f(t,y), y(0) = y0

function y = heun(odefun, y0, t)
  y(:,1) = y0;
  for it = 2:numel(t)
    dt = t(it) - t(it-1);
    fstar = odefun(t(it-1), y(:,it-1));
    ystar = y(:,it-1) + dt * fstar;
    y(:,it) = y(:,it-1) + dt/2 * (fstar + odefun(t(it), ystar));
  endfor
end

function y = forward_euler(odefun, y0, t)
  y(:,1) = y0;
  for it = 2:numel(t)
    dt = t(it) - t(it-1);
    y(:,it) = y(:,it-1) + dt * odefun(t(it-1), y(:,it-1));
  endfor
end

%%% a) — solve with given h
f = @(t,y) % TODO;
y0 = % TODO;
t1 = 0:% TODO h:% TODO T;
y1 = heun(f, y0, t1);
%% report y at the requested t-value

%%% b) — plot
plot(t1, y1)
hold all
plot(t2, y2)
legend('h = ...', 'h = ...')
xlabel t
ylabel y(t)
print -dpng ex_pointb.png

%%% c) — stability (Heun)
%% Apply method to y' = lambda y to get u(n+1) = R(h*lambda) * u(n)
%%   Forward Euler:  R(z) = 1 + z
%%   Heun:           R(z) = 1 + z + z^2/2
%%   RK4:            R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24
%% Stability: |R(z)| < 1 with z = h * (min ∂f/∂y).
%% Solve the inequality for h.

%%% d) — error vs. "exact" with tiny h
T = % TODO;
h_ref = 1e-3;
tex = 0:h_ref:T;
yex = heun(f, y0, tex);
hvect = [0.1, 0.05, 0.025, 0.0125];
for ii = 1:numel(hvect)
  h = hvect(ii);
  t = 0:h:T;
  y = heun(f, y0, t);
  err(ii) = abs(y(end) - yex(end));
endfor
err
%% err results

%%% e) — theory comment
%% Heun is order 2 → halving h should divide error by 2^2 = 4.
%% Check err(1:end-1) ./ err(2:end) ≈ 4.
%% Forward Euler order 1 → factor 2. RK4 order 4 → factor 16.


%%% T_quadrature — COMPOSITE NUMERICAL INTEGRATION

%%% a) — implement composite rule (signature given in exam!)
function I = comp_simpson(a, b, M, f)
  x = linspace(a, b, M+1);
  I = 0;
  for ii = 1:M
    I += (1/6*f(x(ii)) + 2/3*f((x(ii)+x(ii+1))/2) + 1/6*f(x(ii+1))) * (x(ii+1)-x(ii));
  end
end

function I = comp_trap(a, b, M, f)
  x = linspace(a, b, M+1);
  I = 0;
  for ii = 1:M
    I += (f(x(ii)) + f(x(ii+1)))/2 * (x(ii+1) - x(ii));
  end
end

function I = comp_midpoint(a, b, M, f)
  x = linspace(a, b, M+1);
  I = 0;
  for ii = 1:M
    I += f((x(ii)+x(ii+1))/2) * (x(ii+1)-x(ii));
  end
end

%%% b) — evaluate for n = 4, 8, 16, 32
f = @(x) % TODO;
format long e
ii = 1;
for n = [4 8 16 32]
  I(ii) = comp_simpson(% TODO a, % TODO b, n, f);
  ii = ii + 1;
end
%% I values

%%% Question c) — error vs exact
Iex = % TODO;
err = abs(I - Iex)
%% err(1:end-1) ./ err(2:end) should ≈ 16 for Simpson (order 4), 4 for trap (order 2)

%%% d) — loglog plot
loglog(1./[4 8 16 32], err, 1./[4 8 16 32], 1./[4 8 16 32].^4)
legend('error', 'reference order 4 slope', 'location', 'southeast')
axis tight
xlabel('h')
ylabel('error')
print -dpng ex_pointd.png
%% Parallel slopes in loglog → numerical order matches theoretical order.

%%% e) — exact result trick
%% Simpson exact for polynomials of degree ≤ 3.
%% Trapezoidal exact for polynomials of degree ≤ 1.
%% Midpoint exact for polynomials of degree ≤ 1.
%% If f is a polynomial of degree ≤ exactness, error is 0 for any n.


%%% T_interp — POLYNOMIAL INTERPOLATION

%% Lagrangian / Vandermonde approach
x = [% TODO nodes];
y = [% TODO values];
n = length(x) - 1;       % polynomial degree
c = polyfit(x, y, n);    % coefficients, highest power first
xx = linspace(min(x), max(x), 200);
yy = polyval(c, xx);
plot(x, y, 'o', xx, yy, '-')
xlabel x; ylabel y
legend('data', 'interpolant')
print -dpng ex_interp.png

%% Evaluate at a specific point
x_eval = % TODO;
y_eval = polyval(c, x_eval)


%%% T_lsq — LEAST SQUARES BEST FIT

x = [% TODO];
y = [% TODO];
m = % TODO degree of fit;       % m < n (number of nodes - 1)
c = polyfit(x, y, m);           % LS solution
xx = linspace(min(x), max(x), 200);
yy = polyval(c, xx);
plot(x, y, 'o', xx, yy, '-')
xlabel x; ylabel y
legend('data', 'LS fit')
print -dpng ex_lsq.png

%% normal equations
% A = vander(x);                  % full Vandermonde
% A = A(:, end-m:end);            % keep columns for degree m
% c = (A' * A) \ (A' * y(:));     % normal equations
%% (better numerically: c = A \ y(:) — uses QR internally)


%%% T_linsys_direct — DIRECT METHODS FOR Ax = b

A = % TODO;
b = % TODO;

%% General: backslash uses LU with pivoting
x = A \ b;

%% Explicit LU
[L, U, P] = lu(A);
y = L \ (P*b);
x = U \ y;

%% Cholesky (only if A is SPD)
R = chol(A);            % A = R'*R
y = R' \ b;
x = R \ y;

%% Conditioning check
cond(A)
%% cond(A) >> 1 → ill-conditioned → small perturbations in b cause large errors in x


%%% T_linsys_iter — JACOBI / GAUSS-SEIDEL

function [x, niter, res_vec] = jacobi(A, b, x0, tol, kmax)
  D = diag(diag(A));
  N = D;
  M = D - A;          % so A = N - M, iteration: x = N\(M*x + b)
  x = x0; k = 0; res_vec = [];
  r = b - A*x;
  while norm(r) > tol*norm(b) && k < kmax
    x = N \ (M*x + b);
    r = b - A*x;
    res_vec = [res_vec; norm(r)];
    k = k + 1;
  end
  niter = k;
end

function [x, niter, res_vec] = gauss_seidel(A, b, x0, tol, kmax)
  N = tril(A);
  M = N - A;
  x = x0; k = 0; res_vec = [];
  r = b - A*x;
  while norm(r) > tol*norm(b) && k < kmax
    x = N \ (M*x + b);
    r = b - A*x;
    res_vec = [res_vec; norm(r)];
    k = k + 1;
  end
  niter = k;
end

%% Convergence check — spectral radius of the iteration matrix B = N \ M
B_jac = eye(size(A)) - diag(1./diag(A)) * A;       % iteration matrix for Jacobi
rho = max(abs(eig(B_jac)))
%% rho < 1 ⇔ method converges (independent of x0). Smaller rho → faster.


%%% T_diff — NUMERICAL DIFFERENTIATION (rare but easy)

%% Forward / Backward / Centered first derivative at x0 with step h
% forward:  ( f(x0+h) - f(x0) ) / h           ; order 1
% backward: ( f(x0) - f(x0-h) ) / h           ; order 1
% centered: ( f(x0+h) - f(x0-h) ) / (2*h)     ; order 2
% second derivative centered:
%           ( f(x0+h) - 2*f(x0) + f(x0-h) ) / h^2   ; order 2