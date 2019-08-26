clear;
close all;

pkg load symbolic
format long

syms x;                               % x (is a global symbol)
global ITER_LIMIT;                    % Limit of iterations
ITER_LIMIT = 10000;


% ============================== Method 1 ====================================
function [xAprox, iter] = sne_fd_1(expr, x0, tol, graf)
  %{
  Steffensen's Method
  
  Arguments:
    f {string} -- polynomial whose solution must be found
    x0 {float, int} -- initial value to start iterations
    tol {float, int} -- tolerance that indicates the stop condition
    graf {int} -- flag that indicates if a graph must be done

  Returns:
    xAprox {float} -- root approximation
    _iter {int} -- amount of iterations required
  %}

  % -------------------------- Local variables -------------------------------
  global ITER_LIMIT
  iter = 0;                           % Amount of iterations
  f = function_handle(expr);          % Use expression as function
  xAprox = x0;                        % Initialize x_n
  xNext = 0;                          % Initialize x_(n+1)
  error = abs(f(xAprox));             % First errror calculation
  
  % ------------------------- Steffensen's Method ---------------------------
  try
    while (error > tol)
      if (iter > ITER_LIMIT)
        warning("Iteration limit reached");
        return;
      endif
      div = f(xAprox + f(xAprox)) - f(xAprox);
      if (div ~= 0)
        xNext = xAprox - (f(xAprox)^2) / div;
      else
        warning("[Math error]: Division by zero");
        return;
      end
      xAprox = xNext;
      error = abs(f(xAprox));
      iter++;
    endwhile
  catch
    warning("[Math error]: Problem found executing the method");
    return;
  end
endfunction


% ============================== Method 2 ====================================
%
% Yun-Petkovic Method
%
% Métodos iterativos aplicados a la ecuación de Kepler. Page 112.
%
% Arguments:
%   f  {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   x1 {float, int} - interval low value
%   x2 {float, int} - interval high value
%   tol {float, int} - tolerance that indicates the stop condition
%   graf {int} - flag that indicates if a graph must be done
%
% Returns:
%   xAprox {float} - root approximation
%   iter {int} - amount of iterations required
%
function [xAprox, iter] = sne_fd_2(f, x0, x1, x2, tol, graf=1)
  if (typeinfo(f) != "string")
    error("f must be a string");
  endif

  if (graf != 0 && graf != 1)
    error("graf must be 0 or 1");
  endif

  syms fx(x);
  fx(x) = f;
  xn = [x0];
  iter = 0;

  try
    error = [abs(double(fx(xn(end))))];
    hk = 1;
    ak = x1;
    bk = x2;

    while (error(end) > tol)
      xk = xn(end);

      xk_next = xk - double(fx(xk)) * (2 * hk / (double(fx(bk)) - double(fx(ak))));

      xn = [xn xk_next];
      error = [error abs(double(fx(xk_next)))];
      
      hk = xk_next - xk;
      ak = xk_next - hk;
      bk = xk_next + hk;

      iter += 1;
    endwhile

    if (graf == 1)
      k = 0:1:iter;
      plotFunction(k, error, "Yun-Petkovic Method")
    endif
  catch
    error("f has an unknown function.");
  end
  
  xAprox = xn(end);
  return;
endfunction


% ============================== Method 3 ====================================
%
% Jain Method
%
% Metodos iterativos optimos para la resolucion de ecuaciones no lineales. Page 1. Equation 1.
%
% Arguments:
%   f  {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   tol {float, int} - tolerance that indicates the stop condition
%   graf {int} - flag that indicates if a plot must be done
%
% Returns:
%   xAprox {float} - root approximation
%   iter {int} - amount of iterations required
%
function [xAprox, iter] = sne_fd_3(f, x0, x1, x2, tol, graf=1)
  if (typeinfo(f) != "string")
    error("f must be a string");
  endif

  if (graf != 0 && graf != 1)
    error("graf must be 0 or 1");
  endif

  syms fx(x);
  fx(x) = f;
  xn = [x0];
  iter = 0;

  try
    error = [abs(double(fx(xn(end))))];

    while (error(end) > tol)
      xk = xn(end);

      y = double(fx(xk));
      yk = steffensen_method(f, x0, iter);

      xk_next = xk - y**3 / ((double(fx(xk + y)) - y) * (y - double(fx(yk))));

      xn = [xn xk_next];
      error = [error abs(double(fx(xk_next)))];

      iter += 1;
    endwhile

    if (graf == 1)
      k = 0:1:iter;
      plotFunction(k, error, "Jain Method")
    endif
  catch
    error("f has an unknown function.");
  end
  
  xAprox = xn(end);
  return;
endfunction


% ============================== Method 4 ====================================
%





% =========================== Auxiliaries functions =============================

% Steffensen Method
%
% This function is used to calculate some necessary values for other functions
%
% Arguments:
%   exp {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   n {int} - number of iterations
%
% Returns:
%   xAprox {float} - root approximation
%
function xAprox = steffensen_method(exp, x0, n)
  syms fx(x);
  fx(x) = exp;
  xAprox = x0;
  iter = 0;
  
  while (iter <= n)
    yk = double(fx(xAprox));
    df = double(diff(fx)(xAprox));
    df2 = double(diff(diff(fx))(xAprox));
    
    div = (2 * df^2 - yk * df2);
    xAprox = xAprox - (2 * yk * df) / div;
    iter += 1;
  endwhile
  return;
endfunction


% This function is used to plot iterations vs error
%
% Arguments:
%   k {iterable} - an iterable with x axis values
%   error {iterable} - an iterable with y axis values
%   title {string} - plot title
%
% Returns:
%   This function doesn't return
function plotFunction(k, error, plotTitle)
  plot(k, error);
  title(plotTitle);
  xlabel("Iterations k");
  ylabel("Error |f(xk)|");
endfunction
