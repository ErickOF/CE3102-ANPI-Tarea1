clear;
close all;

pkg load symbolic
format long


syms x;                               % x (is a global symbol)
global ITER_LIMIT;                    % Limit of iterations
ITER_LIMIT = 99;


% ============================== Method 1 ====================================
function [xAprox, _iter] = sne_ud_1(expr, x0, tol, graf)
  %{
  Halley's Method
 
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
  _iter = 0;                          % Amount of iterations
  f = function_handle(expr);          % Use expression as function
  xAprox = x0;                            % Initialize x_n
  xNext = 0;                          % Initialize x_(n+1)
  error = abs(f(xAprox));                 % First errror calculation
  
  % --------------------------- Halley's Method ------------------------------
  try
    while (error > tol)
      if (_iter > ITER_LIMIT)
        warning("Iteration limit reached");
        return;
      endif
      fPrime = function_handle(diff(expr));
      fPrime2 = function_handle(diff(diff(expr)));
      div = 2 * fPrime(xAprox)^2 - f(xAprox) * fPrime2(xAprox);
      if (div ~= 0)
        xNext = xAprox - (2* f(xAprox) * fPrime(xAprox)) / div;
      else
        warning("[Math error]: Division by zero");
        return;
      end
      xAprox = xNext;
      error = abs(f(xAprox));
      _iter++;
    endwhile
  catch
    warning("[Math error]: Problem found executing the method");
    return;
  end
endfunction


% ============================== Method 2 ====================================
function [xAprox, _iter] = sne_ud_2(expr, x0, tol, graf)
  %{
  Frontini's and Sormani's Method
  
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
  _iter = 0;                          % Amount of iterations
  f = function_handle(expr);          % Use expression as function
  xAprox = x0;                            % Initialize x_n
  xNext = 0;                          % Initialize x_(n+1)
  error = abs(f(xAprox));                 % First errror calculation
  
  % ------------------- Frontini's and Sormani's Method----=------------------
  try
    while (error > tol)
      if (_iter > ITER_LIMIT)
        warning("Iteration limit reached");
        return;
      endif
      fPrime = function_handle(diff(expr));
      div1 = fPrime(xAprox);
      if (div1 ~= 0)
        div2 = fPrime(xAprox - 0.5 * f(xAprox) / div1);
        if (div2 ~= 0)
          xNext = xAprox - f(xAprox)/div2;
        else
          warning("[Math error]: Division by zero");
        end 
      else
        warning("[Math error]: Division by zero");
        return;
      end
      xAprox = xNext;
      error = abs(f(xAprox));
      _iter++;
    endwhile
  catch
    warning("[Math error]: Problem found executing the method");
    return;
  end
endfunction


% ============================== Method 3 ====================================
%
% Chebyshev Method
%    
% Metodos iterativos aplicados a la ecuación de Kepler. Page 62.
%
% Arguments:
%   f {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   tol {float, int} - tolerance that indicates the stop condition
%   graph {int} - flag that indicates if a graph must be done
%
% Returns:
%   xAprox {float} - root approximation
%   iter {int} - amount of iterations required
%
function [xAprox, iter] = sne_ud_3(f, x0, tol, graf=1)
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
      df = diff(fx);
      df2 = diff(diff(fx));

      Lf = y * double(df2(xk)) / double(df(xk))^2;

      xk_next = xk - (1 + 0.5 * Lf) * y / double(df(xk));

      xn = [xn xk_next];
      error = [error abs(double(fx(xk_next)))];

      iter += 1;
    endwhile

    if (graf == 1)
      k = 0:1:iter;
      plotFunction(k, error, "Chebyshev Method");
    endif
  catch
    error("f has an unknown function.");
  end
  
  xAprox = xn(end);
  return;
endfunction


% ============================== Method 4 ====================================
%
% Newton-Secant Method
%
% Arguments:
%   f {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   tol {float, int} - tolerance that indicates the stop condition
%   graf {int} - flag that indicates if a graf must be done
%
% Returns:
%   xAprox {float} - root approximation
%   iter {int} - amount of iterations required
%
function [xAprox, iter] = sne_ud_4(f, x0, tol, graf=1)
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
      df = diff(fx);
      yk = xk - y / double(df(xk));

      xk_next = xk - (y / (y - double(fx(yk)))) * y / double(df(xk));

      xn = [xn xk_next];
      error = [error abs(double(fx(xk_next)))];

      iter += 1;
    endwhile

    if (graf == 1)
      k = 0:1:iter;
      plotFunction(k, error, "Newton-Secant Method")
    endif
  catch
    error("f has an unknown function.");
  end
  
  xAprox = xn(end);
  return;
endfunction


% ============================== Method 5 ====================================
%
% Danby Burkardt Method
%
% Metodos iterativos aplicados a la ecuación de Kepler. Page 120.
%
% Arguments:
%   f {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   tol {float, int} - tolerance that indicates the stop condition
%   graf {int} - flag that indicates if a graf must be done
%
% Returns:
%   xAprox {float} - root approximation
%   iter {int} - amount of iterations required
%
function [xAprox, iter] = sne_ud_5(f, x0, tol, graf=1)
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

      yk = double(fx(xk));
      df = double(diff(fx)(xk));
      df2 = double(diff(diff(fx))(xk));
      df3 = double(diff(diff(diff(fx)))(xk));
      
      Lf = yk * df2 / df^2;
      Ldf = df * df3 / df2^2;

      xk_next = xk - 1.5 * ((2 - Lf)**2 / (6 - 9 * Lf + 3 * Lf**2 + Lf**2 * Ldf)) * (yk / df);

      xn = [xn xk_next];
      error = [error abs(double(fx(xk_next)))];

      iter += 1;
    endwhile

    if (graf == 1)
      k = 0:1:iter;
      plotFunction(k, error, "Danby Burkardt Method");
    endif
  catch
    error("f has an unknown function.");
  end
  
  xAprox = xn(end);
  return;
endfunction


% ============================== Method 6 ====================================
%



% =========================== Auxiliaries functions =============================

% This function is used to plot iterations vs error

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