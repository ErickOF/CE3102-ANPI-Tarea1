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
    graph = 1
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
    graph = 1
  catch
    warning("[Math error]: Problem found executing the method");
    return;
  end
endfunction


% ============================== Method 3 ====================================
%
% Chebyshev Method
%    
% Metodos iterativos aplicados a la ecuaci�n de Kepler. Page 62.
%
% Arguments:
%   f {string} - polynomial whose solution must be found
%   x0 {float, int} - initial value to start iterations
%   tol {float, int} - tolerance that indicates the stop condition
%   graph {int} - flag that indicates if a graph must be done
%
% Returns:
%   xn {float} - root approximation
%   _iter {int} - amount of iterations required
%
function [xAprox, iter] = sne_ud_3(f, x0, tol, graf=1)
  if (typeinfo(f) != "string")
    error('f must be a string');
  endif

  if (graf != 0 && graf != 1)
    error('graf must be 0 or 1');
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

      xn = [xn xk_next]
      error = [error abs(double(fx(xk_next)))]

      iter += 1;
        
      %if graf == 1:
      %  k = np.linspace(0, _iter, _iter + 1)
      %  plotFunction(k, error, 'Chebyshev Method')
    endwhile
  catch
    error("f has an unknown function.");
  end
  
  xAprox = xn(end);
  return;
endfunction