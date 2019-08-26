pkg load symbolic
format long

syms x;                               % x (is a global symbol)
global ITER_LIMIT;                    % Limit of iterations
ITER_LIMIT = 10000;


% ============================== Method 1 ====================================
function [xn, itera, graph] = sne_fd_1(expr, x0, tol)
  %{
  Steffensen's Method
  Arguments:
      expr {symbolic} -- polynomial whose solution must be found
      x0 {double} -- initial value to start iterations
      tol {double} -- tolerance that indicates the stop condition

  Returns:
      xn {double} -- root approximation
      itera {int} -- amount of iterations required
      graph {int} -- flag that indicates if a graph must be done
  %}

  % -------------------------- Local variables -------------------------------
  global ITER_LIMIT
  itera = 0;                          % Amount of iterations
  graph = 0;                          % Wether the graph will be shown
  f = function_handle(expr);          % Use expression as function
  xn = x0;                            % Initialize x_n
  xNext = 0;                          % Initialize x_(n+1)
  error = abs(f(xn));                 % First errror calculation
  
  % ------------------------- Steffensen's Method ---------------------------
  try
    while (error > tol)
      if (itera > ITER_LIMIT)
        warning("Iteration limit reached");
        return;
      endif
      div = f(xn + f(xn)) - f(xn);
      if (div ~= 0)
        xNext = xn - (f(xn)^2) / div;
      else
        warning("[Math error]: Division by zero");
        return;
      end
      xn = xNext;
      error = abs(f(xn));
      itera++;
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
% Metodos iterativos aplicados a la ecuación de Kepler. Page 62.
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