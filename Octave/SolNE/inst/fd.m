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
    graph = 1
  catch
    warning("[Math error]: Problem found executing the method");
    return;
  end
endfunction