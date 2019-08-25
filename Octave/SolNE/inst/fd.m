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


[xn, itera, graph] = sne_fd_1(x^2-exp(x), 0.5, 0.000001)