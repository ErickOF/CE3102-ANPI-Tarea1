pkg load symbolic
format long

syms x;                               % x (is a global symbol)
global ITER_LIMIT;                    % Limit of iterations
ITER_LIMIT = 99;


% ============================== Method 1 ====================================
function [xn, itera, graph] = sne_ud_1(expr, x0, tol)
  %{
  Halley's Method
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
  
  % --------------------------- Halley's Method ------------------------------
  try
    while (error > tol)
      if (itera > ITER_LIMIT)
        warning("Iteration limit reached");
        return;
      endif
      fPrime = function_handle(diff(expr));
      fPrime2 = function_handle(diff(diff(expr)));
      div = 2 * fPrime(xn)^2 - f(xn) * fPrime2(xn);
      if (div ~= 0)
        xNext = xn - (2* f(xn) * fPrime(xn)) / div;
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


% ============================== Method 2 ====================================
function [xn, itera, graph] = sne_ud_2(expr, x0, tol)
  %{
  Frontini's and Sormani's Method
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
  
  % ------------------- Frontini's and Sormani's Method----=------------------
  try
    while (error > tol)
      if (itera > ITER_LIMIT)
        warning("Iteration limit reached");
        return;
      endif
      fPrime = function_handle(diff(expr));
      div1 = fPrime(xn);
      if (div1 ~= 0)
        div2 = fPrime(xn - 0.5 * f(xn) / div1);
        if (div2 ~= 0)
          xNext = xn - f(xn)/div2;
        else
          warning("[Math error]: Division by zero");
        end 
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