function [tValues, yValues] = adamsBashforth...
  (fun, range, deltaT, initialValues)
  
  if size( initialValues )(2) ~= 1
      error('initialValues must be a column vector.\n');
  endif
  
  order = size(fun(0, initialValues))(1);
    
  iterations = ( range(2)-range(1) )/deltaT;
  iterations = ceil( iterations ) + 1;
    
  tValues = zeros(1, iterations);
  yValues = zeros(order, iterations);
    
  tValues(1) = range(1);
  tValues(2) = tValues(1) + deltaT;
  tValues(3) = tValues(2) + deltaT;
  
  # RK4 for first 3 points
  [f1 f2 f3] = calculateFirstPoints(fun, range, deltaT, initialValues);
  
  yValues(:, 1) = f1;
  yValues(:, 2) = f2;
  yValues(:, 3) = f3;
  
  t_n = tValues(3);
  y_n = yValues(:, 3);
  for iteration=3:iterations
    
    f_n   = fun(t_n, y_n);
    f_nm1 = fun(t_n - deltaT, yValues(:, iteration-1));
    f_nm2 = fun(t_n - 2*deltaT, yValues(:, iteration-2));
    y_np1 = y_n + ((23/12)*f_n -(16/12)*f_nm1 +(5/12)*f_nm2)*deltaT;
    
    tValues(iteration+1) = t_n + deltaT;
    yValues(:, iteration+1) = y_np1;
      
    t_n = t_n + deltaT;
    y_n = y_np1;
  endfor
    
endfunction

function [f1 f2 f3] = calculateFirstPoints(fun, range, deltaT, initialValues)
  newRangeInitialPoint = range(1);
  newRangeFinalPoint = range(1) + deltaT*3;
  newRange = [newRangeInitialPoint newRangeFinalPoint];
  
  [t y] = rungeKutta(fun, newRange, deltaT, initialValues);
  
  f1 = y(:, 1);
  f2 = y(:, 2);
  f3 = y(:, 3);
endfunction
