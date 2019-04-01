function [tValues, yValues] = CrankNicholsorForNonLinear...
    (fun, jacobian, range, deltaT, initialValues)
    
    if size( initialValues )(2) ~= 1
      error('initialValues must be a column vector.\n');
    endif
  
    order = size(fun(0, initialValues))(1);
    
    iterations = ( range(2)-range(1) )/deltaT;
    iterations = ceil( iterations ) + 1;
    
    tValues = zeros(1, iterations);
    yValues = zeros(order, iterations);
    
    tValues(1) = range(1);
    yValues(:, 1) = initialValues;
    
    t_n = tValues(1);
    y_n = yValues(:, 1);
    for iteration=2:iterations
      t_np1 = t_n + deltaT;
      g = @(t_np1, y_np1)[y_np1 - y_n - deltaT*fun(t_np1, y_np1)];
      dg = @(t_np1, y_np1)[eye(order) - deltaT*(...
      0.5*jacobian(t_n, y_n) + 0.5*jacobian(t_np1, y_np1)...
      )];
    
      y_np1 = newtonRaphson(g, dg, t_np1, y_n);
    
      tValues(iteration) = t_np1;
      yValues(:, iteration) = y_np1;
      
      t_n = t_np1;
      y_n = y_np1;
    endfor
  
endfunction
