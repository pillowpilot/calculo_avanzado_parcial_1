function [tOut, yOut] = rungeKutta(equations, range, deltaT, yInitialValues)
  
  iterations = ceil( (range(2)-range(1)) / deltaT ) + 1;
  
  tOut = zeros(1, iterations);
  yOut = zeros(length(yInitialValues), iterations);
  
  t = range(1);
  yValues = yInitialValues;
  
  tOut(1) = t;
  yOut(:, 1) = yValues;
  for iteration = 2:iterations
    
    k1 = calculateK1(equations, t, yValues);
    k2 = calculateK2(equations, t, yValues, deltaT, k1);
    k3 = calculateK3(equations, t, yValues, deltaT, k2);
    k4 = calculateK4(equations, t, yValues, deltaT, k3);
    
    yValues = yValues + (1/6)*(k1 + 2*k2 + 2*k3 + k4)*deltaT;
    t = t + deltaT;
    
    tOut(iteration) = t;
    yOut(:, iteration) = yValues;
  end
  
endfunction

function k = calculateK1(equations, t, yValues)
  k = equations(t, yValues);
endfunction

function k = calculateK2(equations, t, yValues, deltaT, k1)
  t = t + 0.5*deltaT;
  yValues = yValues + 0.5*k1*deltaT;
  k = equations(t, yValues);
endfunction

function k = calculateK3(equations, t, yValues, deltaT, k2)
  t = t + 0.5*deltaT;
  yValues = yValues + 0.5*k2*deltaT;
  k = equations(t, yValues);
endfunction

function k = calculateK4(equations, t, yValues, deltaT, k3)
  t = t + deltaT;
  yValues = yValues + k3*deltaT;
  k = equations(t, yValues);
endfunction