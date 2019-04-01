function root = newtonRaphson(fun, jacobian, t, initialGuess)

  if size( initialGuess )(2) ~= 1
    error('initialGuess must be a column vector.\n');
  endif
  
  if size( fun(t, initialGuess) )(2) ~= 1
    error('fun(t, initialGuess) should output a column vector.\n');
  endif
  
  root = initialGuess;
  do
    correction = jacobian(t, root) \ fun(t, root);
    root = root - correction;
  until norm(correction) < 1e-6
  
endfunction
