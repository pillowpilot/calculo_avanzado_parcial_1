u_max = 10;
K     = 10;
K_s   = 10;
k_d   = 1/10;
k_e   = 1/10;
k_h   = 1/10;

fun = @(t, y)[...
  u_max*(1-y(1)/K)*(y(3)/(K_s+y(3)))*y(1)-k_d*y(1)-k_e*y(1);...
  k_d*y(1)-k_h*y(2);...
  k_e*y(1)+k_h*y(2)-u_max*(1-y(1)/K)*(y(3)/(K_s+y(3)))*y(1)];
  
initialValues = [1; 100; 0];
range = [0 100];

[t, y] = ode45(fun, range, initialValues);

finalValues = y( size(y)(1), :)' # [9.7748; 9.7795; 81.4457]
plot(t, y);