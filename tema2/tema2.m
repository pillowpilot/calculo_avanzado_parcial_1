w_c = 2 * pi; 
m = 1; 
K = m * w_c^2;
initialValues = [0 .001];

# Z = 0 and w -> w_c
Z = 0;

figure();

id = 1;
for alpha = [.01, .05, .1, .3, .5, .7, .9, .99, .999]
  w = alpha * w_c;
  F = @(t) ( 10^(-2) * sin(w*t) );
  fun = @(t, y) [y(2); (-(K/m)*y(1) -(Z/m)*y(2) +(F(t)/m))];
  
  [t, y] = ode45(fun, [0 50], initialValues);
  
  subplot(3, 3, id);
  plot(t, y);
  title(['w = ' num2str(w) ', \alpha = ' num2str(alpha)]);
  
  id = id + 1;
endfor

figure();

alpha = 0.999999;
w = alpha * w_c;
F = @(t) ( 10^(-2) * sin(w*t) );
fun = @(t, y) [y(2); (-(K/m)*y(1) -(Z/m)*y(2) +(F(t)/m))];
[t, y] = ode45(fun, [0 1000], initialValues);
plot(t, y);

# Z != 0 and w = (1-1e-6)w_c 

figure()

id = 1;
for Z = [.01 .05 .1 .3 .5 .7 1 2 3]
  alpha = 1;
  w = alpha * w_c;
  F = @(t) ( 10^(-2) * sin(w*t) );
  fun = @(t, y) [y(2); (-(K/m)*y(1) -(Z/m)*y(2) +(F(t)/m))];
  [t, y] = ode45(fun, [0 40], initialValues);
  
  subplot(3, 3, id);
  plot(t, y);
  title(['Z = ' num2str(Z)]);
  
  id = id + 1;
endfor
