fun = @(t, y)...
  [-2*(y(1)-1)*(y(2)-2)^2+(y(1)-1)*(y(2)-2)-(y(2)-2)-3*(y(1)-1);
  -1*(y(1)-1)+(y(1)-1)^2*(y(2)-2)];

jacobian = @(t, y)...
  [-2*(y(2)-2)^2+(y(2)-2)-3, -4*(y(1)-1)*(y(2)-2)+(y(1)-1)-1; ...
  -1+2*(y(2)-2)*(y(1)-1), (y(1)-1)^2];
  
stable_point = [1; 2];

stable_point_matrix = jacobian(0, stable_point); # [-3 -1; -1 0]
[lambda1 lambda2] = eig(stable_point_matrix);    # [-3.302 0.302] -> unstable

# Como uno de los auto valores es positivo, el punto [1 2] es inestable.
# Pequeñas oscilaciones harán que el sistema tienda al infinito (se aleje del 
# punto de estabilidad).

deltaT = 0.1; # 1/max{eigenvalues} is approx .3
range = [0 2]; 

# Change center and increase radii to explore other regions
center = stable_point;
radii = .5;

# Change total_points to control plot density
total_points = 400;

for iteration=1:total_points
  initialPoint = center + radii*rand(size(stable_point)) - (radii/2);
  [t y] = adamsBashforth(fun, range, deltaT, initialPoint);
  hold on;
  plot(y(1, :), y(2, :), '-', 'Color', 'k');
  plot(initialPoint(1), initialPoint(2), 'o', 'Color', 'k');
endfor

yL = get(gca,'YLim');
xL = get(gca,'XLim');
line(xL, [2 2], 'Color','r');
line([1 1], yL,'Color','r');

hold off;