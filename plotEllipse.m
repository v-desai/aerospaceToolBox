%% plot Ellipse
function plotEllipse(state,covariance)
a = 3*covariance(1,1)^.5;
b = 3*covariance(2,2)^.5;

% [V,D] = eig(covariance);
% if D(2,2) > D(1,1)
%     eigvec = V(:,2);
% else
%     eigvec = V(:,1);
% end
% angle = atan(eigvec(2)/eigvec(1));
angle = 0;

x0 = state(1,1);
y0 = state(2,1);

syms x y
% Equation of an ellipse rotated by an angle
% eqn = ((x-x0)*cos(angle)+(y-y0)*sin(angle))^2/a^2 + ((x-x0)*sin(angle)-(y-y0)*cos(angle))^2/b^2 == 1;
eqn = (x-x0)^2/a^2 + (y-y0)^2/b^2 == 1;

xmin = -1.1*a+x0;
xmax = x0+1.1*a;
ymin = -1.1*b+y0;
ymax = y0+1.1*b;

scatter(x0,y0)
hold on
fimplicit(eqn,[xmin,xmax,ymin,ymax])
title(' ')
xlabel('Radial [km]')
ylabel('In-track [km]')
grid on
% axis equal