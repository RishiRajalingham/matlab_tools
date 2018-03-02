function h = rr_polarPcolor(rho, theta, X)


XX = X;
XX(size(X,1) + 1, size(X,2) + 1) = nan;

% transform data in polar coordinates to Cartesian coordinates.
x = (rho)'*cos(theta);
y = (rho)'*sin(theta);

% plot data on top of grid
h = pcolor(x,y,XX);

end