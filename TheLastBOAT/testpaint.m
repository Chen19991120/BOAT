
r = [0:0.02:1]';
t = 0:1:360;
x = r*sind(t)+5;
y = r*cosd(t);
z = (x.^2+y.^2)/4;
surf(x,y,z)