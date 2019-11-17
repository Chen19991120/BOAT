%先得运行MS0_FT到30度，即theta1改到31
x=-sqrt(H/A):0.01*sqrt(H/A):sqrt(H/A);
z=A*x.^2;zh=0*x+H;
plot(x,z);hold on;plot(x,zh);hold on;
d=tan(pi/6)^2+4*A*c;
tmin=(tan(pi/6)-sqrt(d))/(2*A);
tmax=(tan(pi/6)+sqrt(d))/(2*A);
t=tmin:0.01*(tmax-tmin):tmax;
z=A*t.^2;
zl=tan(pi/6)*t+c;
shadedplot(t,z,zl,'r');
