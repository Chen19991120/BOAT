%先得运行ms91_LT到135度
THETA=135;
x=-sqrt(H/A):0.01*sqrt(H/A):sqrt(H/A);
z=A*x.^2;zh=0*x+H;
plot(x,z);hold on;plot(x,zh);hold on;
d=tan(THETA/180*pi)^2+4*A*c;
tmax=(tan(THETA/180*pi)+sqrt(d))/(2*A);
tmin=(H-c)/tan(THETA/180*pi);
t=tmin:0.01*(tmax-tmin):tmax;
z=A*t.^2;
zl=tan(THETA/180*pi)*t+c;
zh=0*t+H;
shadedplot(t,zh,zl,'r');hold on;
o=tmax:0.01*(sqrt(H/A)-tmax):sqrt(H/A);
oc=A*o.^2;
oh=0*o+H;
shadedplot(o,oc,oh,'g');hold on;