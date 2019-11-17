%90度之前
T=i;
t=-sqrt(H/A):0.01*sqrt(H/A):sqrt(H/A);
x=t;
z=A*x.^2;zh=0*t+H;
y=tan((T)/180*pi)*x+c;%吃水线
yb=tan((T-90)/180*pi)*(x-real(cobx))+real(cobz);%浮力作用线
ym=tan((T-90)/180*pi)*(x-real(comx))+real(comz);%重力作用线
plot(x,z,'r');hold on;grid on;%抛物线
plot(t,zh,'r');hold on;%顶
plot(x,y,'b');hold on;%吃水线
plot(x,yb,'g');hold on;%浮力作用线
plot(x,ym,'y');hold on;%重力作用线
text(real(comx),real(comz),'o','color','k');hold on;%标注质心
text(real(cobx),real(cobz),'o','color','k');hold on;%标注浮心
axis equal;