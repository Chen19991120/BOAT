%90��֮ǰ
T=i;
t=-sqrt(H/A):0.01*sqrt(H/A):sqrt(H/A);
x=t;
z=A*x.^2;zh=0*t+H;
y=tan((T)/180*pi)*x+c;%��ˮ��
yb=tan((T-90)/180*pi)*(x-real(cobx))+real(cobz);%����������
ym=tan((T-90)/180*pi)*(x-real(comx))+real(comz);%����������
plot(x,z,'r');hold on;grid on;%������
plot(t,zh,'r');hold on;%��
plot(x,y,'b');hold on;%��ˮ��
plot(x,yb,'g');hold on;%����������
plot(x,ym,'y');hold on;%����������
text(real(comx),real(comz),'o','color','k');hold on;%��ע����
text(real(cobx),real(cobz),'o','color','k');hold on;%��ע����
axis equal;