%先得运行MSFT_179到175度，即theta1改到176
theta=175;
x=-sqrt(H/A):0.01*sqrt(H/A):sqrt(H/A);
z=A*x.^2;zh=0*x+H;
plot(x,z);hold on;plot(x,zh);hold on;%初始界限
d=tan(theta*pi/180)^2+4*A*c;
tmin=(tan(theta*pi/180)-sqrt(d))/(2*A);
tmax=(tan(theta*pi/180)+sqrt(d))/(2*A);
tL=-sqrt(H/A);tR=-tL;
%t=tmin:0.01*(tmax-tmin):tmax;
t=tL:0.01*(tmin-tL):tmin;
z=A*t.^2;%曲面
zl=tan(theta*pi/180)*t+c;%吃水线
ch=H+0*t;%甲板
shadedplot(t,ch,z,'b');hold on;%左半边
M=tmin:0.01*(tmax-tmin):tmax;
hM=0*M+H;%甲板
lM=tan(theta*pi/180)*M+c;%吃水线
shadedplot(M,hM,lM,'g');hold on;%中间
R=tmax:0.01*(tR-tmax):tR;
zr=A*R.^2;%吃水线
cr=H+0*R;%甲板
shadedplot(R,cr,zr,'r');hold on;%右半边
