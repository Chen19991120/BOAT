%�ȵ�����MSFT_179��175�ȣ���theta1�ĵ�176
theta=175;
x=-sqrt(H/A):0.01*sqrt(H/A):sqrt(H/A);
z=A*x.^2;zh=0*x+H;
plot(x,z);hold on;plot(x,zh);hold on;%��ʼ����
d=tan(theta*pi/180)^2+4*A*c;
tmin=(tan(theta*pi/180)-sqrt(d))/(2*A);
tmax=(tan(theta*pi/180)+sqrt(d))/(2*A);
tL=-sqrt(H/A);tR=-tL;
%t=tmin:0.01*(tmax-tmin):tmax;
t=tL:0.01*(tmin-tL):tmin;
z=A*t.^2;%����
zl=tan(theta*pi/180)*t+c;%��ˮ��
ch=H+0*t;%�װ�
shadedplot(t,ch,z,'b');hold on;%����
M=tmin:0.01*(tmax-tmin):tmax;
hM=0*M+H;%�װ�
lM=tan(theta*pi/180)*M+c;%��ˮ��
shadedplot(M,hM,lM,'g');hold on;%�м�
R=tmax:0.01*(tR-tmax):tR;
zr=A*R.^2;%��ˮ��
cr=H+0*R;%�װ�
shadedplot(R,cr,zr,'r');hold on;%�Ұ��
