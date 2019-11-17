%���91-theta2��ʱ����ԭ���صļ���
cof=2/sqrt(B);
theta2=i;
i=91;
while(i<theta2)
    theta=i/180*pi;
    fprintf("%d��:",i);
    c=getC(theta,H,A,B,vp);%��ˮ�߽ؾ�
    d=tan(theta)^2+4*A*c;
    xmin=-sqrt(H/A);%���޵�x
    xL=(H-c)/tan(theta);%�󽻵�x
    xR=(tan(theta)+sqrt(d))/(2*A);%�ҽ���x
    xmax=sqrt(H/A);%�Ҽ��޵�x
    cobx=getX(c,xmin,xL,xR,xmax,theta,H,A,B);%����x
    cobz=getZ(c,xmin,xL,xR,xmax,theta,H,A,B);%����z
    cobc=cobz+cobx/tan(theta);%���������߽ؾ�
    L=getL(cobx,cobz,comz,theta);%����
    Ms=L*vp*1000*9.8;%��ԭ����
    COMZ0=[COMZ,comz];COMZ=COMZ0;%����
    COBC0=[COBC,cobc];COBC=COBC0;%����
    MS0=[MS,Ms];
    MS=MS0;
    T0=[T,i];
    T=T0;
    Zero0=[Zero,0];
    Zero=Zero0;
    fprintf("c is %f,cobx is %f ,cobz is %f,cobc is %f,L is %f ,Ms is %f\n",c,cobx,cobz,cobc,L,Ms);
    i=i+1;
end
function c=getC(theta,H,A,B,vp)%�̶��Ƕ���ؾ�
    cof=4/(3*sqrt(B));
    i=1;
    cmin=H+tan(theta)*sqrt(H/A);%���Ͻ���
    cmax=H-tan(theta)*sqrt(H/A);%���Ͻ���
    c=double((cmin+cmax)/2);
    while(i<100)
        d=tan(theta)^2+4*A*c;
        xmin=-sqrt(H/A);
        xL=(H-c)/tan(theta);
        xR=(tan(theta)+sqrt(d))/(2*A);
        xmax=sqrt(H/A);
        fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        v1=quadv(fun1,xL,xR);v1=v1*cof;%ˮ����
        fun2=@(x) (H-A*x.^2).^1.5;
        v2=quadv(fun2,xmin,xL);v2=v2*cof;%ˮ����
        v=quadv(fun2,xmin,xmax);v=v*cof;%�����
        volumn=v-v1-v2;
    if(volumn<vp)%���С��̫����
       if(c<cmax)
           cmax=c;
       end
    end
    if(volumn>vp)%�����̫����
        if(c>cmin)
            cmin=c;
        end
    end
    if(volumn==vp)
        break;
    end
        c=double((cmin+cmax)/2);
        i=i+1;
    end
end
function cobx=getX(c,xmin,xL,xR,xmax,theta,H,A,B)%����x������
    cof=4/(3*sqrt(B));
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v1=quadv(fun1,xL,xR);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;
    v2=quadv(fun2,xmin,xL);v2=v2*cof;
    v=quadv(fun2,xmin,xmax);v=v*cof;
    volumn=v-v1-v2;
    %СvΪ�������V�ǹ���x���������
    Fun1=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    V1=quadv(Fun1,xL,xR);V1=V1*cof;
    Fun2=@(x) x.*(H-A*x.^2).^1.5;
    V2=quadv(Fun2,xmin,xL);V2=V2*cof;
    %fprintf("volumn is %f, Vx1 is %f,",volumn,V1);%�������
    cobx=-(V1+V2)/volumn;
end
function cobz=getZ(c,xmin,xL,xR,xmax,theta,H,A,B)%����z������
    cof=4/(3*sqrt(B));Cof=1/sqrt(B);
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;%�ұ�
    v1=quadv(fun1,xL,xR);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;%���
    v2=quadv(fun2,xmin,xL);v2=v2*cof;
    v=quadv(fun2,xmin,xmax);v=v*cof;
    volumn=v-v1-v2;
    %СvΪ�������V�ǹ���z���������
    Fun1=@(x) ((tan(theta)*x+c).^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    V1=quadv(Fun1,xL,xR);V1=V1*Cof;
    Fun2=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x.^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
    V2=quadv(Fun2,xmin,xL);V2=V2*Cof;
    V=quadv(Fun2,xmin,xmax);V=V*Cof;
    Vz=V-V1-V2;
    cobz=Vz/volumn;
end
function L=getL(cobx,cobz,comz,theta)%�����ġ����ġ������ԭ����
    %L=(tan(theta)*(comz-cobz)-cobx)/sqrt(1+tan(theta)^2);
    %L=((cobz-comz)-tan(theta-pi/2)*cobx)*sqrt(1+tan(theta)^2);
    %L=(tan(theta)*(comz-cobz)+0-cobx)/sqrt(1+tan(theta)^2);%���ǲ��Թ�ʽ
    comx=0;L=(cobz-comz)*sin(theta)+(cobx-comx)*cos(theta);
end