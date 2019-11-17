%[theta1-89��]ʱ����c
i=theta1;
while(i<90)
    t=i/180*pi;
    c=getC(t,H,A,B,vp);%��ˮ�߽ؾ�
    d=tan(t)^2+4*A*c;
    xL=(tan(t)-sqrt(d))/(2*A);%�󽻵�x
    xM=(H-c)/tan(t);%�м��x�����ҽ���x
    xR=sqrt(H/A);%�Ҽ��޵�x
    cobx=getX(c,xL,xM,xR,t,H,A,B);%����x
    cobz=getZ(c,xL,xM,xR,t,H,A,B);%����z
    cobc=cobz+cobx/tan(theta);%���������߽ؾ�
    L=getL(cobx,cobz,comz,t);%����
    Ms=L*vp*1000*9.8;%��ԭ����
    COMZ0=[COMZ,comz];COMZ=COMZ0;%����
    COBC0=[COBC,cobc];COBC=COBC0;%����
    MS0=[MS,Ms];
    MS=MS0;
    T0=[T,i];
    T=T0;
    Zero0=[Zero,0];
    Zero=Zero0;
    fprintf("%d��:c is %f,cobx is %f,cobz is %f,cobc is %f,L is %f,Ms is %f\n",i,c,cobx,cobz,cobc,L,Ms);
    i=i+1;
end

function c=getC(theta,H,A,B,vp)%ʹ�ýǶȵõ�c��ֵ
    i=1;
    cmin=-tan(theta)^2/(4*A);%c����
    cmax=H;cof=4/(3*sqrt(B));%c����
    c=double((cmin+cmax)/2);%���ַ�������ֵ
    while(i<100)
        d=tan(theta)^2+4*A*c;
        xL=(tan(theta)-sqrt(d))/(2*A);%�󽻵�x
        xM=(H-c)/tan(theta);%�м��x�����ҽ���
        xR=sqrt(H/A);%�Ҽ��޵�x
        fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        v1=quadv(fun1,xL,xM);v1=v1*cof;
        fun2=@(x) (H-A*x.^2).^1.5;
        v2=quadv(fun2,xM,xR);v2=v2*cof;
        volumn=v1+v2;
    if(volumn>vp)%�����̫����
       if(c<cmax)
           cmax=c;
       end
    end
    if(volumn<vp)%���С��������
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
function cobx=getX(c,xL,xM,xR,theta,H,A,B)%����x������
    cof=4/(3*sqrt(B));
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v1=quadv(fun1,xL,xM);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;
    v2=quadv(fun2,xM,xR);v2=v2*cof;
    
    Fun1=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    V1=quadv(Fun1,xL,xM);V1=V1*cof;
    Fun2=@(x) x.*(H-A*x.^2).^1.5;
    V2=quadv(Fun2,xM,xR);V2=V2*cof;
    cobx=(V1+V2)/(v1+v2);%Сv�����������v�����x����
end
function cobz=getZ(c,xL,xM,xR,theta,H,A,B)%����z������
    COF=1/sqrt(B);cof=4/(3*sqrt(B));
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v1=quadv(fun1,xL,xM);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;
    v2=quadv(fun2,xM,xR);v2=v2*cof;
    %fprintf("vp is %f,",v1+v2);
    Fun1=@(x) ((tan(theta)*x+c).^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(tan(theta)*x+c-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    V1=quadv(Fun1,xL,xM);V1=V1*COF;
    Fun2=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
    V2=quadv(Fun2,xM,xR);V2=V2*COF;
    %fprintf("v1 is %f,v2 is %f,V1 is %f,V2 is %f",v1,v2,V1,V2);
    cobz=(V1+V2)/(v1+v2);%Сv�����������v�����x����
end
function L=getL(cobx,cobz,comz,theta)%�����ġ����ġ������ԭ����
    L=-(tan(theta)*(comz-cobz)-cobx)/sqrt(1+tan(theta)^2);
end