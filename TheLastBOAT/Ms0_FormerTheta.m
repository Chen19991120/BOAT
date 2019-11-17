%���0-theta1ʱ����ԭ���صļ���
theta1=i;
i=1;
while(i<theta1)
    theta=i/180*pi;
    fprintf("%d��:",i);
    c=getC(theta,H,A,B,vp);%��ˮ�߽ؾ�
    d=tan(theta)^2+4*A*c;
    xmin=(tan(theta)-sqrt(d))/(2*A);%�󽻵�x
    xmax=(tan(theta)+sqrt(d))/(2*A);%�ҽ���x
    cobx=getX(c,xmin,xmax,theta,A,B);%����x
    cobz=getZ(c,xmin,xmax,theta,A,B);%����z
    cobc=cobz+cobx/tan(theta);%���������߽ؾ�
    L=getL(cobx,cobz,comz,theta);%����
    COMZ0=[COMZ,comz];COMZ=COMZ0;%����
    COBC0=[COBC,cobc];COBC=COBC0;%����
    Ms=L*vp*1000*9.8;%��ԭ����
    MS0=[MS,Ms];MS=MS0;%������ͼ
    T0=[T,i];T=T0;%�¶���ͼ
    Zero0=[Zero,0];Zero=Zero0;%����
    fprintf("c is %f,cobx is %f,cobz is %f,cobc is %f,L is %f,Ms is %f\n",c,cobx,cobz,cobc,L,Ms);
    i=i+1;
end
function c=getC(theta,H,A,B,vp)%�̶��Ƕ���ؾ�
    cof=4/(3*sqrt(B));%��ϵ��
    i=1;
    cmin=0;%�ؾ�����
    cmax=H-tan(theta)*sqrt(H/A);%�ؾ�����
    c=double((cmin+cmax)/2);%���ַ��������
    while(i<100) 
        d=tan(theta)^2+4*A*c;
        fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        xmin=(tan(theta)-sqrt(d))/(2*A);
        xmax=(tan(theta)+sqrt(d))/(2*A);
        volumn=quadv(fun,xmin,xmax);
        volumn=volumn*cof;
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
function cobx=getX(c,xmin,xmax,theta,A,B)%����x������
    cof=4/(3*sqrt(B));
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(fun,xmin,xmax);v=v*cof;
    Fun=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    V=quadv(Fun,xmin,xmax);V=V*cof;
    
    cobx=V/v;
end
function cobz=getZ(c,xmin,xmax,theta,A,B)%����z������
    cof=4/(3*sqrt(B));Cof=1/sqrt(B);
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(fun,xmin,xmax);v=v*cof;
    
    Fun=@(x) ((tan(theta)*x+c).^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(tan(theta)*x+c-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    V=quadv(Fun,xmin,xmax);V=V*Cof;
    cobz=V/v;
end
function L=getL(cobx,cobz,comz,theta)%�����ġ����ġ������ԭ����
    %L=-(tan(theta)*(comz-cobz)-cobx)/sqrt(1+tan(theta)^2);%equal to��zb-zm��sin��theta��+xbcos��theta��
    L=sin(theta)*(cobz-comz)+cos(theta)*(cobx);%�Ѳ�������������ʽ����һ��
end