%���118��-179ʱ����ԭ���صļ���
%A��B��comz���޸�
cof=2/sqrt(B);
i=theta2;
while(i<180)
    theta=i/180*pi;
    fprintf("%d��:",i);
    c=getC(theta,H,A,B,vp);%��ˮ�߽ؾ�
    d=tan(theta)^2+4*A*c;
    xL=(tan(theta)-sqrt(d))/(2*A);%�󽻵�x
    xR=(tan(theta)+sqrt(d))/(2*A);%�ҽ���x
    xmin=-sqrt(H/A);xmax=sqrt(H/A);%���Ҽ���ֵ��x
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
    cmin=0;
    cmax=H+tan(theta)*sqrt(H/A);%���Ͻ���
    c=double((cmin+cmax)/2);
    while(i<100)
        d=tan(theta)^2+4*A*c;
        xL=(tan(theta)-sqrt(d))/(2*A);
        xR=(tan(theta)+sqrt(d))/(2*A);
        xmin=-sqrt(H/A);xmax=sqrt(H/A);
        fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        v=quadv(fun,xL,xR);v=v*cof;
        Fun=@(x) (H-A*x.^2).^1.5;
        V=quadv(Fun,xmin,xmax);V=V*cof;
        volumn=V-v;
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
    %fprintf("c is %f,volumn is %f,",c,volumn);
end
function cobx=getX(c,xmin,xL,xR,xmax,theta,H,A,B)%����x������
    cof=4/(3*sqrt(B));
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(fun,xL,xR);v=v*cof;
    Fun=@(x) (H-A*x.^2).^1.5;
    V=quadv(Fun,xmin,xmax);V=V*cof;
    volumn=V-v;
    %СvΪ�������V�ǹ���x���������
    funx=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(funx,xL,xR);v=v*cof;
    Funx=@(x) x.*(H-A*x.^2).^1.5;
    V=quadv(Funx,xmin,xmax);V=V*cof;
    vx=V-v;
    %fprintf("volumn is %f,Vx is %f,",volumn,vx);%�������
    cobx=(vx)/volumn;
end
function cobz=getZ(c,xmin,xL,xR,xmax,theta,H,A,B)%����z������
    cof=4/(3*sqrt(B));Cof=1/sqrt(B);
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;%ˮ��
    v=quadv(fun,xL,xR);v=v*cof;
    Fun=@(x) (H-A*x.^2).^1.5;
    V=quadv(Fun,xmin,xmax);V=V*cof;
    volumn=V-v;
    %СvΪ�������V�ǹ���x���������
    funz=@(x) ((tan(theta)*x+c)^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(tan(theta)*x+c-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    v=quadv(funz,xL,xR);v=v*Cof;
    Funz=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x.^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
    V=quadv(Funz,xmin,xmax);V=V*Cof;
    vz=V-v;
    %fprintf("volumn is %f,vz is %f,",volumn,vz);%�������
    cobz=(vz)/volumn;
end
function L=getL(cobx,cobz,comz,theta)%�����ġ����ġ������ԭ����
    comx=0;L=(cobz-comz)*sin(theta)+(cobx-comx)*cos(theta);
    %L=(tan(theta)*(comz-cobz)+(0-cobx))/sqrt(1+tan(theta)^2);%���µ�Ч���Ѳ���
end