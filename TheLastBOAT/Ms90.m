%90�㸴ԭ����
c=get90c(H,A,B,vp);
xmin=c;
xmax=sqrt(H/A);
COF=1/sqrt(B);cof=4/(3*sqrt(B));
funx1=@(x) x.*(H-A*x.^2).^1.5;%x����
vx1=quadv(funx1,xmin,xmax);vx1=vx1*cof;
funx2=@(x) (H-A*x.^2).^1.5;%���з�ĸ���������
v2=quadv(funx2,xmin,xmax);v2=v2*cof;
cobx=vx1/v2;%����x

funz1=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
vz1=quadv(funz1,xmin,xmax);vz1=vz1*COF;
cobz=vz1/v2;%����z
cobc=cobz;%���������߽ؾ�
L=cobz-comz;%����
Ms=L*vp*1000*9.8;%��ԭ����
i=90;
COMZ0=[COMZ,comz];COMZ=COMZ0;%����
COBC0=[COBC,cobc];COBC=COBC0;%����
MS0=[MS,Ms];
    MS=MS0;
    T0=[T,i];
    T=T0;
    Zero0=[Zero,0];
    Zero=Zero0;
fprintf("%d��:c is %f,cobx is %f,cobz is %f,cobc is %f,L is %f,Ms is %f\n",i,c,cobx,cobz,cobc,L,Ms);
function c=get90c(H,A,B,vp)
    cof=4/(3*sqrt(B));
    i=1;
    cmin=0;
    cmax=sqrt(H/A);
    c=double((cmin+cmax)/2);
    while(i<=100)
        if(mod(i,100)==0)
        %    fprintf("%dth: c is %f\n",i,c);%�������
        end
        xmin=c;
        xmax=sqrt(H/A);
        fun=@(x) (H-A*x.^2).^1.5;
        volumn=quadv(fun,xmin,xmax);
        volumn=volumn*cof;
    if(volumn<vp)%�����̫����
       if(c<cmax)
           cmax=c;
       end
    end
    if(volumn>vp)%���С��������
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