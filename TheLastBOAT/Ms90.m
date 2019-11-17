%90°复原力矩
c=get90c(H,A,B,vp);
xmin=c;
xmax=sqrt(H/A);
COF=1/sqrt(B);cof=4/(3*sqrt(B));
funx1=@(x) x.*(H-A*x.^2).^1.5;%x分子
vx1=quadv(funx1,xmin,xmax);vx1=vx1*cof;
funx2=@(x) (H-A*x.^2).^1.5;%所有分母，体积积分
v2=quadv(funx2,xmin,xmax);v2=v2*cof;
cobx=vx1/v2;%浮心x

funz1=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
vz1=quadv(funz1,xmin,xmax);vz1=vz1*COF;
cobz=vz1/v2;%浮心z
cobc=cobz;%浮心作用线截距
L=cobz-comz;%力臂
Ms=L*vp*1000*9.8;%复原力矩
i=90;
COMZ0=[COMZ,comz];COMZ=COMZ0;%重心
COBC0=[COBC,cobc];COBC=COBC0;%稳心
MS0=[MS,Ms];
    MS=MS0;
    T0=[T,i];
    T=T0;
    Zero0=[Zero,0];
    Zero=Zero0;
fprintf("%d°:c is %f,cobx is %f,cobz is %f,cobc is %f,L is %f,Ms is %f\n",i,c,cobx,cobz,cobc,L,Ms);
function c=get90c(H,A,B,vp)
    cof=4/(3*sqrt(B));
    i=1;
    cmin=0;
    cmax=sqrt(H/A);
    c=double((cmin+cmax)/2);
    while(i<=100)
        if(mod(i,100)==0)
        %    fprintf("%dth: c is %f\n",i,c);%测试输出
        end
        xmin=c;
        xmax=sqrt(H/A);
        fun=@(x) (H-A*x.^2).^1.5;
        volumn=quadv(fun,xmin,xmax);
        volumn=volumn*cof;
    if(volumn<vp)%体积大，太高了
       if(c<cmax)
           cmax=c;
       end
    end
    if(volumn>vp)%体积小，不够高
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