%[theta1-89°]时候求c
i=theta1;
while(i<90)
    t=i/180*pi;
    c=getC(t,H,A,B,vp);%吃水线截距
    d=tan(t)^2+4*A*c;
    xL=(tan(t)-sqrt(d))/(2*A);%左交点x
    xM=(H-c)/tan(t);%中间点x，即右交点x
    xR=sqrt(H/A);%右极限点x
    cobx=getX(c,xL,xM,xR,t,H,A,B);%浮心x
    cobz=getZ(c,xL,xM,xR,t,H,A,B);%浮心z
    cobc=cobz+cobx/tan(theta);%浮心作用线截距
    L=getL(cobx,cobz,comz,t);%力臂
    Ms=L*vp*1000*9.8;%复原力矩
    COMZ0=[COMZ,comz];COMZ=COMZ0;%重心
    COBC0=[COBC,cobc];COBC=COBC0;%稳心
    MS0=[MS,Ms];
    MS=MS0;
    T0=[T,i];
    T=T0;
    Zero0=[Zero,0];
    Zero=Zero0;
    fprintf("%d°:c is %f,cobx is %f,cobz is %f,cobc is %f,L is %f,Ms is %f\n",i,c,cobx,cobz,cobc,L,Ms);
    i=i+1;
end

function c=getC(theta,H,A,B,vp)%使用角度得到c的值
    i=1;
    cmin=-tan(theta)^2/(4*A);%c下限
    cmax=H;cof=4/(3*sqrt(B));%c上限
    c=double((cmin+cmax)/2);%二分法迭代求值
    while(i<100)
        d=tan(theta)^2+4*A*c;
        xL=(tan(theta)-sqrt(d))/(2*A);%左交点x
        xM=(H-c)/tan(theta);%中间点x，即右交点
        xR=sqrt(H/A);%右极限点x
        fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        v1=quadv(fun1,xL,xM);v1=v1*cof;
        fun2=@(x) (H-A*x.^2).^1.5;
        v2=quadv(fun2,xM,xR);v2=v2*cof;
        volumn=v1+v2;
    if(volumn>vp)%体积大，太高了
       if(c<cmax)
           cmax=c;
       end
    end
    if(volumn<vp)%体积小，不够高
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
function cobx=getX(c,xL,xM,xR,theta,H,A,B)%浮心x轴坐标
    cof=4/(3*sqrt(B));
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v1=quadv(fun1,xL,xM);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;
    v2=quadv(fun2,xM,xR);v2=v2*cof;
    
    Fun1=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    V1=quadv(Fun1,xL,xM);V1=V1*cof;
    Fun2=@(x) x.*(H-A*x.^2).^1.5;
    V2=quadv(Fun2,xM,xR);V2=V2*cof;
    cobx=(V1+V2)/(v1+v2);%小v代表体积，大v代表带x积分
end
function cobz=getZ(c,xL,xM,xR,theta,H,A,B)%浮心z轴坐标
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
    cobz=(V1+V2)/(v1+v2);%小v代表体积，大v代表带x积分
end
function L=getL(cobx,cobz,comz,theta)%由重心、浮心、倾角求复原力臂
    L=-(tan(theta)*(comz-cobz)-cobx)/sqrt(1+tan(theta)^2);
end