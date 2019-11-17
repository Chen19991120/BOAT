%求出91-theta2°时，复原力矩的计算
cof=2/sqrt(B);
theta2=i;
i=91;
while(i<theta2)
    theta=i/180*pi;
    fprintf("%d°:",i);
    c=getC(theta,H,A,B,vp);%吃水线截距
    d=tan(theta)^2+4*A*c;
    xmin=-sqrt(H/A);%左极限点x
    xL=(H-c)/tan(theta);%左交点x
    xR=(tan(theta)+sqrt(d))/(2*A);%右交点x
    xmax=sqrt(H/A);%右极限点x
    cobx=getX(c,xmin,xL,xR,xmax,theta,H,A,B);%浮心x
    cobz=getZ(c,xmin,xL,xR,xmax,theta,H,A,B);%浮心z
    cobc=cobz+cobx/tan(theta);%浮心作用线截距
    L=getL(cobx,cobz,comz,theta);%力臂
    Ms=L*vp*1000*9.8;%复原力矩
    COMZ0=[COMZ,comz];COMZ=COMZ0;%重心
    COBC0=[COBC,cobc];COBC=COBC0;%稳心
    MS0=[MS,Ms];
    MS=MS0;
    T0=[T,i];
    T=T0;
    Zero0=[Zero,0];
    Zero=Zero0;
    fprintf("c is %f,cobx is %f ,cobz is %f,cobc is %f,L is %f ,Ms is %f\n",c,cobx,cobz,cobc,L,Ms);
    i=i+1;
end
function c=getC(theta,H,A,B,vp)%固定角度求截距
    cof=4/(3*sqrt(B));
    i=1;
    cmin=H+tan(theta)*sqrt(H/A);%左上交点
    cmax=H-tan(theta)*sqrt(H/A);%右上交点
    c=double((cmin+cmax)/2);
    while(i<100)
        d=tan(theta)^2+4*A*c;
        xmin=-sqrt(H/A);
        xL=(H-c)/tan(theta);
        xR=(tan(theta)+sqrt(d))/(2*A);
        xmax=sqrt(H/A);
        fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        v1=quadv(fun1,xL,xR);v1=v1*cof;%水上右
        fun2=@(x) (H-A*x.^2).^1.5;
        v2=quadv(fun2,xmin,xL);v2=v2*cof;%水上左
        v=quadv(fun2,xmin,xmax);v=v*cof;%总体积
        volumn=v-v1-v2;
    if(volumn<vp)%体积小，太高了
       if(c<cmax)
           cmax=c;
       end
    end
    if(volumn>vp)%体积大，太低了
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
function cobx=getX(c,xmin,xL,xR,xmax,theta,H,A,B)%浮心x轴坐标
    cof=4/(3*sqrt(B));
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v1=quadv(fun1,xL,xR);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;
    v2=quadv(fun2,xmin,xL);v2=v2*cof;
    v=quadv(fun2,xmin,xmax);v=v*cof;
    volumn=v-v1-v2;
    %小v为体积，大V是关于x的体积积分
    Fun1=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    V1=quadv(Fun1,xL,xR);V1=V1*cof;
    Fun2=@(x) x.*(H-A*x.^2).^1.5;
    V2=quadv(Fun2,xmin,xL);V2=V2*cof;
    %fprintf("volumn is %f, Vx1 is %f,",volumn,V1);%测试输出
    cobx=-(V1+V2)/volumn;
end
function cobz=getZ(c,xmin,xL,xR,xmax,theta,H,A,B)%浮心z轴坐标
    cof=4/(3*sqrt(B));Cof=1/sqrt(B);
    fun1=@(x) (tan(theta)*x+c-A*x.^2).^1.5;%右边
    v1=quadv(fun1,xL,xR);v1=v1*cof;
    fun2=@(x) (H-A*x.^2).^1.5;%左边
    v2=quadv(fun2,xmin,xL);v2=v2*cof;
    v=quadv(fun2,xmin,xmax);v=v*cof;
    volumn=v-v1-v2;
    %小v为体积，大V是关于z的体积积分
    Fun1=@(x) ((tan(theta)*x+c).^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    V1=quadv(Fun1,xL,xR);V1=V1*Cof;
    Fun2=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x.^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
    V2=quadv(Fun2,xmin,xL);V2=V2*Cof;
    V=quadv(Fun2,xmin,xmax);V=V*Cof;
    Vz=V-V1-V2;
    cobz=Vz/volumn;
end
function L=getL(cobx,cobz,comz,theta)%由重心、浮心、倾角求复原力臂
    %L=(tan(theta)*(comz-cobz)-cobx)/sqrt(1+tan(theta)^2);
    %L=((cobz-comz)-tan(theta-pi/2)*cobx)*sqrt(1+tan(theta)^2);
    %L=(tan(theta)*(comz-cobz)+0-cobx)/sqrt(1+tan(theta)^2);%都是测试公式
    comx=0;L=(cobz-comz)*sin(theta)+(cobx-comx)*cos(theta);
end