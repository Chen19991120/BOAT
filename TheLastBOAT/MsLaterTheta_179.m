%求出118°-179时，复原力矩的计算
%A、B、comz已修改
cof=2/sqrt(B);
i=theta2;
while(i<180)
    theta=i/180*pi;
    fprintf("%d°:",i);
    c=getC(theta,H,A,B,vp);%吃水线截距
    d=tan(theta)^2+4*A*c;
    xL=(tan(theta)-sqrt(d))/(2*A);%左交点x
    xR=(tan(theta)+sqrt(d))/(2*A);%右交点x
    xmin=-sqrt(H/A);xmax=sqrt(H/A);%左右极限值点x
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
    cmin=0;
    cmax=H+tan(theta)*sqrt(H/A);%左上交点
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
    %fprintf("c is %f,volumn is %f,",c,volumn);
end
function cobx=getX(c,xmin,xL,xR,xmax,theta,H,A,B)%浮心x轴坐标
    cof=4/(3*sqrt(B));
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(fun,xL,xR);v=v*cof;
    Fun=@(x) (H-A*x.^2).^1.5;
    V=quadv(Fun,xmin,xmax);V=V*cof;
    volumn=V-v;
    %小v为体积，大V是关于x的体积积分
    funx=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(funx,xL,xR);v=v*cof;
    Funx=@(x) x.*(H-A*x.^2).^1.5;
    V=quadv(Funx,xmin,xmax);V=V*cof;
    vx=V-v;
    %fprintf("volumn is %f,Vx is %f,",volumn,vx);%测试输出
    cobx=(vx)/volumn;
end
function cobz=getZ(c,xmin,xL,xR,xmax,theta,H,A,B)%浮心z轴坐标
    cof=4/(3*sqrt(B));Cof=1/sqrt(B);
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;%水上
    v=quadv(fun,xL,xR);v=v*cof;
    Fun=@(x) (H-A*x.^2).^1.5;
    V=quadv(Fun,xmin,xmax);V=V*cof;
    volumn=V-v;
    %小v为体积，大V是关于x的体积积分
    funz=@(x) ((tan(theta)*x+c)^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(tan(theta)*x+c-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    v=quadv(funz,xL,xR);v=v*Cof;
    Funz=@(x) (H^2-A^2*x.^4).*sqrt(H-A*x.^2)-2*A*x.^2.*(H-A*x.^2).^1.5/3-(H-A*x.^2).^2.5/5;
    V=quadv(Funz,xmin,xmax);V=V*Cof;
    vz=V-v;
    %fprintf("volumn is %f,vz is %f,",volumn,vz);%测试输出
    cobz=(vz)/volumn;
end
function L=getL(cobx,cobz,comz,theta)%由重心、浮心、倾角求复原力臂
    comx=0;L=(cobz-comz)*sin(theta)+(cobx-comx)*cos(theta);
    %L=(tan(theta)*(comz-cobz)+(0-cobx))/sqrt(1+tan(theta)^2);%上下等效，已测试
end