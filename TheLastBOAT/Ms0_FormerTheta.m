%求出0-theta1时，复原力矩的计算
theta1=i;
i=1;
while(i<theta1)
    theta=i/180*pi;
    fprintf("%d°:",i);
    c=getC(theta,H,A,B,vp);%吃水线截距
    d=tan(theta)^2+4*A*c;
    xmin=(tan(theta)-sqrt(d))/(2*A);%左交点x
    xmax=(tan(theta)+sqrt(d))/(2*A);%右交点x
    cobx=getX(c,xmin,xmax,theta,A,B);%浮心x
    cobz=getZ(c,xmin,xmax,theta,A,B);%浮心z
    cobc=cobz+cobx/tan(theta);%浮心作用线截距
    L=getL(cobx,cobz,comz,theta);%力臂
    COMZ0=[COMZ,comz];COMZ=COMZ0;%重心
    COBC0=[COBC,cobc];COBC=COBC0;%稳心
    Ms=L*vp*1000*9.8;%复原力矩
    MS0=[MS,Ms];MS=MS0;%力矩作图
    T0=[T,i];T=T0;%温度作图
    Zero0=[Zero,0];Zero=Zero0;%横线
    fprintf("c is %f,cobx is %f,cobz is %f,cobc is %f,L is %f,Ms is %f\n",c,cobx,cobz,cobc,L,Ms);
    i=i+1;
end
function c=getC(theta,H,A,B,vp)%固定角度求截距
    cof=4/(3*sqrt(B));%常系数
    i=1;
    cmin=0;%截距下限
    cmax=H-tan(theta)*sqrt(H/A);%截距上限
    c=double((cmin+cmax)/2);%二分法迭代求解
    while(i<100) 
        d=tan(theta)^2+4*A*c;
        fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
        xmin=(tan(theta)-sqrt(d))/(2*A);
        xmax=(tan(theta)+sqrt(d))/(2*A);
        volumn=quadv(fun,xmin,xmax);
        volumn=volumn*cof;
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
function cobx=getX(c,xmin,xmax,theta,A,B)%浮心x轴坐标
    cof=4/(3*sqrt(B));
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(fun,xmin,xmax);v=v*cof;
    Fun=@(x) x.*(tan(theta)*x+c-A*x.^2).^1.5;
    V=quadv(Fun,xmin,xmax);V=V*cof;
    
    cobx=V/v;
end
function cobz=getZ(c,xmin,xmax,theta,A,B)%浮心z轴坐标
    cof=4/(3*sqrt(B));Cof=1/sqrt(B);
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    v=quadv(fun,xmin,xmax);v=v*cof;
    
    Fun=@(x) ((tan(theta)*x+c).^2-A^2*x.^4).*sqrt(tan(theta)*x+c-A*x.^2)-2*A*x.^2.*(tan(theta)*x+c-A*x.^2).^1.5/3-(tan(theta)*x+c-A*x.^2).^2.5/5;
    V=quadv(Fun,xmin,xmax);V=V*Cof;
    cobz=V/v;
end
function L=getL(cobx,cobz,comz,theta)%由重心、浮心、倾角求复原力臂
    %L=-(tan(theta)*(comz-cobz)-cobx)/sqrt(1+tan(theta)^2);%equal to（zb-zm）sin（theta）+xbcos（theta）
    L=sin(theta)*(cobz-comz)+cos(theta)*(cobx);%已测试上下两个公式数据一样
end