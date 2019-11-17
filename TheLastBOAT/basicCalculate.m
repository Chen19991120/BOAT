%基础数据的运算，H、A、B、L、V总、浮心、吃水线
%comz=0.078175856225677;
clear all;
format long;
H=0.16;kA=1.65;A=4/(kA^2*H);kb=3;B=4/(kb^2*H);m=1.26;%my
comx=0;
COMZ=[];%重心
COBC=[];%稳心
MS=[];%力矩数据
T=[];%倾斜角度数数据
Zero=[];%水平线数据
%k=(sqrt(5)-1)/2;%长高宽的公比
%ka=宽/高，kb是长/高
%the above is my data
L=0.084;
mz=0.89247;zz=B*0.042^2;zz=zz+0.04;%半径与龙骨高
mg=0.1;zg=zz+0.02+0.25;%重物重心加0.02m半径+0.25m高
mc=m-mz-mg;%zc=(H^2.5/sqrt(B))*(1-1/5);%仅对龙骨那个面进行积分求质心
zc=2*H/3;
xmin=0;
cof=8/(3*sqrt(B));
fun=@(x)(H-A*x.^2).^1.5;
v=quadv(fun,0,sqrt(H/A));
v=v*cof;%总体积
comz=(mc*zc+mg*zg+mz*zz)/(mc+mg+mz);%质心z
%吃水线
vp=m/1000;i=1;
cmin=0;cmax=H;
c=double((cmin+cmax));
while(i<1000)
    func=@(x) (c-A*x.^2).^1.5;
    xmax=sqrt(c/A);
    volumn=quadv(func,0,xmax);volumn=volumn*cof;
    %fprintf("%dth:c is %f\n",i,c);
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
%已求出吃水线，可以求浮心啦
xmax=sqrt(c/A);
funb=@(x) (c^2-A^2*x.^4).*sqrt(c-A*x.^2)-2*A*x.^2.*(c-A*x.^2).^1.5/3-(c-A*x.^2).^2.5/5;
cobz=quadv(funb,0,xmax);cofb=2*1000/(m*sqrt(B));
cobz=cobz*cofb;
%comz=2*H/3;%my
%comz=85/1000;
fprintf("A is %f,B is %f,H is %f\n",A,B,H);
fprintf("Vall is %f,waterline is %f\n",v,c);
fprintf("质心M(0,0,%f)\n",comz);
fprintf("浮心B(0,0,%f)\n",cobz);