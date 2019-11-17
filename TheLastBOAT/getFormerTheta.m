%求出0-theta时最大的theta使得，z=tan(theta)(x-sqrt(H/A))+H与曲面围成的体积大于v排
%已修正A、B
%clear all;
i=1;
cof=4/(3*sqrt(B));%函数的常系数
while(i<90)
    theta=i/180*pi;%theta从0到90度递增，求出theta的变化点
    %fprintf("%dth calculate:",i);%原本打印用来测试的
    c=H-tan(theta)*sqrt(H/A);
    d=tan(theta)^2+4*A*c;%求根公式判别式b^2-4Ac
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;%匿名函数
    xmin=(tan(theta)-sqrt(d))/(2*A);
    xmax=sqrt(H/A);
    volumn=quadv(fun,xmin,xmax);%数值积分
    volumn=volumn*cof;
    %fprintf("volumn is %f\n",volumn);%原本打印用来测试的
    if(volumn<vp)%随theta变大，volumn变小，最后比所需v排小的时候退出，求得theta1
        break;
    end
    i=i+1;
end
