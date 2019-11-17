%求出（90-theta）时最大的theta使得，z=tan(theta)(x+sqrt(H/A))+H与曲面围成的体积小于v排
%此时v排<m/ρ水
%已修正A、B
i=91;
cof=4/(3*sqrt(B));%常系数
while(i<180)
    theta=i/180*pi;
    %fprintf("%dth calculate:",i);%测试输出
    c=H+tan(theta)*sqrt(H/A);
    d=tan(theta)^2+4*A*c;
    %fprintf("d is %f,",d);%测试输出
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    xmin=-sqrt(H/A);
    xmax=(tan(theta)+sqrt(d))/(2*A);
    volumn=quadv(fun,xmin,xmax);
    volumn=volumn*cof;
    Fun=@(x) (H-A*x.^2).^1.5;
    Cof=8/(3*sqrt(B));%常系数
    V=quadv(Fun,0,-xmin);
    V=V*Cof;
    %fprintf("V is %f,v水上 is %f,",V,volumn);%测试输出
    volumn=V-volumn;
    %fprintf("volumn is %f\n",volumn);%测试输出
    if(volumn<vp)
        break;
    end
    i=i+1;
end
