%�����90-theta��ʱ����thetaʹ�ã�z=tan(theta)(x+sqrt(H/A))+H������Χ�ɵ����С��v��
%��ʱv��<m/��ˮ
%������A��B
i=91;
cof=4/(3*sqrt(B));%��ϵ��
while(i<180)
    theta=i/180*pi;
    %fprintf("%dth calculate:",i);%�������
    c=H+tan(theta)*sqrt(H/A);
    d=tan(theta)^2+4*A*c;
    %fprintf("d is %f,",d);%�������
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;
    xmin=-sqrt(H/A);
    xmax=(tan(theta)+sqrt(d))/(2*A);
    volumn=quadv(fun,xmin,xmax);
    volumn=volumn*cof;
    Fun=@(x) (H-A*x.^2).^1.5;
    Cof=8/(3*sqrt(B));%��ϵ��
    V=quadv(Fun,0,-xmin);
    V=V*Cof;
    %fprintf("V is %f,vˮ�� is %f,",V,volumn);%�������
    volumn=V-volumn;
    %fprintf("volumn is %f\n",volumn);%�������
    if(volumn<vp)
        break;
    end
    i=i+1;
end
