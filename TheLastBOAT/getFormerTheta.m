%���0-thetaʱ����thetaʹ�ã�z=tan(theta)(x-sqrt(H/A))+H������Χ�ɵ��������v��
%������A��B
%clear all;
i=1;
cof=4/(3*sqrt(B));%�����ĳ�ϵ��
while(i<90)
    theta=i/180*pi;%theta��0��90�ȵ��������theta�ı仯��
    %fprintf("%dth calculate:",i);%ԭ����ӡ�������Ե�
    c=H-tan(theta)*sqrt(H/A);
    d=tan(theta)^2+4*A*c;%�����ʽ�б�ʽb^2-4Ac
    fun=@(x) (tan(theta)*x+c-A*x.^2).^1.5;%��������
    xmin=(tan(theta)-sqrt(d))/(2*A);
    xmax=sqrt(H/A);
    volumn=quadv(fun,xmin,xmax);%��ֵ����
    volumn=volumn*cof;
    %fprintf("volumn is %f\n",volumn);%ԭ����ӡ�������Ե�
    if(volumn<vp)%��theta���volumn��С����������v��С��ʱ���˳������theta1
        break;
    end
    i=i+1;
end
