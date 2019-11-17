function c=getWaterLine(m,H,A,B)
vp=m/1000;cof=8/(3*sqrt(B));
i=1;
cmin=0;cmax=H;
c=double((cmin+cmax));
while(i<=100)
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
end