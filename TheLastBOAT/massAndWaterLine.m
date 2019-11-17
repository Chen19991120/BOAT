MZmax=1;
i=0;
mFe=[];
C=[];
while(i<MZmax*1000)
    MZ=i/1000;
    mFe0=[mFe,MZ];
    mFe=mFe0;
    
    M=mg+mc+MZ;
    c=getWaterLine(M,H,A,B);
    C0=[C,c];
    C=C0;
    i=i+1;
end
plot(mFe,C);
grid on;
