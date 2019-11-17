
%x3d = linspace(-0.25,.25,100);
%y3d = linspace(-.12,.12,100);
%plotHull(x3d,y3d,H,2*sqrt(H/B),2*sqrt(H/A));hold on;
plot(T,real(COBC));hold on;
plot(T,COMZ);hold on;