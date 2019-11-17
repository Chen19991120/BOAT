function res = plotHull(x,y,D,L,B)
    [X,Y] = meshgrid(x,y);
    res = ((D/(B/2).^2).*Y.^2) + ((D/(L/2).^2).*X.^2);
    indicies = find(res > D);
    res(indicies) = NaN;
    surf(X,Y,res)
    colormap default
    xlabel('x')
    ylabel('y')
    zlabel('z')
    shading('faceted')
    axis('equal')
end