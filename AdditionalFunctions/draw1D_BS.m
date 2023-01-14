function draw1D_BS(ER,Ey,Hx,dz)

z = linspace(0,length(ER)*dz,length(ER));
imagesc([0 length(ER)*dz], [1.5 -1.5] , ER);
colormap(flipud(gray));
hold on
plot(z,Ey,'w','linewidth',4)
plot(z,Ey,'b','linewidth',2)
plot(z,Hx,'w','linewidth',4)
plot(z,Hx,'r','linewidth',2)
hold off

end