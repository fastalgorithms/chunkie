function fieldplot(u,chnkr,xg,yg,xylim,ngr,fontsize)
F1=zeros(ngr);
for k=1:ngr
  F1(1:ngr,k)=u((k-1)*ngr+(1:ngr));
end
fh =  findobj('type','figure');
nf = length(fh);
figure(nf+1)
imagesc(xg,yg,real(F1));      
colormap(jet)
axis xy
axis equal
colorbar
% figure(nf+2)
% imagesc(xg,yg,imag(F1));      
% colormap(jet)
% axis xy
% axis equal
% colorbar

hold on
plot(chnkr)
% np=length(z)/2;
% xy=1.1*xylim;
% zext=[xy(1);0;z(1:np);1;xy(2);xy(2)+1i*xy(3);xy(1)+1i*xy(3);xy(1)];
% fill(real(zext),imag(zext),'w','EdgeColor','w')  
% zext=[xy(2);1;z(np+1:2*np);0;xy(1);xy(1)+1i*xy(4);xy(2)+1i*xy(4);xy(2)];
% fill(real(zext),imag(zext),'w','EdgeColor','w')  
%title('Field $u({\bf x})$','Interpreter','LaTeX')
xlabel('$x_1$','Interpreter','LaTeX','FontSize',fontsize)
ylabel('$x_2$','Interpreter','LaTeX','FontSize',fontsize)
axis(xylim)
axis equal