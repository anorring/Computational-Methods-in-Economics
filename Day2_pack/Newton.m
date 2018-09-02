clear;

syms x

eq = 50+x - x^2;

d_eq = diff(eq);

f = matlabFunction(eq);
df = matlabFunction(d_eq);

x_old = 2;

x_grid = (0:0.1:19)';

plot(x_grid,f(x_grid),'LineWidth',1.8)
xlabel('x','FontSize',12,'fontname','times')
ylabel('f(x)=50+x - x^2','FontSize',12,'fontname','times')
hold on

for i = 1:5

    x_new = x_old-f(x_old)/df(x_old);

    pause;
    plot([x_old],[f(x_old)],'r.','markersize',26)
    pause
    plot([x_old,x_new],[f(x_old)+df(x_old)*(x_old-x_old),f(x_old)+df(x_old)*(x_new-x_old)],'k','LineWidth',1.8)
    plot([x_new],[f(x_old)+df(x_old)*(x_new-x_old)],'r>','LineWidth',1.8)
    pause
    plot([x_new,x_new],[f(x_old)+df(x_old)*(x_new-x_old),f(x_new)],'k--','LineWidth',1.8)

    x_old = x_new;

end

hold off

x_old = 2;

plot([0],[f(x_old)],'r.','markersize',26);
xlabel('Iteration','FontSize',12,'fontname','times')
ylabel('Function value','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
axis([-1 11 -300 50])
hold on

for i = 1:10

    x_new = x_old-f(x_old)/df(x_old);

    x_old = x_new;
    
    pause;
    plot([i],[f(x_old)],'r.','markersize',26)
        

end



