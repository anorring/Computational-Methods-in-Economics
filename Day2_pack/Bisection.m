clear;

% Function will be f(x)=exp(x)-exp(2.2087). For this function we know that
% f(0)<0 and f(4)>0, so that will do. Thus we will pick a=0 and b=4.

syms x

eq = exp(x)-exp(2.2087);

f = matlabFunction(eq);

a = 0;

b = 4;

x_grid = (0:0.01:4.5)';

plot([0 0],[a b],'LineWidth',1.8,'color','k')
xlabel('Iteration','FontSize',12,'fontname','times')
ylabel('Interval','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
axis([-1 11 -1 5])
hold on
plot(-1:11,2.2087*ones(13,1),'--')


for i = 1:10

    
    x = (a+b)/2;

    plot([i-1],[x],'b.','markersize',26)
    
    ff = f(x);
    
    if ff*f(b)<=0
        a = x;
    else
        b = x;
    end
    
%     [a,b]
    pause;
    plot([i i],[a b],'LineWidth',1.8,'color','k')
    
end

plot([i],(a+b)/2,'b.','markersize',26)

% But we can also solve this using Newton's method

pause;

x_grid = (0:0.01:4.5)';

d_eq = diff(eq);
df = matlabFunction(d_eq);

x_old = 4;

plot([0],[x_old],'r.','markersize',26)

for i = 1:10

    x_new = x_old-f(x_old)/df(x_old);

    x_old = x_new;
    
    plot([i],[x_old],'r.','markersize',26)
    pause;

end

