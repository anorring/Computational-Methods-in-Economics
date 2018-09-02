clear all;
close all;

% Parameter values

delta = 0.98;
eta = 1/3;

% Declare symbolic variables and functions

syms y x z

fs = % Insert code for the objective function here

dfs = % Insert code for the derivative

% Convert functions to executable/numerical

f = matlabFunction(fs,'vars',[y x z]);
df = matlabFunction(dfs,'vars',[y x z]);

z = % insert code here for z

% Building the grid for x

xg = % Insert code to generate a grid for x which has dimensions 100x1, centered at z, with a range of +/-20%


% Initial guess for y

yg = 0.02-0.3*xg+1.3*z;
yg0 = yg;	% Store for plotting later

% Iteration (loop)

metric = 1;				% Distance from 0

while metric>1e-6
    
    fg = f(yg,xg,z);
    dfg = df(yg,xg,z);
    yg = % Insert code for Newton's method here
    
    metric = max(abs(fg))
    
end

% Plot the results

% Main functions
figure;
p1 = plot(xg,yg,'linewidth',1.6);		% Solution
hold on
p2 =plot(xg,yg0,'linewidth',1.6);		% Initial guess
p3 =plot(xg,xg,'k','linewidth',1.6);	% 45 degree line
hold off
% Legend and labels
legend1 = legend([p1,p2,p3],'Solution','Initial guess','45^0 line');
set(legend1,'fontname','times','Location','southeast','FontSize',12)
set(gca,'FontSize',12,'fontname','times')
xlabel('Grid for x','FontSize',12,'fontname','times')
ylabel('Solution, y','FontSize',12,'fontname','times')