clear;

% Parameter values

delta = 0.98;
eta = 1/3;

% Declare symbolic variables and functions

syms y x z

fs =  % Insert code here for the objective function

dfys = % Insert code her for the symbolic derivatives
dfzs = % Insert code here for the symbolic derivatives

% Convert functions to executable/numerical

f = matlabFunction(fs,'vars',[y x z]);
dfy = matlabFunction(dfys,'vars',[y x z]);
dfz = matlabFunction(dfzs,'vars',[y x z]);

% Building the grid for x (centered around z)

z = % insert code here
xg = % Insert code for the grid for x (cetered around z, +- 20%)

% Initial guess for y and z=g(y)

yg = 0.02-0.3*xg+1.3*z;
yg0 = yg;                   % Store for plotting later
ygo = yg;                   % Initial guess for z=g(y) (values)
zg = % Insert code for the initial guess for g(y) (approximated function)
dyg = gradient(yg);         % Gradient of y = g(g(y))
dzg = gradient(zg);         % Gradient of z = g(y)
dzy = dzg./dyg;             % d g(y) / d y or dz/dy

% Iteration (loop)

metric2 = 1;    % Distance from 0 (second loop)


while metric2>1e-6  % Second loop
    
    metric1 = 1;    % Distance from 0 (first loop)
    
    while metric1>1e-6  % First loop    
    
        fg = f(yg,xg,zg);
        dfg = % Insert code for the total derivative
        yg = % Insert code for Newton's algorithm

        metric1 = max(abs(fg))      % Update distance
        
        zg = interp1(xg,ygo,yg);    % Update the z points, keeping g(y) constant
        dyg = gradient(yg);
        dzg = gradient(zg);
        dzy = % Insert code for dz/dy
    
    end % End of first loop
    
    zg = interp1(xg,yg,yg);     % Update z=g(y)
    ygo = yg;                   % Update points for interpolation
    
    fg = f(ygo,xg,zg);          % Update y=g(g(y))
    
    metric2 = max(abs(fg))      % Update distance
    
end % End of second loop

% Plot the results

% Main functions
figure;
p4 = plot(xg,eta*delta*xg.^(eta),'o');  % Exact solution
hold on
p1 = plot(xg,yg,'linewidth',1.6);       % Solution
p2 = plot(xg,yg0,'linewidth',1.6);      % Initial guess
p3 = plot(xg,xg,'k','linewidth',1.6);   % 45 degree line
hold off
% Legend and labels
legend1 = legend([p1,p2,p3,p4],'Solution','Initial guess','45^0 line','Exact solution');
set(legend1,'fontname','times','Location','southeast','FontSize',12)
set(gca,'FontSize',12,'fontname','times')
xlabel('Grid for x','FontSize',12,'fontname','times')
ylabel('Solution, y','FontSize',12,'fontname','times')