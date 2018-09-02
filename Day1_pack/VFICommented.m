% Clear out old stuff in the workspace/ old figures.

clear all;
close all;

% Parameters

alpha = 0.3;            % The capital share of output.
beta = (1.04)^(-1/4);   % The discount factor.
delta = 0.02;              % The rate of depreciation.
sigma = 1;              % The coefficient of relative risk aversion. 
N = 1000;                % # of gridpoints.

% Declare the utility function - this is CRRA utility, which has a
% different form depending on the value of sigma 
% hint - remember l'hopital's rule?

if sigma == 1
    u = @(c) % insert your code here
else
    u = @(c) % insert your code here
end

% State space (grid)

kss = % insert your code here for the steady state level of capital.
kgrid = linspace(0.9*kss,1.1*kss,N)';                      % The grid, an N x 1 vector, is defined around the steady state.
kmatrix = repmat(kgrid,1,N);                               % An, N x N, matrix with kgrid as columns, N times.
kmatrixprime = kmatrix';                                   % The same but transposed.

% Productivity values

zg = 1.02;
zb = 0.98;

% And the associated transition matrix

P = % insert your code here so that pr(z'=g|z=g) = 0.7, pr(z'=b|z=g)=0.3, pr(z'=g|z=b)=0.3, pr(z'=b|z=b)=0.7

% Declare consumption matrices for each grid point and for each possible
% choice on the grid. These are two N-by-N matrices

cmatrix_g = % insert your code here for consumption in the good state
cmatrix_b = % insert your code here for consumption in the bad state

% Remove negative values

cmatrix_g(cmatrix_g<0) = 0;  % the inner bracket returns an array with values 0 or 1 depending on whether the condition is met
cmatrix_b(cmatrix_b<0) = 0;  

% Utility matrix

utility_g = u(cmatrix_g)';
utility_b = u(cmatrix_b)';

% The initial value functions will be two N-by-1 vectors of zeros.

v_g0 = zeros(N,1);
v_b0 = zeros(N,1);

% We are now ready to go ahead using value function iteration.

metric = 1;    % start the metric off at a distance of 1
its = 0;       % initialise the iteration counter

while metric>1e-6
    
    V_g = repmat(v_g0,1,N);
    V_b = repmat(v_b0,1,N);

    v_g1 = (max(%insert your code here)))';
    
    v_b1 = (max(%insert your  code here)))';    
    
    metric = max(max(abs([v_g1-v_g0 v_b1-v_b0]))); %the code will break here if the the metric is small enough
    
    v_g0 = v_g1;   % update the value functions
    v_b0 = v_b1;
    its = its + 1  % no semicolon here means matlab spits the results out to the screen - so we can see the iterations
end

% Done!

% Let's find the policy function

[v_g1,kprime_gIX] = (max(utility_g+beta*(P(1,1)*V_g+P(1,2)*V_b))); % kprime_gIX is the index of the max

[v_b1,kprime_bIX] = (max(utility_b+beta*(P(2,1)*V_g+P(2,2)*V_b))); 

kprime_g = kgrid(kprime_gIX); % find the capital level which max utility
kprime_b = kgrid(kprime_bIX);

% Plot the decision rules in the (k,k') plane

if delta == 1 && sigma == 1 % for question vi)
    % Main functions
    figure;
    % plot the good state capital policy, numerical solution
    p1 = plot(%your code,%your code,'linewidth',1.6,'color','b');
    hold on
    % plot the bad state capital policy, numerical solution
    plot(%your code,% your code,'linewidth',1.6,'color','b');
    title('Policy function(s)');
    % plot the good state capital policy, analytic solution
    p2 = plot(%your code,'linewidth',1.6,'color','r');
    % plot the bad state capital policy, analytic solution
    plot(%your code,'linewidth',1.6,'color','r')
    plot(kgrid,kgrid,'k')
    % Legends and labels
    legend1 = legend([p1,p2],'Numerical solution','Exact solution');
    set(legend1,'fontname','times','Location','best','FontSize',12)
    set(gca,'FontSize',12,'fontname','times')
    xlabel('Capital today,k','FontSize',12,'fontname','times')
    ylabel('Capital tomorrow, k''','FontSize',12,'fontname','times')
    hold off
else
    % Main functions
    figure;
    % plot the good state capital policy
    p1 = plot(%your code, %your code,'linewidth',1.6,'color','b');
    hold on
    % plot the bad state capital policy
    p2 = plot(%your code,your code,'linewidth',1.6,'color','r');
    title('Policy function(s)');
    plot(kgrid,kgrid,'k')
    % Legends and labels
    legend1 = legend([p1,p2],'Good state','Bad state');
    set(legend1,'fontname','times','Location','best','FontSize',12)
    set(gca,'FontSize',12,'fontname','times')
    xlabel('Capital today,k','FontSize',12,'fontname','times')
    ylabel('Capital tomorrow, k''','FontSize',12,'fontname','times')
    hold off
end

% Let's simulate the model

T = 10000;   % Length of simulation
% insert your code here to set the random seed
e = % insert your code here to create a Tx1 'random' vector
ksim = zeros(T+1,1);
zsim = ksim;
z = [zg;zb];
kprime = [kprime_gIX',kprime_bIX'];
ksim(1) = find(kgrid==interp1(kgrid,kgrid,kss,'nearest'));
zsim(1) = 1;

% make 10,000 draws from a uniform distribution on [0; 1]
for t = 1:T
    if % insert your code here
        % insert your code here
    else
        % insert your code here
    end
    ksim(t+1) = kprime(ksim(t),zsim(t));
end

k = kgrid(ksim);
zz = z(zsim(1:end-1));
y = z(zsim(1:end-1)).*k(1:end-1).^(alpha);
in = k(2:end)-(1-delta)*k(1:end-1);
c = y-in;

% we only want to plot the last 100 periods
zplot = % insert your code here
yplot = % insert your code here
iplot = % insert your code here
cplot = % insert your code here

% Let's plot these results

figure;
subplot(2,2,1);
plot(%your code,'linewidth',1.6,'color','k');
ylabel('Productivity, z','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,2);
plot(%your code,'linewidth',1.6,'color','k');
ylabel('Output, y','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,3);
plot(% your code,'linewidth',1.6,'color','k');
ylabel('Investment, i','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,4);
plot(%your code,'linewidth',1.6,'color','k');
ylabel('Consumption, c','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')

% And let's calculate some moments

muk = % your code, mean of k
mui = % your code, mean of i
muy = % your code, mean of y
muc = % your code, mean of c

x = log([zz y in c]);

cor = % insert your code here to calculate the correlation matrix
% remember to calculate the RELATIVE standard deviation of the logarithm of each variable to output.


disp('(Log) correlations between z, y, i, and c, respectively:') % this prints it to your command window
disp(cor)
disp('Diagonal entries are standard deviations relative to output')

% Now find the unconditional distribution

Tg = zeros(N,N);
Tb = zeros(N,N);

for i = 1:N
    Tg(i,kprime_gIX(i)) = 1; %pick out only the states that maximise value function
    Tb(i,kprime_bIX(i)) = 1;
end

Tm = [P(1,1)*Tg,P(1,2)*Tg;P(2,1)*Tb,P(2,2)*Tb];

Tm = sparse(Tm); % squeezes out the zero elements

opts.disp = 0; % don't display some things
[V,D] = % insert your code here to find the long run distribution as the eigenvector associated with
%a unit eigenvalue
V = V/sum(V); % normalise to one
dist = V(1:N)+V(N+1:end);

% Let's plot it and compare it to the simulated distribution.

% Create a histogram for the simulated distribution

h = zeros(N,1);
h1000 = zeros(N,1);

for i = 1:N
    h(i) = sum(k==kgrid(i));
    h1000(i) = sum(k(1:1000)==kgrid(i));
end

% Main functions
figure; 
p3 = plot(kgrid,h1000./(1000),'LineWidth',1.8);
hold on
p1 = plot(kgrid,h./(T+1),'LineWidth',1.8);
p2 = plot(kgrid,dist,'LineWidth',1.8);
% Legends and labels
xlabel('Capital stock','FontSize',12,'fontname','times')
ylabel('Frequency/Probability','FontSize',12,'fontname','times')
title('Unconditional distribution of capital','FontSize',12,'fontname','times')
legend1 = legend([p2,p3,p1],'Theoretical distribution','Simulated frequency (T=1,000)','Simulated frequency (T=10,000)');
set(legend1,'fontname','times','Location','best','FontSize',12)
set(gca,'FontSize',12,'fontname','times')





