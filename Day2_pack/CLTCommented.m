clear;

% Specify a weird looking distribution

xbar = % insert your code here

a = % insert your code here

pdf = % insert your code here

% Done, we can plot it.

xgrid = linspace(0,1,1000);
figure; 
% insert your code here - don't forget to label & title your figure!






% Calculate CDF
cdf = % insert your code here

% Now calculate the mean and the variance of the distribution. I'm too lazy
% to do this by hand, so I use numerical integration instead.

% Mean:
mpdf = @(x) x.*pdf(x);
mu = quad(mpdf,0,1);

% Variance:
vpdf = @(x) pdf(x).*(x-mu).^2;
sigma2 = quad(vpdf,0,1); 

% The inverse of the CDF. This I did by hand. But you could use
% interpolation techniques if you would like to.
inv_cdf = % insert code here

% Plot the cdf and its inverse
figure;
subplot(1,2,1)
% insert code here to plot CDF
xlabel('Support, x','FontSize',12,'fontname','times')
ylabel('Cumulative distribution','FontSize',12,'fontname','times')
title('Cumulative Distribution Function','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(1,2,2)
% insert code here to plot inverse CDF
xlabel('Cumulative distribution','FontSize',12,'fontname','times')
ylabel('Support, x','FontSize',12,'fontname','times')
title('Inverse Cumulative Distribution Function','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')


% Now let's check out the central limit theorem.
 

%Insert code here to generate Sn 
 
 
 
 
 
 
 
 
 
% Let's plot these results too and see how it compares to a N(mu,sigma^2).
figure; 
h = histogram(Sn+mu,'Normalization','pdf');
hold on
x = linspace(min(Sn)+mu,max(Sn)+mu,ceil(sqrt(M)));
sigma = sqrt(sigma2);
f = exp(-(x-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(x,f,'LineWidth',1.8)
xlabel('Support, x','FontSize',12,'fontname','times')
ylabel('Density','FontSize',12,'fontname','times')
title('Histogram vs. Probability Density Function','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')

% Check variance
sigma2_hat = % insert code

% And check normality.
[t,p] = jbtest(Sn);

