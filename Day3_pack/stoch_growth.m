clear;

beta = 1.03^(-1/4);
gamma = 1;
alpha = 0.3;
delta = 0.02;
sigma = 0.23;
rho = 0;
N = 100;
kss = ((1/beta-1+delta)/alpha)^(1/(alpha-1));
kgr = linspace(0.4*kss,1.7*kss,N)';

P = [(1+rho)/2 (1-rho)/2;(1-rho)/2 (1+rho)/2];

syms km k kpg kpb s

C = exp((2-s)*sigma+(1-s)*sigma)*km.^(alpha)+(1-delta)*km-k;
Cg = exp(sigma)*k.^(alpha)+(1-delta)*k-kpg;
Cb = exp(-sigma)*k.^(alpha)+(1-delta)*k-kpb;
Rg = 1+exp(sigma)*alpha*k^(alpha-1)-delta;
Rb = 1+exp(-sigma)*alpha*k^(alpha-1)-delta;

fs = C.^(-gamma)-beta*( ((2-s)*P(1,1)-(1-s)*P(2,1))*(Rg*Cg.^(-gamma))...
    +((2-s)*P(1,2)-(1-s)*P(2,2))*(Rb*Cb.^(-gamma)) );

dfks = diff(fs,k);
dfkpgs = diff(fs,kpg);
dfkpbs = diff(fs,kpb);

f = matlabFunction(fs,'vars',[km k kpg kpb s]);

dfk = matlabFunction(dfks,'vars',[km k kpg kpb s]);

dfkpg = matlabFunction(dfkpgs,'vars',[km k kpg kpb s]);
dfkpb = matlabFunction(dfkpbs,'vars',[km k kpg kpb s]);

kg = 0.99*kgr+0.01*kss;
kb = kg*0.99;

k = [kg kb];

kpg_policy = griddedInterpolant(kgr,kg);
kpb_policy = griddedInterpolant(kgr,kb);

kpg = [kpg_policy(kg), kpg_policy(kb)];
kpb = [kpb_policy(kg), kpb_policy(kb)];

dkpg = (gradient(kpg')')./(gradient(k')');
dkpb = (gradient(kpb')')./(gradient(k')');

metric2 = 1;
tic
while metric2>1e-8
    
    metric1 = 1;

    while metric1>1e-8

        dfgr = dfk(kgr,k(:,1),kpg(:,1),kpb(:,1),1)+dfkpg(kgr,k(:,1),kpg(:,1),kpb(:,1),1).*dkpg(:,1) ...
            +dfkpb(kgr,k(:,1),kpg(:,1),kpb(:,1),1).*dkpb(:,1);
        
        dfbr = dfk(kgr,k(:,2),kpg(:,2),kpb(:,2),2)+dfkpg(kgr,k(:,2),kpg(:,2),kpb(:,2),2).*dkpg(:,2) ...
            +dfkpb(kgr,k(:,2),kpg(:,2),kpb(:,2),2).*dkpb(:,2);

        eg = f(kgr,kg,kpg(:,1),kpb(:,1),1);
        eb = f(kgr,kb,kpg(:,2),kpb(:,2),2);

        kg = kg-eg./dfgr;
        kb = kb-eb./dfbr;
        
        k = [kg kb];

        kpg = [kpg_policy(kg), kpg_policy(kb)];
        kpb = [kpb_policy(kg), kpb_policy(kb)];

        dkpg = (gradient(kpg')')./(gradient(k')');
        dkpb = (gradient(kpb')')./(gradient(k')');
        
        metric1 = max(max(abs([eg eb])));

    end
    
    kpg_policy = griddedInterpolant(kgr,kg);
    kpb_policy = griddedInterpolant(kgr,kb);  
    
    eg = f(kgr,kg,kpg(:,1),kpb(:,1),1);
    eb = f(kgr,kb,kpg(:,2),kpb(:,2),2);
    
    metric2 = max(max(abs([eg eb])));

end
toc

figure;
p1 = plot(log(kgr),log(kg)-log((1-delta))-log(kgr),'linewidth',1.6);
hold on
p2 = plot(log(kgr),log(kb)-log((1-delta))-log(kgr),'linewidth',1.6);
hold off
legend1 = legend([p1,p2],'\sigma=0.23','\sigma=-0.23');
set(legend1,'fontname','times','Location','northeast','FontSize',12)
set(gca,'FontSize',12,'fontname','times')
xlabel('(Log of) Capital','FontSize',12,'fontname','times')
ylabel('(Log of) Investment','FontSize',12,'fontname','times')

kg = min(kg,max(kgr));
kg = max(kg,min(kgr));

kb = min(kb,max(kgr));
kb = max(kb,min(kgr));

Cg = exp(sigma)*kgr.^(alpha)+(1-delta)*kgr-kg;
Cb = exp(-sigma)*kgr.^(alpha)+(1-delta)*kgr-kb;

C = [Cg;Cb];

Yg = exp(sigma)*kgr.^(alpha);
Yb = exp(-sigma)*kgr.^(alpha);

Y = [Yg;Yb];

Ing = kg-(1-delta)*kgr;
Inb = kb-(1-delta)*kgr;

In = [Ing;Inb];

F = griddedInterpolant(kgr,kgr,'next');
kgd = F(kg);
kbd = F(kb);

Tg = zeros(N,N);
Tb = zeros(N,N);

for i = 1:N
    
    ixg = find(kgr==kgd(i));
    Tg(i,ixg) = 1-(kgr(ixg)-kg(i))/(kgr(ixg)-kgr(ixg-1));
    Tg(i,ixg-1) = 1-Tg(i,ixg);
    
    ixb = find(kgr==kbd(i));
    if ixb == 1
        Tb(i,ixb) = 1;
    else
    Tb(i,ixb) = 1-(kgr(ixb)-kb(i))/(kgr(ixb)-kgr(ixb-1));
    Tb(i,ixb-1) = 1-Tb(i,ixb);
    end
    
end

T = [P(1,1)*Tg,P(1,2)*Tg;
    P(2,1)*Tb,P(2,2)*Tb];

T = sparse(T);

opts.disp = 0;
[V,D] = eigs(T',1,1,opts);
V = V/sum(V);
dist = V(1:N)+V(N+1:end);

Kgr = [kgr;kgr];

Vimp = V;
Vimp(:,2) = zeros(N*2,1);
Vimp(N+1:end,2) = dist;

Cimp(:,1) = Vimp(:,1)'*C;
Cimp(:,2) = Vimp(:,2)'*C;

Yimp(:,1) = Vimp(:,1)'*Y;
Yimp(:,2) = Vimp(:,2)'*Y;

Inimp(:,1) = Vimp(:,1)'*In;
Inimp(:,2) = Vimp(:,2)'*In;

zgr = [exp(sigma)*ones(N,1);exp(-sigma)*ones(N,1)];

Zimp(:,1) = Vimp(:,1)'*zgr;
Zimp(:,2) = Vimp(:,2)'*zgr;

Kimp(:,1) = Vimp(:,1)'*Kgr;
Kimp(:,2) = Vimp(:,2)'*Kgr;


TT = 15;

for t = 3:TT
    Vimp(:,t) = T'*Vimp(:,t-1);
    Cimp(:,t) = Vimp(:,t)'*C;
    Yimp(:,t) = Vimp(:,t)'*Y;
    Inimp(:,t) = Vimp(:,t)'*In;
    Zimp(:,t) = Vimp(:,t)'*zgr;
    Kimp(:,t) = Vimp(:,t)'*Kgr;
end
    
figure;
subplot(2,2,1);
plot(log(Zimp)-log(Zimp(1)),'linewidth',1.6);
ylabel('Productivity, z','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,2);
plot(log(Yimp)-log(Yimp(1)),'linewidth',1.6);
ylabel('Output, y','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,3);
plot(Inimp,'linewidth',1.6);
ylabel('Investment, i','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')
subplot(2,2,4);
plot(log(Cimp)-log(Cimp(1)),'linewidth',1.6);
ylabel('Consumption, c','FontSize',12,'fontname','times')
xlabel('Time (quarters)','FontSize',12,'fontname','times')
set(gca,'FontSize',12,'fontname','times')

