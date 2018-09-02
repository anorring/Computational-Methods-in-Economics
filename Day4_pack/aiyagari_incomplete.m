clear;

gamma = 2;
mu = 0.5;
alpha = 1/3;
phi = 0;
beta = 1.03^(-1/12);
delta = 0.02;
r = 1/beta-1e-4-1;
rh = 1/beta-1e-5-1;
rl = -0.01;
N = 100;

n = 0.94;
T = [(n-(1-n)*1/3)/n,1-(n-(1-n)*1/3)/n;1/3,1-1/3];
P = T;

bgrid = linspace(0,log(500-phi+1),N)';
bp = exp(bgrid)+phi-1;

bg1 = bp;
bb1 = bp;
bg = bp;
bb = bp;

mean_b = 1;
    
bg1 = bp;
bb1 = bp;
bg = bp;
bb = bp;
    
r = (rh+rl)/2;

w = (1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

metric = 1;

while metric > 1e-6

    gg = griddedInterpolant(bg1,bp,'linear');
    gb = griddedInterpolant(bb1,bp,'linear');

    Eg = beta*(1+r)*( T(1,1)*(bp*(1+r)+w-max(gg(bp),phi)).^(-gamma)+T(1,2)*(bp*(1+r)+mu*w-max(gb(bp),phi)).^(-gamma) );
    Eb = beta*(1+r)*( T(2,1)*(bp*(1+r)+w-max(gg(bp),phi)).^(-gamma)+T(2,2)*(bp*(1+r)+mu*w-max(gb(bp),phi)).^(-gamma) );

    bg1 = (Eg.^(-1/gamma)-w+bp)./(1+r);
    bb1 = (Eb.^(-1/gamma)-mu*w+bp)./(1+r);

    metric = max(max(abs([bg1-bg bb1-bb])));

    bg = bg1;
    bb = bb1;

end

bg = max(gg(bp),phi);
bb = max(gb(bp),phi);

F = griddedInterpolant(bp,bp,'next');
bgd = F(max(gg(bp),phi));
bbd = F(max(gb(bp),phi));

Tg = zeros(N,N);
Tb = zeros(N,N);

for i = 1:N
    
    ixg = find(bp==bgd(i));
    Tg(i,ixg) = 1-(bp(ixg)-bg(i))/(bp(ixg)-bp(ixg-1));
    Tg(i,ixg-1) = 1-Tg(i,ixg);
    
    ixb = find(bp==bbd(i));
    if ixb == 1
        Tb(i,ixb) = 1;
    else
    Tb(i,ixb) = 1-(bp(ixb)-bb(i))/(bp(ixb)-bp(ixb-1));
    Tb(i,ixb-1) = 1-Tb(i,ixb);
    end
    
end

M = [P(1,1)*Tg,P(1,2)*Tg;
    P(2,1)*Tb,P(2,2)*Tb];

M = sparse(M);

opts.disp = 0;
[V,D] = eigs(M',1,1,opts);
V = V/sum(V);
dist = V(1:N)+V(N+1:end);

mean_b = bp'*dist;

