% This code partially uses the ps3 sample code provided by Prof. Lanteri.

clear; clc;
%% Parameters %%%%%%%
par.beta = 0.96;
par.r = 0.035;
par.phi = -0.5;
par.w = 1;
par.tau = 1/72;

par.p = 0.95;
par.q = 0.1;

par.u = @(x) log(x);

grid.A = linspace(-0.5,9.5,1000)'; % Possible levels of bond, current period 
grid.C = exp(linspace(log(0.1),log(10),1000))'; % Possible levels of consumption

%%%% Value Function Iteration %%%%%%%%%%

vfie.V = zeros(length(grid.A),1);

vfie.auxC = repmat(grid.C',length(grid.A),1);
vfie.auxA = repmat(grid.A,1,length(grid.C));

vfie.auxB = ((1-par.tau)+(1+par.r).*vfie.auxA-vfie.auxC); % bond next period
vfie.auxB(vfie.auxB<par.phi) = NaN; % filter the value hit the borrowing limit


vfiu.V = zeros(length(grid.A),1);

vfiu.auxC = repmat(grid.C',length(grid.A),1);
vfiu.auxA = repmat(grid.A,1,length(grid.C));

vfiu.auxB = ((1/4)+(1+par.r).*vfiu.auxA-vfiu.auxC); 
vfiu.auxB(vfiu.auxB<par.phi) = NaN;


vfi.err = 1; % Initialize error
vfi.tol = 1e-6; % tolerance
vfi.count = 0; % counter for repeatition


while vfi.err > vfi.tol  & vfi.count < 1000
    vfie.interpV = griddedInterpolant(grid.A,vfie.V);
    vfiu.interpV = griddedInterpolant(grid.A,vfiu.V);

    vfie.newV = par.u(vfie.auxC) + par.beta.*(par.p.*vfie.interpV(vfie.auxB)+(1-par.p).*vfiu.interpV(vfie.auxB));
    vfiu.newV = par.u(vfiu.auxC) + par.beta.*((1-par.q).*vfie.interpV(vfiu.auxB)+par.q.*vfiu.interpV(vfiu.auxB));

    [vfie.newV, vfie.index] = max(vfie.newV,[],2);
    [vfiu.newV, vfiu.index] = max(vfiu.newV,[],2);

    vfie.err = max(abs(vfie.V-vfie.newV));
    vfiu.err = max(abs(vfiu.V-vfiu.newV));
    vfi.err = max(vfie.err,vfiu.err);

    vfie.V = vfie.newV;
    vfiu.V = vfiu.newV;
    vfi.count = vfi.count + 1;
end

%% Plots of policy functions

figure(1)
subplot(4,1,1)
oute.V = vfie.V;
plot(grid.A,oute.V)
xlabel('a_{t}')
ylabel('V_{t}')
title('Value Function with state e')
subplot(4,1,2)
oute.polC = grid.C(vfie.index); % consumption policy for x=e
plot(grid.A,oute.polC)
xlabel('a_{t}')
ylabel('c_{t}')
title('Policy Function of Consumption with state e')
subplot(4,1,3)
outu.V = vfiu.V;
plot(grid.A,outu.V)
xlabel('a_{t}')
ylabel('V_{t}')
title('Value Function with state u')
subplot(4,1,4)
outu.polC = grid.C(vfiu.index); % consumption policy for x=u
plot(grid.A,outu.polC)
xlabel('a_{t}')
ylabel('c_{t}')
title('Policy Function of Consumption with state u')
saveas(gcf,'f1_pf.jpeg');

%% Simulation %%%%%%%%%

sim.T = 1000;
sim.M = [par.p (1-par.q); (1-par.p) par.q];
sim.mc = dtmc(sim.M); % simulate the markov chain
sim.x = simulate(sim.mc,sim.T) % the employment state
% 1 is employed, 2 is not employed.
sim.x = (-1).* (sim.x-2);
% now 1 is employed, 0 is not employed.

sim.interpCe = griddedInterpolant(grid.A,oute.polC);
sim.interpCu = griddedInterpolant(grid.A,outu.polC);
sim.a(1) = 0; % initial bond 
sim.y(1) = sim.x(1) * (1-par.tau) + (1-sim.x(1)) * (1/4);
sim.c(1) = sim.x(1) * sim.interpCe(sim.a(1)) + (1-sim.x(1)) * sim.interpCu(sim.a(1));

for i = 2:sim.T
    sim.a(i) = sim.y(i-1) + (1+par.r)*sim.a(i-1) - sim.c(i-1); % by budget constraint
    sim.y(i) = sim.x(i) * (1-par.tau) + (1-sim.x(i)) * (1/4);
    sim.c(i) = sim.x(i) * sim.interpCe(sim.a(i)) + (1-sim.x(i)) * sim.interpCu(sim.a(i));
end

% Plot of simulation 
figure(2)
subplot(3,1,1)
plot(1:sim.T,sim.y)
xlabel('t')
ylabel('y_{t}')
title('Simulation of Employment')
subplot(3,1,2)
plot(1:sim.T,sim.a)
xlabel('t')
ylabel('a_{t}')
title('Simulation of Asset')
subplot(3,1,3)
plot(1:sim.T,sim.c)
xlabel('t')
ylabel('c_{t}')
title('Simulation of Consumption')
saveas(gcf,'f2_sim.jpeg');


%% Distribution %%%%%%%%

% Bond policy function
oute.polB = (1-par.tau)+(1+par.r).* grid.A - oute.polC;
outu.polB = (1/4)+(1+par.r).* grid.A - outu.polC;

% indicator function: Find the closest level of bond
for i = 1:1000
    [~,ldiste.index(i)] = min(abs(grid.A - oute.polB(i)));
    [~,ldistu.index(i)] = min(abs(grid.A - outu.polB(i)));
end

ldiste.l = (1/2000).*ones(1000,1); % initializing, employed
ldistu.l = (1/2000).*ones(1000,1); % initializing, unemployed

ldist.err = 1;
ldist.tol = 1e-6;
ldist.count = 0;

while ldist.err > ldist.tol & ldist.count < 1000
    for i = 1:1000
        ldiste.newl(i) = par.p .* (ldiste.index==i) * ldiste.l + (1-par.q) .* (ldistu.index==i) * ldistu.l;
        ldistu.newl(i) = (1-par.p) .* (ldiste.index==i) * ldiste.l + par.q .* (ldistu.index==i) * ldistu.l;
    end 

    ldiste.err = max(abs(ldiste.newl'-ldiste.l));
    ldistu.err = max(abs(ldistu.newl'-ldistu.l));
    ldist.err = max(ldiste.err,ldistu.err);

    ldiste.l = ldiste.newl';
    ldistu.l = ldistu.newl';

    ldist.count = ldist.count + 1;
end

% plot the distribution

figure(3)
subplot(2,1,1)
plot(grid.A,ldiste.l)
xlabel('a')
ylabel('\lambda(a,e)')
title('Distribution of \lambda with Employment State')
subplot(2,1,2)
plot(grid.A,ldistu.l)
xlabel('a')
ylabel('\lambda(a,u)')
title('Distribution of \lambda with Unemployment State')
saveas(gcf,'f3_dist.jpeg');

% check market clearing
totalbond = ldiste.l'*grid.A+ldistu.l'*grid.A;
% The result is 0.2279

%%% Market clearing interest rate

% Try the same procedure with
par.r = 0.025;
% and the total bond is -0.06. This gives a lower bound of r
% Similarly we can reduce the upper bound to 0.03

clear; clc;

par.beta = 0.96;
par.phi = -0.5;
par.w = 1;
par.tau = 1/72;

par.p = 0.95;
par.q = 0.1;

par.u = @(x) log(x);

par.rmax = 0.03;
par.rmin = 0.025;

r.err = 1;
r.tol = 1e-6;

r.count = 0;

grid.A = linspace(-0.5,9.5,1000)'; % Possible levels of bond, current period 
grid.C = exp(linspace(log(0.1),log(10),1000))'; % Possible levels of consumption

while r.err > r.tol & r.count < 10 % WARNING: this iteration is SLOW!
    par.r = 0.5*par.rmax + 0.5*par.rmin;

    vfie.V = zeros(length(grid.A),1);

    vfie.auxC = repmat(grid.C',length(grid.A),1);
    vfie.auxA = repmat(grid.A,1,length(grid.C));
    vfie.auxB = ((1-par.tau)+(1+par.r).*vfie.auxA-vfie.auxC); % bond next period
    vfie.auxB(vfie.auxB<par.phi) = NaN; % filter the value hit the borrowing limit


    vfiu.V = zeros(length(grid.A),1);

    vfiu.auxC = repmat(grid.C',length(grid.A),1);
    vfiu.auxA = repmat(grid.A,1,length(grid.C));
    vfiu.auxB = ((1/4)+(1+par.r).*vfiu.auxA-vfiu.auxC); 
    vfiu.auxB(vfiu.auxB<par.phi) = NaN;

    vfi.err = 1; % Initialize error
    vfi.tol = 1e-6; % tolerance
    vfi.count = 0; % counter for repeatition

    while vfi.err > vfi.tol  & vfi.count < 1000
        vfie.interpV = griddedInterpolant(grid.A,vfie.V);
        vfiu.interpV = griddedInterpolant(grid.A,vfiu.V);

        vfie.newV = par.u(vfie.auxC) + par.beta.*(par.p.*vfie.interpV(vfie.auxB)+(1-par.p).*vfiu.interpV(vfie.auxB));
        vfiu.newV = par.u(vfiu.auxC) + par.beta.*((1-par.q).*vfie.interpV(vfiu.auxB)+par.q.*vfiu.interpV(vfiu.auxB));

        [vfie.newV, vfie.index] = max(vfie.newV,[],2);
        [vfiu.newV, vfiu.index] = max(vfiu.newV,[],2);

        vfie.err = max(abs(vfie.V-vfie.newV));
        vfiu.err = max(abs(vfiu.V-vfiu.newV));
        vfi.err = max(vfie.err,vfiu.err);

        vfie.V = vfie.newV;
        vfiu.V = vfiu.newV;
        vfi.count = vfi.count + 1;
    end
    oute.polC = grid.C(vfie.index);
    outu.polC = grid.C(vfiu.index);
    oute.polB = (1-par.tau)+(1+par.r).* grid.A - oute.polC;
    outu.polB = (1/4)+(1+par.r).* grid.A - outu.polC;

    % indicator function: Find the closest level of bond
    for i = 1:1000
        [~,ldiste.index(i)] = min(abs(grid.A - oute.polB(i)));
        [~,ldistu.index(i)] = min(abs(grid.A - outu.polB(i)));
    end

    ldiste.l = (1/2000).*ones(1000,1); % initializing, employed
    ldistu.l = (1/2000).*ones(1000,1); % initializing, unemployed

    ldist.err = 1;
    ldist.tol = 1e-6;
    ldist.count = 0;

    while ldist.err > ldist.tol & ldist.count < 1000
        for i = 1:1000
            ldiste.newl(i) = par.p .* (ldiste.index==i) * ldiste.l + (1-par.q) .* (ldistu.index==i) * ldistu.l;
            ldistu.newl(i) = (1-par.p) .* (ldiste.index==i) * ldiste.l + par.q .* (ldistu.index==i) * ldistu.l;
        end 

        ldiste.err = max(abs(ldiste.newl'-ldiste.l));
        ldistu.err = max(abs(ldistu.newl'-ldistu.l));
        ldist.err = max(ldiste.err,ldistu.err);

        ldiste.l = ldiste.newl';
        ldistu.l = ldistu.newl';

        ldist.count = ldist.count + 1;
    end

    totalbond = ldiste.l'*grid.A+ldistu.l'*grid.A;

    r.err = abs(totalbond)

    if totalbond>0
        par.rmax = par.r;
        par.rmin = par.rmin;
    end
    if totalbond<0
        par.rmax = par.rmax;
        par.rmin = par.r;
    end

    r.count = r.count + 1
end

figure(4)
subplot(2,1,1)
plot(grid.A,ldiste.l)
xlabel('a')
ylabel('\lambda(a,e)')
title('Distribution of \lambda with Employment State')
subplot(2,1,2)
plot(grid.A,ldistu.l)
xlabel('a')
ylabel('\lambda(a,u)')
title('Distribution of \lambda with Unemployment State')
saveas(gcf,'f3_dist.jpeg');