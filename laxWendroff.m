% Lax-Wendroff Method For Test Case -2 
gamma  = 1.4;
xmin=-10;
xmax=10;
N=50;
dx=(xmax-xmin)/N;
% Initialization at t(time)=0
for i =  1:N                 
    if i<=N/2
        rho(i) = 1.0;    %density
        P(i) = 100000;      %pressure
        u(i) = 0.0;      %velocity
    else
        rho(i) = 0.010;
        P(i) = 1000;
        u(i) = 0.0;
    end
end

% u-VECTOR
E = P/(gamma-1)+0.5*rho.*u.^2; %energy
U1 = rho ;
U2 = rho.*u;       
U3 = E;
lambda = 8.02 * 10^(-4);

tmax = 0.1;
dt = 4.01 * 10^(-4);
steps = tmax/dt;

for n = 1:steps 
    F1 = rho .*u; 
    F2 = rho .* u .^ 2; 
    F3 = u .* (E + P);

    for i = 2 : N-1
        U1(i) = U1(i) - lambda/2*(F1(i+1) - F1(i-1)) - lambda*lambda/2*(F1(i+1) - 2*F1(i) + F1(i-1));
        U2(i) = U2(i) - lambda/2*(F2(i+1) - F2(i-1)) - lambda*lambda/2*(F2(i+1) - 2*F2(i) + F2(i-1));
        U3(i) = U3(i) - lambda/2*(F3(i+1) - F3(i-1)) - lambda*lambda/2*(F3(i+1) - 2*F3(i) + F3(i-1));
    end 

    % Boundary Conditions 
    E(N) = E(N-1);
    U2(N) = U2(N-1);
    U3(N) = U3(N-1);
    E(1) = E(2);
    U2(1) = U2(2);
    U3(1) = U3(2);

    rho = U1;
    u = U2 ./ rho ;
    P = U3 - 0.5*rho.*u.^2 ;
    P = P .* (gamma-1); 
end    
% Plotting Scripts Go Here
pressure = P;
density = rho ;
velocity = u; 
sound = (1.4 * (pressure./density)).^(0.5);
mach = velocity ./ sound;
entropy = (8.314)/(1.4 - 1) * (log (pressure ./ (density .^ 1.4)));
figure(1)
subplot(231)
plot(x, pressure, "--or")
xlabel('X','fontSize',10);
ylabel('pressure','fontSize',10);

subplot(232)
plot(x, velocity, "--or")
xlabel('X','fontSize',10);
ylabel('velocity','fontSize',10);

subplot(233)
plot(x, sound, "--or")
xlabel('X','fontSize',10);
ylabel('speed of sound','fontSize',10);

subplot(234)
plot(x, density, "--or")
xlabel('X','fontSize',10);
ylabel('density','fontSize',10);

subplot(235)
plot(x, entropy, "--or")
xlabel('X','fontSize',10);
ylabel('entropy','fontSize',10);

subplot(236)
plot(x, mach, "--or")
xlabel('X','fontSize',10);
ylabel('mach number','fontSize',10);