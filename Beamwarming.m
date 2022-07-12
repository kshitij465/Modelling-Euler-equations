% Beam Warming using Van Leer (Test Case-1)
gamma  = 1.4;
lambda = 1.069*10^(-3) ;
x_min=-10;
x_max=10;
N=50;
d_x=(x_max-x_min)/N;
for i =  1:N                  % Initialization at t(time)=0
    if i<=N/2
        rho(i) = 1.0;    %density
        P(i) = 100000;      %pressure
        u(i) = 0.0;      %velocity
    else
        rho(i) = 0.125;
        P(i) = 10000;
        u(i) = 0.0;
    end
end
 E = P/(gamma-1)+0.5*rho.*u.^2; %energy
 U1 = rho ;
 U2 = rho.*u;       
 U3 = E;
t=0;
t_end = 0.01;
dt  = 4.276 * 10^(-4); 
while t<t_end
        
        a = (gamma*P./rho).^0.5;   % speed of sound
        M = u ./ a ;
        for i = 1:N 
            if M(i) <= -1
                FP1(i) = 0;
                FP2(i) = 0;
                FP3(i) = 0;

                FN1(i) = (rho(i)*a(i))*M(i);
                FN2(i) = (rho(i)*a(i)*a(i))*((gamma)*M(i)*M(i) + 1)/gamma ;
                FN3(i) = (rho(i)*a(i)*a(i)*a(i)*M(i))*(0.5*M(i)*M(i) + 1/(gamma - 1));
            elseif M(i) < 1 && M(i) > -1
                FP1(i) = (rho(i)*a(i)/4)*(M(i) + 1)^2;
                FP2(i) = (rho(i)*a(i)/4)*(M(i) + 1)^2*((gamma-1)*u(i) + 2*a(i))/gamma ;
                FP3(i) = (rho(i)*a(i)/4)*(M(i) + 1)^2*((gamma-1)*u(i) + 2*a(i))^2/(2*(gamma+1)*(gamma-1));

                FN1(i) = -(rho(i)*a(i)/4)*(M(i) - 1)^2;
                FN2(i) = -(rho(i)*a(i)/4)*(M(i) - 1)^2*((gamma-1)*u(i) - 2*a(i))/gamma ;
                FN3(i) = -(rho(i)*a(i)/4)*(M(i) - 1)^2*((gamma-1)*u(i) - 2*a(i))^2/(2*(gamma+1)*(gamma-1));
            else 
                FP1(i) = (rho(i)*a(i))*M(i);
                FP2(i) = (rho(i)*a(i)*a(i))*((gamma)*M(i)*M(i) + 1)/gamma ;
                FP3(i) = (rho(i)*a(i)*a(i)*a(i)*M(i))*(0.5*M(i)*M(i) + 1/(gamma - 1));

                FN1(i) = 0;
                FN2(i) = 0;
                FN3(i) = 0;
            end
        end
        for i = 1:N 
            Ub1(i) = U1(i) - lambda*(FP1(i) - FP1(i-1) - FN1(i-1) + FN1(i));
            Ub2(i) = U2(i) - lambda*(FP2(i) - FP2(i-1) - FN2(i-1) + FN2(i));
            Ub3(i) = U3(i) - lambda*(FP3(i) - FP3(i-1) - FN3(i-1) + FN3(i));
        end 
        for i = 3 : N-2
            U1(i) = 0.5*(Ub1(i) + U1(i)) -lambda/2*(FP1(i)-FP1(i-1)) - lambda/2*(FN1(i+1)-FN1(i)) - lambda/2*(FP1(i)-2*FP1(i-1) + FP1(i-2)) + lambda/2*(FN1(i+2)-2*FN1(i+1) + FN1(i)); 
            U2(i) = 0.5*(Ub2(i) + U2(i)) -lambda/2*(FP2(i)-FP2(i-1)) - lambda/2*(FN2(i+1)-FN2(i)) - lambda/2*(FP2(i)-2*FP2(i-1) + FP2(i-2)) + lambda/2*(FN2(i+2)-2*FN2(i+1) + FN2(i)); 
            U3(i) = 0.5*(Ub3(i) + U3(i)) -lambda/2*(FP3(i)-FP3(i-1)) - lambda/2*(FN3(i+1)-FN3(i)) - lambda/2*(FP3(i)-2*FP3(i-1) + FP3(i-2)) + lambda/2*(FN3(i+2)-2*FN3(i+1) + FN3(i)); 
        end 
        
    %   Boundary Conditions
        U1(N) = U1(N-1);
        U2(N) = U2(N-1);
        U3(N) = U3(N-1);
        U1(1) = U1(2);
        U2(1) = U2(2);
        U3(1) = U3(2);
       
        rho = U1;
        u = U2 ./ rho ;
        P = U3 - 0.5*rho.*u.^2 ;
        P = P .* (gamma-1); 
        
        t = t+dt;          % time increment                    
end

pressure = P;
density = rho ;
velocity = u; 
sound = (1.4 * (pressure./density)).^(0.5);
mach = velocity ./ sound;
entropy =(8.314)/(1.4 - 1) * (log (pressure ./ (density .^ 1.4)));
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