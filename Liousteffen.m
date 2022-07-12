% Liou Stefan method for test case -1
gamma  = 1.4;
x_min=-10;
x_max=15;
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
 E = P/(gamma-1)+0.5*rho.*u.^2; %
 U2 = rho.*u;       
 U3 = rho.*E;
t=0;
t_end = 0.01;
dt  = 4.276 * 10^(-4); 
while t<t_end
        
        a = (gamma*P./rho).^0.5;   % speed of sound
        M = u ./ a ;
        for i = 1:N 
            if M(i) <= -1
                Mp(i) =  0;
                Mn(i) = M(i);
            elseif M(i) < 1 && M(i) > -1
                Mp(i) = (M+1)^2/4 ;
                Mn(i) = (M-1)^2/4;
            else 
                Mp(i) = M(i);
                Mn(i) = 0;
            end
        end 
        H = .5*u.^2+a.^2/(gamma-1);         % H = (energy+pressure)/density
        for i = 1:N 
            FP1(i) = max(0, Mp(i) + Mn(i))*rho(i)*a(i);
            FP2(i) = max(0, Mp(i) + Mn(i))*rho(i)*a(i)*u(i);
            FP3(i) = max(0, Mp(i) + Mn(i))*rho(i)*a(i)*H(i);

            FN1(i) = min(0, Mp(i) + Mn(i))*rho(i)*a(i);
            FN2(i) = min(0, Mp(i) + Mn(i))*rho(i)*a(i)*u(i);
            FN3(i) = min(0, Mp(i) + Mn(i))*rho(i)*a(i)*H(i);
            if M(i) <= -1
                FP2(i) = FP2(i);
                FN2(i) = FN(i) + 1;
            elseif M(i) < 1 && M(i) > -1
                FP2(i) = FP2(i) + 0.5*(1 + M(i));
                FN2(i) = FN2(i) + 0.5*(1 - M(i));
            else 
                FP2(i) = FP2(i) + 1;
                FN2(i) = FN2(i);
            end
        end
        for i = 1:N-1                     %intercell numerical flux
            Fhp1(i) = FP1(i)+FN1(i+1);
            Fhn1(i+1) = FP1(i)+FN1(i+1);
            Fhp2(i) = FP2(i)+FN2(i+1);
            Fhn2(i+1) = FP2(i)+FN2(i+1);
            Fhp3(i) = FP3(i)+FN3(i+1);
            Fhn3(i+1) = FP3(i)+FN3(i+1);
        end
  
       for i=2:N-1
        rhon(i) = rho(i)-dt*(Fhp1(i)-Fhn1(i));        
        U2(i) = rho(i).u(i)-dt(Fhp2(i)-Fhn2(i));        
        U3(i) = rho(i).E(i)-dt(Fhp3(i)-Fhn3(i));    
       end
       
      % Boundary Conditions
       rhon(N) = rhon(N-1);
       U2(N) = U2(N-1);
       U3(N) = U3(N-1);
       rhon(1) = rhon(2);
       U2(1) = U2(2);
       U3(1) = U3(2);
       
       u = U2./rhon;                       % velocity at t+dt
       E  = U3./rhon;                      % energy at t+dt
       P = (gamma-1)*(E-0.5*rhon.*u.^2);   %pressure at t+dt
       
       rho = rhon; % new density       
       t = t+dt;          % time increment                    
end

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