%% Chapter 3 - CO2 Uptake model, uncorrected

%% some housekeeping

close all;
clear all;

fontname = 'Arial'; %OUP recommended font 
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

%% Define constants - ALL FOR COL-0 HERE
Lz          = 91.3e-6;    % Tube height (m)
phi_c0      = 0.0164;   % CO2 conc at lower boundary (mol/m^3)
phi_cLz     = phi_c0;   % CO2 conc at upper boundary (mol/m^3)

D_c         = 0.139e-4; % CO2 diffusion constant (m^2/s)
G_c         = 0.00042;     % CORRECTED conductance coeff. G(m/s)
Lx          = 0.001;     % Leaf width (m)
d           = 3.296e-5;    % Minimum palisade cell width (m)?
rhovec      = [10e-7; 10e-6; 10e-5; 10e-4];   % Tube diameter vector (m)
zvec        = (0:Lz/100:Lz)';                 % Height vector (m)
Am = Lx^2;

%% CO2 concentration vs air channel depth 
figure(1)
line_color = ['b' 'g' 'm' 'r'];
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k_c = sqrt(4*G_c/(D_c*rho));
    phivec_c = (phi_cLz*sinh(k_c.*zvec)+phi_c0*sinh(k_c.*(Lz-zvec)))/(sinh(k_c*Lz));
    plot(phivec_c,zvec,'Color', line_color(rr),'LineWidth',2); grid on; hold on;
    axis([0.007 phi_c0 0 Lz]); axis square;
end
xlabel('\phi_{CO_2}(z), mol/m^{3}','FontSize',18);
ylabel('z (m)','FontSize',18);
llgd=legend('\rho=10e^{-7}', '\rho=10e^{-6}','\rho=10e^{-5}','\rho=10e^{-4}');
llgd.FontSize = 14;

%% CO2 uptake vs diameters, for single air channel
figure(2);
srhovec = linspace(10e-9,10e-2,50);
CO2singlevec = zeros(length(rhovec),1);
for rr=1:length(srhovec)
    rho = srhovec(rr);
    k = sqrt(4*G_c/(D_c*rho));
    kk = (k/sqrt(rho))*Lz;
    CO2singlevec(rr) =(pi*G_c*(rho^1.5))*(phi_cLz+phi_c0)*(coth(kk)-csch(kk));
end
semilogx(srhovec,CO2singlevec,'k','LineWidth',3); grid on;
axis([srhovec(1) srhovec(end) CO2singlevec(1) CO2singlevec(end-1)]); axis square;
xlabel('Air channel diameter \rho (m)','FontSize',18);
ylabel('U_{CO_2} (\rho) (mols^{-1})','FontSize',18);%% CO2 uptake vs diameters, for 2D with no separation
ax=gca; ax.YAxis.Exponent = -6;

%% CO2 uptake vs diameters, for 2D with no sep
figure(3);
srhovec = linspace(10e-9,10e-2,100);
CO2singlevec = zeros(length(rhovec),1);
for rr=1:length(srhovec)
    rho = srhovec(rr);
    k = sqrt(4*G_c/(D_c*rho));
    kk = (k/sqrt(rho))*Lz;
    CO2singlevec(rr) =((pi*G_c*(rho^1.5))*(phi_cLz+phi_c0)*(coth(kk)-csch(kk)))/(rho^2);
end
loglog(srhovec,CO2singlevec,'k','LineWidth',3); grid on;
xlabel('Air channel diameter \rho (m)','FontSize',18);
ylabel('U^{leaf}_{CO_2} (\rho) (molm^{-2}s^{-1})','FontSize',18);

%% CO2 uptake vs diameters, for 2D with cell thickness d
figure(4);
rhovec = logspace(-7,-2);
nCO2tubevec = zeros(length(rhovec),1);
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k = sqrt(4*G_c/(D_c*rho));
    nCO2tubevec(rr) = (pi*G_c*(rho^1.5)/((sqrt(rho))*k*((rho+d)^2)*sinh(k*Lz)))...
                *(phi_cLz*cosh(k*Lz)-phi_c0-phi_cLz+phi_c0*cosh(k*Lz));
end

plot(rhovec, nCO2tubevec);
loglog(rhovec,nCO2tubevec,'k','LineWidth',3,'MarkerSize',10); grid on;
axis([rhovec(1) rhovec(end) 1e-6 1.5e-4]); axis square;
xlabel('Air channel diameter \rho (m)','FontSize',18);
ylabel('U^{leaf}_{CO_2} (\rho) (molm^{-2}s^{-1})','FontSize',18);
 hold on;
 %xline(1.8e-5, '-g'); xline(8e-5, '-g');
plot([2.48e-5 2.48e-5], [1e-6 1.5e-4], 'Color', [0.00000  0.55294  0.97647], 'LineWidth', 2); %Col0
plot([3.45e-5 3.45e-5], [1e-6 1.5e-4], 'Color', [0.51765,  0.00000,  0.80392], 'LineWidth', 2); %Arp3
plot([3.70e-5 3.70e-5], [1e-6 1.5e-4], 'Color', [0.00000,  0.62353,  0.50588], 'LineWidth', 2); %qua2
plot([2.48e-5 2.48e-5], [1e-6 1.5e-4], 'Color', [1.00000,  0.35294,  0.68627], 'LineWidth', 2); % epf2oe
plot([2.87e-5 2.87e-5], [1e-6 1.5e-4], 'Color', [1.00000,  0.43137,  0.22745], 'LineWidth', 2); % focl
plot([2.20e-5 2.20e-5], [1e-6 1.5e-4], 'Color', [0.00000,  0.98824,  0.81176], 'LineWidth', 2); % atml1
plot([1.07e-4 1.07e-4], [1e-6 1.5e-4], 'Color', [0.88627,  0.00392,  0.20392], 'LineWidth', 2); % re6
 
 ldg = legend('Model output', 'Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6');
 lgd.FontSize = 24;