%% Water chapter
close all;
clear all;

%% some housekeeping

fontname = 'Arial';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);


%% Define constants - FOR COL-0 ONLY HERE
Lz          = 91.3e-6;    % Tube height (m)
phi_c0      = 0.0164;   % CO2 conc at lower boundary (mol/m^3)
phi_cLz     = phi_c0;   % CO2 conc at upper boundary (mol/m^3)
phi_h0      = 0.49;      % H2O conc at lower boundary (mol/m^3)
phi_hLz     = phi_h0;   % H2O conc at upper boundary (mol/m^3)
phi_cell    = 55300;      % H2O conc in leaf tissue (mol/m^3)
phi_hsat    = 1.44; %saturation H2o concentration

D_c         = 0.139e-4; % CO2 diffusion constant (m^2/s)
D_h         = 0.282e-4; % H2O diffusion constant (m^2/s)
R           = 8.314; % gas constant
T           = 298.15; %temp in Kelvin
G_c         = 2.50e-3;     % CO2 wall flux constant (m/s)
kAdj        = 3.753e-8; % leaf hydraulic conductance adjusted to 25 deg C
G_h         = kAdj*R*T;        % H2O conductance coeff (m/s)
Lx          = 0.001;     % Leaf width (m)
d           = 35e-6;    % Minimum palisade cell width (m)?
rhovec      = [10e-7; 10e-6; 10e-5; 10e-4; 10e-3];   % Tube diameter vector (m)
zvec        = (0:Lz/100:Lz)';                 % Height vector (m)


%% Plotting concentrations
figure(11)
% %CO2
% subplot(1,2,1)
% rhovec      = [10e-7; 10e-6; 10e-5; 10e-4];   % Tube diameter vector (m)
% line_color = ['b' 'g' 'm' 'r'];
% for rr=1:length(rhovec)
%     rho = rhovec(rr);
%     k_c = sqrt(4*G_c/(D_c*rho));
%     phivec_c = (phi_cLz*sinh(k_c.*zvec)+phi_c0*sinh(k_c.*(Lz-zvec)))/(sinh(k_c*Lz));
%     plot(phivec_c,zvec,'Color', line_color(rr),'LineWidth',2); grid on; hold on;
%     axis([0.007 phi_c0 0 Lz]); %axis square;
% end
% xlabel('\phi_{CO_2}(z), mol/m^{3}','FontSize',14);
% ylabel('z (m)','FontSize',14);
% llgd=legend('\rho=10e^{-7}', '\rho=10e^{-6}','\rho=10e^{-5}','\rho=10e^{-4}');
% llgd.FontSize = 14;
% %H20
% subplot(1,2,2)
line_color = ['b' 'g' 'm' 'r'];
rhovec      = [10e-7; 10e-6; 10e-5; 10e-4];   % Tube diameter vector (m)
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k_h = sqrt(4*G_h/(D_h*rho));
    phivec_h = ((phi_hLz-phi_cell)*sinh(k_h.*zvec)...
                +(phi_h0-phi_cell)*sinh(k_h.*(Lz-zvec))...
                +phi_cell*sinh(k_h*Lz))/sinh(k_h*Lz);
    plot(phivec_h,zvec,'Color', line_color(rr),'LineWidth',2); grid on; hold on;
    axis([0 10 0 Lz]); %axis square;
end
plot(phi_hsat*ones(length(zvec),1),zvec,'k','LineWidth',2);
xlabel('\phi_{H_{2}O}(z), mol/m^{3}','FontSize',18);
ylabel('Air channel depth (m)','FontSize',18);
llgd=legend('\rho=10e^{-7}', '\rho=10e^{-6}','\rho=10e^{-5}','\rho=10e^{-4}', 'Saturation');
llgd.FontSize = 18;

%% Plotting water loss  
figure(12); right_color = [0 244 7]/255;
% water loss
rhovec = logspace(-7,-2);
nH2Otubevec = zeros(length(rhovec),1);
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k = sqrt(4*G_h/(D_h*rho));
    nH2Otubevec(rr) = (1/(rho+d)^2)*(pi*G_h*rho*phi_cell*Lz...
                        +((pi*G_h*(rho^1.5))/(sqrt(4*G_h/D_h)))*...
                        (phi_cell-phi_hLz-phi_h0)*(coth(k*Lz)-csch(k*Lz)));    
end
loglog(rhovec,nH2Otubevec,'k','LineWidth',3,'MarkerSize',10, 'color',right_color); grid on;
axis([rhovec(1) rhovec(end) 1e-1 2e1]);
axis square;
ylabel('\Gamma_{leaf},Water loss rate (\rho) (molm^{-2}s^{-1})','FontSize',18);
xlabel('Air channel diameter \rho (m)','FontSize',18);
%% Plotting leaf's dilemma
fig = figure(13);
left_color = [0 0 0]/255;
right_color = [0 244 7]/255;
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% CO2 uptake
rhovec = logspace(-7,-2);
nCO2tubevec = zeros(length(rhovec),1);
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k = sqrt(4*G_c/(D_c*rho));
    nCO2tubevec(rr) = (pi*G_c*(rho^1.5)/((sqrt(rho))*k*((rho+d)^2)*sinh(k*Lz)))...
                *(phi_cLz*cosh(k*Lz)-phi_c0-phi_cLz+phi_c0*cosh(k*Lz));
end
plot(rhovec, nCO2tubevec);
loglog(rhovec,nCO2tubevec,'k','LineWidth',3,'MarkerSize',10); grid on;hold on;
xlabel('Air channel diameter \rho (m)','FontSize',18);
ylabel('U^{leaf}_{CO_2} (\rho) (molm^{-2}s^{-1})','FontSize',18);
plot([2.48e-5 2.48e-5], [1e-7 1.5e-4], 'Color', [0.00000  0.55294  0.97647], 'LineWidth', 2); %Col0
plot([3.45e-5 3.45e-5], [1e-7 1.5e-4], 'Color', [0.51765,  0.00000,  0.80392], 'LineWidth', 2); %Arp3
plot([3.70e-5 3.70e-5], [1e-7 1.5e-4], 'Color', [0.00000,  0.62353,  0.50588], 'LineWidth', 2); %qua2
plot([2.48e-5 2.48e-5], [1e-7 1.5e-4],':', 'Color', [1.00000,  0.35294,  0.68627], 'LineWidth', 2); % epf2oe
plot([2.87e-5 2.87e-5], [1e-7 1.5e-4], 'Color', [1.00000,  0.43137,  0.22745], 'LineWidth', 2); % focl
plot([2.20e-5 2.20e-5], [1e-7 1.5e-4], 'Color', [0.00000,  0.98824,  0.81176], 'LineWidth', 2); % atml1
plot([1.07e-4 1.07e-4], [1e-7 1.5e-4], 'Color', [0.88627,  0.00392,  0.20392], 'LineWidth', 2); % re6
% water loss
rhovec = logspace(-7,-2);
nH2Otubevec = zeros(length(rhovec),1);
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k = sqrt(4*G_h/(D_h*rho));
    nH2Otubevec(rr) = (1/(rho+d)^2)*(pi*G_h*rho*phi_cell*Lz...
                        +((pi*G_h*(rho^1.5))/(sqrt(4*G_h/D_h)))*...
                        (phi_cell-phi_hLz-phi_h0)*(coth(k*Lz)-csch(k*Lz)));    
end
yyaxis right
loglog(rhovec,nH2Otubevec,'k','LineWidth',3,'MarkerSize',10, 'color',right_color);
set(gca, 'YDir','reverse')
%axis([rhovec(1) rhovec(end) 1e-9 1e-1]); axis square;
ylabel('\Gamma_{leaf},Water loss rate (\rho) (molm^{-2}s^{-1})','FontSize',18);

 ldg = legend('CO_{2} Uptake rate', 'Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6','Water loss rate');
 lgd.FontSize = 24;
 

%% Overlay on plot of CO2 uptake of n tubes