%% Sensitivity Analysis prior to fitting/predictive modelling
close all;
clear all;

%% Some housekeeping 

fontname = 'Arial';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

%% Define constants - ALL FOR COL-0 HERE
Lz          = 91.3e-6;    % Tube height (m)
phi_c0      = 0.0164;   % CO2 conc at lower boundary (mol/m^3)
phi_cLz     = phi_c0;   % CO2 conc at upper boundary (mol/m^3)

D_c         = 0.139e-4; % CO2 diffusion constant (m^2/s)
G_c         = 0.00344;     % Uncorrected conductance coeff. G(m/s)
Lx          = 0.001;     % Leaf width (m)
d           = 3.296e-5;    % Minimum palisade cell width (m)?
%rhovec      = [10e-7; 10e-6; 10e-5; 10e-4];   % Tube diameter vector (m)
zvec        = (0:Lz/100:Lz)';                 % Height vector (m)
Am = Lx^2;


%% Assumption verifications
figure(5);

phi_cLz2    = phi_c0/2;
phi_cLz3    = phi_c0/4;
GGvec = logspace(-5, -1);
Gtubevec = zeros(length(GGvec),1);
rho = 23e-6;
for gg=1:length(GGvec)
    G_c = GGvec(gg);
    k = sqrt(4*G_c/(D_c*rho));
    Gtubevec(gg) = (pi*G_c*(rho^1.5)/((sqrt(rho))*k*((rho+d)^2)*sinh(k*Lz)))...
                *(phi_cLz*cosh(k*Lz)-phi_c0-phi_cLz+phi_c0*cosh(k*Lz));
            
    Gtubevec2(gg) = (pi*G_c*(rho^1.5)/((sqrt(rho))*k*((rho+d)^2)*sinh(k*Lz)))...
                *(phi_cLz2*cosh(k*Lz)-phi_c0-phi_cLz2+phi_c0*cosh(k*Lz));
            
    Gtubevec3(gg) = (pi*G_c*(rho^1.5)/((sqrt(rho))*k*((rho+d)^2)*sinh(k*Lz)))...
                *(phi_cLz3*cosh(k*Lz)-phi_c0-phi_cLz3+phi_c0*cosh(k*Lz));
end

h(1) = loglog(GGvec,Gtubevec,'k','LineWidth',2,'MarkerSize',10); hold on; grid minor;
h(2) = loglog(GGvec,Gtubevec2,'b','LineWidth',2,'MarkerSize',10);
h(3) = loglog(GGvec,Gtubevec3,'m','LineWidth',2,'MarkerSize',10);

xlabel('Conductance coefficient G_{CO_2} (ms^{-1})','FontSize',18);
ylabel('U^{leaf}_{CO_2}(G) (molm^{-2}s^{-1})','FontSize',18);
X = [GGvec(:,21) GGvec(:,22)];
Y = [Gtubevec(21,:) Gtubevec(22,:)];
h1 = area(X,Y,'LineStyle',':');
h1.FaceColor = [0 1 0];

XX = [GGvec(:,32) GGvec(:,33)];
YY = [Gtubevec(32,:) Gtubevec(33,:)];
h2 = area(XX,YY,'LineStyle',':');
h2.FaceColor = [1 0 0];

llgd=legend('100% CO_2','50% CO_2','25% CO_2', 'Desired G','Uncorrected G');
llgd.FontSize = 14;