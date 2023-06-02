%% Predictive modelling - parameter estimation
close all;
clear all;

fontname = 'Arial';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);

%% Define constants 
phi_c0      = 0.0164;   % CO2 conc at lower boundary (mol/m^3)
phi_cLz     = phi_c0;   % CO2 conc at upper boundary (mol/m^3)

D_c         = 0.139e-4; % CO2 diffusion constant (m^2/s)
G_c         = 0.00344;     % Uncorrected conductance coeff. G(m/s)
Lx          = 0.001;     % Leaf width (m)
%zvec        = (0:Lz/100:Lz)';                 % Height vector (m)
Am = Lx^2;
% for gm conversion
T           = convtemp(25,'C','K');
R           = 8.3144;
syms G;

%% Col-0
rhovec     = [22.35e-6 19.55e-6 22.27e-6 26.4e-6 33.43e-6];
rhoCol0 = rhovec;
dvec       = [35.56e-6 29.45e-6 29.74e-6 32.65e-6 37.39e-6];
A400vec    = [13.38e-6 14.16e-6 13.78e-6 15.04e-6 15.48e-6];
A400Col0 = A400vec;
Lzvec      = [82.5e-6 88e-6 112.75e-6 88e-6 85.25e-6];
LzvecCol0 = Lzvec;
real_gmvec = [1.41e-6 1.32e-6 1.33e-6 1.34e-6 1.53e-6];
gmCol0 = real_gmvec;
% calculate Eta_Col0
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:20
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_Col0(:,i) = mode(Gvec);
    % Now reverse eng gm
    %without pressure from original formula (Woodward) - sharkey tool accounts for it
    gm = (G_Col0(:,i))/(R*T);
    Eta_Col0(:,i) = gm/real_gm;
end

%% important for later
Eta_control = mean(Eta_Col0);
alpha_Col0 = Eta_Col0/Eta_control;

% lz vs alpha plot procedurally here:
figure(6); grid on; hold on;
s = scatter(Lzvec',alpha_Col0');
s.MarkerEdgeColor = [0.00000  0.55294  0.97647];
s.MarkerFaceColor = [0.00000  0.55294  0.97647];

% Correction of G and A for Col0
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_Col0(:,exclude_idx) = (mean(Eta_Col0(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_Col0(i);
    Ahat_Col0(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% arp3
rhovec     = [38.36e-6 30.02e-6 35.08e-6];
rhoarp3 = rhovec;
dvec       = [37.85e-6 30.07e-6 35.5e-6];
A400vec    = [13.77e-6 13.51e-6 13.42e-6];
A400arp3 = A400vec;
Lzvec      = [107.25e-6 121e-6 88e-6];
Lzvecarp3 = Lzvec;
real_gmvec = [1.22e-6 1.08e-6 1.15e-6];
gmarp3 = real_gmvec;

% calculate Eta
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:20
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_arp3(:,i) = mode(Gvec);
    % Now reverse eng gm
    gm = (G_arp3(:,i))/(R*T);
    Eta_arp3(:,i) = gm/real_gm;
end
alpha_arp3 = Eta_arp3/Eta_control;

% lz vs alpha plot procedurally here:
s = scatter(Lzvec',alpha_arp3');
s.MarkerEdgeColor = [0.51765,  0.00000,  0.80392];
s.MarkerFaceColor = [0.51765,  0.00000,  0.80392];

% Correction of G and A for arp3
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_arp3(:,exclude_idx) = (mean(Eta_arp3(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_arp3(i);
    Ahat_arp3(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% qua2
rhovec     = [36.87e-6 35.71e-6 36.83e-6 38.56e-6];
rhoqua2 = rhovec;
dvec       = [40.15e-6 37.87e-6 42.33e-6 40.68e-6];
A400vec    = [16.54e-6 17.63e-6 18.72e-6 15.85e-6];
A400qua2 = A400vec;
Lzvec      = [96.25e-6 93.5e-6 85.25e-6 82.5e-6];
Lzvecqua2 = Lzvec;
real_gmvec = [1.96e-6 1.79e-6 1.73e-6 1.79e-6];
gmqua2 = real_gmvec;

% calculate Eta
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:40
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_qua2(:,i) = mode(Gvec);
    % Now reverse eng gm
    gm = (G_qua2(:,i))/(R*T);
    Eta_qua2(:,i) = gm/real_gm;
end
alpha_qua2 = Eta_qua2/Eta_control;

% lz vs alpha plot procedurally here:
s = scatter(Lzvec',alpha_qua2');
s.MarkerEdgeColor = [0.00000,  0.62353,  0.50588];
s.MarkerFaceColor = [0.00000,  0.62353,  0.50588];

% Correction of G and A for qua2
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_qua2(:,exclude_idx) = (mean(Eta_qua2(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_qua2(i);
    Ahat_qua2(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% epf2oe
rhovec     = [25.55e-6 20.49e-6 28.47e-6];
rhoepf = rhovec;
dvec       = [45.35e-6 25.34e-6 27.53e-6];
A400vec    = [13.16e-6 12.55e-6 15.45e-6];
A400epf = A400vec;
Lzvec      = [96.25e-6 121e-6 99e-6];
Lzvecepf = Lzvec;
real_gmvec = [1.42e-6 1.31e-6 1.65e-6];
gmepf = real_gmvec;

% calculate Eta
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:20
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_epf(:,i) = mode(Gvec);
    % Now reverse eng gm
    gm = (G_epf(:,i))/(R*T);
    Eta_epf(:,i) = gm/real_gm;
end
alpha_epf = Eta_epf/Eta_control;

% lz vs alpha plot procedurally here:
s = scatter(Lzvec',alpha_epf');
s.MarkerEdgeColor = [1.00000,  0.35294,  0.68627];
s.MarkerFaceColor = [1.00000,  0.35294,  0.68627];

% Correction of G and A for epf
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_epf(:,exclude_idx) = (mean(Eta_epf(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_epf(i);
    Ahat_epf(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% focl
rhovec     = [29.76e-6 25.77e-6 29.83e-6 29.66e-6 28.46e-6];
rhofocl = rhovec;
dvec       = [39.26e-6 35.05e-6 42.93e-6 37.32e-6 38.47e-6];
A400vec    = [12.37e-6 9.07e-6 8.83e-6 11.44e-6 12.92e-6];
A400focl = A400vec;
Lzvec      = [143e-6 145.75e-6 123.75e-6 143e-6 123.75e-6];
Lzvecfocl = Lzvec;
real_gmvec = [1.25e-6 1.08e-6 1.12e-6 1.23e-6 1.24e-6];
gmfocl = real_gmvec;

% calculate Eta
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:40
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_focl(:,i) = mode(Gvec);
    % Now reverse eng gm
    gm = (G_focl(:,i))/(R*T);
    Eta_focl(:,i) = gm/real_gm;
end
alpha_focl= Eta_focl/Eta_control;

% lz vs alpha plot procedurally here:
s = scatter(Lzvec',alpha_focl');
s.MarkerEdgeColor = [1.00000,  0.43137,  0.22745];
s.MarkerFaceColor = [1.00000,  0.43137,  0.22745];

% Correction of G and A for focl
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_focl(:,exclude_idx) = (mean(Eta_focl(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_focl(i);
    Ahat_focl(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% ATML1
rhovec     = [22.79e-6 22.5e-6 22.43e-6 20.11e-6];
rhoL1 = rhovec;
dvec       = [34.11e-6 24.5e-6 34.48e-6 36.53e-6];
A400vec    = [14.68e-6 16.18e-6 14.84e-6 15.15e-6];
A400L1 = A400vec;
Lzvec      = [156.75e-6 176e-6 167.75e-6 134.75e-6];
LzvecL1 = Lzvec;
real_gmvec = [1.57e-6 1.61e-6 1.77e-6 1.79e-6];
gmL1 = real_gmvec;

% calculate Eta
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:20
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_L1(:,i) = mode(Gvec);
    % Now reverse eng gm
    gm = (G_L1(:,i))/(R*T);
    Eta_L1(:,i) = gm/real_gm;
end
alpha_L1 = Eta_L1/Eta_control;
% lz vs alpha plot procedurally here:
s = scatter(Lzvec',alpha_L1');
s.MarkerEdgeColor = [0.00000,  0.98824,  0.81176];
s.MarkerFaceColor = [0.00000,  0.98824,  0.81176];

% Correction of G and A for L1
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_L1(:,exclude_idx) = (mean(Eta_L1(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_L1(i);
    Ahat_L1(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% re6
rhovec     = [170.5e-6 176e-6 181.5e-6];
rhore6 = rhovec;
dvec       = [61.21e-6 48.03e-6 72.61e-6];
A400vec    = [10.25e-6 12.26e-6 11.21e-6];
A400re6 = A400vec;
Lzvec      = [170.5e-6 176e-6 181.5e-6];
Lzvecre6 = Lzvec;
real_gmvec = [1.18e-6 1.45e-6 1.42e-6];
gmre6 = real_gmvec;

% calculate Eta
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    A400 = A400vec(i);
    Lz = Lzvec(i);
    real_gm = real_gmvec(i);
    for n=1:40
     vpasolve(((pi*G*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*G/(D_c*rho))*...
         ((rho+d)^2)*sinh(sqrt(4*G/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*G/(D_c*rho))...
         *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*G/(D_c*rho))*Lz)))/Am == A400, G, 'Random',true);
     double(ans);
     Gvec(:,n) = ans;
    end
    G_re6(:,i) = mode(Gvec);
    % Now reverse eng gm
    gm = (G_re6(:,i))/(R*T);
    Eta_re6(:,i) = gm/real_gm;
end
alpha_re6 = Eta_re6/Eta_control;

% lz vs alpha plot procedurally here:
s = scatter(Lzvec',alpha_re6');
s.MarkerEdgeColor = [0.88627,  0.00392,  0.20392];
s.MarkerFaceColor = [0.88627,  0.00392,  0.20392];
ldg = legend('Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6');
lgd.FontSize = 24;
xlabel('Air channel depth L_{z} (m)','FontSize',18);
ylabel('Deviation factor \alpha','FontSize',18);

% Correction of G and A for re6
for exclude_idx=1:size(rhovec,2)
    desired_idx = [1:(exclude_idx -1), (exclude_idx+1):size(rhovec,2)]; 
    Ghat_re6(:,exclude_idx) = (mean(Eta_re6(:,desired_idx)))*real_gmvec(:,exclude_idx)*R*T;
end
for i=1:length(rhovec)
    rho = rhovec(i);
    d = dvec(i);
    Lz = Lzvec(i);
    GG = Ghat_re6(i);
    Ahat_re6(i) = ((pi*GG*(Lx^2)*(rho^1.5)/((sqrt(rho))*sqrt(4*GG/(D_c*rho))*...
     ((rho+d)^2)*sinh(sqrt(4*GG/(D_c*rho))*Lz)))*(phi_cLz*cosh(sqrt(4*GG/(D_c*rho))...
     *Lz)-phi_c0-phi_cLz+phi_c0*cosh(sqrt(4*GG/(D_c*rho))*Lz)))/Am;
end
%% Plotting
figure(7); grid on; hold on; axis square;

s = scatter(Ahat_Col0',A400Col0');
s.MarkerEdgeColor = [0.00000  0.55294  0.97647];
s.MarkerFaceColor = [0.00000  0.55294  0.97647];

s = scatter(Ahat_arp3',A400arp3');
s.MarkerEdgeColor = [0.51765,  0.00000,  0.80392];
s.MarkerFaceColor = [0.51765,  0.00000,  0.80392];

s = scatter(Ahat_qua2',A400qua2');
s.MarkerEdgeColor = [0.00000,  0.62353,  0.50588];
s.MarkerFaceColor = [0.00000,  0.62353,  0.50588];

s = scatter(Ahat_epf',A400epf');
s.MarkerEdgeColor = [1.00000,  0.35294,  0.68627];
s.MarkerFaceColor = [1.00000,  0.35294,  0.68627];

s = scatter(Ahat_focl',A400focl');
s.MarkerEdgeColor = [1.00000,  0.43137,  0.22745];
s.MarkerFaceColor = [1.00000,  0.43137,  0.22745];

s = scatter(Ahat_L1',A400L1');
s.MarkerEdgeColor = [0.00000,  0.98824,  0.81176];
s.MarkerFaceColor = [0.00000,  0.98824,  0.81176];

s = scatter(Ahat_re6',A400re6');
s.MarkerEdgeColor = [0.88627,  0.00392,  0.20392];
s.MarkerFaceColor = [0.88627,  0.00392,  0.20392];

xlabel('U^{leaf}_{CO_2} (molm^{-2}s^{-1}), Corrected','FontSize',14);
ylabel('A_{400}, Measured CO_{2} Assimilation Rate (molm^{-2}s^{-1})','FontSize',14);
refline(1,0);
axis([0.5e-5 3.5e-5 0.6e-5 2.25e-5])
axis normal;
ldg = legend( 'Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6', 'Reference line');
lgd.FontSize = 24;

%% Plot to recheck main results from Chapter 3

% CO2 concentration, revisited
% Define constants - ALL FOR COL-0 HERE
Lz          = 91.3e-6;    % Tube height (m)
phi_c0      = 0.0164;   % CO2 conc at lower boundary (mol/m^3)
phi_cLz     = phi_c0;   % CO2 conc at upper boundary (mol/m^3)

D_c         = 0.139e-4; % CO2 diffusion constant (m^2/s)
G_c         = mean(Ghat_Col0);     % corrected conductance coeff. G(m/s)
Lx          = 0.001;     % Leaf width (m)
d           = 3.2958e-05;    % Minimum palisade cell width (m)?
rhovec      = [2.4800e-05];   % Tube diameter vector (m)
% mean rho 2.4800e-05
zvec        = (0:Lz/100:Lz)';                 % Height vector (m)
Am = Lx^2;

%% CO2 uptake vs diameters, for 2D with cell thickness d
figure(8);
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
axis([rhovec(1) rhovec(end) 1e-6 2.7e-5]); 
xlabel('Air channel diameter \rho (m)','FontSize',18);
ylabel('U^{leaf}_{CO_2} (\rho) (molm^{-2}s^{-1}), Corrected','FontSize',18);
 hold on;
 %xline(1.8e-5, '-g'); xline(8e-5, '-g');
plot(mean(rhoCol0), mean(A400Col0), 'p', 'MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [0.00000  0.55294  0.97647]); %Col0
plot(mean(rhoarp3), mean(A400arp3),'p', 'MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [0.51765,  0.00000,  0.80392]); %Arp3
plot(mean(rhoqua2), mean(A400qua2), 'p','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [0.00000,  0.62353,  0.50588]); %qua2
plot(mean(rhoepf), mean(A400epf), 'p','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [1.00000,  0.35294,  0.68627]); % epf2oe
plot(mean(rhofocl), mean(A400focl), 'p','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [1.00000,  0.43137,  0.22745]); % focl
plot(mean(rhoL1), mean(A400L1), 'p','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [0.00000,  0.98824,  0.81176]); % atml1
plot(mean(rhore6), mean(A400re6), 'p','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor', [0.88627,  0.00392,  0.20392]); % re6
 
ldg = legend('Model output', 'Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6');
lgd.FontSize = 24;
axis normal;

figure(9)
rhovec      = [2.4800e-05]; 
G_c         = 0.00042; 
for rr=1:length(rhovec)
    rho = rhovec(rr);
    k_c = sqrt(4*G_c/(D_c*rho));
    phivec_c = (phi_cLz*sinh(k_c.*zvec)+phi_c0*sinh(k_c.*(Lz-zvec)))/(sinh(k_c*Lz));
    plot(phivec_c,zvec,'LineWidth',2); grid on; hold on;
    axis([0.0162 phi_c0 0 Lz]); axis square;
end
xlabel('\phi_{CO_2}(z), mol/m^{3}','FontSize',18);
ylabel('z (m)','FontSize',18);

%% gm vs Lz
figure(10); hold on; grid on;
subplot(1,2,1)
hold on; grid on;
s = scatter(gmCol0',LzvecCol0');
s.MarkerEdgeColor = [0.00000  0.55294  0.97647];
s.MarkerFaceColor = [0.00000  0.55294  0.97647];

s = scatter(gmarp3',Lzvecarp3');
s.MarkerEdgeColor = [0.51765,  0.00000,  0.80392];
s.MarkerFaceColor = [0.51765,  0.00000,  0.80392];

s = scatter(gmqua2',Lzvecqua2');
s.MarkerEdgeColor = [0.00000,  0.62353,  0.50588];
s.MarkerFaceColor = [0.00000,  0.62353,  0.50588];

s = scatter(gmepf',Lzvecepf');
s.MarkerEdgeColor = [1.00000,  0.35294,  0.68627];
s.MarkerFaceColor = [1.00000,  0.35294,  0.68627];

s = scatter(gmfocl',Lzvecfocl');
s.MarkerEdgeColor = [1.00000,  0.43137,  0.22745];
s.MarkerFaceColor = [1.00000,  0.43137,  0.22745];

s = scatter(gmL1',LzvecL1');
s.MarkerEdgeColor = [0.00000,  0.98824,  0.81176];
s.MarkerFaceColor = [0.00000,  0.98824,  0.81176];

s = scatter(gmre6',Lzvecre6');
s.MarkerEdgeColor = [0.88627,  0.00392,  0.20392];
s.MarkerFaceColor = [0.88627,  0.00392,  0.20392];

xlabel('Mesophyll conductance, g_m (µmolm^{-2}s^{-1}Pa^{-1})','FontSize',14);
ylabel('Air channel depth L_z, (m)','FontSize',14);
ldg = legend('Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6');
lgd.FontSize = 24;
axis square;

subplot(1,2,2)
hold on; grid on;
s = scatter(gmCol0',rhoCol0');
s.MarkerEdgeColor = [0.00000  0.55294  0.97647];
s.MarkerFaceColor = [0.00000  0.55294  0.97647];

s = scatter(gmarp3',rhoarp3');
s.MarkerEdgeColor = [0.51765,  0.00000,  0.80392];
s.MarkerFaceColor = [0.51765,  0.00000,  0.80392];

s = scatter(gmqua2',rhoqua2');
s.MarkerEdgeColor = [0.00000,  0.62353,  0.50588];
s.MarkerFaceColor = [0.00000,  0.62353,  0.50588];

s = scatter(gmepf',rhoepf');
s.MarkerEdgeColor = [1.00000,  0.35294,  0.68627];
s.MarkerFaceColor = [1.00000,  0.35294,  0.68627];

s = scatter(gmfocl',rhofocl');
s.MarkerEdgeColor = [1.00000,  0.43137,  0.22745];
s.MarkerFaceColor = [1.00000,  0.43137,  0.22745];

s = scatter(gmL1',rhoL1');
s.MarkerEdgeColor = [0.00000,  0.98824,  0.81176];
s.MarkerFaceColor = [0.00000,  0.98824,  0.81176];

s = scatter(gmre6',rhore6');
s.MarkerEdgeColor = [0.88627,  0.00392,  0.20392];
s.MarkerFaceColor = [0.88627,  0.00392,  0.20392];

xlabel('Mesophyll conductance, g_m (µmolm^{-2}s^{-1}Pa^{-1})','FontSize',14);
ylabel('Air channel diameter, (m)','FontSize',14);
ldg = legend('Col-0', 'arp3', 'qua2', 'EPF2-OE','focl1-1','ATML1_{pro}:KRP1', 're6');
lgd.FontSize = 24;
axis square;

%% Error calculations

errCol0 = abs(A400Col0./Ahat_Col0);
MerrCol0 = mean(errCol0)
errarp3 = mean(abs(A400arp3./Ahat_arp3));
errqua2 = mean(abs(A400qua2./Ahat_qua2));
errepf= mean(abs(A400epf./Ahat_epf));
errfocl = mean(abs(A400focl./Ahat_focl));
errL1 = mean(abs(A400L1./Ahat_L1));
errre6 = mean(abs(A400re6./Ahat_re6));


% mean alphas
Malpha_Col0= mean(alpha_Col0);
Malpha_arp3= mean(alpha_arp3);
Malpha_qua2= mean(alpha_qua2);
Malpha_epf= mean(alpha_epf);
Malpha_focl= mean(alpha_focl);
Malpha_L1= mean(alpha_L1);
Malpha_re6= mean(alpha_re6);

MAhat_Col0= mean(Ahat_Col0);
MAhat_arp3= mean(Ahat_arp3);
MAhat_qua2= mean(Ahat_qua2);
MAhat_epf= mean(Ahat_epf);
MAhat_focl= mean(Ahat_focl);
MAhat_L1= mean(Ahat_L1);
MAhat_re6= mean(Ahat_re6);

