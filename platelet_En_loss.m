%% Function for calculating platelet sedimentation
function buoyancy_flux=platelet_En_loss(K,supercooling,size_cutoff,...
    bins,b)

c_i=2000;
Tm=273.15;
S_sw=34;
L=334774;
delta_T=1.853*(S_sw/28);
%supercooling=0.05;
H_s=(c_i*(Tm-delta_T));

%% Enthalphy Calculation for solidification in water column
Phi=((c_i*(Tm-delta_T-supercooling)+L)-H_s)/L;
new_ice=(1-Phi);

%% Population distibution of platelet sizes
%size_cutoff=0.1;
%bins=100;
%b=100;
R=linspace(0,size_cutoff,bins);
population_distribution=exp(-b*R);
a=new_ice/(sum(population_distribution)*(size_cutoff/bins));
mass_distribution=a*population_distribution;

% Plot of platelet size vs. mass distribution
%figure
%plot(R,mass_distribution)

%% Velocity distribution as a function of platelet size
% Stokes problem
g=10;
rho_br=1000;
rho_i=800;
delta_rho=rho_br-rho_i;
aspect_ratio=1/100;
drag_coeff=2;
%velocity=sqrt((g*delta_rho*aspect_ratio*R)/(drag_coeff*rho_br));
velocity=K*sqrt(R);

% plot of platelet size vs. rise velocity
%figure
%plot(R,velocity)

%% Platelet size vs. mass flux
%figure
%plot(R,mass_distribution.*velocity)

%% Total buoyancy mass flux in kg/s moving upward
buoyancy_flux=sum((mass_distribution.*velocity)*(size_cutoff/bins));




