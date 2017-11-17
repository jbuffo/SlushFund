
%% Vectors for each time chunk (1000000) of iteration for plotting
Tplot=[];
Splot=[];
phiplot=[];
wplot=[];
Stotplot=[];
plateletplot=[];
SinIplot=[];

%% Time loops
time=[0:1000:10000];
for k=1:length(time)-1;

%% Input Parameters 
Ttop=250;                %% Ice shelf top temp
Tbottom=271.25;             %% '       ' bottom temp (-1.95C)
Tm=273.15;               %% Melt temperature of pure ice
k_i=2;                   %% Diffusion Coefficient for Heat (ice)
k_br=.6;                 %% '                            ' (brine)
k_s=2*10^-9;                   %% Diffusion Coefficient for Salt
dt=50;                    %% Time step
dz=.01;                   %% Spatial step
H=1.5;                   %% Thickness of ice shelf
c_br=3985;               %% Specific heat of seawater
c_i=2000;                %% '              ' ice
L=334774;                %% Latent heat of fusion ice<->water
S_top=100;               %% Salinity at top of shelf
S_bottom=34;             %% '         ' bottome of shelf
rho_i=917;               %% Density of Ice
rho_br=1028;             %% Density of Brine
t_start=time(k);
t_end=time(k+1);
tf=t_end-t_start;                %% Final time
Tol=.1;                  %% Error Tolerance (includes T&Phi&S)
m=2;                     %% Cementation exponent
g=10;                    %% Gravity
k_ibr=100;               %% Partition Coefficient
mu=1.88*10^-3;           %% Viscosity

Depth=[0:-1:-(H/dz)];
Precipitated_Salt=0*[0:dz:H]';

%% Sedimentation Parameters
phi_crit=0.7;            %% Critical porosity for sedimentation
a_s=1/50;                %% Platelet aspect ration
c_d=2;                   %% Drag coeff.
supercooling=0.04;       %% Supercooling in water column
size_cutoff=0.1;         %% Platelet size cutoff
bins=100;                %% Number of bins for platelet size dist.
b=100;                   %% Coefficient for population dist.
K=sqrt((g*(rho_br-rho_i)*a_s)/(c_d*rho_br)); %% Variational parameter for Stokes velocity calculation

%K_sweep=[0.075,0.091104,0.1,0.15];
%for m=2;
%    K=K_sweep(m);
buoyancy_flux=platelet_En_loss(K,supercooling,size_cutoff,bins,b);

%% Freeze over date since real world BC starts @ April ~April 15th (June 3=49,Jul 25=101,Sept 11=149, Sept 18=156, Apr 21=6)
freeze_date=(86400*149)/dt;

%% Real world boundary condition forcing
t=[-2000000:dt:31557600];
TtopBC=-9*sin((2*pi/31557600)*t)-18.15+273.15;
for i=1:length(TtopBC)
    if TtopBC(i)<250
        TtopBC(i)=250;
    %elseif TtopBC(i)>259.61
    %    TtopBC(i)=259.0;
    else
        TtopBC(i)=TtopBC(i);
    end
end
TbottomBC=0.9*heaviside(t-15778800)+271.25+0.05*sin((2*pi/89400)*t)+0.025*sin((2*pi/86400)*t);
TbottomBCmod=0.9*heaviside(t-15778800)+271.25+0.05*sin((2*pi/89400)*t)+0.025*sin((2*pi/86400)*t);
SC=0*[t_start:dt:31557600];
for i=1:length(TbottomBC);
    if TbottomBC(i)<271.25;
        SC(i)=(271.25-TbottomBC(i));
        TbottomBCmod(i)=271.25;
    else
        TbottomBCmod(i)=TbottomBC(i);
    end
end
if t_start>=15778800;
TbottomBCmod(1:end)=272.15;
else
TbottomBCmod=TbottomBCmod;
end

%% Initial Condition Vectors for Phi and T and w and S and Platelet Fraction

if t_start==0
    T_initial=0*[0:dz:H]+Tbottom;
    S_initial=0*[0:dz:H]+S_bottom;
    phi_initial=0*[0:dz:H]+1;
    w_initial=0*[0:dz:H];
    Platelet_track=0*[0:dz:H];
    S_ice=0*[0:dz:H];
else
    T_initial=IC(:,1)';
    S_initial=IC(:,2)';
    phi_initial=IC(:,3)';
    w_initial=IC(:,4)';
    Platelet_track=IC(:,5)';
    S_ice=IC(:,6)';
end

%% Vectors for plotting
Temperature=[];
Thold=[];
Liquid_Fraction=[];
phihold=[];
Salinity=[];
Shold=[];
Darcy_Velocity=[];
whold=[];
Temp_mov=[];
Platelet_frac=[];
phold=[];
Salt_in_Ice=[S_ice'];
sinihold=[];

%% Looping Over Time
for n=1:tf/dt
    if n==1;
        %[NewTemp,NewPhi,NewS,Neww]=...
        %    ESM_project_T_phi_S_advection_bars2(k_i,k_br,...
        %    k_s,dt,dz,H,c_i,c_br,L,Tol,phi_initial,T_initial,...
        %    S_initial,S_bottom,rho_i,rho_br,Tm,m,w_initial,n,Ttop,Tbottom);
        [NewTemp,NewPhi,NewS,Neww,IceSalt]=one_D_adv_ARD_w_fix(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,Tol,phi_initial,...
        T_initial,S_initial,S_bottom,rho_i,...
        rho_br,Tm,m,w_initial,n,TtopBC(n+freeze_date),Tbottom,k_ibr,S_ice);
        NewTemp(1)=TtopBC(n+freeze_date);
        NewTemp(length(NewTemp))=Tbottom;
        %NewS(1)=0;
        %NewS(length(NewS))=S_bottom;
        Temperature=NewTemp;
        Thold=[Thold NewTemp];
        Liquid_Fraction=NewPhi;
        phihold=[phihold NewPhi];
        Salinity=NewS;
        Shold=[Shold NewS];
        Darcy_Velocity=Neww;
        whold=[whold Neww];
        Salt_in_Ice=Salt_in_Ice+IceSalt';
        sinihold=[sinihold Salt_in_Ice+IceSalt'];
    else
        [NewTemp,NewPhi,NewS,Neww,IceSalt]=one_D_adv_ARD_w_fix(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,Tol,Liquid_Fraction',...
        Temperature',Salinity',S_bottom,rho_i,...
        rho_br,Tm,m,Darcy_Velocity',n,TtopBC(n+freeze_date),Tbottom,k_ibr,Salt_in_Ice');
        NewTemp(1)=TtopBC(n+freeze_date);
        NewTemp(length(NewTemp))=Tbottom;
        %NewS(1)=0;
        %NewS(length(NewS))=S_bottom;
        count=0;
        for i=1:H/dz;
            if NewPhi(i)<1;
                count=count+1;
            else
                break
            end
        end
        solid_frac_flux=buoyancy_flux*dt;
        if NewPhi(count)>phi_crit;
            NewPhi(count)=NewPhi(count)-solid_frac_flux;
            Platelet_track(count)=Platelet_track(count)+solid_frac_flux;
        else
            NewPhi(count+1)=NewPhi(count+1)-solid_frac_flux;
            Platelet_track(count+1)=Platelet_track(count+1)+solid_frac_flux;
        end
        
        %for i=1:(H/dz)+1;
        %    if NewPhi(i)>0.9;
        %        NewTemp(i)=Tm-(1.853*(34/28))-0.05;
        %        NewS(i)=34;
        %    else
        %        NewTemp(i)=NewTemp(i);
        %    end
        %end
        Temperature=NewTemp;
        Thold=[Thold NewTemp];
        Liquid_Fraction=NewPhi;
        phihold=[phihold NewPhi];
        Salinity=NewS;
        Shold=[Shold NewS];
        Darcy_Velocity=Neww;
        whold=[whold Neww];
        Platelet_frac=Platelet_track';
        phold=[phold Platelet_track'];
        Salt_in_Ice=Salt_in_Ice+IceSalt';
        sinihold=[sinihold Salt_in_Ice+IceSalt'];

        if n==500 || n==1000 || n==1500 || n==2000 || n==4000 || n==6000 || n==8000 || n==10000 || n==12000 ...
            || n==14000 || n==16000 || n==18000;
            disp(n)
            datestr(now)
        else
        end
        %% Movie creation (comment out unless making movies)
        %set(gcf,'visible','off')
        %fig=figure;
        %    subplot(1,5,1);
        %    image(NewTemp,'CDataMapping','scaled');
        %    title('Temperature');
        %    colorbar;
        %    subplot(1,5,2);
        %    image(NewPhi,'CDataMapping','scaled');
        %    title('Liquid Fraction');
        %    colorbar;
        %    subplot(1,5,3);
        %    image(NewS,'CDataMapping','scaled');
        %    title('Salinity');
        %    colorbar;
        %    subplot(1,5,4);
        %    image(Neww,'CDataMapping','scaled');
        %    title('Brine Velocity');
        %    colorbar;
        %    subplot(1,5,5);
        %    image(Salt_in_Ice(:,n-1)+IceSalt,'CDataMapping','scaled');
        %    title('Salt in Ice');
        %    colorbar;
        %temp=getframe(fig);
        %Temp_mov=[Temp_mov temp];
        %% End of movie making code
        
    end
end

Tplot=[Tplot Temperature];
Splot=[Splot Salinity];
phiplot=[phiplot Liquid_Fraction];
wplot=[wplot Darcy_Velocity];
Stotplot=[Stotplot Liquid_Fraction.*Salinity];
plateletplot=[plateletplot Platelet_frac];
SinIplot=[SinIplot Salt_in_Ice];

%% Output Images
%figure
%subplot(2,4,1)
%image(Temperature,'CDataMapping','scaled')
%title('Temperature')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,2)
%image(Liquid_Fraction,'CDataMapping','scaled')
%title('Liquid Fraction')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,3)
%image(Salinity,'CDataMapping','scaled')
%title('Salinity')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,4)
%image(Darcy_Velocity(1:end,10:end),'CDataMapping','scaled')
%title('Darcy Velocity')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,5)
%image(Liquid_Fraction.*Salinity,'CDataMapping','scaled')
%title('Total Salt (ppt)')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,6)
%image(Platelet_frac,'CDataMapping','scaled')
%title('Platelet Ice Fraction')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,7)
%plot(Platelet_frac(:,end),Depth,0*Depth+0.1,Depth)
%title('Platelet Fraction vs. Depth')
%xlabel('Platelet Fraction')
%ylabel('Depth (cm)')
%subplot(2,4,8)
%text(0.0,0.0,['Supercooled ' num2str(supercooling) 'K'])
%text(0.0,0.1,['Ttop= ' num2str(Ttop) ' Kelvin'])
%text(0.0,0.2,['Tbottom= ' num2str(Tbottom) ' Kelvin'])
%text(0.0,0.3,['tf= ' num2str(t_end/86400) ' days'])
%text(0.0,0.4,['K= ' num2str(K)])
%text(0.0,0.5,['dt= ' num2str(dt) ' sec'])
%axis off

%implay(Temp_mov)
%movie2avi(Temp_mov,'ice_shelf.avi')

%save(strcat('ice_core_w_fix_Emod_250_full_density_plot_sc_',...
%    num2str(supercooling),'_freeze_day_',num2str(freeze_date),'_',num2str(t_start),'_to_',num2str(t_end),'_sec_K_',num2str(K),'_dt_',...
%    num2str(dt),'Ra_c_1_exp_prob_dist_phi_c_',num2str(phi_crit),'_sizecutoff_',...
%    num2str(size_cutoff),'_c_d_',num2str(c_d),'_a_s_',num2str(a_s),'.mat'),'Temperature',...
%    'Salinity','Liquid_Fraction','Darcy_Velocity','Platelet_frac','Salt_in_Ice','Thold','Shold',...
%    'phihold','whold','phold','sinihold')

save(strcat('ice_core_W3_with_sc_',num2str(t_start),'_to_',num2str(t_end),'.mat'),'Temperature',...
    'Salinity','Liquid_Fraction','Darcy_Velocity','Platelet_frac','Salt_in_Ice','Thold','Shold',...
    'phihold','whold','phold','sinihold')

%% Initial Conditions for next time loop
IC=[Temperature Salinity Liquid_Fraction Darcy_Velocity Platelet_frac Salt_in_Ice];
end

figure
subplot(2,4,1)
plot(Tplot,Depth)
title('Temperature vs. Depth over time')
xlabel('Temperature (K)')
ylabel('Depth (cm)')
subplot(2,4,2)
plot(phiplot,Depth)
title('Liquid Fraction vs. Depth over time')
xlabel('Liquid Fraction')
ylabel('Depth (cm)')
subplot(2,4,3)
plot(Splot,Depth)
title('Salinity vs. Depth over time')
xlabel('Salinity (psu)')
ylabel('Depth (cm)')
subplot(2,4,4)
plot(wplot,Depth)
title('Darcy Velocity vs. Depth over time')
xlabel('Darcy Velocity (m/s)')
ylabel('Depth (cm)')
subplot(2,4,5)
plot(SinIplot,Depth)
title('Salt in Ice vs. Depth over time')
xlabel('Salt in Ice (ppt)')
ylabel('Depth (cm)')
subplot(2,4,6)
plot(phiplot.*Splot+(1-phiplot).*SinIplot,Depth)
title('Total Salt vs. Depth over time')
xlabel('Total Salt (psu)')
ylabel('Depth (cm)')
subplot(2,4,7)
plot(plateletplot,Depth,0*Depth+0.1,Depth)
title('Platelet Ice Fraction vs. Depth over time')
xlabel('Platelet Ice Fraction')
ylabel('Depth (cm)')
subplot(2,4,8)
text(0.0,0.0,['Supercooled ' num2str(supercooling) 'K'])
%text(0.0,0.1,['Ttop= ' num2str(Ttop) ' Kelvin'])
%text(0.0,0.2,['Tbottom= ' num2str(Tbottom) ' Kelvin'])
text(0.0,0.3,['tf= ' num2str(t_end/86400) ' days'])
text(0.0,0.4,['K= ' num2str(K)])
text(0.0,0.5,['dt= ' num2str(dt) ' sec'])
axis off
