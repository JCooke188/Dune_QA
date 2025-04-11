% Created by Justin Cooke
% Purpose of this script is to view the energy spectra as the flow marches
% into the dune field

%% Start

clc;
clear;
close all;

%% Load Data

myDir = dir('./Data/Sept13/x*');

Ndir = length(myDir);

ixyz = 0;
iu = 0;

for i = 1:Ndir
    

    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    myPath = strcat(myFolder,'/',myName);
    
    data = load(myPath);
    
    if contains(myPath,'README')
        ixyz = ixyz + 1;
        x(ixyz) = round(data(1,2));
        z_local{ixyz} = data(:,end);
    elseif contains(myPath,'comp(u,0)')
        iu = iu + 1;
        u{iu} = data(:,4:end);
    end
        
end

clear my* Ndir i* data

myDir = dir('../Amplitude Modulation/DuneField/WSS/WallStress*');

surf_N = length(myDir);

for i = 1:surf_N
   
    myName = myDir(i).name;
    myFolder = myDir(i).folder;
    
    loadMe = strcat(myFolder,'/',myName);
    
    if contains(loadMe,'README')
        surf_xyz = load(loadMe);
        surf_xyz = surf_xyz(:,2:end);
    else
        surf_tau = load(loadMe);
        surf_tau = surf_tau(:,4:end);
    end
    
    
end

clear myName myFolder loadMe



%% Create a standard global z

nu = 0.32e-4;
delta_abl = 300;
u_tau_dune = 0.14;
u_tau_flat = 0.12;
delta_nu_dune = nu/u_tau_dune;
delta_nu_flat = nu/u_tau_flat;

z_global = linspace(1.5,400,100);
zplus_dunes = z_global./delta_nu_dune;
zplus_flat = z_global./delta_nu_flat;
zdelta_global = z_global./delta_abl;

delta_ibl = [0.1,0.06356,0.1665,0.1665,0.1373,0.1372,0.2968,0.4362,...
    0.4362,0.5289].*delta_abl;

%% Create Mean and Fluctuating Quantities 

N_u = length(u);

for i = 1:N_u

    tempu = u{i};
    U{i} = mean(tempu);
    up{i} = tempu - U{i};
    urms{i} = rms(up{i});
    
end
    
clear temp*

%% FFT Stuff

[Nt,Nz] = size(u{1});
dt = 0.65;
Lt = dt*Nt;
dw  = 2*pi/Lt;
n_w  = -Nt/2 : 1 : Nt/2 - 1;
omega = n_w.*dw;

t = linspace(0,Lt,Nt);
T = t.*(U{end}(end)/delta_abl);
total_T = T(end);

%% Do FFT

% 1 = zdelta 0.01 = 3m 
% 5 = zdelta 0.05 = 15m
% 8 = zdelta 0.1 = 30m 
% 23 = zdelta 0.3 = 90m
% 38 = zdelta 0.5 = 150m

zselection = [1,5,8,23,38];
NzSelect = length(zselection);

for i = 1:N_u
    for j = 1:NzSelect
    
      up_temp = up{i}(:,zselection(j));

      u_hat_w = (1/Nt) * fftshift(fft(up_temp));
      psi_w = abs(u_hat_w).^2;
      E_hat_w{i,j} = psi_w ./ (dw);
         
    % Pre-multiplied Energy Spectra
      E_w{i,j} = E_hat_w{i,j}.*omega';
    end 
      
end

clear urms_temp psi_w 

%% Calculate maximum wavenumber allowable

kmax = [4.189,2.0944,1.0472,0.5236,0.2618];
gridSize = [0.75,1.50,3.0,6.0];

for i = 1:N_u
    for j = 1:NzSelect
        Ubar_infty{i}(j) = U{i}(zselection(j));
        w_cutoff{i}(j) = Ubar_infty{i}(j)*kmax(j)/(2*pi());
        ii = find(abs(omega) <= w_cutoff{i}(j),1);
        myCutoff{i}(j) = ii;
    end
end



%% Plot Spectrum at diff Z elevations and confirm slope in inertial range

dx_extra = logspace(0,-1,4);

dx_1 = dx_extra.^-2;
dx_2 = dx_extra.^(-5/3);

close all

set(0,'defaultTextInterpreter','latex');

for j = 1:NzSelect-1
figure(j);
    aTitle = strcat ('$z = $',num2str(zselection(j))); 
    for i = 1:N_u
        nexttile;
        loglog(omega((end/2):end-myCutoff{i}(j)),...
            smoothdata(E_hat_w{i,j}((end/2):end-myCutoff{i}(j)))); hold on;
        loglog(dx_extra,dx_2/800,'k-','LineWidth',2);
        myTitle = strcat ('$\hat{x} = $',num2str(i-2));
        title(myTitle);
        xlabel('$\omega$');
        ylabel('$E(\omega)$');
    end
end

%% IBL Freq.

% omega_ibl = Uc / delta_i

% Find Uc

N_ibl = length(delta_ibl);
U_infty_ibl = zeros(9,1);
omega_ibl = zeros(9,1);

for i = 1:N_ibl
    ii = find(z_global >= delta_ibl(i),1);
    
    U_infty_ibl(i) = U{i+1}(ii);
    omega_ibl(i) = U_infty_ibl(i) / delta_ibl(i);   
    
end

omega_asl = U{1}(8) / 30;

%% Now plot at fixed x as a function of z

close all;

zColors = ["#2a4858","#0c9bba","#44d6a9","#fafa70"];
newColors = ["#2a4858","#2436db","#fa007f","#ff9a00"];



for j = 2:N_u-1
 figure(j);
%nexttile;
for i = 1:NzSelect-1

    loglog(omega((end/2):end-myCutoff{j}(i)),...
        smoothdata(E_hat_w{j,i}((end/2):end-myCutoff{j}(i))),...
        'Color',newColors(i),'LineWidth',2); hold on;
    set(gca,'FontName','SansSerif','FontSize',14);
    xlabel('$\omega$','FontName','SansSerif','FontSize',20);
    ylabel('$E(\omega)$','FontName','SansSerif','FontSize',20);
end

xline(omega_ibl(j-1),'LineWidth',1.5);
thisTitle = strcat('$\hat{x} =  $',num2str(round(x(j)) - 1850),' m');
title(thisTitle,'FontName','SansSerif','FontSize',20);
xlim([10^(-3) 10^1]);
ylim([10^(-6) 3*10^0]);
% legend('z/\delta = 0.01','z/\delta = 0.05','z/\delta = 0.1',...
%     'z/\delta = 0.3','Location','NorthEast');
end

% Alkali Flat Only

figure();
for i = 1:4
loglog(omega((end/2):end-myCutoff{1}(i)),...
    smoothdata(E_hat_w{1,i}((end/2):end-myCutoff{1}(i))),...
    'Color',newColors(i),'LineWidth',2); hold on;
end
xline(omega_asl,'LineWidth',1.5);
xlabel('$\omega$','FontName','SansSerif','FontSize',14);
ylabel('$E(\omega)$','FontName','SansSerif','FontSize',14);
title('Alkali Flat','FontName','SansSerif','FontSize',18);


%% Truncated Plot

% close all;

tiledlayout(2,2);
nexttile;
for i = 1:4
loglog(omega((end/2):end-myCutoff{1}(i)),...
    smoothdata(E_hat_w{1,i}((end/2):end-myCutoff{1}(i))),...
    'Color',newColors(i),'LineWidth',2); hold on;

end
set(gca,'FontName','SansSerif','FontSize',18);
xline(omega_asl,'LineWidth',1.5);
xlabel('$\omega$','FontName','SansSerif','FontSize',24);
ylabel('$E(\omega)$','FontName','SansSerif','FontSize',24);
title('Alkali Flat','FontName','SansSerif','FontSize',18);
legend('z/\delta = 0.01','z/\delta = 0.06','z/\delta = 0.10',...
            'z/\delta = 0.30','\omega_i = U_i/\delta_i',...
            'Location','SouthWest');
xlim([10^(-3) 10^1]);
ylim([10^(-6) 3*10^0]);

grabTheseFigs = [2,6,10];        
        
for j = 1:length(grabTheseFigs)
nexttile;
    k = grabTheseFigs(j);
    for i = 1:NzSelect-1
        loglog(omega((end/2):end-myCutoff{k}(i)),...
            smoothdata(E_hat_w{k,i}((end/2):end-myCutoff{k}(i))),...
            'Color',newColors(i),'LineWidth',2); hold on;
        set(gca,'FontName','SansSerif','FontSize',14);
        xlabel('$\omega$','FontName','SansSerif','FontSize',20);
        ylabel('$E(\omega)$','FontName','SansSerif','FontSize',20);
    end
    set(gca,'FontName','SansSerif','FontSize',18);
    xline(omega_ibl(k),'LineWidth',1.5);
    if j == 1
        makeATitle = strcat('$\hat{x}_1 = $',num2str(round(x(k)) - 1850),' m'); 
    elseif j == 2
        makeATitle = strcat('$\hat{x}_5 = $',num2str(round(x(k)) - 1850),' m');
    else
        makeATitle = strcat('$\hat{x}_9 = $',num2str(round(x(k)) - 1850),' m'); 
    end
    title(makeATitle,'Interpreter','Latex','FontSize',18);
    xlim([10^(-3) 10^1]);
    ylim([10^(-6) 3*10^0]);
end


%% Plotting at each x every z
% Goal is to see at what x/delta does the magnitude begin to overlap

% First, calculate the spectrum
for i = 1:N_u
    for j = 1:Nz
    
      up_temp = up{i}(:,j);

      allz_u_hat_w = (1/Nt) * fftshift(fft(up_temp));
      psi_w = abs(allz_u_hat_w).^2;
      allz_E_hat_w{i,j} = psi_w ./ (dw);
         
    % Pre-multiplied Energy Spectra
      allz_E_w{i,j} = allz_E_hat_w{i,j}.*omega';
    end 
      
end

clear urms_temp psi_w 

%% IBL Height from correlation

xhat = x(2:end) - 1850;
dibl_corr = 0.29.*(xhat).^(0.71);
dibl_corr(1:2) = 30;


%% Now plot at each x 

close all;


for i = 1:N_u-1
    
    red_count  = 0;
    blue_count = 0;

    figure(i)
    for j = 1:Nz
        if z_global(j) <= dibl_corr(i)
            red_count = red_count + 1;
            loglog(omega(end/2:end),smoothdata(allz_E_hat_w{i+1,j}(end/2:end)),...
                'r-',"LineWidth",1); hold on
        elseif z_global(j) >= dibl_corr(i)
            blue_count = blue_count + 1;
            loglog(omega(end/2:end),smoothdata(allz_E_hat_w{i+1,j}(end/2:end)),...
                'b-',"LineWidth",1); hold on
        end

        
    end

    store_red(i) = red_count;
    store_blue(i) = blue_count;

end



%% End