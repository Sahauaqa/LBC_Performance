clc
clear
close all
%% This is a code you can, it tells you what is the expected BC performance based on your inputs
% The assumptions are: ST_seawater = 1.9* ST_NaCl, using (Romer, 2013, JCP_B)
% Inputs for this code:
% Your requirements: dT [Kelvin]; total Area [m^2]; C0 and Cyield [ppm]; Recovery Rate [%]
% Additional inputs: dC in a Soret cell, Linear temperature percentage, Concentration separation ratio

% Outputs for this code:
% The number of cells required in M (row) and N (column) for the BC,
% Heat flux [kW/m^2] and total Pressure drop [Pa]
% Yield [L/day]

% !!! ======== For a specific separation ratio, it selects the optimal [M,N] combination for minimal amount of cells ======== !!!
% !!! ======== For the selected M and N, it then select the optimal sepration ratio that gives max yield ======== !!!

dT = 40; %[Kelvin]
totalA = 5; %[m^2]
C0 = 70000; %[ppm]
targetC = 35000; %[ppm]
targetRR = 50; %[%]
aspectRrange = 100:50:500; %aspect ratio (L/W) of the LBC you want, can be a numberic array
R_T = 65;% [%] Linear temperature percentage
R_C = [0.5:0.05:1]; % Concentration separation ratio, can be a numberic array
resultfname = "LBC_Resutls_Test.txt";

% Input the steady state concentration difference FdC between top and bottom boundaries in a Soret cell
% when 100% linear T difference = FdT and
% feedwater C = FC0
% (F is just because we were using Fortran to generate these values originally, no need to wonder why there is an 'F')
% Values here are based on seawater performance 
FdC = 4000; % [ppm] % This FdC is what's between the boundaries!
FdT = 42; % [K] % Can be different from C0 for the LBC, dC will be scaled
FC0 = 35000; % [ppm] % Can be different from C0 for the LBC, dC will be scaled

% From this point onwards, be careful when changing the code, it may casue
% issues
C0 = C0/1E6; % [Weight fraction]
targetC = targetC/1E6; % [Weight fraction]
targetRR = targetRR/100;  % [Weight fraction]
FdC = FdC/1E6; % [Weight fraction]
FC0 = FC0/1E6; % [Weight fraction]
FdCdev = FdC*0.95; % Fully developed defined as when dC>=0.95*dC(t inf)
% dC can be scaled based on:
% > dC is propotional to dT
% > dC is propotional to C0(1-C0), where C0 the initial C at the cell inlet
dCdev = FdCdev/(FC0*(1-FC0))*(C0*(1-C0)) / FdT*dT;

%% Form a [R_C = dC/dC_dev] -> [t_cha/tau] lookup table
% Reading in transient thermodiffusion data calculated with Fortran custom code, working fluid is NaCl
% ST_NaCl, using (Romer, 2013, JCP_B); D_NaCl, using (Caldwell, 1973)
fname=sprintf('transientC.txt');
M1 = readmatrix(fname);
time = M1(:,1);
y = M1(:,2);
Cr = M1(:,6);
height = 1; % Channel height [mm]

nth = find(time>0,1)-1;
T = length(time)/nth;
for i = 1 : T
    time_new(:,i) = time((i-1)*nth+1:i*nth);
    y_new(:,i) = y((i-1)*nth+1:i*nth);
    Cr_new(:,i) = Cr((i-1)*nth+1:i*nth);
end
% Concentration
FdCtrans = Cr_new(1,:)-Cr_new(end,:); % C difference between boundaries
FdCmax = Cr_new(1,end)-Cr_new(end,end); % max possible dC
dCRatio = FdCtrans./FdCmax; % fully-developed dC defined as when dC_dev/dC(t inf) = 0.95
% Time
dev_i = find(dCRatio>0.95, 1); % Index of when dC is fully developed
Ftau = time_new(1,dev_i); % tau, C profile development time
FdCdev = FdCtrans(dev_i); % dC_dev
%Look up table for R_C -> t_cha/tau
tRatio = time_new(1,:)/Ftau; %t_cha/tau
% FdCdev = Cr_new(1,dev_i)-Cr_new(end,dev_i);
R_Crange = FdCtrans./FdCdev; % R_C



%% Form a [R_T = h_linearT / total h] -> [dC(R_T) / dC(R_T = 100%)] lookup table
% If using symmetric T profile, use txt file 'SymmetricT.txt'
fname2=sprintf('SymmetricT.txt');
% If using asymmetric T profile, use txt file 'AsymmetricT.txt'
%fname2=sprintf('AsymmetricT.txt');
M2 = readmatrix(fname2);
R_Trange = M2(:,1); % [R_T = h_linearT / total h]
dCRatio_RTrange = M2(:,2); % [dC(R_T) / dC(R_T = 100%)] for different R_T
RatioLow = interp1(R_Trange, dCRatio_RTrange, R_T, 'linear', 'extrap'); % Linear interpolation, for top low C part, [dC(R_T) / dC(R_T = 100%)]
RatioHigh = RatioLow; % For btm high C part, [dC(R_T) / dC(R_T = 100%)]

%% Burgers calculation based on the read data
Mlong_rec = zeros(length(R_C), 1);
Nfull_rec = zeros(length(R_C), 1);
celltime = zeros(length(R_C), 1);
totaltime = zeros(length(R_C), 1);

for k = 1:length(R_C)
    rc = R_C(k);
    
    % Now we want to calculate M and N,
    % For a particular N, there is a equilibrium M that if you increase M further, there will be little further drop in C
    Nfull = 10; % Nfull is no of full cells in parallel to the flow direction, Nfull us equivalnet to j in paper
    Mlong = Nfull*100; % Mlong equivalnet to i in paper. From experiences, Meq < Nfull*100
    % You can have other strategies to generate a random Mlong, but note Mlong > Meq(Nfull)
    [Cyield, ~, ~] = BC_yieldC(C0, targetRR, rc, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
    step = 5; % Random number that is not too slow (step=1) or too rough

    while Cyield > targetC %Increase N until Cyield is low enough
        Nfull = Nfull+step;
        Mlong = Nfull*100;
        [Cyield, ~, ~] = BC_yieldC(C0, targetRR, rc, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
    end

    % Now we first reduce M to Meq
    [Cyield, CexcessM, Meq] = BC_yieldC(C0, targetRR, rc, Mlong, Nfull, dCdev, RatioLow, RatioHigh);

    Mlong = Meq;

    % Now we try to adjust M, N so that we get the same C_yield with fewer rows (thus higher LBC flow rate)
    totalRow_old = Mlong;
    totalRow_new = totalRow_old-1;
    while totalRow_new < totalRow_old
        totalRow_old = Mlong;
        Nfull = Nfull+1; % Increase N until we get max possible N
        Mlong = Mlong-1; % Decrease M gradually
        [Cyield, ~, ~] = BC_yieldC(C0, targetRR, rc, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
        % If for N = fixed number, we try to find a mininal M that satisfy Cyield <= target C
        while Cyield < targetC
            Mlong = Mlong-1;
            [Cyield, ~, ~] = BC_yieldC(C0, targetRR, rc, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
        end
        Mlong = Mlong+1;
        totalRow_new = Mlong; % Update total cell number
    end
    Nfull_rec(k) = Nfull-1;
    Mlong_rec(k) = totalRow_new;
    R_Ci = find(R_Crange>rc,1);
    celltime(k) = Ftau*tRatio(R_Ci); % [s] Time in individual cell
    totaltime(k) = celltime(k)*(Mlong*2); % [s], Total time it takes
end

figNamestr = sprintf("BC design when desalinating from %d to %d ppm", C0*1E6, targetC*1E6);
figure('Name', figNamestr)

% Show the initial condition in the 1st panel
subplot(2,3,1)
plot([0,100,100],[100,100,0])
txty = 90;
text(5,txty,"The operation conditions")
text(5,txty-10,sprintf(">The temperautre difference is %.0f Kelvin", dT))
text(5,txty-20,sprintf(">Total area of this Buregers Cascade is %.0f m^2", totalA))
text(5,txty-30,sprintf(">The input concnetration is %.0f ppm and ", C0*1E6))
text(5,txty-35,sprintf("the output is %.0f ppm, with a recovery rate of %.0f %%", targetC*1E6, targetRR*1E2))
text(5,txty-45,sprintf(">The aspect ratio is between %.0f and %.0f,", aspectRrange(1), aspectRrange(end)))
text(5,txty-50,sprintf("with a step of %.0f", aspectRrange(2)-aspectRrange(1)))
set(gca,'XTick',[], 'YTick', [])
set(gca,"FontSize",12)
pbaspect([1 1 1])


% Plot normalised separation time
subplot(2,3,2)
i1 = find(R_C==1);
nomralisedCelltime = celltime./celltime(i1);
nomralisedLBCtime = totaltime./totaltime(i1);
yyaxis left
plot(R_C,nomralisedCelltime)
hold on
xlabel("Separation ratio: dC_{cell} / dC_{fully developed}")
ylabel("Normalised separation time in each channel")
yyaxis right
plot(R_C,nomralisedLBCtime)
ylabel("Normalised separation time in LBC")
hold off
grid on
pbaspect([1 1 1])
set(gca,"FontSize",12)


% Volumetric flow rate (L/day])
% Flow speed (v) = l ./ cell_time
% Cross sectional area (cs) = Nfull * w * h
% Total flow rate: v*cs = h * l * w * Nfull/ cell_time
% =  h * (l * Mlong*2) * (w * Nfull) / (cell_time * Mlong*2) = h* totalA / totalT
Mlong1 = Mlong_rec(i1);
Nfull1 = Nfull_rec(i1);
Q_yield = 1/1000 * totalA ./ totaltime * targetRR; % Unit: [m^3/s]
Q_yield = Q_yield*1000*3600*24; % Unit : [L/m^3]
subplot(2,3,3)
yyaxis left
plot(R_C,Q_yield)
xlabel("Separation ratio: dC_{cell} / dC_{fully developed}")
ylabel("Freshwater volumetric flow rate if height is 1 mm (L/day)")
grid on
pbaspect([1 1 1])
set(gca,"FontSize",12)
optimal_index = find(Q_yield == max(Q_yield));
optimalSR = R_C(optimal_index);
optmialcelltime = celltime(optimal_index);
text(optimalSR, 0.8*max(Q_yield), ["Optimal separation ratio is" optimalSR]);
% Plot the number of cells required for different R_c
yyaxis right
plot(R_C, (Mlong_rec.*Nfull_rec)/(Mlong1*Nfull1))
ylabel("Nomralised number of cells requried")
grid on
pbaspect([1 1 1])
set(gca,"FontSize",12)

% Plot the contour plots if fully separated in each cell
subplot(2,3,4)
title("C contour if using optimal sepration ratio")
Mlong = Mlong_rec(optimal_index);
Nfull = Nfull_rec(optimal_index);
plotCcontour(Mlong, Nfull, optimalSR, C0, targetRR, dCdev, RatioLow, RatioHigh)


Hrange = 1:5; % [mm]
Nfull = Nfull_rec(optimal_index);
Mlong = Mlong_rec(optimal_index);
[lgdtxt, Pdrop_rec, flux, optimalYield] = parameters(totalA, targetRR, height, optmialcelltime, Hrange, aspectRrange, dT, Nfull, Mlong);
flux = flux/R_T;
%Pressure drop as a function of cell parameters
subplot(2,3,5)
plot(Hrange,Pdrop_rec/1E3,"o","linewidth", 2)
ylabel("Pressure drop, \Delta{\itP} (kPa)");
ytickformat("%.1f")
xlabel("Cell height, {\ith} (mm)");
leg = legend(lgdtxt);
title(leg,"BC Dimension");
grid on
pbaspect([1 1 1])


% Yield and heat flux when changing height
subplot(2,3,6)
yyaxis left
plot(Hrange, flux,"-","linewidth", 2)
ylabel("Heat flux, {\itq} (kW m^{-2})");
hold on
yyaxis right
plot(Hrange, optimalYield,"--", "linewidth", 2)
xlabel("Cell height, {\ith} (mm)");
ylabel("Freshwater yield (L day^{-1})");
xlim([1,5])
xtickformat("%.1f")
ytickformat("%.2f")
set(gca,"FontSize",12)
grid on
set(gca,"fontsize", 12)
pbaspect([1 1 1])
hold off

resultID = fopen(resultfname,"a");
fprintf(resultID,'Symmetric with %.d pc linear T, ratio High %.4f, ratio Low, %.4f\n',R_T, RatioHigh, RatioLow);
fprintf(resultID,'dC/dC_{fullydev}: \n');
fprintf(resultID,'%.2f, ',R_C);
fprintf(resultID,'\n');
fprintf(resultID,'Normalised separation time, \n');
fprintf(resultID,'%.3f, ',nomralisedLBCtime);
fprintf(resultID,'\n');
fprintf(resultID,'Total number of cells, \n');
fprintf(resultID,'%.1f, ',Mlong_rec.*Nfull_rec);
fprintf(resultID,'\n');
fprintf(resultID,'Freshwater yield (L/day) when h = 1mm, \n');
fprintf(resultID,'%.1f, ',Q_yield);
fprintf(resultID,'\n');
fprintf(resultID,'Cell height (mm), \n');
fprintf(resultID,'%.d, ',Hrange);
fprintf(resultID,'\n');
fprintf(resultID,'Heat Flux (kW/m^2), \n');
fprintf(resultID,'%.1f, ',flux);
fprintf(resultID,'\n');
fprintf(resultID,'Optimal freshwater yield (L/day) at optimal SR%.2f, \n',optimalSR);
fprintf(resultID,'%.1f, ',optimalYield);
fprintf(resultID,'\n');
fprintf(resultID,'\n');
fclose(resultID);


%% ========================================================================
function [dimensiontxt, Pdrop, Flux, optimalYield] = parameters(totalA, RR, height, SR_t, Hrange, aspectRrange, delT, Nfull, Mlong)
%delta T between the cold side of the Burgers cascade and the bluk water for cooling [K]
delT_conv = 20;

% Cell parameters
M = Mlong*2;
N = Nfull;
Wrange = sqrt(totalA./aspectRrange)./N;
Lrange = sqrt(totalA.*aspectRrange)./M;

%Convert height unit to m
Hrange = Hrange/1000;
height = height/1000;

% Thermophysical properties of the fluid
% Water @ 40deg C
miu = 0.6527*10^(-3); %dynamic viscosity [Pa.s]
rho = 992.2; %density [kg/m]
k = 0.641; %thermal conudctivity [W/m.K]

%individual cell parameters
Pdrop = zeros(length(Hrange),length(Wrange));
Flux = zeros(length(Hrange),1);
optimalYield = zeros(length(Hrange),1);
for hh = 1:length(Hrange)
    for ww = 1:length(aspectRrange)
        H = Hrange(hh);
        W = Wrange(ww);
        L = Lrange(ww); %Cell length = totalArea/totalWidth/M  [m]
        t_cell = SR_t*(H^2)/(height^2); %Time taken for C profile to fully develop in a cell [s], for 1mm height, it takes SR_h s, and t is proportional to h^2
        Velocity = L/t_cell; %Max velocity for C profile to fully develop [m/s]
        FR = (Velocity)*(H)*(W)*N; %Flow rate in the entire Burgers cascade [m^3/s]

        %Pressure
        raspect = H/W; %aspect ratio [none] = H[m]/W[m]
        Dh = 2*(H)*(W/(H+W)); %Hydrodynamic diameter [m] = 2*H[m]*W[m]/(H[m]*W[m])
        Re = (rho*(Velocity)*Dh)/miu; %Reylaigh number = rho[kg/m^3]*V[m/s]*Dh[m] / miu[Pa.s]
        f = 64/Re; %Friction factor, C = 64, for values of C, refer to Table 8.3 on p450
        delP_cell = f*rho*(Velocity)^2*(L)/(2*Dh); %Pressure drop in a cell [Pa] = f*rho[kg/m^3]*V[m/s]^2*L[m]/(2*Dh[m])
        delP_burgers = delP_cell*N*M; %Pressure drop in the entire Burgers cascade [Pa]
        Wpump = FR*delP_burgers; %Pump work [W] =  FR[m^3/s] * delP [Pa]

        %Matrix
        t_burgers = 1/(RR*FR*1000*3600); % [hr/L] Time [hour] taken to produce 1L of fresh(er) water for a certain recoery rate [hr/L]
        yield = RR*FR*1000*3600*24; %'Fresh water yield [L/day]
        Ee = Wpump*t_burgers/1000; %Electrical energy consumed to produce 1L of fresh water [kWh/L]
        Lburgers = (L)*M; %Total length of the burgers [m]
        Wburgers = (W)*N; %Total width of the burgers [m]

        Pdrop(hh,ww) = delP_burgers;
        optimalYield(hh,1) = yield;
    end
    %Heat flux
    q = (k*delT/(H))/1000; %kW/m^2 = (k[W/m.K]*delT[K]/H[m])/1000
    h_conv = (q*1000)/delT_conv; %heat transfer coefficient [W/K] = (q[kW/m^2]*1000)/delT[K];
    Flux(hh,1) = q/1;
end

Wrange = Wrange*1000;
Lrange = Lrange.*1000;
dimensiontxt = strings(length(aspectRrange),1); %String for pressure drop plot legend
for ww = length(aspectRrange):-1:1
    dimensiontxt(ww) = sprintf('r = %.0f with %.1f by %.1f cell',aspectRrange(ww), Wrange(ww), Lrange(ww));
end
end



function plotCcontour(Mlong, Nfull, sr, C0, targetRR, dCdev, RatioLow, RatioHigh)
[~, C] = BC_yieldC(C0, targetRR, sr, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
N = Nfull*2+1; % In the Saiki paper, 2 cells (shifted) is 1 long cell
M =  Mlong*2; % 1 long cell = 2 rectangular cells

% Replace NaN with the average of neighbouring items so that contour plot can be plotted
for m =1:M
    for n=1:N
        if C(m,n) == 0
            if m == 1
                C(m,n) = C0;
            elseif m == M
                C(m,n) = (C(m,n-1)+C(m,n+1))/2;
            elseif n == 1
                C(m,n) = (C(m-1,n)+ C(m+1,n))/2;
            elseif n == N
                C(m,n) = (C(m-1,n)+ C(m+1,n))/2;
            else
                Csurr = [C(m,n-1),C(m-1,n),C(m+1,n),C(m,n+1)];
                C(m,n) = mean(Csurr);
            end
            Velocity(m,n) = NaN;
        end
    end
end
xlim([1,N/2])
ylim([1,M/2])
[X,Y]=meshgrid(1:M,1:N);
% contourf(X,Y,C')
S = image((1:N)/2,(M:-2:0)/2,C*1000,'CDataMapping','scaled');
%set(gca,'CLim',[29,31]);
%image(C,'CDataMapping','scaled')
hold on
c = colorbar('Ticks',0:10:150);
c.Ruler.TickLabelFormat = '%.0f';
c.Label.String = '{\itC} (10^{3} ppm)';
c.Label.FontSize = 12;
%S.EdgeColor='none'; %shading interp %Smooth the color
contour(Y/2,M/2-X/2,C'*1000,[5,10,20,30,34,36,40,50,60,80,100],'LineColor','k','ShowText','on') %Add the black lines
ylabel("{\itm} for M_{long}")
xlabel("{\itn} for N_{full}")
set(gca, "Fontsize",12)
hold off
end





function [Cyield, C, Meq] = BC_yieldC(C0,targetRR, SR, Mlong, Nfull, dCdev, RatioLow, RatioHigh)
% Assumptions:
% Perfect bifurcation at the horizontal midplane
% dC in each cell can be scaled based on:
% > dC is propotional to dT
% > dC is propotional to C0(1-C0), where C0 the initial C at the cell inlet
dC = dCdev*SR; % achieved dC is smaller due to not long enough time
N = Nfull*2+1; % Matlab cannot do 1/2 indices
M = Mlong*2; % Matlab cannot do 1/2 indices
C = zeros(M,N);
Velocity = zeros(M,N);

% Calculating the C in each cell
% First row, C0 everywhere:
for n=2:2:N
    C(1,n) = C0;
    Velocity(1,n) = 1;
end
% Following row, always be based on the previous row
for m = 2:M
    % mode(m,2)==1 is equivalent to i=n+1/2 in manuscript, it has N full cells
    % and no half cells
    if mod(m,2) == 1
        for n=2:2:N
            if n==2
                dC2 = RatioLow*SR*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1) + Velocity(m-1,n+1)/2;
                C(m,n) = (C(m-1,n-1)*Velocity(m-1,n-1) + (C(m-1,n+1)-dC2/4.32)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            elseif n == 2*Nfull
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1);
                C(m,n) = ((C(m-1,n-1)+dC1/4.32)*Velocity(m-1,n-1)/2 + C(m-1,n+1)*Velocity(m-1,n+1))/Velocity(m,n);
            else
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1)/2;
                C(m,n) = ((C(m-1,n-1)+dC1/4.32)*Velocity(m-1,n-1)/2+(C(m-1,n+1)-dC2/4.32)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            end
        end
        % mode(m,2)==0 is equivalent to i=n in manuscript, it has N-1 full cells
        % and 2 half cells
    elseif mod(m,2) == 0
        for n=1:2:N
            if n==1
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n+1)/2;
                C(m,n) = C(m-1,n+1)-dC2/4.32;
            elseif n == N
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2;
                C(m,n) = C(m-1,n-1)+dC1/4.32;
            else
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1)/2;
                C(m,n) = ((C(m-1,n-1)+dC1/4.32)*Velocity(m-1,n-1)/2+(C(m-1,n+1)-dC2/4.32)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            end
        end
    end
end

% Calculate C_yield based on recovery rate
SamplingCell = fix(N*targetRR);
LowerM = 0;
LowerV = 0;
for n = 1:SamplingCell
    LowerM = LowerM + C(M,n)*Velocity(M,n);
    LowerV = LowerV + Velocity(M,n);
end
LowerC = LowerM/LowerV;
Cyield= LowerC;

% Fing Meq that satisfies C_drop(j)>= 95% C_drop(equlibrium)
LowerC = zeros(1,M);
LowerRR = zeros(1,M);

for m = 2:2:M %each row
    LowerM = 0;
    LowerV = 0;
    %for n = N:-1:(N-SamplingCell)
    for n = 1:SamplingCell
        LowerM = LowerM + C(m,n)*Velocity(m,n);
        LowerV = LowerV + Velocity(m,n);
    end
    LowerC(m) = LowerM/LowerV;
    LowerC(m-1) = LowerC(m);
    LowerRR(m) = LowerV/Nfull;
    LowerRR(m-1) = LowerRR(m);
end

Cdrop = C0-LowerC;
Meq = find(Cdrop>Cdrop(end)*0.95,1);
Meq = Meq/2+0.5; %Convert to full-length cell

% % If purpose is to increase concentration, e.g. salination
% SamplingCell = fix(N*targetRR);
% higherM = 0;
% higherV = 0;
% for n = N:-1:(N-SamplingCell+1)
%     higherM = higherM + C(M,n)*Velocity(M,n);
%     higherV = higherV + Velocity(M,n);
% end
% higherC = higherM/higherV;
% Cyield= higherC;
end


