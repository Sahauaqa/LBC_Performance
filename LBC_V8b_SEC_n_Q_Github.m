clc
clear
close all

%% This is a code you can play with, it tells you what is the expected BC performance based on your inputs
% Inputs for this code:
% Your requirements: dT [Kelvin]; T_mean [degC]; 
% Fluid properties: S_T of the fluid [10^{-3} K^{-1}]

% Outputs for this code:
% SEC (KWh/m3 of yield) and unit area flow rate (m3 of yield/m2)

% The parameters this program vary to minimise SEC:
% R_C: linear concentration range
% RR: recovery rate
% M, N: Row and column of channels within the LBC


dT = 80; %[Kelvin]
Tmean = 50; %[degC]
ST = 10; %[10^{-3} K^{-1}], input the ST at ||25 degC|| please!
D = 1.35; %[10^{-9} m^2s^{-1}], input the D at ||25 degC|| please!
resultfname = "SECnQ_NaOH.txt";


% function SECnQin1(dT, Tmean, ST, D, resultfname)
C0_rg = [10:10:200]*1e3; %Feed water concentration, [ppt]
Cy_rg = [11:10:251]*1e3; %Yield concentration [ppt]

C0_rg = C0_rg/1E6; % [Weight fraction]
Cy_rg = Cy_rg /1E6; % [Weight fraction]
sr = 1; %Splite ratio, 1 is perfect split, 0 if complete mixing

for C0 = C0_rg
    flag = 0;

    for Cy = Cy_rg
        RRrg = 20:10:80; % Recovery rate range [%]
        RRrg = RRrg/100;  % [Weight fraction]

        aspectR = 50000; %aspect ratio (L/W) of the LBC you want, it won't affect SEC or Q/A
        totalA = 100; % [m2] total area occupied by the LBC saline channels
        RT = 60;% [%] Linear temperature percentage
        RCrg = 0.5 : 0.1 : 1.0; % Concentration separation ratio, can be a numberic array

        % Base condtion
        % The steady state concentration difference FdC between top and bottom boundaries in a Soret cell
        % when 100% linear T difference = FdT
        % feedwater C = FC0
        % feed solution is NaCl/H2O
        FdC = 3000; % [ppm] % This FdC is what's between the boundaries for NaCl @ 25 degC
        FdC = FdC/1e6; % [wt frac]
        FdT = 42; % [K] % Can be different from C0 for the LBC, dC will be scaled
        FC0 = 35000; % [ppm] % Can be different from C0 for the LBC, dC will be scaled
        FC0 = FC0/1e6; % [wt frac]
        Ftau = 47.5; % tau, C profile development time, Based on tau = h^2/(pi^2 D), D = 2.132*10^-9 m^2/s
        % Add species correction
        FdC = FdC*ST/1.203; %FdC(25degC, X) = FdC(25degC, NaCl)*ST(25degC, X)/ST(25degC, NaCl)
        Ftau = Ftau/(D/1.498); %Ftau(25degC, X) = Ftau(25degC, NaCl)* (D(25degC, X)/D(25degC, NaCl))^{-1}
        %Add temperautre correction, assuming ST(T) and D(T) curve is the same for X and NaCl
        FdC = FdC*(5*(1-exp((9.87-Tmean)/55))) / 1.203; %FdC(Tmean, X) = FdC(25degC, X)*ST(Tmean, NaCl)/ST(25degC, NaCl)
        Ftau = Ftau/((0.44+0.0423*Tmean) / 1.498); %Ftau(Tmean, X) = Ftau(25degC, X)* (D(Tmean, NaCl)/D(25degC, NaCl))^{-1}

        FdCdev = FdC*0.95; % Fully developed defined as when dC>=0.95*dC(t inf)

        % Form a [R_C = dC/dC_dev] -> [t_cha/tau] lookup table
        % Reading in transient thermodiffusion data calculated with Fortran custom code, working fluid is NaCl
        % ST_NaCl, using (Romer, 2013, JCP_B); D_NaCl, using (Caldwell, 1973)
        fname=sprintf('transientC.txt');
        M1 = readmatrix(fname);
        time = M1(:,1);
        y = M1(:,2);
        Cr = M1(:,3);
        height = 1; % Channel height [mm]

        nth = find(time>0,1)-1;
        T = length(time)/nth;
        for i = 1 : T
            time_new(:,i) = time((i-1)*nth+1:i*nth);
            y_new(:,i) = y((i-1)*nth+1:i*nth);
            Cr_new(:,i) = Cr((i-1)*nth+1:i*nth);
        end
        % Concentration
        dCtrans = Cr_new(1,:)-Cr_new(end,:); % C difference between boundaries
        dCmax = Cr_new(1,end)-Cr_new(end,end); % max possible dC
        dCRatio = dCtrans./dCmax; % fully-developed dC defined as when dC_dev/dC(t inf) = 0.95
        % Time
        dev_i = find(dCRatio>0.95, 1); % Index of when dC is fully developed
        dCdev_RC = dCtrans(dev_i); % dC_dev
        %Look up table for R_C -> t_cha/tau
        txttau = time_new(1,dev_i); % tau as in the txt file
        tRatio = time_new(1,:)/txttau; %t_cha/tau
        R_Crange = dCtrans./dCdev_RC; % R_C

        %% Form a [R_T = h_linearT / total h] -> [dC(R_T) / dC(R_T = 100%)] lookup table
        % If using symmetric T profile, use txt file 'SymmetricT.txt'
        fname2=sprintf('SymmetricT.txt');
        % If using asymmetric T profile, use txt file 'AsymmetricT.txt'
        %fname2=sprintf('AsymmetricT.txt');
        M2 = readmatrix(fname2);
        R_Trange = M2(:,1); % [R_T = h_linearT / total h]
        dCRatio_RTrange = M2(:,2); % [dC(R_T) / dC(R_T = 100%)] for different R_T
        RatioLow = interp1(R_Trange, dCRatio_RTrange, RT, 'linear', 'extrap'); % Linear interpolation, for top low C part, [dC(R_T) / dC(R_T = 100%)]
        RatioHigh = RatioLow; % For btm high C part, [dC(R_T) / dC(R_T = 100%)]

        % dC can be scaled based on:
        % > dC is propotional to dT
        % > dC is propotional to C0(1-C0), where C0 the initial C at the cell inlet
        dCdev = sr*FdCdev / (FC0*(1-FC0)) * (C0*(1-C0)) / FdT*dT;

        %% Burgers calculation based on the read data
        Mlong_rec = zeros(length(RRrg), length(RCrg));
        Nfull_rec = zeros(length(RRrg), length(RCrg));
        celltime = zeros(length(RRrg), length(RCrg));
        totaltime = zeros(length(RRrg), length(RCrg));
        mockQ = zeros(length(RRrg), length(RCrg));


        % !!! ======== For a specific separation ratio, it selects the optimal [M,N] combination for minimal amount of cells ======== !!!
        for RRi = 1:length(RRrg)
            for RCj = 1:length(RCrg)

                targetRR = RRrg(RRi); %Recovery rate for this loop
                RC = RCrg(RCj); %Concentration spearation ratio for this loop

                % Now we want to calculate M and N,
                % For a particular N, there is a equilibrium M that if you increase M further, there will be little further drop in C
                Nfull = 5; % Nfull is no of full cells in parallel to the flow direction, Nfull us equivalnet to j in paper
                Mlong = Nfull*100; % Mlong equivalnet to i in paper. From experiences, Meq < Nfull*100
                % You can have other strategies to generate a random Mlong, but note Mlong > Meq(Nfull)

                step = 5; % Random number that is not too slow (step=1) or too rough
                if Cy>C0 %TDS

                    Cymax = C0/targetRR-1000/1e6; % Max possible C of the yield if all salts concentrated into yield
                    Cymax_old = 1000; % Max possible C of the yield if all salts concentrated into yield, for the previous RR
                    if RRi>1
                        Cymax_old = C0/RRrg(RRi-1)-1000/1e6;
                    end
                    % Avoids situation where salination at certian recovery
                    % rate is not possible at there is not enough salt in
                    % the solution
                    if RRi == 1 && Cy>Cymax
                        Cy = Cymax;
                    elseif RRi>1 && Cy>Cymax
                        flag = 1;
                        break
                    elseif RRi>1 && Cy>Cymax_old
                        flag = 1;
                        break
                    end

                    [Cyield, ~, ~] = BC_yieldC(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                    while Cyield < Cy %Increase N until Cyield is high enough
                        Nfull = Nfull+step;
                        Mlong = Nfull*100;
                        [Cyield, ~, ~] = BC_yieldC(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                    end

                    % Now we first reduce M to Meq
                    [~, ~, Meq] = BC_yieldC(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                    Mlong = Meq;
                    % Now we try to adjust M, N so that we get the same C_yield with fewer rows (thus higher LBC flow rate)
                    Mold = Mlong;
                    Mnew = Mold-1;
                    while Mnew < Mold && Mlong>2
                        Mold = Mlong;
                        Nfull = Nfull+1; % Increase N until we get max possible N
                        Mlong = Mlong-1; % Decrease M gradually
                        [Cyield, ~, ~] = BC_yieldC(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                        % If for N = fixed number, we try to find a mininal M that satisfy Cyield >= target C
                        while Cyield > Cy && Mlong>1
                            Mlong = Mlong-1;
                            [Cyield, ~, ~] = BC_yieldC(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                        end
                        Mlong = Mlong+1;
                        Mnew = Mlong; % Update total cell number
                    end

                elseif Cy<C0 %TDD
                    [Cyield, ~, ~] = BC_yieldC_TDD(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                    while Cyield > Cy %Increase N until Cyield is low enough
                        Nfull = Nfull+step;
                        Mlong = Nfull*100;
                        [Cyield, ~, ~] = BC_yieldC_TDD(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                    end

                    % Now we first reduce M to Meq
                    [~, ~, Meq] = BC_yieldC_TDD(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                    Mlong = Meq;
                    % Now we try to adjust M, N so that we get the same C_yield with fewer rows (thus higher LBC flow rate)
                    Mold = Mlong;
                    Mnew = Mold-1;
                    while Mnew < Mold && Mlong>2
                        Mold = Mlong;
                        Nfull = Nfull+1; % Increase N until we get max possible N
                        Mlong = Mlong-1; % Decrease M gradually
                        [Cyield, ~, ~] = BC_yieldC_TDD(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                        % If for N = fixed number, we try to find a mininal M that satisfy Cyield <= target C
                        while Cyield < Cy && Mlong>1
                            Mlong = Mlong-1;
                            [Cyield, ~, ~] = BC_yieldC_TDD(C0, targetRR, RC, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
                        end
                        Mlong = Mlong+1;
                        Mnew = Mlong; % Update total cell number
                    end

                else% No concentration change
                    Nfull = 1;
                    Mnew = 1;
                end

                Nfull_rec(RRi,RCj) = Nfull-1;
                Mlong_rec(RRi,RCj) = Mnew;
                R_Ci = find(R_Crange>RC,1);
                celltime(RRi,RCj) = Ftau*tRatio(R_Ci); % [s] Time in individual cell
                totaltime(RRi,RCj) = celltime(RRi)*(Mnew*2); % [s], Total time it takes
                mockQ(RRi,RCj) = 1 ./ totaltime(RRi) * targetRR; % Unit: [m^3 s^{-1} per m^2 footprint area]
                %Q_yield(RRi,RCj) = Q_yield(RRi)*1000*3600*24; % Unit : [L day^{-1} per m^2 footprint area]

            end
            if flag ==1
                break
            end
        end

        if Cy ~= C0
            % [optimal_i,optimal_j]  = find(totaltime == min(min(totaltime))); % Index of the optimal recovery rate, obtained when Q_yield is maximised
            [optimal_i,optimal_j]  = find(mockQ == max(max(mockQ))); % Index of the optimal recovery rate, obtained when Q_yield is maximised

            optimalRR = RRrg(optimal_i(1));
            optimalRC = RCrg(optimal_j(1));
            optmialcelltime = celltime(optimal_i(1), optimal_j(1));
            optimalN = Nfull_rec(optimal_i(1), optimal_j(1));
            optimalM = Mlong_rec(optimal_i(1), optimal_j(1));

            Hrange = 1:2; % [mm]
            [~, Pdrop_rec, flux, optimalYield] = parameters(totalA, optimalRR, height, optmialcelltime, Hrange, aspectR, dT, RT, optimalN, optimalM);

            %Shading effects
            if Cy >= C0  %TDS
                flux = flux*0.4;
            else %TDD
                flux = flux*0.5;
            end

            calYield = (optimalYield(1)/1000); % Yield for calculation, [L/day] to [m^3/day] per unit footprint area
            optimalqQ = flux(1)*24/calYield; % SEC = flux [kW/m^2] *24 [hr/day] / [m^3/day/m^2] = SEC [kWh/m^3]
            optimalExergy = optimalqQ * (dT/(Tmean+273.15)); % X = Q* (Carnot Cycle Efficiency)

        else
            optimalRR = 1;
            optimalRC = 0;
            optimalqQ = 0;
            optimalExergy = 0;
            optimalYield = 10;
        end

        resultID = fopen(resultfname,"a");
        fprintf(resultID,'C0 [ppt], %d, Cyield [ppt], %d, Optimal RT [pc], %.1f, Optimal R_C, %d, SECwShd [kWh/m^3], %d, X p Yield [kWh/m^3], %d, Q [m^3/day m^2], %.2f', C0*1000, Cy*1000, optimalRR*100, optimalRC, optimalqQ, optimalExergy(1), optimalYield(1)/1000);
        fprintf(resultID,'\n');

    end
end

%Plotting and saving the data in plots
fclose(resultID);
M1 = readmatrix(resultfname);

%Contour plot of SEC, x is C0, y is Cout
x = M1(:,2); % C0 in ppt
y = M1(:,4); % Cy in ppt

% I = M1(:,6); %Column 6 is optimal recovery rate
% I = M1(:,8); %Column 8 is optimal concentration separaiton ratio (R_C)
SEC = M1(:,10); %Column 10 is SEC
% I = M1(:,12); %Column 12 is exergy per yield
FlowRate = M1(:,14); %Column 14 is volumetric flow rate per unit area

X = reshape(abs(x), [], 1);
Y = reshape(abs(y), [], 1);
L = reshape(abs(SEC), [], 1);
% Create the interpolation
n_ptsX = 190;
n_ptsY = 240;
F = TriScatteredInterp(X,Y,L);
[yi,xi] = ndgrid(max(Y(:)):(-abs(max(Y(:))))/n_ptsY:min(Y(:)), ...
    min(X(:)):abs(max(X(:)))/n_ptsX:max(X(:)));
vi = F(xi,yi);

figure(1)
set(gcf,'units','points','position',[100,300,300,200]) %Plot size
% S = image(vi,'CDataMapping','scaled');
contourf(xi(1,:),yi(:,1),vi,'ShowText','on') %Add the black lines
c = colorbar();
c.Ruler.TickLabelFormat = '%.0f';
c.Label.String = 'SEC (kWh m^{-3})';
c.Label.FontSize = 12;
%S.EdgeColor='none'; %shading interp %Smooth the color
ytickformat("%.1f")
xtickformat('%.1f')

ylabel("{\itC}_{yield} (10^{3} ppm)")
xlabel("{\itC}_{0} (10^{3} ppm)")
set(gca, "Fontsize",10)

exportgraphics(gcf,'SEC.pdf','ContentType','vector')
SECdatafile='SEC.xls';
writematrix(vi, SECdatafile);

%Contour plot of flow rate per area, x is C0, y is Cout
L2 = reshape(abs(FlowRate), [], 1);
F2 = TriScatteredInterp(X,Y,L2);
vi2 = F2(xi,yi);

figure(2)
set(gcf,'units','points','position',[600,300,300,200]) %Plot size
% S = image(vi,'CDataMapping','scaled');
contourf(xi(1,:),yi(:,1),vi2,'ShowText','on') %Add the black lines
c = colorbar();
c.Ruler.TickLabelFormat = '%.0f';
c.Label.String = 'Q (m^{3} m^{-2})';
c.Label.FontSize = 12;
%S.EdgeColor='none'; %shading interp %Smooth the color
ytickformat("%.1f")
xtickformat('%.1f')

ylabel("{\itC}_{yield} (10^{3} ppm)")
xlabel("{\itC}_{0} (10^{3} ppm)")
set(gca, "Fontsize",10)

exportgraphics(gcf,'FlowRate p Area.pdf','ContentType','vector')
FRdatafile='FR.xls';
writematrix(vi2,FRdatafile);
% end
fclose("all")

%% ========================================================================
function [dimensiontxt, Pdrop, HeatRate, YieldFlowRate] = parameters(totalA, RR, height, celltime, Hrange, aspectRrg, delT, R_T, Nfull, Mlong)
%delta T between the cold side of the Burgers cascade and the bluk water for cooling [K]
delT_conv = 20;
Recri = 1000; % Critical Reyleigh number

% Cell parameters
M = Mlong*2;
N = Nfull;
Wrange = sqrt(totalA./aspectRrg)./N; %[m]
Lrange = sqrt(totalA.*aspectRrg)./M; %[m]
WBurgers_rg = sqrt(totalA./aspectRrg); %[m]
LBurgers_rg = sqrt(totalA.*aspectRrg); %[m]

%Convert height unit to m
Hrange = Hrange/1000;
height = height/1000;

% Thermophysical properties of the fluid
% Water @ 100deg C
miu = 0.283*10^(-3); %dynamic viscosity [Pa.s]
rho = 958.4; %density [kg/m]
k = 0.680; %thermal conudctivity [W/m.K]

%individual cell parameters
Pdrop = zeros(length(Hrange),length(Wrange));
HeatRate = zeros(length(Hrange),1);
YieldFlowRate = zeros(length(Hrange),1);
for hh = 1:length(Hrange)
    for ww = 1:length(aspectRrg)
        H = Hrange(hh);
        W = Wrange(ww);
        L = Lrange(ww); %Cell length = totalArea/totalWidth/M  [m]
        t_cell = celltime*(H^2)/(height^2); %Time taken for C profile to fully develop in a cell [s], for 1mm height, it takes SR_h s, and t is proportional to h^2
        Velocity = L/t_cell; %Max velocity for C profile to fully develop [m/s]

        % Velocity2 = 1000*miu/(rho*H); % Max velocity for laminar flow
        %
        % display(Velocity)
        % display(Velocity2)
        %
        % if Velocity>Velocity2
        %     Velocity = Velocity2;
        % end

        FR = (Velocity)*(H)*(W)*N/totalA; %Flow rate in the entire Burgers cascade [m^3/s per unit footprint area]

        %Pressure
        raspect = H/W; %aspect ratio [none] = H[m]/W[m]
        Dh = 2*(H)*(W/(H+W)); %Hydrodynamic diameter [m] = 2*H[m]*W[m]/(H[m]*W[m])
        Re = (rho*(Velocity)*Dh)/miu; %Reylaigh number = rho[kg/m^3]*V[m/s]*Dh[m] / miu[Pa.s]
        if Re>Recri
            Velocity = Recri*miu/(rho*Dh);%Re = rho*u*L_c/miu, Re_cri<=1000, thus u<Re_cri*miu/(rho*L_c)
        end
        f = 64/Re; %Friction factor, C = 64, for values of C, refer to Table 8.3 on p450
        delP_cell = f*rho*(Velocity)^2*(L)/(2*Dh); %Pressure drop in a cell [Pa] = f*rho[kg/m^3]*V[m/s]^2*L[m]/(2*Dh[m])
        delP_burgers = delP_cell*N*M; %Pressure drop in the entire Burgers cascade [Pa]
        Wpump = FR*delP_burgers; %Pump work [W] =  FR[m^3/s] * delP [Pa]

        %Matrix
        t_burgers = 1/(RR*FR*1000*3600); % [hr/L] Time [hour] taken to produce 1L of fresh(er) water for a certain recoery rate [hr/L]
        yield = RR*FR*1000*3600*24; %Yield [L/day]
        Ee = Wpump*t_burgers/1000; %Electrical energy consumed to produce 1L of fresh water [kWh/L]
        Lburgers = (L)*M; %Total length of the burgers [m]
        Wburgers = (W)*N; %Total width of the burgers [m]

        Pdrop(hh,ww) = delP_burgers;
        YieldFlowRate(hh,1) = yield;
    end
    %Heat flux
    q = (k*delT/(H))/1000; %kW/m^2 = (k[W/m.K]*delT[K]/H[m])/1000
    q = q/(R_T/100); % Scale the flux based on reduced linear T range
    %h_conv = (q*1000)/delT_conv; %heat transfer coefficient [W/K] = (q[kW/m^2]*1000)/delT[K];
    HeatRate(hh,1) = q;
end

Wrange = Wrange*1000;
Lrange = Lrange.*1000;
dimensiontxt = strings(length(aspectRrg),1); %String for pressure drop plot legend
for ww = length(aspectRrg):-1:1
    dimensiontxt(ww) = sprintf('r = %.0f with %.1f (%.1f) by %.1f (%.1f) LBC',aspectRrg(ww), WBurgers_rg(ww), Wrange(ww), LBurgers_rg(ww), Lrange(ww));
end
end

% Make a sound at the end of the script
Fs = 8192; % Sampling frequency (Hz)
duration = 0.5; % Duration of the tone (seconds)
frequency = 1000; % Frequency of the tone (Hz)
t = 0:1/Fs:duration; % Time vector
y = sin(2*pi*frequency*t); % Generate a sine wave tone
sound(y, Fs); % Play the tone

function plotCcontour(Mlong, Nfull, R_C, C0, targetRR, dCdev, RatioLow, RatioHigh)
[~, C] = BC_yieldC(C0, targetRR, R_C, Mlong, Nfull, dCdev, RatioLow, RatioHigh);
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
c = colorbar('Ticks',0:50:250);
c.Ruler.TickLabelFormat = '%.0f';
c.Label.String = '{\itC} (10^{3} ppm)';
c.Label.FontSize = 12;
%S.EdgeColor='none'; %shading interp %Smooth the color
contour(Y/2,M/2-X/2,C'*1000,[10:20:250],'LineColor','k','ShowText','on') %Add the black lines
ylabel("{\itm} for M_{long}")
xlabel("{\itn} for N_{full}")
set(gca, "Fontsize",12)
hold off
end





function [Cyield, C, Meq] = BC_yieldC(C0,targetRR, R_C, Mlong, Nfull, dCdev, RatioLow, RatioHigh)
% Assumptions:
% Perfect bifurcation at the horizontal midplane
% dC in each cell can be scaled based on:
% > dC is propotional to dT
% > dC is propotional to C0(1-C0), where C0 the initial C at the cell inlet

dC = dCdev*R_C; % achieved dC is smaller due to not long enough time

if (Mlong<2)
    Mlong = 2;
end
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
                dC2 = RatioLow*R_C*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1) + Velocity(m-1,n+1)/2;
                C(m,n) = (C(m-1,n-1)*Velocity(m-1,n-1) + (C(m-1,n+1)-dC2/4)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            elseif n == 2*Nfull
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1);
                C(m,n) = ((C(m-1,n-1)+dC1/4)*Velocity(m-1,n-1)/2 + C(m-1,n+1)*Velocity(m-1,n+1))/Velocity(m,n);
            else
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1)/2;
                C(m,n) = ((C(m-1,n-1)+dC1/4)*Velocity(m-1,n-1)/2+(C(m-1,n+1)-dC2/4)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            end
        end
        % mode(m,2)==0 is equivalent to i=n in manuscript, it has N-1 full cells
        % and 2 half cells
    elseif mod(m,2) == 0
        for n=1:2:N
            if n==1
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n+1)/2;
                C(m,n) = C(m-1,n+1)-dC2/4;
            elseif n == N
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2;
                C(m,n) = C(m-1,n-1)+dC1/4;
            else
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1)/2;
                C(m,n) = ((C(m-1,n-1)+dC1/4)*Velocity(m-1,n-1)/2+(C(m-1,n+1)-dC2/4)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            end
        end
    end
end

% Calculate C_yield based on recovery rate
SamplingCell = fix(N*targetRR);
if SamplingCell < 1
    SamplingCell = 1;
end
higherM = 0;
higherV = 0;
% display(M)
for n = N:-1:(N-SamplingCell+1)
    higherM = higherM + C(M,n)*Velocity(M,n);
    higherV = higherV + Velocity(M,n);
end
higherC = higherM/higherV;
Cyield= higherC;

% Fing Meq that satisfies C_drop(j)>= 95% C_drop(equlibrium)
higherC = zeros(1,M);
higherRR = zeros(1,M);

for m = 2:2:M %each row
    higherM = 0;
    higherV = 0;
    %for n = N:-1:(N-SamplingCell)
    for n = N:-1:(N-SamplingCell)
        higherM = higherM + C(m,n)*Velocity(m,n);
        higherV = higherV + Velocity(m,n);
    end
    higherC(m) = higherM/higherV;
    higherC(m-1) = higherC(m);
    higherRR(m) = higherV/Nfull;
    higherRR(m-1) = higherRR(m);
end

Cincre = higherC-C0;
Meq = find(Cincre>Cincre(end)*0.95,1);
Meq = Meq/2+0.5; %Convert to full-length cell
% display(Cyield)

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

function [Cyield, C, Meq] = BC_yieldC_TDD(C0,targetRR, SR, Mlong, Nfull, dCdev, RatioLow, RatioHigh)
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
                C(m,n) = (C(m-1,n-1)*Velocity(m-1,n-1) + (C(m-1,n+1)-dC2/4)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            elseif n == 2*Nfull
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1);
                C(m,n) = ((C(m-1,n-1)+dC1/4)*Velocity(m-1,n-1)/2 + C(m-1,n+1)*Velocity(m-1,n+1))/Velocity(m,n);
            else
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1)/2;
                C(m,n) = ((C(m-1,n-1)+dC1/4)*Velocity(m-1,n-1)/2+(C(m-1,n+1)-dC2/4)*Velocity(m-1,n+1)/2)/Velocity(m,n);
            end
        end
        % mode(m,2)==0 is equivalent to i=n in manuscript, it has N-1 full cells
        % and 2 half cells
    elseif mod(m,2) == 0
        for n=1:2:N
            if n==1
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n+1)/2;
                C(m,n) = C(m-1,n+1)-dC2/4;
            elseif n == N
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2;
                C(m,n) = C(m-1,n-1)+dC1/4;
            else
                dC1 = RatioHigh*dC*C(m-1,n-1)*(1-C(m-1,n-1))/(C0*(1-C0));
                dC2 = RatioLow*dC*C(m-1,n+1)*(1-C(m-1,n+1))/(C0*(1-C0));
                Velocity(m,n) = Velocity(m-1,n-1)/2 + Velocity(m-1,n+1)/2;
                C(m,n) = ((C(m-1,n-1)+dC1/4)*Velocity(m-1,n-1)/2+(C(m-1,n+1)-dC2/4)*Velocity(m-1,n+1)/2)/Velocity(m,n);
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


