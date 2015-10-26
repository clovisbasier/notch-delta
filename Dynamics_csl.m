%% Reset
clear
clc
clf

%% (Plot 1 - S,Dtrans)

% sp1 = subplot(2,2,1);
% 
% D0 = 200;                %free cis-Delta
% S0 = 0;                  %intracellular domain of Notch
% 
% tf = 60;                 %t final
% dt = 0.001;              %step
% tfs = tf/dt - 1;        %tf ajusted for step
% 
% betaN = 1;               %production rate of Notch
% betaD = 0;               %production rate of Delta
% gamma = 0.1;             %combined degradation and dilution rate
% gammaS = 0.1;            %rate of decay of S
% kc = 0.2;                %strenght of cis-inhibition
% kt = 2;                  %strenght of transactivation
% Ntrans = 0;              %Notch in neighboring cell
% 
% for Dtrans = [0.063 0.084 0.11 0.15 0.20 0.26 0.35 0.46 0.62 0.82 1.1 1.4]
%      
%     N0 = betaN./(gamma + D0./kc + Dtrans./kt);                             %free Notch
%     
%     y = elow(N0,D0,S0,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans);
%     S = y(3,:);
%     x = y(4,:);
%     pS = plot(x,S);               %plot S
%     hold on
% end
% hold off
% 
% title('')
% legend('0.063 ug.ml^-^1','0.084 ug.ml^-^1','0.11 ug.ml^-^1','0.15 ug.ml^-^1',...
%     '0.20 ug.ml^-^1','0.26 ug.ml^-^1','0.35 ug.ml^-^1','0.46 ug.ml^-^1',...
%     '0.62 ug.ml^-^1','0.82 ug.ml^-^1','1.1 ug.ml^-^1','1.4 ug.ml^-^1'...
%     ,'Location','northwest')
% xlabel('Time (h)')
% ylabel('Reporter (10^4 a.u.)')

%% Plot 2 - S,Dtrans,t

% sp2 = subplot(1,2,2);
% 
% D0 = 200;                %free cis-Delta
% S0 = 0;                  %intracellular domain of Notch
% 
% tf = 60;                 %t final
% dt = 0.001;              %step
% %tfs = tf/dt - 1;        %tf ajusted for step
% 
% betaN = 1;               %production rate of Notch
% betaD = 0;               %production rate of Delta
% gamma = 0.1;             %combined degradation and dilution rate
% gammaS = 0.1;            %rate of decay of S
% kc = 0.2;                %strenght of cis-inhibition
% kt = 2;                  %strenght of transactivation
% Ntrans = 0;              %Notch in neighboring cell
% 
% 
% for Dtrans = [0.063 0.084 0.11 0.15 0.20 0.26 0.35 0.46 0.62 0.82 1.1 1.4] %trans-Delta = concentration of neighboring Delta = Dplate
%     
%     N0 = betaN./(gamma + D0./kc + Dtrans./kt);                             %free Notch
%     
%     y = elow(N0,D0,S0,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans); %run Euler aproximation
%     S = y(3,:);
%     x = y(4,:);
%     
%     treshold = 0.25;                                                       %determination of treshold point S(treshold)
%     k = size(S,2);                                                         
%     while S(k) > treshold && k > 0
%         k = k-1;
%     end
%     
%     scatter3(ones(1,size(S,2)-k).*Dtrans,x(1,k+1:size(S,2)),...            %plot S part above treshold in green
%         S(1,k+1:size(S,2)),0.1,'.','green');
%     hold on
%     scatter3(ones(1,k).*Dtrans,x(1,1:k),S(1,1:k),0.1,'.','black');         %plot S part under treshold in black
%     hold on
% end
% hold off
% 
% 
% title('Notch response to both {\itcis-} and {\ittrans-}Delta')
% xlabel('Dplate (ug.ml^-^1)')
% ylabel('Time (h)')
% ylim([0 tf])
% zlabel('Reporter (10^4 a.u.)')
% 
% sp2_S = S;
% clear betaD betaN D0 dt Dtrans gamma gammaS k kc kt N0 Ntrans S S0 tf treshold x y

%% (Plot 3 - receiver/sender)
% 
% sp3 = subplot(2,2,3);
% 
% tf = 15;                %t final
% dt = 0.001;             %step
% tfs = tf/dt;            %tf ajusted for step
% 
% betaN1 = 10;             %production rate of Notch cell 1
% betaD1 = 8;              %production rate of Delta cell 1
% betaN2 = 9;              %production rate of Notch cell 2
% betaD2 = 10;             %production rate of Delta cell 2
% gamma = 0.1;             %combined degradation and dilution rate
% gammaS = 0.1;            %rate of decay of S
% kc = 0.2;                %strenght of cis-inhibition
% kt = 2;                  %strenght of transactivation
% 
% D = 0;                %free cis-Delta
% S = 0;                %intracellular domain of Notch
% N = 0;                %free Notch
% 
% N1 = [N zeros(1,tfs)];   %preallocating space
% D1 = [D zeros(1,tfs)];   %preallocating space
% S1 = [S zeros(1,tfs)];   %preallocating space
% 
% 
% N2 = [N zeros(1,tfs)];   %preallocating space
% D2 = [D zeros(1,tfs)];   %preallocating space
% S2 = [S zeros(1,tfs)];   %preallocating space
% 
% 
% for k = 1:tfs
%     y = elowSS(N1(k),D1(k),S,dt,betaN1,betaD1,gamma,gammaS,kc,kt,N2(k),D2(k)); %cell1
%     N1(k+1) = y(1);
%     D1(k+1) = y(2);
%     S1(k+1) = y(3);
%     
%     y = elowSS(N2(k),D2(k),S,dt,betaN2,betaD2,gamma,gammaS,kc,kt,N1(k),D1(k)); %cell2
%     N2(k+1) = y(1);
%     D2(k+1) = y(2);
%     S2(k+1) = y(3);
% end
% 
% S1S2 = zeros(tf,2);   %preallocating space
% for k = 1:tf
%     S1S2(k,1) = S1(k+1000*(k-1));
%     S1S2(k,2) = S2(k+1000*(k-1));
% end
% 
% x = 0:dt:tf-dt;                %x axis scale
% 
% colormap('winter')
% imagesc(S1S2)
% 
% title('Receiver - sender (1D, 2 cells)')
% xlabel('Cells (1D)')
% ylabel('Time (h)')
% c = colorbar('westoutside');
% c.Label.String = 'Reporter (10^4 a.u.)';

%% Plot  - 21 cells 1D

sp4 = subplot(4,2,1);

betaDp = [20 zeros(1,10)];        %betaD in 1D (linear decay)
for k = 1:10
    betaDp(k+1) = betaDp(k) - 2;
end
betaD = zeros(1,21);%preallocating space
for k = 1:10
    betaD(k) = betaDp(12-k);
    betaD(k+10) = betaDp(k);
end
betaD(21) = betaDp(11);

betaN = ones(1,21).*10;             %betaN in 1D (constant)

tf = 15;                %t final
dt = 0.001;             %step
tfs = tf/dt - 1;        %tf ajusted for step

gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation

D = zeros(tf/dt,21);     %free cis-Delta
S = zeros(tf/dt,21);     %intracellular domain of Notch
N = zeros(tf/dt,21);     %free Notch

for k = 1:tfs
    for l = 1:21
        %cells interact with left (l-1) and right (l+1) neighbors
        if l == 1                       %exception : left
            y = elowSS(N(k,l),D(k,l),S(k,l),dt,betaN(l),betaD(l),gamma,gammaS,...
                kc,kt,N(k,l+1),D(k,l+1));           
            N(k+1,l) = y(1);
            D(k+1,l) = y(2);
            S(k+1,l) = y(3);
        elseif l == 21                  %exception : right
            y = elowSS(N(k,l),D(k,l),S(k,l),dt,betaN(l),betaD(l),gamma,gammaS,...
                kc,kt,N(k,l-1),D(k,l-1));           
            N(k+1,l) = y(1);                        
            D(k+1,l) = y(2);
            S(k+1,l) = y(3);
        else
            y = elowSS(N(k,l),D(k,l),S(k,l),dt,betaN(l),betaD(l),gamma,gammaS,...
                kc,kt,N(k,l-1) + N(k,l+1),D(k,l-1) + D(k,l+1));
            N(k+1,l) = y(1);
            D(k+1,l) = y(2);
            S(k+1,l) = y(3);
        end
    end
end

S21 = zeros(tf,21);%preallocating space
for k = 1:tf
    for l = 1:21      %computing time steps (each hour)
        S21(k,l) = S(k+1000*(k-1),l);
    end
end


plot(1:21,betaN)                %plot betaN for each cell
hold on
plot(1:21,betaD)                %plot betaD for each cell

xlim([1 21])

title('1D simulation, 21 cells')
ylabel('Production rate (a.u.)')
legend('betaN0','betaD0')

sp5 = subplot(4,2,3);

colormap('winter')
imagesc(S21)                    %plot evolution S in time
xlabel('Cells (1D)')
ylabel('Time (h)')

sp5_D = D;
sp5_S = S;
sp5_N = N;
sp5_S21 = S21;
clear betaD betaDp betaN D dt gamma gammaS k kc kt I N S S21 tf tfs y

%% Plot  - 21*4 cells 2D

sp6 = subplot(2,2,3);

betaDp = [20 zeros(1,10)];          %betaD in 1D (linear decay)
for k = 1:10
    betaDp(k+1) = betaDp(k) - 2;
end
betaD = zeros(1,21);%preallocating space
for k = 1:10
    betaD(k) = betaDp(12-k);
    betaD(k+10) = betaDp(k);
end
betaD(21) = betaDp(11);

betaN = ones(1,21).*10;             %betaN in 1D (constant)

for k = 1:4                         %same betaD and betaN for cell on a column
    betaD(k,:) = betaD(1,:);
    betaN(k,:) = betaN(1,:);
end

tf = 15;                %t final
dt = 0.001;             %step
tfs = tf/dt - 1;        %tf ajusted for step

gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation

D = zeros(4,21,tf/dt);   %free cis-Delta
S = zeros(4,21,tf/dt);   %intracellular domain of Notch
N = zeros(4,21,tf/dt);   %free Notch

for k = 1:tfs
    for l = 1:21
        for m = 1:4
            %cells interact with left (l-1), right (l+1), top (m-1), bottom (m+1) neighbors
            if m == 1 && l == 1             %exception : top left
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m+1,l,k) + N(m,l+1,k),...
                    D(m+1,l,k) + D(m,l+1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 1 && l == 21        %exception : top right
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m+1,l,k) + N(m,l-1,k),...
                    D(m+1,l,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 1                   %exception : top
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m+1,l,k) + N(m,l+1,k) + N(m,l-1,k),...
                    D(m+1,l,k) + D(m,l+1,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 4 && l == 1         %exception : bottom left
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m,l+1,k),...
                    D(m-1,l,k) + D(m,l+1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 4 && l == 21        %exception : bottom right
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 4                   %exception : bottom
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m,l+1,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m,l+1,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif l == 1                   %exception : left
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m+1,l,k) + N(m,l+1,k),...
                    D(m-1,l,k) + D(m+1,l,k) + D(m,l+1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif l == 21                  %exception : right
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m+1,l,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m+1,l,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            else
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m+1,l,k) + N(m,l+1,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m+1,l,k) + D(m,l+1,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            end
        end
    end
end

S21_4 = S(:,:,tfs+1);

colormap('winter')
imagesc(S21_4)                    %plot evolution S final
title('2D simulation, 21*4 cells')
xlabel('Cells')
ylabel('Cells')

sp6_D = D;
sp6_S = S;
sp6_N = N;
sp6_S21_4 = S21_4;

clear betaD betaDp betaN D dt gamma gammaS k kc kt I N S S21 tf tfs y m l

%% Plot  - 21*21 cells 2D

sp7 = subplot(2,2,4);


betaD = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;
         1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1;
         1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1;
         1 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 2 1;
         1 2 3 4 5 5 5 5 5 5 5 5 5 5 5 5 5 4 3 2 1;
         1 2 3 4 5 6 6 6 6 6 6 6 6 6 6 6 5 4 3 2 1;
         1 2 3 4 5 6 7 7 7 7 7 7 7 7 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 8 8 8 8 8 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 9 9 9 9 9 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 9 10 10 10 9 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 9 10 11 10 9 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 9 10 10 10 9 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 9 9 9 9 9 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 8 8 8 8 8 8 8 7 6 5 4 3 2 1;
         1 2 3 4 5 6 7 7 7 7 7 7 7 7 7 6 5 4 3 2 1;
         1 2 3 4 5 6 6 6 6 6 6 6 6 6 6 6 5 4 3 2 1;
         1 2 3 4 5 5 5 5 5 5 5 5 5 5 5 5 5 4 3 2 1;
         1 2 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 2 1;
         1 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 1;
         1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1;
         1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
betaN = ones(21,21).*5;

tf = 15;                %t final
dt = 0.001;             %step
tfs = tf/dt - 1;        %tf ajusted for step

gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation

D = zeros(21,21,tf/dt);   %free cis-Delta
S = zeros(21,21,tf/dt);   %intracellular domain of Notch
N = zeros(21,21,tf/dt);   %free Notch

for k = 1:tfs
    for l = 1:21
        for m = 1:21
            %cells interact with left (l-1), right (l+1), top (m-1), bottom (m+1) neighbors
            if m == 1 && l == 1             %exception : top left
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m+1,l,k) + N(m,l+1,k),...
                    D(m+1,l,k) + D(m,l+1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 1 && l == 21        %exception : top right
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m+1,l,k) + N(m,l-1,k),...
                    D(m+1,l,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 1                   %exception : top
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m+1,l,k) + N(m,l+1,k) + N(m,l-1,k),...
                    D(m+1,l,k) + D(m,l+1,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 21 && l == 1         %exception : bottom left
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m,l+1,k),...
                    D(m-1,l,k) + D(m,l+1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 21 && l == 21        %exception : bottom right
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif m == 21                   %exception : bottom
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m,l+1,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m,l+1,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif l == 1                   %exception : left
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m+1,l,k) + N(m,l+1,k),...
                    D(m-1,l,k) + D(m+1,l,k) + D(m,l+1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            elseif l == 21                  %exception : right
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m+1,l,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m+1,l,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            else
                y = elowSS(N(m,l,k),D(m,l,k),S(m,l,k),...
                    dt,betaN(m,l),betaD(m,l),gamma,gammaS,kc,kt,...
                    N(m-1,l,k) + N(m+1,l,k) + N(m,l+1,k) + N(m,l-1,k),...
                    D(m-1,l,k) + D(m+1,l,k) + D(m,l+1,k) + D(m,l-1,k));           
                N(m,l,k+1) = y(1);
                D(m,l,k+1) = y(2);
                S(m,l,k+1) = y(3);
            end
        end
    end
end
S21_4 = S(:,:,tfs+1);

betaN_D = zeros(21,21);
for l = 1:21
    for m = 1:21
        if m == 1 && l == 1             %exception : top left
            betaN_D(m,l) = (betaN(m,l) + betaD(m+1,l)/4 + betaD(m,l+1)/4)...
                          /(betaD(m,l) + betaN(m+1,l)/4 + betaN(m,l+1)/4);          
        elseif m == 1 && l == 21        %exception : top right
            betaN_D(m,l) = (betaN(m,l) + betaD(m+1,l)/4 + betaD(m,l-1)/4)...
                          /(betaD(m,l) + betaN(m+1,l)/4 + betaN(m,l-1)/4);
        elseif m == 1                   %exception : top
            betaN_D(m,l) = (betaN(m,l) + betaD(m+1,l)/4 + betaD(m,l+1)/4 + betaD(m,l-1)/4)...
                          /(betaD(m,l) + betaN(m+1,l)/4 + betaN(m,l+1)/4 + betaN(m,l-1)/4);
        elseif m == 21 && l == 1         %exception : bottom left
            betaN_D(m,l) = (betaN(m,l) + betaD(m-1,l)/4 + betaD(m,l+1)/4)...
                          /(betaD(m,l) + betaN(m-1,l)/4 + betaN(m,l+1)/4);
        elseif m == 21 && l == 21        %exception : bottom right
            betaN_D(m,l) = (betaN(m,l) + betaD(m-1,l)/4 + betaD(m,l-1)/4)...
                          /(betaD(m,l) + betaN(m-1,l)/4 + betaN(m,l-1)/4);
        elseif m == 21                   %exception : bottom
            betaN_D(m,l) = (betaN(m,l) + betaD(m-1,l)/4 + betaD(m,l+1)/4 + betaD(m,l-1)/4)...
                          /(betaD(m,l) + betaN(m-1,l)/4 + betaN(m,l+1)/4 + betaN(m,l-1)/4);
        elseif l == 1                   %exception : left
            betaN_D(m,l) = (betaN(m,l) + betaD(m-1,l)/4 + betaD(m+1,l)/4 + betaD(m,l+1)/4)...
                          /(betaD(m,l) + betaN(m-1,l)/4 + betaN(m+1,l)/4 + betaN(m,l+1)/4);
        elseif l == 21                  %exception : right
            betaN_D(m,l) = (betaN(m,l) + betaD(m+1,l)/4 + betaD(m-1,l)/4 + betaD(m,l-1)/4)...
                          /(betaD(m,l) + betaN(m+1,l)/4 + betaN(m-1,l)/4 + betaN(m,l-1)/4);
        else
            betaN_D(m,l) = (betaN(m,l) + betaN(m+1,l) + betaN(m-1,l) + betaN(m,l+1) + betaN(m,l-1))...
                          /(betaD(m,l) + betaD(m+1,l) + betaD(m-1,l) + betaD(m,l+1) + betaD(m,l-1));
        end
    end
end

colormap('winter')
imagesc(S21_4)                    %plot evolution S final
xlabel('Cells')
ylabel('Cells')
title('2D simulation, 21*21 cells')
c = colorbar;
c.Label.String = 'Reporter (10^4 a.u.)';

% sp8 = subplot(4,2,4);
% colormap('winter')
% imagesc(betaN_D)
% xlabel('Cells')
% ylabel('Cells')
% title('betaN / betaD (neighborhood)')
% colorbar;

sp9 = subplot(2,2,2);
colormap(sp9,'cool')
imagesc(betaD-betaN)
ylabel('Cells')
title('betaD - betaN')
colorbar();

sp7_D = D;
sp7_S = S;
sp7_N = N;
sp7_S21_4 = S21_4;
sp8_betaN_D = betaN_D;
sp9_betaN = betaN;
sp9_betaD = betaD;

clear D dt gamma gammaS k kc kt I N S S21 tf tfs y m l

