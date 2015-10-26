%% Reset
clear
clc
clf

%% Plot 1 - S,Dtrans

sp1 = subplot(2,2,1);

D0 = 200;                %free cis-Delta
S0 = 0;                  %intracellular domain of Notch

tf = 60;                 %t final
dt = 0.001;              %step
%tfs = tf/dt - 1;        %tf ajusted for step

betaN = 1;               %production rate of Notch
betaD = 0;               %production rate of Delta
gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation
Ntrans = 0;              %Notch in neighboring cell

for Dtrans = [0.063 0.084 0.11 0.15 0.20 0.26 0.35 0.46 0.62 0.82 1.1 1.4]
     
    N0 = betaN./(gamma + D0./kc + Dtrans./kt);                             %free Notch
    
    y = elow(N0,D0,S0,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans);
    S = y(3,:);
    x = y(4,:);
    pS = plot(x,S);               %plot S
    hold on
end
hold off

title('')
legend('0.063 ug.ml^-^1','0.084 ug.ml^-^1','0.11 ug.ml^-^1','0.15 ug.ml^-^1',...
    '0.20 ug.ml^-^1','0.26 ug.ml^-^1','0.35 ug.ml^-^1','0.46 ug.ml^-^1',...
    '0.62 ug.ml^-^1','0.82 ug.ml^-^1','1.1 ug.ml^-^1','1.4 ug.ml^-^1'...
    ,'Location','northwest')
xlabel('Time (h)')
ylabel('Reporter (10^4 a.u.)')

%% Plot 2 - S,Dtrans,t

sp2 = subplot(1,2,2);

D0 = 200;                %free cis-Delta
S0 = 0;                  %intracellular domain of Notch

tf = 60;                 %t final
dt = 0.001;              %step
%tfs = tf/dt - 1;        %tf ajusted for step

betaN = 1;               %production rate of Notch
betaD = 0;               %production rate of Delta
gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation
Ntrans = 0;              %Notch in neighboring cell


for Dtrans = [0.063 0.084 0.11 0.15 0.20 0.26 0.35 0.46 0.62 0.82 1.1 1.4] %trans-Delta = concentration of neighboring Delta = Dplate
    
    N0 = betaN./(gamma + D0./kc + Dtrans./kt);                             %free Notch
    
    y = elow(N0,D0,S0,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans); %run Euler aproximation
    S = y(3,:);
    x = y(4,:);
    
    treshold = 0.25;                                                       %determination of treshold point S(treshold)
    k = size(S,2);                                                         
    while S(k) > treshold && k > 0
        k = k-1;
    end
    
    scatter3(ones(1,size(S,2)-k).*Dtrans,x(1,k+1:size(S,2)),...            %plot S part above treshold in green
        S(1,k+1:size(S,2)),0.1,'.','green');
    hold on
    scatter3(ones(1,k).*Dtrans,x(1,1:k),S(1,1:k),0.1,'.','black');         %plot S part under treshold in black
    hold on
end
hold off


title('Notch response to both {\itcis-} and {\ittrans-}Delta')
xlabel('Dplate (ug.ml^-^1)')
ylabel('Time (h)')
ylim([0 tf])
zlabel('Reporter (10^4 a.u.)')

%% Plot 3 - receiver/sender

sp3 = subplot(2,2,3);

tf = 15;                %t final
dt = 0.001;             %step
tfs = tf/dt;            %tf ajusted for step

betaN1 = 10;             %production rate of Notch cell 1
betaD1 = 8;              %production rate of Delta cell 1
betaN2 = 9;              %production rate of Notch cell 2
betaD2 = 10;             %production rate of Delta cell 2
gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation

D = 0;                %free cis-Delta
S = 0;                %intracellular domain of Notch
N = 0;                %free Notch

N1 = [N zeros(1,tfs)];   %preallocating space
D1 = [D zeros(1,tfs)];   %preallocating space

N2 = [N zeros(1,tfs)];   %preallocating space
D2 = [D zeros(1,tfs)];   %preallocating space

for k = 1:tfs
    y = elowSS(N1(k),D1(k),S,dt,betaN1,betaD1,gamma,gammaS,kc,kt,N2(k),D2(k)); %cell1
    N1(k+1) = y(1);
    D1(k+1) = y(2);
    S1(k+1) = y(3);
    
    y = elowSS(N2(k),D2(k),S,dt,betaN2,betaD2,gamma,gammaS,kc,kt,N1(k),D1(k)); %cell2
    N2(k+1) = y(1);
    D2(k+1) = y(2);
    S2(k+1) = y(3);
end

S1S2 = zeros(tf,2);   %preallocating space
for k = 1:tf
    S1S2(k,1) = S1(k+1000*(k-1));
    S1S2(k,2) = S2(k+1000*(k-1));
end

x = 0:dt:tf-dt;                %x axis scale

colormap('winter')
imagesc(S1S2)

title('Receiver - sender (1D, 2 cells)')
xlabel('Cells (1D)')
ylabel('Time (h)')
c = colorbar('westoutside');
c.Label.String = 'Reporter (10^4 a.u.)';


%% Plot  - 21 cells 1D

sp4 = subplot(4,2,1);

betaDp = [20 zeros(1,10)];        %betaD in 1D (linear decay)
for k = 1:10
    betaDp(k+1) = betaDp(k) - 2;
end
betaD = zeros(1,21);                %preallocating space
for k = 1:10
    betaD(k) = betaDp(12-k);
    betaD(k+10) = betaDp(k);
end
betaD(21) = betaDp(11);

betaN = ones(1,21).*10;             %betaN in 1D (constant)

tf = 15;                 %t final
dt = 0.001;             %step
tfs = tf/dt - 1;        %tf ajusted for step

gamma = 0.1;             %combined degradation and dilution rate
gammaS = 0.1;            %rate of decay of S
kc = 0.2;                %strenght of cis-inhibition
kt = 2;                  %strenght of transactivation
Ntrans = 0;              %Notch in neighboring cell

D = zeros(tf/dt,21);          %free cis-Delta
S = zeros(tf/dt,21);          %intracellular domain of Notch
N = zeros(tf/dt,21);          %free Notch

for k = 1:tfs
    for l = 1:21
        if l == 1
            y = elowSS(N(k,l),D(k,l),S(k,l),dt,betaN(l),betaD(l),gamma,gammaS,...
                kc,kt,N(k,l+1),D(k,l+1));
            N(k+1,l) = y(1);
            D(k+1,l) = y(2);
            S(k+1,l) = y(3);
        elseif l == 21
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

S21 = zeros(tf,21);   %preallocating space
for k = 1:tf
    for l = 1:21
        S21(k,l) = S(k+1000*(k-1),l);
    end
end


plot(1:21,betaN)
hold on
plot(1:21,betaD)

xlim([1 21])

title('1D simulation, 21 cells')
ylabel('Production rate (a.u.)')
legend('betaN0','betaD0')

sp5 = subplot(4,2,3);

colormap('winter')
imagesc(S21)
xlabel('Cells (1D)')
ylabel('Time (h)')

