%% Reset
clear
clc
clf

%% Parameters

%N = 200;                       %free Notch
%D = 200;                       %free cis-Delta
%S = 0;                         %intracellular domain of Notch

%tf = 60;                        %t final
%dt = 0.001;                    %step
%tfs = tf/dt - 1;               %tf ajusted for step

%x = 0:dt:tf-dt;                %x axis scale

%betaN = 1;                    %production rate of Notch
%betaD = 1;                    %production rate of Delta
%gamma = 1;                     %combined degradation and dilution rate
%gammaS = 1;                    %rate of decay of S
%kc = 2;                        %strenght of cis-inhibition
%kt = 1;                        %strenght of transactivation
%Ntrans = 0;                    %Notch in neighboring cell
%Dtrans = 1.45;                 %trans-Delta = concentration of neighboring Delta = Dplate

%% Plot 1 - N D S,t

%elow(N,D,S,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans)
y = elow(0,0,0,10,0.001,1,0,0.1,0.1,0.2,2,0,0.62);

N = y(1,:);
D = y(2,:);
S = y(3,:);
x = y(4,:);

sp1 = subplot(2,1,1);
pN = plot(x,N);               %plot N
hold on
pD = plot(x,D);               %plot D
pS = plot(x,S);               %plot S
hold off

title('')
legend('free Notch','free cis-Delta','Reporter')
xlabel('Time (h)')

pN.LineWidth = 2;
pD.LineWidth = 2;
pS.LineWidth = 2;

%% Plot 2 - S,Dtrans

sp2 = subplot(2,2,3);
for Dtrans = [0.063 0.084 0.11 0.15 0.20 0.26 0.35 0.46 0.62 0.82 1.1 1.4]
    %elow(N,D,S,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans)
    y = elow(0,0,0,10,0.001,1,0,0.1,0.1,0.2,2,0,Dtrans);
    S = y(3,:);
    x = y(4,:);
    pS = plot(x,S);               %plot S
    hold on
end
hold off

title('')
legend('0.063 ug.ml^-^1','0.084 ug.ml^-^1','0.11 ug.ml^-^1','0.15 ug.ml^-^1',...
    '0.20 ug.ml^-^1','0.26 ug.ml^-^1','0.35 ug.ml^-^1','0.46 ug.ml^-^1',...
    '0.62 ug.ml^-^1','0.82 ug.ml^-^1','1.1 ug.ml^-^1','1.4 ug.ml^-^1')
xlabel('Time (h)')
ylabel('Reporter (10^4 a.u.)')

pN.LineWidth = 2;
pD.LineWidth = 2;
pS.LineWidth = 2;

%% Plot 3 - S,Dtrans,t

sp4 = subplot(2,2,4);

for Dtrans = [0.063 0.084 0.11 0.15 0.20 0.26 0.35 0.46 0.62 0.82 1.1 1.4]
    %elow(N,D,S,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans)
    y = elow(0,200,0,10,0.001,1,0,0.1,0.1,0.2,2,0,Dtrans);
    S = y(3,:);
    x = y(4,:);
    
    k = 1;
    while S(k) < 0.05 && k < size(S,2)
        k = k + 1;
    end
    
    scatter3(ones(1,k).*Dtrans,x(1,1:k),S(1,1:k),0.1,'.','black');    %plot S
    hold on
    scatter3(ones(1,size(S,2)-k).*Dtrans,x(1,1:size(S,2)-k),S(1,1:size(S,2)-k),0.1,'.','green');
end
hold off


title('Notch response to both {\itcis-} and {\ittrans-}Delta')
xlabel('Dplate (ug.ml^-^1)')
ylabel('Time (h)')
ylim([0 5])
zlabel('Reporter (10^4 a.u.)')
