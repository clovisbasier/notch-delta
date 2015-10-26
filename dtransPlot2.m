function dtransPlot2()

% Plot 2 - S,Dtrans,t

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

sp2_S = S;
clear betaD betaN D0 dt Dtrans gamma gammaS k kc kt N0 Ntrans S S0 tf treshold x y


end

