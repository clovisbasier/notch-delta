function dtransPlot1

    % (Plot 1 - S,Dtrans)

    sp1 = subplot(2,2,1);

    D0 = 200;                %free cis-Delta
    S0 = 0;                  %intracellular domain of Notch

    tf = 60;                 %t final
    dt = 0.001;              %step
    tfs = tf/dt - 1;        %tf ajusted for step

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
 
end