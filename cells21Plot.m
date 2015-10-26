function cells21Plot
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
end