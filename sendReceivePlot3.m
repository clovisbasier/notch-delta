function sendReceivePlot3()

    %% (Plot 3 - receiver/sender)

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
    S1 = [S zeros(1,tfs)];   %preallocating space


    N2 = [N zeros(1,tfs)];   %preallocating space
    D2 = [D zeros(1,tfs)];   %preallocating space
    S2 = [S zeros(1,tfs)];   %preallocating space


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

end