function cells21x21Plot()
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

end