function cells21x4Plot()

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
end