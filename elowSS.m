function y = elowSS(N,D,S,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans)


    %% Equation (Euler approximation)

        %prod Notch - degradation - cis-inhibiton - transactivation of Notch
        dN = ( betaN - gamma.*N - D.*N./kc - Dtrans.*N./kt ).*dt;
        %prod cis-Detla - degradation - cis-inhibiton - transactivation of Delta        
        dD = ( betaD - gamma.*D - D.*N./kc - D.*Ntrans./kt ).*dt;
        %transactivation - degradation
        dS = ( Dtrans.*N./kt - gammaS.*S ).*dt;

        N2 = N + dN;     %approximation affine N
        D2 = D + dD;     %approximation affine D
        S2 = S + dS;     %approximation affine S
        
    y = [N2 D2 S2];