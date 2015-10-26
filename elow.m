function y = elow(N,D,S,tf,dt,betaN,betaD,gamma,gammaS,kc,kt,Ntrans,Dtrans)

    tfs = tf/dt - 1;               %tf ajusted for step

    %% Equation (Euler approximation)

    N = [N zeros(1,tfs)];   %preallocating space
    D = [D zeros(1,tfs)];   %preallocating space
    S = [S zeros(1,tfs)];   %preallocating space

    for k = 1:tfs
        %prod Notch - degradation - cis-inhibiton - transactivation of Notch
        dN = ( betaN - gamma.*N(k) - D(k).*N(k)./kc - Dtrans.*N(k)./kt ).*dt;
        %prod cis-Detla - degradation - cis-inhibiton - transactivation of Delta        
        dD = ( betaD - gamma.*D(k) - D(k).*N(k)./kc - D(k).*Ntrans./kt ).*dt;
        %transactivation - degradation
        dS = ( Dtrans.*N(k)./kt - gammaS.*S(k) ).*dt;

        N(k+1) = N(k) + dN;     %approximation affine N
        D(k+1) = D(k) + dD;     %approximation affine D
        S(k+1) = S(k) + dS;     %approximation affine S
    end
    
    x = 0:dt:tf-dt;                %x axis scale
    
    y = [N;D;S;x];