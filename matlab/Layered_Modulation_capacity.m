clear all;
close all ;

alpha = 1
gamma = 0:0.01:1 ;

beta = gamma.*alpha ;

SNR_dB = 0:3:30 ;
SNR = 10.^(SNR_dB./10);

theta = 0:0.01:1 ;
theta = theta .* (pi/4) ;

beta0 = 0.25  ;

N_theta = length( theta ) ;

N_1 = 4 ;
L1 = zeros( N_1, 1 ) ;
N_2 = 4 ;
L2 = zeros( N_2, 1 ) ;
N_s = N_1 * N_2 ;
s = zeros( N_s, N_theta ) ;

        A1 = alpha ;

        L1(1) = +1 + j ;
        L1(2) = -1 + j ;
        L1(3) = -1 - j ;
        L1(4) = +1 - j ;
        
        L11 = L1.*A1 ;


RUNS = 5000 ;

A2 = sqrt(beta0) ;
for n = 1 : length(theta)
    
    x(n) = cos(pi/4-theta(n)) * sqrt(2) ;
    y(n) = sin(pi/4-theta(n)) * sqrt(2) ;
    
    L2(1, n) = +x(n) + j*y(n) ;
    L2(2, n) = -x(n) + j*y(n) ;
    L2(3, n) = -x(n) - j*y(n) ;
    L2(4, n) = +x(n) - j*y(n) ;
    
    for p = 0 : (N_1 - 1)
        for q = 1 : N_2
            k = p*N_1 + q ;
            s(k,n) = L1(p+1).*A1 + L2(q, n).*A2 ;
        end
    end     
end

    L20(1, :) = L1(1) * A2 * ones(1,N_theta) ; 
    L20(2, :) = L1(2) * A2 * ones(1,N_theta) ; 
    L20(3, :) = L1(3) * A2 * ones(1,N_theta) ;
    L20(4, :) = L1(4) * A2 * ones(1,N_theta) ; 


figure(1)
plot( x , y ) ;
grid;

figure(2)
plot( reshape( real(s), N_s*length(theta), 1 ), reshape( imag(s), N_s*length(theta), 1 ), '.' , real(L1), imag(L1),'o', real(s(:,1)), imag(s(:,1)), 'd'  );
grid;



for m = 1 : length(gamma)
        
    L22 = L2.*A1.*sqrt( gamma(m) ) ;
        
    D_Euclid( m, : ) = min( abs( L11(2).*A1 + L22(4,:) - ( L11(1) + L22( 2,:) ) ), abs( L11(2) + L22(4,:) - ( L11(1) + L22( 3,:) ) ) ) ;
    
end

figure(1000)
mesh( theta./(pi/4)*45, gamma, D_Euclid )

for m = 1 : length(SNR)
    SNR( m ) 
    for n = 1 : length(theta)
        % A2 = sqrt(beta0) ;
        sigma = A1/sqrt( SNR( m ) ) ;
 
        for p = 1 : N_s
            
            for t = 1 : RUNS
             
                w = ( randn(1) + j*randn(1) ) * sigma ;
                                
                for q = 1 : N_s
                    d2 = ( abs( s(p,n) + w - s(q,n) ) )^2 - ( abs(w) )^2 ;
                    QQ( t, q ) = exp( - d2 / ( 2*sigma^2 ) ) ;                  
                end
                
            end            
            
            E_Q(p) = mean( log2( sum( QQ, 2 ) ) ) ;                      
            
        end
        
               
        Capacity(m,n) = log2( N_s ) - mean( E_Q ) ; 
    
    end
end


figure(3)
mesh( theta./(pi/4)*45, SNR, Capacity ) ;

figure(4)
plot( SNR_dB, Capacity(:,1), SNR_dB, max( Capacity, [], 2  ), '-d', 1:10, log2(1+10.^((1:10)./10)*(1+beta0^2)),'--'  )
xlabel('SNR (dB)')
ylabel('Spetral Efficiency')
grid;