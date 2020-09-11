                                                                  
                                             
clear all;
close all ;

P = 900 ;
alpha = 1 ;
gamma = 0.0005 : 0.01 : 0.99999 ;      %% The power splitting in terms of minimum Euclid distance
gamma_dB = -10:-1:-30 ;
% beta    = 10.^(gamma_dB./10)*alpha ;
beta = gamma*alpha ;
beta0   = 0.5*alpha ;


g1 = 0.20 ;
g2 = 1.00 ;
sigma = 1;

for m = 1 : length( beta )
    
    P1 = P * beta( m ) ;
    P2 = P * ( 1 - beta( m ) ) ;

    R1(1, m) = log2( 1 + P1*g1/(sigma+P2*g1) ) ;
    R2(1, m) = log2( 1 + P2*g2/sigma ) ;        

    R1(2, m) = log2( 1 + P1*g1/(sigma+P2*g1) ) ;
    R2(2, m) = log2( 1 + P*g2/sigma ) ;        
        
    R1(3, m) = beta(m)*log2( 1 + P1*g1/beta(m)/sigma ) ;
    R2(3, m) = (1-beta(m))*log2( 1 + P2*g2/(1-beta(m))/sigma ) ;        
    
end


SNR_dB  = 0:2:30 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0:0.5:0.5 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;


N_1 = 16 ;                   %% QPSK/16QAM
N_2 = 4 ;                   %% QPSK
N_s = N_1 * N_2 ;           %% The total constellation size => 16

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;

P_1 = mean( abs( s1 ).^2 ) 
P_2 = mean( abs( qammod( [0:1:N_2-1]' , N_2, 0 ) ).^2 ) 

s = zeros( N_1*N_2, L_theta ) ;

RUNS = 400 ;
clock
tic

s2 = qammod( [0:1:N_2-1]' , N_2, 0 ) ;    

for n = 1 : length( beta )
    
    P1 = P * beta( n ) ;
    P2 = P * ( 1 - beta( n ) ) ;
    
    for p = 1 : N_2
        w2 = ( randn(N_2,RUNS) + j*randn(N_2,RUNS) ) * sigma ;
        E_QPSK1(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2*g1/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2*g1/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
        E_QPSK2(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2*g2/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2*g2/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
    end
    
    R_QPSK1( n ) = log2( N_2 ) - mean( E_QPSK1 ) ; 
    R_QPSK2( n ) = log2( N_2 ) - mean( E_QPSK2 ) ;    
    
    for m = 1 : length( theta )
        
       
        for k = 1 : N_2
            s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( P1/P_1 ) + s2( k )*sqrt( P2/P_2 )*exp(j*theta(m)) ;
        end
        
        for p = 1 : N_s
            w = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
            E_Q1(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g1) + w - s0*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
            E_Q2(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g2) + w - s0*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );            
        end
        
        R10( m, n ) = log2( N_s ) - mean( E_Q1 ) ;
        R20( m, n ) = log2( N_s ) - mean( E_Q2 ) ; 

        R11( m, n ) = R10( m, n ) - R_QPSK1( n ) ;        
        R21( m, n ) = R20( m, n ) - R_QPSK2( n ) ;      
           
    end
end


for m = 1 : length( beta )
    
    R10opt( m ) = max( R10(:,m) ) ;
    R11opt( m ) = R10opt( m ) - R_QPSK1( m ) ;

    R20opt( m ) = max( R20(:,m) ) ;
    R21opt( m ) = R20opt( m ) - R_QPSK2( m ) ;
    
end

save(  strcat( 'layered_modualtion_QPSK_QPSK', num2str(now), '.mat' ) ) ;
clock
toc

figure( 10 ) ;
plot( R2(1,:), R1(1,:) , R2(2,:), R1(2,:) , R2(3,:), R1(3,:) , R_QPSK2, R11(1,:) ,'--', R_QPSK2, R20(1,:), '--' , R_QPSK2 , R_QPSK1, '-.' ) ;
grid ;
xlabel('Enhancement-Layer Rate R2') ;
ylabel('Base-Layer Rate R1; Total Rate R') ;


figure( 20 ) ;
plot( R1(1,:), R2(1,:), R1(2,:), R2(2,:), R1(3,:), R2(3,:), R11(1,:), R_QPSK2, '--', R11, R20(1,:), '--' , R11opt(1,:) , R20opt(1,:), '-.' , R11opt(1,:) , R_QPSK2, '-.'  ) ;
grid ;
xlabel('Base-Layer Rate R1') ;
ylabel('Enhancement-Layer Rate R2; Total Rate R') ;


figure( 30 ) ;
plot( beta , R10(1,:), beta, R_QPSK1,'--', beta, R11(1,:), '-.', beta, R20(1,:) , '-+' , beta, R_QPSK2, '-s', beta, R21(1,:), '-d' ) ;
grid ;
xlabel('Power-Splitting \beta') ;
ylabel('Achievable Rate R') ;


return


for m = 1:1:L_theta
    s2( :, m ) = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    for n = 1:1:N_2
        s( (1+(n-1)*N_1):(n*N_1) , m ) = s1 + s2( n, m ).*sqrt(P_1/P_2*beta0) ;
    end
end


for m = 1 : length( theta )
    
    s2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    
    for n = 1 : length( beta )
        
        for k = 1 : N_2
            s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( (1-beta(n))*P/P_1 ) + s2( k ).*sqrt( beta(n)*P/P_2 ) ;
        end
        
        for r = 1 : length(SNR)
            
            sigma = sqrt( P_1/SNR( r ) ) ;
            
            for p = 1 : N_s
                
                w = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
                
                E_Q(p) = mean( log2( sum( exp( - ( abs( s0(p) + w - s0*ones(1,RUNS) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
                
            end
            
            Capacity( r, m, n) = log2( N_s ) - mean( E_Q ) ;             
            
        end
    
    end
end


for m = 1 : length( theta )
    
    for n = 1 : length( beta )
        
        for r = 1 : length(SNR)
            
            Cap( n, r, m ) = Capacity( r, m, n) ;             
            
        end
    
    end
end


for m = 1 : length(beta)
    for n = 1 : length(SNR)
        C_max( n, m ) = max( Capacity( n, :, m)  ) ;
    end
end


save(  strcat( 'layered_modualtion_QPSK_QPSK', num2str(now), '.mat' ) ) ;
clock
toc


figure(10)
mesh( theta./pi*180, SNR_dB, Capacity( :, :, 1) ) ;
grid;


figure(20)
%plot( beta, Capacity( 1, 1, :),beta, Capacity( 6, 1, :), beta, Capacity( 11, 1, :), beta, Capacity( 16, 1, :) )
plot( beta', Cap( :, 1, 1), beta', Cap( :, 6, 1 ), beta', Cap( :, 11, 1 ), beta', Cap( :, 16, 1) )


figure(100)
plot( SNR_dB, Capacity( :, 1, 10), SNR_dB, C_max( :, 10 ) )
grid ;