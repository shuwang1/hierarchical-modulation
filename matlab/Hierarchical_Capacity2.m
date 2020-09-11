                                                                  
                                             
clear all;
close all ;

P = 800 ;
alpha = 1 ;
gamma = 0.001 : 0.01 : 0.99999 ;      %% The power splitting in terms of minimum Euclid distance
gamma_dB = -10:-1:-30 ;
% beta    = 10.^(gamma_dB./10)*alpha ;
beta = gamma*alpha ;
beta0   = 0.5*alpha ;


g1 = 0.12 ;
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

theta   = 0:0.5:0.9 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;

N_1 = 16 ;                  %% QPSK/16QAM
N_2 = 4 ;                   %% QPSK
N_s = N_1 * N_2 ;           %% The total constellation size => 16

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;

P_1 = mean( abs( s1 ).^2 ) 
P_2 = mean( abs( qammod( [0:1:N_2-1]' , N_2, 0 ) ).^2 ) 

RUNS = 500 ;

clock
tic

for m = 1 : length( theta )
    
    s2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    
    for n = 1 : length( beta )
        
        P1 = P * beta( n ) ;
        P2 = P * ( 1 - beta( n ) ) ;
        
        for k = 1 : N_2
            s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( P1/P_1 ) + s2( k ).*sqrt( P2/P_2 ) ;
        end
        
        for p = 1 : N_s
            w = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
            E_Q1(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g1) + w - s0*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
            E_Q2(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g2) + w - s0*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );            
        end
        
        R10_16QAM_QPSK( m, n ) = log2( N_s ) - mean( E_Q1 ) ;
        R20_16QAM_QPSK( m, n ) = log2( N_s ) - mean( E_Q2 ) ; 

        for p = 1 : N_2
            w2 = ( randn(N_2,RUNS) + j*randn(N_2,RUNS) ) * sigma ;
            E_Q12(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2*g1/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2*g1/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
            E_Q22(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2*g2/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2*g2/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
        end
        
        R12_QPSK1( m, n ) = log2( N_2 ) - mean( E_Q12 ) ; 
        R11_16QAM( m, n ) = R10_16QAM_QPSK( m, n ) - R12_QPSK1( m, n ) ;

        R22_QPSK1( m, n ) = log2( N_2 ) - mean( E_Q22 ) ;
        R21_16QAM( m, n ) = R20_16QAM_QPSK( m, n ) - R22_QPSK1( m, n ) ;      
           
    end
end

save(  strcat( 'layered_modualtion_16QAM_QPSK', num2str(now), '.mat' ) ) ;

figure( 20 ) ;
plot( R1(1,:) , R2(1,:), R1(2,:) , R2(2,:), '-.', R1(3,:) , R2(3,:), '--', R11_16QAM(1,:) , R20_16QAM_QPSK(1,:), '--', R11_16QAM(1,:) , R22_QPSK1(1,:), '-.' ) ;
grid ;
xlabel('R1, the rate for bad user, bps/symbol') ;
ylabel('R2, the rate for good user, bps/symbol') ;

figure( 30 ) ;
plot( beta , R10_16QAM_QPSK(1,:), beta, R12_QPSK1(1,:),'--', beta, R11_16QAM(1,:), '-.', beta, R20_16QAM_QPSK(1,:) , '-+' , beta, R22_QPSK1(1,:), '-s', beta, R21_16QAM(1,:), '-d' ) ;
grid ;
xlabel('Power-Splitting \gamma') ;
ylabel('Achievable Rate R') ;

pause %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear s0 s2 w E_Q1 E_Q2

N_1 = 4 ;                  %% QPSK/16QAM
N_2 = 4 ;                   %% QPSK
N_s = N_1 * N_2 ;           %% The total constellation size => 16

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;

P_1 = mean( abs( s1 ).^2 ) 
P_2 = mean( abs( qammod( [0:1:N_2-1]' , N_2, 0 ) ).^2 ) 

for m = 1 : length( theta )
    
    s2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    
    for n = 1 : length( beta )
        
        P1 = P * beta( n ) ;
        P2 = P * ( 1 - beta( n ) ) ;
        
        for k = 1 : N_2
            s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( P1/P_1 ) + s2( k ).*sqrt( P2/P_2 ) ;
        end
        
        for p = 1 : N_s
            w = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
            E_Q1(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g1) + w - s0*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
            E_Q2(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g2) + w - s0*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );            
        end
        
        R10_QPSK_QPSK( m, n ) = log2( N_s ) - mean( E_Q1 ) ;
        R20_QPSK_QPSK( m, n ) = log2( N_s ) - mean( E_Q2 ) ; 

        for p = 1 : N_2
            w2 = ( randn(N_2,RUNS) + j*randn(N_2,RUNS) ) * sigma ;
            E_Q12(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2*g1/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2*g1/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
            E_Q22(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2*g2/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2*g2/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
        end
        
        R12_QPSK2( m, n ) = log2( N_2 ) - mean( E_Q12 ) ; 
        R11_QPSK( m, n ) = R10_QPSK_QPSK( m, n ) - R12_QPSK2( m, n ) ;

        R22_QPSK2( m, n ) = log2( N_2 ) - mean( E_Q22 ) ;
        R21_QPSK( m, n ) = R20_QPSK_QPSK( m, n ) - R22_QPSK2( m, n ) ;      
           
    end
end

save(  strcat( 'layered_modualtion_16QAM_QPSK_QPSK', num2str(now), '.mat' ) ) ;
clock
toc

figure( 50 ) ;
plot( R1(1,:) , R2(1,:), R1(2,:) , R2(2,:), '-.', R1(3,:) , R2(3,:), '--', R11_16QAM(1,:) , R20_16QAM_QPSK(1,:), '--' , R11_QPSK(1,:) , R20_QPSK_QPSK(1,:), '--' , R11_16QAM(1,:) , R22_QPSK1(1,:), '--' , R11_QPSK(1,:) , R22_QPSK2(1,:), '--') ;
grid ;
xlabel('R1, the rate for bad user, bps/symbol') ;
ylabel('R2, the rate for good user, bps/symbol') ;

figure( 70 ) ;
plot( beta , R10_QPSK_QPSK(1,:), beta, R12_QPSK2(1,:),'--', beta, R11_QPSK(1,:), '-.', beta, R20_QPSK_QPSK(1,:) , '-+' , beta, R22_QPSK2(1,:), '-s', beta, R21_QPSK(1,:), '-d' ) ;
grid ;
xlabel('Power-Splitting \gamma') ;
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
        
        P1 = P * beta( n ) ;
        P2 = P * ( 1 - beta( n ) ) ;
       
        for k = 1 : N_2
            s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( P1/P_1 ) + s2( k ).*sqrt( P2/P_2 ) ;
        end
        
        for r = 1 : length(SNR)
            
            sigma = sqrt( P/SNR( r ) ) ;
            
            for p = 1 : N_s
                w = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
                E_Q(p) = mean( log2( sum( exp( - ( abs( s0(p) + w - s0*ones(1,RUNS) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
            end            
            Capacity( r, m, n) = log2( N_s ) - mean( E_Q ) ;             

            for p = 1 : N_2
                w2 = ( randn(N_2,RUNS) + j*randn(N_2,RUNS) ) * sigma ;
                E_Q2(p) = mean( log2( sum( exp( - ( abs( s2(p).*sqrt( P2/P_2 ) + w2 - s2*ones(1,RUNS).*sqrt( P2/P_2 ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
            end
            Capacity2( r, m, n) = log2( N_s ) - mean( E_Q2 ) ; 
            
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


figure(10)
mesh( theta./pi*180, SNR_dB, Capacity( :, :, 1) ) ;
grid;


figure(20)
%plot( beta, Capacity( 1, 1, :),beta, Capacity( 6, 1, :), beta, Capacity( 11, 1, :), beta, Capacity( 16, 1, :) )
plot( beta', Cap( :, 1, 1), beta', Cap( :, 6, 1 ), beta', Cap( :, 11, 1 ), beta', Cap( :, 16, 1) )


figure(100)
plot( SNR_dB, Capacity( :, 1, 10), SNR_dB, C_max( :, 10 ) )
grid ;