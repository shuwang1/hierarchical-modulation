                                             
clear all;
close all ;

P = 200 ;
alpha = 1 ;
gamma = 0.0000 : 0.1 : 1 ;      %% The power splitting in terms of minimum Euclid distance
gamma_dB = -10:-1:-30 ;
% beta    = 10.^(gamma_dB./10)*alpha ;
beta = gamma*alpha ;
beta0   = 0.5*alpha ;


g1 = 0.25 ;
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

% theta   = 0:0.1:1 ;
theta   = 0 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;


N_1 = 4 ;                   %% QPSK/16QAM
N_2 = 4 ;                   %% QPSK
N_s = N_1 * N_2 ;           %% The total constellation size => 16

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;

P_1 = mean( abs( s1 ).^2 ) 
P_2 = mean( abs( qammod( [0:1:N_2-1]' , N_2, 0 ) ).^2 ) 

s = zeros( N_1*N_2, L_theta ) ;

RUNS = 500 ;
clock
tic

s_16QAM = qammod( [0:1:16-1]' , 16, 0 ) ;
s_QPSK  = qammod( [0:1:4-1]' , 4, 0 ) ;
s2 = qammod( [0:1:N_2-1]' , N_2, 0 ) ;    

sigma = sigma/sqrt(2) ;
for n = 1 : length( beta )
    
    beta(n)
    
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
            s0_16QAM( (1+(k-1)*16):(k*16), 1 ) = s_16QAM*sqrt( P1/P_1 ) + s2( k )*sqrt( P2/P_2 )*exp(j*theta(m)) ;
            s0_QPSK( (1+(k-1)*4):(k*4), 1 )    = s_QPSK*sqrt( P1/P_1 ) + s2( k )*sqrt( P2/P_2 )*exp(j*theta(m)) ;
        end
        
        N_s = N_2*16 ;
        for p = 1 : N_s
            w_16QAM = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
            E1_16QAM(p) = mean( log2( sum( exp( - ( abs( s0_16QAM(p)*sqrt(g1) + w_16QAM - s0_16QAM*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w_16QAM).^2 ) ./ (2*sigma^2) ) ) ) );
            E2_16QAM(p) = mean( log2( sum( exp( - ( abs( s0_16QAM(p)*sqrt(g2) + w_16QAM - s0_16QAM*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w_16QAM).^2 ) ./ (2*sigma^2) ) ) ) );            
        end
        
        R10_16QAM( m, n ) = log2( N_s ) - mean( E1_16QAM ) ;
        R20_16QAM( m, n ) = log2( N_s ) - mean( E2_16QAM ) ; 

        R11_16QAM( m, n ) = R10_16QAM( m, n ) - R_QPSK1( n ) ;        
        R21_16QAM( m, n ) = R20_16QAM( m, n ) - R_QPSK2( n ) ;      
  
        N_s = N_2*4 ;
        for p = 1 : N_s
            w_QPSK = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
            E1_QPSK(p) = mean( log2( sum( exp( - ( abs( s0_QPSK(p)*sqrt(g1) + w_QPSK - s0_QPSK*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w_QPSK).^2 ) ./ (2*sigma^2) ) ) ) );
            E2_QPSK(p) = mean( log2( sum( exp( - ( abs( s0_QPSK(p)*sqrt(g2) + w_QPSK - s0_QPSK*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w_QPSK).^2 ) ./ (2*sigma^2) ) ) ) );            
        end
        
        R10_QPSK( m, n ) = log2( N_s ) - mean( E1_QPSK ) ;
        R20_QPSK( m, n ) = log2( N_s ) - mean( E2_QPSK ) ; 

        R11_QPSK( m, n ) = R10_QPSK( m, n ) - R_QPSK1( n ) ;        
        R21_QPSK( m, n ) = R20_QPSK( m, n ) - R_QPSK2( n ) ;      
        
    end
end


for m = 1 : length( beta )
    
    R10_16QAM_opt( m ) = max( R10_16QAM(:,m) ) ;
    R11_16QAM_opt( m ) = R10_16QAM_opt( m ) - R_QPSK1( m ) ;

    R20_16QAM_opt( m ) = max( R20_16QAM(:,m) ) ;
    R21_16QAM_opt( m ) = R20_16QAM_opt( m ) - R_QPSK2( m ) ;

    R10_QPSK_opt( m ) = max( R10_QPSK(:,m) ) ;
    R11_QPSK_opt( m ) = R10_QPSK_opt( m ) - R_QPSK1( m ) ;

    R20_QPSK_opt( m ) = max( R20_QPSK(:,m) ) ;
    R21_QPSK_opt( m ) = R20_QPSK_opt( m ) - R_QPSK2( m ) ;
    
end

save(  strcat( 'layered_modualtion_QPSK_QPSK', num2str(now), '.mat' ) ) ;
clock
toc

figure( 10 ) ;
plot( R1(1,:), R2(1,:), R1(2,:), R2(2,:), R1(3,:), R2(3,:), R11_16QAM(1,:), R_QPSK2, '--', R11_16QAM(1,:), R20_16QAM(1,:), '--' , R11_QPSK(1,:), R_QPSK2, '-.', R11_QPSK(1,:), R20_QPSK(1,:), '-.' ) ;
grid ;
xlabel('Base-Layer Rate R1') ;
ylabel('Enhancement-Layer Rate R2; Total Rate R') ;


figure( 20 ) ;
plot( R1(1,:), R2(1,:), R1(2,:), R2(2,:), R1(3,:), R2(3,:), R11_16QAM(1,:), R_QPSK2, '--', R11_16QAM(1,:), R20_16QAM(1,:), '--' , R11_16QAM_opt(1,:) , R20_16QAM_opt(1,:), '-.' , R11_16QAM_opt(1,:) , R_QPSK2, '-.'  ) ;
grid ;
xlabel('Base-Layer Rate R1') ;
ylabel('Enhancement-Layer Rate R2; Total Rate R') ;


figure( 30 ) ;
plot( R1(1,:), R2(1,:), R1(2,:), R2(2,:), R1(3,:), R2(3,:), R11_QPSK(1,:), R_QPSK2, '--', R11_QPSK(1,:), R20_QPSK(1,:), '--' , R11_QPSK_opt(1,:) , R20_QPSK_opt(1,:), '-.' , R11_QPSK_opt(1,:) , R_QPSK2, '-.'  ) ;
grid ;
xlabel('Base-Layer Rate R1') ;
ylabel('Enhancement-Layer Rate R2; Total Rate R') ;


figure( 40 ) ;
plot( beta , R20_16QAM(1,:),  '-d', beta, R20_16QAM_opt , '-+', beta, R21_16QAM(1,:), '-.', beta, R21_16QAM_opt(1,:), beta, R_QPSK2  ) ;
grid ;
xlabel('Power-Splitting \beta') ;
ylabel('Achievable Rate R') ;

figure( 50 ) ;
plot( beta , R20_QPSK(1,:), beta, R_QPSK1,'--', beta, R21_QPSK(1,:), '-.', beta, R20_16QAM_opt , '-+' , beta, R_QPSK2, '-s', beta, R21_16QAM_opt(1,:), '-d' ) ;
grid ;
xlabel('Power-Splitting \eta') ;
ylabel('Achievable Rate R') ;
legend('Total Rate')

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