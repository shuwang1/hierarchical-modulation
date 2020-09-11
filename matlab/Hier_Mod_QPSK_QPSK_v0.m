clear all;
close all ;

alpha = 10 ;
gamma = 0.0 : 0.01 : 1 ;      %% The power splitting in terms of minimum Euclid distance
gamma_dB = -10:-1:-30 ;
% beta    = 10.^(gamma_dB./10)*alpha ;
beta = gamma*alpha ;
beta0   = 0.5*alpha ;

SNR_dB  = 0:2:30 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;


N_1 = 4 ;                  %% QPSK
L1  = zeros( N_1, 1 ) ;
N_2 = 4 ;                   %% QPSK
L2  = zeros( N_2, 1 ) ;
N_s = N_1 * N_2 ;           %% The total constellation size => 16

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;

P_1 = mean( abs( s1 ).^2 ) 
P_2 = mean( abs( qammod( [0:1:N_2-1]' , N_2, 0 ) ).^2 ) 

s = zeros( N_1*N_2, L_theta ) ;
for m = 1:1:L_theta
    s2( :, m ) = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    for n = 1:1:N_2
        s( (1+(n-1)*N_1):(n*N_1) , m ) = s1 + s2( n, m ).*sqrt(P_1/P_2*beta0) ;
    end
end

clock
tic
RUNS = 500 ;

PP = zeros( N_s, RUNS ) ;

for m = 1 : length( theta )
    
    s2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    
    for n = 1 : length( beta )
        
        for k = 1 : N_2
            s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( (alpha-beta(n))/P_1 ) + s2( k ).*sqrt( beta(n)/P_2 ) ;
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


