clear all
close  all
clock
tic
P = 100
RUNS = 2000 ;
theta   = 0:0.01:1 ;
theta   = theta .* (pi/4) ;

N_1 = 16 ;
N_2 = 4 ;

g1 = 1.00 ;
g2 = 1.00 ;

s_16QAM = qammod( [0:1:N_1-1]' , N_1, 0 ) ;
P_16QAM = mean( abs(s_16QAM).^2 ) ;

s_QPSK = qammod( [0:1:N_2-1]' , N_2, 0 ) ;
P_QPSK = mean( abs(s_QPSK).^2 ) ;

%SNR0_dB = -log10(5)*10:0.01:30 ;
SNR0_dB = -20:1:20 ;
SNR0 = 10.^( SNR0_dB./10 );

w0 = randn(1,RUNS) + j*randn(1,RUNS) ;

for m = 1 : length( SNR0 )
    
    sigma = sqrt( P / SNR0( m ) / 2 ) ;

    for p = 1 : N_1
        w1 = ones( N_1, 1 ) * w0 * sigma ;
        E_16QAM(p) = mean( log2( sum( exp( - ( abs( s_16QAM(p).*sqrt( P/P_16QAM ) + w1 - s_16QAM*ones(1,RUNS).*sqrt( P/P_16QAM ) ).^2 - abs(w1).^2 ) ./ (2*sigma^2) ) ) ) );
    end
    R_16QAM0( m ) = log2( N_1 ) - mean( E_16QAM ) ; 
    
    for p = 1 : N_2
        w2 = ones( N_2, 1 ) * w0 ;
        A0 = sqrt( SNR0(m)*2/P_QPSK ) ;
        E_QPSK(p) = mean( log2( sum( exp( - ( abs( s_QPSK(p).*A0 + w2 - s_QPSK*ones(1,RUNS).*A0 ).^2 - abs(w2).^2 ) ./ 2 ) ) ) );
    end
    R_QPSK0( m ) = log2( N_2 ) - mean( E_QPSK ) ; 
end

SNR_dB  = -10:1:20 ;
SNR     = 10.^(SNR_dB./10) ;

beta = 0.8 ;

N_1 = 4 ;                  %% QPSK/16QAM
N_2 = 4 ;                  %% QPSK
s1 = qammod( [0:1:N_1-1]' , N_1, 0 ) ;
s2 = qammod( [0:1:N_2-1]' , N_2, 0 ) ;   
P_1 = P_QPSK ;
P_2 = P_QPSK ;
N_s = N_1 * N_2 ;
for k = 1 : length( SNR )
    sigma = sqrt( P / SNR( k ) / 2 ) ;

    for n = 1 : length( beta )

        P1 = P * beta( n ) ;
        P2 = P * ( 1 - beta( n ) ) ;
        
        SNR1 = SNR(k)*beta(n) ;
        SNR2 = SNR(k)*( 1-beta(n) );

        R_QPSK1( n, k ) = R_QPSK0( round( (log10(SNR1)*10+20)/1 ) + 1) ;
        R_QPSK2( n, k ) = R_QPSK0( round( (log10(SNR2)*10+20)/1 ) + 1) ;

        for m = 1 : length( theta )

            for p = 1 : N_2
                s0( (1+(p-1)*N_1):(p*N_1), 1 ) = s1*sqrt( P1/P_1 ) + s2( p )*sqrt( P2/P_2 )*exp(j*theta(m)) ;
            end

            for p = 1 : N_s
                w = ones( N_s, 1 ) * w0 * sigma ;
                E_Q1(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g1) + w - s0*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
                E_Q2(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g2) + w - s0*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );            
            end

            R10( m, n, k ) = log2( N_s ) - mean( E_Q1 ) ;
            R20( m, n, k ) = log2( N_s ) - mean( E_Q2 ) ; 

            R11( m, n, k ) = R10( m, n, k ) - R_QPSK1( n, k ) ;        
            R21( m, n, k ) = R20( m, n, k ) - R_QPSK2( n, k ) ;      

        end
        
        R10opt( k ) = max( R10( :, n, k) ) ;
        R11opt( k ) = R10opt( k ) - R_QPSK1( n, k ) ;

        R20opt( n, k ) = max( R20(:, n, k ) ) ;
        R21opt( n, k ) = R20opt( n, k ) - R_QPSK2( n, k ) ;
        
    end
    
end

save(  strcat( 'layered_mod_Effective_SINR_', num2str(now), '.mat' ) ) ;
clock
toc

SNR0_dB = -10:1:20 ;
SNR0 = 10.^( SNR0_dB./10 );

for m = 1 : length( SNR0 )
    
    R_16QAM( m ) = R_16QAM0( 10 + m ) ;    
    R_QPSK( m ) = R_QPSK0( 10 + m ) ; 
    
end


delta1 = round( log10(5)*10 / 1 + 1 ) 
delta2 = round( log10(5/4)*10 / 1 + 1 ) 
delta12 = round( log10(5)*10 / 1 - log10(5/4)*10/1 )


figure(100)
plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R_16QAM(delta1:end), SNR0_dB(delta1:end), R_QPSK(delta1:end), SNR0_dB(delta1:end), R_16QAM(delta1:end) - R_QPSK(1:(end-delta1 + 1)),'--', SNR0_dB(delta1:end), R_QPSK(delta12:(end-delta2 )), '-.', SNR0_dB(delta1:end) ,R_QPSK(1:(end-delta1 + 1)),'-.' ) ;
xlim([0 18]);
ylim([0 4.5]) ;
grid ;
legend('Shannon Bound','C_{16QAM}(\gamma)','C_{QPSK}(\gamma)','C_{16QAM}-C_{QPSK}(0.2*\gamma)','C_{QPSK}(0.8*\gamma)','C_{QPSK}(0.2*\gamma)') ;


figure( 1000 ) ;
subplot(1,2,1);
plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R_16QAM(delta1:end), SNR0_dB(delta1:end), R_QPSK(delta1:end), SNR0_dB(delta1:end), R_16QAM(delta1:end) - R_QPSK(1:(end-delta1 + 1)),'--', SNR0_dB(delta1:end), R_QPSK(delta12:(end-delta2 )), '-.', SNR0_dB(delta1:end) ,R_QPSK(1:(end-delta1 + 1)),'-.' ) ;
xlim([0 18]);
ylim([0 4.5]) ;
grid ;
legend('Shannon Bound','C_{16QAM}(\gamma)','C_{QPSK}(\gamma)','C_{16QAM}-C_{QPSK}(0.2*\gamma)','C_{QPSK}(0.8*\gamma)','C_{QPSK}(0.2*\gamma)') ;
subplot(1,2,2)
plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R10opt( delta1:end ), SNR0_dB(delta1:end), R_QPSK(delta1:end), SNR0_dB(delta1:end), R10opt( delta1:end ) - R_QPSK(1:(end-delta1 + 1)),'--', SNR0_dB(delta1:end), R_QPSK(delta12:(end-delta2 )), '-.', SNR0_dB(delta1:end) ,R_QPSK(1:(end-delta1 + 1)),'-.' ) ;
xlim([0 18]);
ylim([0 4.5]) ;
grid ;
legend('Shannon Bound','C_{16QAM}(\gamma)','C_{QPSK}(\gamma)','C_{16QAM}-C_{QPSK}(0.2*\gamma)','C_{QPSK}(0.8*\gamma)','C_{QPSK}(0.2*\gamma)') ;


figure(1005)
plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R_16QAM(delta1:end), SNR0_dB(delta1:end),  R10opt( delta1:end )) ;
xlim([0 18]);
ylim([0 4.5]) ;


return

figure(1010 ) ;
plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R_QPSK(delta1:end), SNR0_dB(delta1:end), R_QPSK( (delta12+1):(end-delta2 + 1 )), '-.', SNR0_dB(delta1:end) ,R_QPSK(1:(end-delta1 + 1)),'--' ) ;
xlim([0 30]);
ylim([0 6]) ;
grid ;
legend('Shannon Bound', 'C_{QPSK}(\gamma)','C_{QPSK}(0.8*\gamma)','C_{QPSK}(0.2*\gamma)') ;
