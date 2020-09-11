clear all;
close all ;

clock

beta  = 1:-0.1:0.00001 ; % 0.8: -0.8 : 0.000001 ; 

theta   = 0:0.05:1 ;
theta   = theta .* (pi/4) ;

P = 100 ;

N_1 = 16 ;                  %% QPSK/16QAM
N_2 = 4 ;                   %% QPSK

s_16QAM = qammod( [0:1:N_1-1]' , N_1, 0 ) ;
P_16QAM = mean( abs(s_16QAM).^2 ) ;

s_QPSK = qammod( [0:1:N_2-1]' , N_2, 0 ) ;
P_QPSK = mean( abs(s_QPSK).^2 ) ;


SNR0_dB = -5 : 0.1 : 30 ;
SNR0    = 10.^(SNR0_dB./10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert from EbNo to SNR.
% Note: Because No = 2*noiseVariance^2, we must add 3 dB
% to get SNR. For details, see Proakis book listed in
% "Selected Bibliography for Performance Evaluation."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EbNo_dB0_QPSK  = SNR0_dB  - 10*log10( 2 )  ;
EbNo_dB0_16QAM = SNR0_dB  - 10*log10( 4 )  ;

ber_QPSK            = berawgn( EbNo_dB0_QPSK,  'qam' , 4  ) ;
ber_16QAM           = berawgn( EbNo_dB0_16QAM, 'qam' , 16 ) ;

ber_QPSK_calc       = qfunc( sqrt(SNR0) ) ;
ber_16QAM_calc      = 0.75*qfunc(sqrt(SNR0/5)) + 0.5*qfunc(3*sqrt(SNR0/5)) - 0.25*qfunc(sqrt(5*SNR0)) ;

runs = ceil( ber_16QAM.^(-1) ).*1000 ;

for m = 1 : length(SNR0_dB) 
    
    RUNS = min( runs( m ), 5.0e6 ) ;
    snr_dB = SNR0_dB( m ) ;    
    
%     b_QPSK  = randint( RUNS, 1, 4 ) ;            
%     s_QPSK = qammod( b_QPSK, 4, 0 ) ;
%     r_QPSK = awgn( s_QPSK, snr_dB, 'measured' ) ;
%     d_QPSK = qamdemod( r_QPSK , 4 , 0 ) ;   
%     [ num_QPSK( m ), ber_QPSK_simu( m ) ] = biterr( b_QPSK, d_QPSK ) ;  
    
    b_16QAM = randint( RUNS, 1, 16 ) ;            
    s_16QAM = qammod( b_16QAM, 16, 0 ) ;
    r_16QAM = awgn( s_16QAM, snr_dB, 'measured' ) ;
    d_16QAM = qamdemod( r_16QAM , 16 , 0 ) ;   
    [ num_16QAM( m ), ber_16QAM_simu( m ) ] = biterr( b_16QAM, d_16QAM ) ;  
    
end
save(  strcat( 'BER_QPSK_16QAM_', num2str(now), '.mat' ) ) ;

% figure( 1000 )
% semilogy( EbNo_dB0_QPSK, ber_QPSK, EbNo_dB0_QPSK, ber_QPSK_calc, '-.' , EbNo_dB0_QPSK, ber_QPSK_simu, '--' ) ;
% grid ;
% ylim([10^(-7) 1])
% xlabel( 'E_b/N_o' ) ;
% ylabel( 'BER' ) ;
% legend( 'QPSK matlab', 'QPSK calc', 'QPSK simu' ) ;

figure( 1003 )
semilogy( EbNo_dB0_16QAM, ber_16QAM, EbNo_dB0_16QAM, ber_16QAM_calc , '-.', EbNo_dB0_16QAM,  ber_16QAM_simu, '--' ) ;
grid ;
ylim([10^(-7) 1])
xlabel( 'E_b/N_o' ) ;
ylabel( 'BER' ) ;
legend( '16QAM matlab', '16QAM calc', '16QAM simu' ) ;


% figure( 1010 )
% semilogy( SNR0_dB, ber_QPSK, SNR0_dB, ber_QPSK_calc, '--' , SNR0_dB , ber_16QAM, SNR0_dB, ber_16QAM_calc , '-.' ) ;
% grid ;
% ylim([10^(-10) 1])
% xlabel( 'E_s/N_o' ) ;
% ylabel( 'BER' ) ;

clock

SNR_dB  = 0:1:30 ;
SNR     = 10.^(SNR_dB./10);

EbNo_dB_QPSK  = SNR_dB  - 10*log10( 2 ) ;
EbNo_dB_16QAM = SNR_dB  - 10*log10( 4 ) ;

for m = 1 : length(SNR_dB)
    
    snr_dB = SNR_dB( m ) ;
    RUNS = min( runs( m ), 5.0e6 ) ;
    
    for n = 1 : length(beta)
        
        for k = 1 : length(theta)
            
            b_16QAM = randint( RUNS, 1, 16 ) ;            
            b_enhan = randint( RUNS, 1, 4  ) ;
            
            s_16QAM = qammod( b_16QAM, 16, 0 ) ;
            s_enhan = qammod( b_enhan,  4, 0 )*exp( j*theta(k) ) ;
            
            r_16QAM = awgn( s_16QAM, snr_dB, 'measured' ) + s_enhan*sqrt( (1-beta(n))/beta(n)*P_16QAM/P_QPSK )/2;
            
            d_16QAM = qamdemod( r_16QAM , 16 , 0 ) ;   

            [ num_16QAM( k, n, m ), BER_16QAM( m, n, k ) ] = biterr( b_16QAM, d_16QAM ) ;  
                        
            BER_theta( 1, k ) = BER_16QAM( m, n, k ) ;
            
            [ diff, index ] = min( abs( ber_16QAM_simu - BER_16QAM( m, n, k ) ) ) ;
            ESNR_16QAM( m, n, k ) = SNR0( index ) ;
            ESNR_dB_16QAM( m, n, k ) = SNR0_dB( index ) ;           
            
        end
        
        [ BER_16QAM_min( m, n ), index ] = min( BER_theta( 1, : ) ) ;       
        theta_16QAM( m , n ) = theta( index ) ; 
        
        ESNR_16QAM_min( m, n ) = ESNR_16QAM( m, n, index ) ;
        ESNR_dB_16QAM_min( m, n ) = ESNR_dB_16QAM( m, n, index ) ;
    end
end

save(  strcat( 'ESNR_16QAM_', num2str(now), '.mat' ) ) ;

figure( 2000 )
semilogy( EbNo_dB0_QPSK, ber_QPSK, EbNo_dB0_QPSK, ber_QPSK_calc, '--' , ...
    EbNo_dB0_16QAM, ber_16QAM, EbNo_dB0_16QAM, ber_QPSK_calc, '--' , ...
    EbNo_dB_16QAM, BER_16QAM(:, 1, 1) , '-.', EbNo_dB_16QAM, BER_16QAM(:, 2, 1) , '-.', EbNo_dB_16QAM, BER_16QAM(:, 3, 1) , '-.' ) ;
grid ;
ylim([10^(-7) 1])
xlabel( 'E_b/N_o' ) ;
ylabel( 'BER' ) ;

figure( 2005 )
semilogy( EbNo_dB0_QPSK, ber_QPSK, EbNo_dB0_QPSK, ber_QPSK_calc, '--' , ...
    EbNo_dB0_16QAM, ber_16QAM, EbNo_dB0_16QAM, ber_16QAM_calc, '--' , EbNo_dB_16QAM - 3, BER_16QAM(:, 1, 1) , '-.'  ) ;
grid ;
ylim([10^(-7) 1])
xlabel( 'E_b/N_o' ) ;
ylabel( 'BER' ) ;
legend( 'QPSK matlab', 'QPSK calc', 'QPSK simu', '16QAM matlab', '16QAM calc', '16QAM simu' ) ;

% figure( 2010 )
% semilogy( SNR0_dB, ber_QPSK.*2, SNR0_dB, ber_QPSK_calc.*2, '--' , SNR0_dB, SER_QPSK( :, 1, 1 ), '-o', SNR0_dB, SER_QPSK( :, 2, 1 ), '--') ;
% grid ;
% ylim([10^(-7) 1])
% xlabel( 'E_s/N_o' ) ;
% ylabel( 'SER' ) ;

% figure( 2020 )
% plot( SNR0_dB, SNR0_dB, SNR_dB , ESNR_dB_QPSK( :, 1, 1 ), '--' ) ;
% hold on ;
% figure( 2020 )
% plot( SNR_dB , ESNR_dB_QPSK( :, 2, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 2), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 3, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 3), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 4, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 4), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 5, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 5), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 6, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 6), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 7, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 7), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 8, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 8), '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 9, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 9), '-.' ) ;
% xlim([0 20])
% grid ;
% xlabel( 'E_s/N_o' ) ;
% ylabel( 'Effective E_s/N_o' ) ;

figure( 2025 )
plot( SNR_dB, SNR_dB, '-', SNR_dB , ESNR_dB_16QAM( :, 1, 1 ), '--' ) ;
hold on ;
figure( 2025 )
plot( SNR_dB , ESNR_dB_16QAM( :, 2, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 2), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 3, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 3), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 4, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 4), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 5, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 5), '-.' ) ;
%     SNR_dB , ESNR_dB_16QAM( :, 6, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 6), '-.' , ...
%     SNR_dB , ESNR_dB_16QAM( :, 7, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 7), '-.' , ...
%     SNR_dB , ESNR_dB_16QAM( :, 8, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 8), '-.' , ...
%     SNR_dB , ESNR_dB_16QAM( :, 9, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 9), '-.' ) ;
xlim([0 20]) ;
grid ;
xlabel( 'E_s/N_o' ) ;
% ylabel( 'Effective E_s/N_o' ) ;

figure( 2027 )
plot( SNR0_dB, SNR0_dB, SNR_dB , ESNR_dB_16QAM( :, 1, 1 ), '--' , SNR_dB , ESNR_dB_16QAM( :, 2, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 2), '-.' ) ;
xlim([0 30])
ylim([0 30])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;

% figure( 2030 )
% plot(SNR_dB, 10.^( (ESNR_dB_QPSK( :, 1, 1 ) - SNR_dB')/10 )     , ...
%     SNR_dB , 10.^( (ESNR_dB_QPSK( :, 2, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^( (ESNR_dB_QPSK_min(:, 2)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_QPSK( :, 3, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 3)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^( (ESNR_dB_QPSK( :, 4, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^( (ESNR_dB_QPSK_min(:, 4)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_QPSK( :, 5, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 5)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_QPSK( :, 6, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 6)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_QPSK( :, 7, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 7)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_QPSK( :, 8, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 8)  - (SNR_dB )')/10 )  , '-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_QPSK( :, 9, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 9)  - (SNR_dB )')/10 )  , '-.' ) ;
% xlim([0 20])
% grid ;
% xlabel( 'E_s/N_o' ) ;
% ylabel( 'Effective E_s/N_o' ) ;

% figure( 2040 )
% plot( SNR_dB  , ESNR_dB_QPSK( :, 1, 1 ) - SNR_dB', ...
%     SNR_dB , ESNR_dB_QPSK( :, 2, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 2) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 3, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 3) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 4, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 4) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 5, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 5) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 6, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 6) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 7, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 7) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 8, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 8) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_QPSK( :, 9, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 9) - (SNR_dB )', '-.' ) ;
% xlim([0 20])
% grid ;
% xlabel( 'E_s/N_o' ) ;
% ylabel( 'Modulation Efficiency' ) ;

figure( 2050 )
plot(SNR_dB, 10.^( (ESNR_dB_16QAM( :, 1, 1 ) - SNR_dB')/10 )     , ...
    SNR_dB , 10.^( (ESNR_dB_16QAM( :, 2, 1 ) - (SNR_dB )')/10 )  , '--' ,  SNR_dB , 10.^( (ESNR_dB_16QAM_min(:, 2)  - (SNR_dB )')/10 )   , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 3, 1 ) - (SNR_dB )')/10 )  , '--' ,  SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 3)  - (SNR_dB )')/10 )   , '-.', ...
    SNR_dB , 10.^( (ESNR_dB_16QAM( :, 4, 1 ) - (SNR_dB )')/10 )  , '--' ,  SNR_dB , 10.^( (ESNR_dB_16QAM_min(:, 4)  - (SNR_dB )')/10 )   , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 5, 1 ) - (SNR_dB )')/10 )  , '--' ,  SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 5)  - (SNR_dB )')/10 )   , '-.' ) ;
%     SNR_dB , 10.^(( ESNR_dB_16QAM( :, 6, 1 ) - (SNR_dB )')/10 )  , '--' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 6)  - (SNR_dB )')/10 )   ,'-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM( :, 7, 1 ) - (SNR_dB )')/10 )  , '--' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 7)  - (SNR_dB )')/10 )   ,'-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM( :, 8, 1 ) - (SNR_dB )')/10 )  , '--' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 8)  - (SNR_dB )')/10 )   ,'-.' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM( :, 9, 1 ) - (SNR_dB )')/10 )  , '--' , ...
%     SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 9)  - (SNR_dB )')/10 )   ,'-.' ) ;
xlim([0 30])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;

figure( 2055 )
plot(SNR_dB, 10.^( (ESNR_dB_16QAM( :, 1, 1 ) - SNR_dB')/10 ) ,  SNR_dB , 10.^( (ESNR_dB_16QAM( :, 2, 1 ) - (SNR_dB )')/10 )  , '--' ,  SNR_dB , 10.^( (ESNR_dB_16QAM_min(:, 2)  - (SNR_dB )')/10 )   , '-.' ,  SNR_dB , 10.^(( ESNR_dB_16QAM( :, 3, 1 ) - (SNR_dB )')/10 )  , '--' ) ;
xlim([0 30])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;

figure( 2060 )
plot( SNR_dB  , ESNR_dB_16QAM( :, 1, 1 ) - SNR_dB', ...
    SNR_dB , ESNR_dB_16QAM( :, 2, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 2) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 3, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 3) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 4, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 4) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 5, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 5) - (SNR_dB )', '-.' ) ;
%     SNR_dB , ESNR_dB_16QAM( :, 6, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 6) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_16QAM( :, 7, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 7) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_16QAM( :, 8, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 8) - (SNR_dB )', '-.' , ...
%     SNR_dB , ESNR_dB_16QAM( :, 9, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 9) - (SNR_dB )', '-.' ) ;
xlim([0 30])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Modulation Efficiency' ) ;

return

g1 = 0.25 ;
g2 = 1.00 ;



RUNS = 1000 ;
clock
tic

%SNR0_dB = -log10(5)*10:0.01:30 ;
SNR0_dB = -50:0.01:50 ;
SNR0 = 10.^( SNR0_dB./10 );

for m = 1 : length( SNR0 )
    
    sigma = sqrt( P / SNR0( m ) / 2 ) ;

    for p = 1 : N_1
        w1 = ( randn(N_1,RUNS) + j*randn(N_1,RUNS) ) * sigma ;
        E_16QAM(p) = mean( log2( sum( exp( - ( abs( s_16QAM(p).*sqrt( P/P_16QAM ) + w1 - s_16QAM*ones(1,RUNS).*sqrt( P/P_16QAM ) ).^2 - abs(w1).^2 ) ./ (2*sigma^2) ) ) ) );
    end
    R_16QAM( m ) = log2( N_1 ) - mean( E_16QAM ) ; 
    
    for p = 1 : N_2
        w2 = ( randn(N_2,RUNS) + j*randn(N_2,RUNS) ) * sigma ;
        E_QPSK(p) = mean( log2( sum( exp( - ( abs( s_QPSK(p).*sqrt( P/P_QPSK ) + w2 - s_QPSK*ones(1,RUNS).*sqrt( P/P_QPSK ) ).^2 - abs(w2).^2 ) ./ (2*sigma^2) ) ) ) );
    end
    R_QPSK( m ) = log2( N_2 ) - mean( E_QPSK ) ; 
end
% 
% delta1 = round( log10(5)*10 / 0.01 + 1 )
% delta2 = round( log10(5/4)*10 / 0.01 + 1 )
% delta12 = round( log10(5)*10 / 0.01 - log10(5/4)*10/0.01 ) 
% 
% clock
% figure( 1000 ) ;
% plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R_16QAM(delta1:end), SNR0_dB(delta1:end), R_QPSK(delta1:end), SNR0_dB(delta1:end), R_16QAM(delta1:end) - R_QPSK(1:(end-delta1 + 1)),'--', SNR0_dB(delta1:end), R_QPSK(delta12:(end-delta2 )), '-.', SNR0_dB(delta1:end) ,R_QPSK(1:(end-delta1 + 1)),'-.' ) ;
% xlim([0 30]);
% ylim([0 6]) ;
% grid ;
% legend('Shannon Bound','C_{16QAM}(SNR)','C_{QPSK}(SNR)','C_{16QAM}-C_{QPSK}(0.2*SNR)','C_{QPSK}(0.8*SNR)','C_{QPSK}(0.2*SNR)') ;
% 
% figure(1010 ) ;
% plot( SNR0_dB(delta1:end), log2(1+SNR0(delta1:end) ), SNR0_dB(delta1:end), R_QPSK(delta1:end), SNR0_dB(delta1:end), R_QPSK( (delta12+1):(end-delta2 + 1 )), '-.', SNR0_dB(delta1:end) ,R_QPSK(1:(end-delta1 + 1)),'--' ) ;
% xlim([0 30]);
% ylim([0 6]) ;
% grid ;
% legend('Shannon Bound', 'C_{QPSK}(SNR)','C_{QPSK}(0.8*SNR)','C_{QPSK}(0.2*SNR)') ;
% return

SNR_dB  = 0:2:30 ;
SNR     = 10.^(SNR_dB./10) ;
RUNS    = 500 ;

N_1 = 4 ;                  %% QPSK/16QAM
N_2 = 4 ;                  %% QPSK
s1 = qammod( [0:1:N_1-1]' , N_1, 0 ) ;
s2 = qammod( [0:1:N_2-1]' , N_2, 0 ) ;   
P_1 = P_QPSK ;
P_2 = P_QPSK ;
for k = 1 : length( SNR )
    sigma = sqrt( P / SNR( k ) / 2 ) ;

    for n = 1 : length( beta )

        P1 = P * beta( n ) ;
        P2 = P * ( 1 - beta( n ) ) ;
        
        SNR1 = SNR(k)*beta(n) ;
        SNR2 = SNR(k)*( 1-beta(n) );

        R_QPSK1( n, k ) = R_QPSK( round( (log10(SNR1+0.00001)*10+50)/0.01 ) + 1) ;
        R_QPSK2( n, k ) = R_QPSK( round( (log10(SNR2+0.00001)*10+50)/0.01 ) + 1) ;

        for m = 1 : length( theta )

            for k = 1 : N_2
                s0( (1+(k-1)*N_1):(k*N_1), 1 ) = s1*sqrt( P1/P_1 ) + s2( k )*sqrt( P2/P_2 )*exp(j*theta(m)) ;
            end

            for p = 1 : N_s
                w = ( randn(N_s,RUNS) + j*randn(N_s,RUNS) ) * sigma ;
                E_Q1(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g1) + w - s0*ones(1,RUNS)*sqrt(g1) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );
                E_Q2(p) = mean( log2( sum( exp( - ( abs( s0(p)*sqrt(g2) + w - s0*ones(1,RUNS)*sqrt(g2) ).^2 - abs(w).^2 ) ./ (2*sigma^2) ) ) ) );            
            end

            R10( m, n, k ) = log2( N_s ) - mean( E_Q1 ) ;
            R20( m, n, k ) = log2( N_s ) - mean( E_Q2 ) ; 

            R11( m, n, k ) = R10( m, n, k ) - R_QPSK1( n, k ) ;        
            R21( m, n, k ) = R20( m, n, k ) - R_QPSK2( n, k ) ;      

        end
    end
end

save(  strcat( 'layered_mod_Effective_SINR_', num2str(now), '.mat' ) ) ;
clock
toc

figure( 1000 )
plot( SNR_dB, log2(1+SNR), SNR_dB, R_16QAM, SNR_dB, R_QPSK ) ;
ylim([0 6])
grid

return
sigma = sqrt(2)/2 ;

L_theta = length( theta ) ;


N_1 = 16 ;                  %% QPSK/16QAM
N_2 = 4 ;                   %% QPSK
N_s = N_1 * N_2 ;           %% The total constellation size => 16

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;
s2 = qammod( [0:1:N_2-1]' , N_2, 0 ) ;    
s = zeros( N_1*N_2, L_theta ) ;

P_1 = mean( abs( s1 ).^2 ) 
P_2 = mean( abs( s2 ).^2 ) 

clock
tic

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

RUNS = 500 ;

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
plot(  R2(1,:), R1(1,:) , R2(2,:), R1(2,:) ,  R2(3,:), R1(3,:) ,  R_QPSK2, R11(1,:) , '--', R_QPSK2, R10(1,:), '--', R_QPSK2, R20(1,:), '--', R_QPSK2, R_QPSK1, '-o' ) ;
grid ;
xlabel('Spectral Efficiency in Bad Channel (|h_1|^2=-10dB), bps/symbol') ;
ylabel('Spectral Efficiency in Good Channel (|h_2|^2=0.0dB), bps/symbol') ;

figure( 20 ) ;
plot( R1(1,:) , R2(1,:), R1(2,:) , R2(2,:), R1(3,:) , R2(3,:), R11(1,:) , R21(1,:), '--' , R11(1,:) , R_QPSK1, '--', R11(1,:) , R20(1,:), '--' ,  R11opt, R21opt, '-.' , R11opt, R_QPSK1, '-.', R11opt, R20opt, '-.') ;
grid ;
xlabel('Spectral Efficiency in Bad Channel (|h_1|^2=-10dB), bps/symbol') ;
ylabel('Spectral Efficiency in Good Channel (|h_2|^2=0.0dB), bps/symbol') ;

figure( 30 ) ;
plot( R1(1,:) , R2(1,:), R1(2,:) , R2(2,:), R1(3,:) , R2(3,:), R11(1,:) , R20(1,:), '--' , R11(1,:) , R_QPSK1, '--', R11opt, R20opt, '-.' , R11opt, R_QPSK1, '-.'  ) ;
grid ;
xlabel('Base-Layer Rate R1') ;
ylabel('Enhancement-Layer Rate R2; Total Rate R') ;

figure( 40 ) ;
plot( beta , R10(1,:), beta, R_QPSK1,'--', beta, R11(1,:), '-.', beta, R20(1,:) , '-+' , beta, R_QPSK2, '-s', beta, R21(1,:), '-d' ) ;
grid ;
xlabel('Power-Splitting \beta') ;
ylabel('Achievable Rate R') ;

figure( 50 ) ;
plot( beta, R11(1,:), beta, R20(1,:) , '--' , beta, R_QPSK2, '-.', beta, R21(1,:), '-o' ) ;
grid ;
xlabel('Power-Splitting \eta') ;
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