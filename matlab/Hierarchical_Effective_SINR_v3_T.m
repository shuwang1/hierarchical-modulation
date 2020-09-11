SNR0_dB = -20:0.0001:20 ;
SNR0    = 10.^(SNR0_dB./10);

EbNo_QPSK  = SNR0/2 ;
EbNo_16QAM = SNR0/4 ;

% Convert from EbNo to SNR.
% Note: Because No = 2*noiseVariance^2, we must add 3 dB
% to get SNR. For details, see Proakis book listed in
% "Selected Bibliography for Performance Evaluation."

EbNo_dB0_QPSK  = SNR0_dB  - 10*log10( 2 ) ;
EbNo_dB0_16QAM = SNR0_dB  - 10*log10( 4 ) ;

ber_QPSK            = berawgn( EbNo_dB0_QPSK,  'qam' , 4  ) ;
ber_16QAM           = berawgn( EbNo_dB0_16QAM, 'qam' , 16 ) ;
ber_QPSK_calc       = qfunc( sqrt(SNR0) ) ;
ber_16QAM_calc      = 0.75*qfunc(sqrt(SNR0/5)) + 0.5*qfunc(3*sqrt(SNR0/5)) - 0.25*qfunc(sqrt(5*SNR0)) ;

figure( 1000 )
semilogy( EbNo_dB0_QPSK, ber_QPSK, EbNo_dB0_QPSK, ber_QPSK_calc, '--' , EbNo_dB0_16QAM , ber_16QAM, EbNo_dB0_16QAM, ber_16QAM_calc , '-.' ) ;
grid ;
ylim([10^(-7) 1])
xlabel( 'E_b/N_o' ) ;
ylabel( 'BER' ) ;

figure( 1010 )
semilogy( SNR0_dB, ber_QPSK, SNR0_dB, ber_QPSK_calc, '--' , SNR0_dB , ber_16QAM, SNR0_dB, ber_16QAM_calc , '-.' ) ;
grid ;
ylim([10^(-10) 1])
xlabel( 'E_s/N_o' ) ;
ylabel( 'BER' ) ;


for m = 1 : length(SNR_dB)
    
    snr_dB = SNR_dB( m ) - 3 ;
    
    for n = 1 : length(beta)
        for k = 1 : length(theta)
            
            [ diff, index ] = min( abs( ber_16QAM_calc - BER_16QAM( m, n, k ) ) ) ;
            ESNR_16QAM( m, n, k ) = SNR0( index ) ;
            ESNR_dB_16QAM( m, n, k ) = SNR0_dB( index ) ;           
            
        end
        
        [ BER_16QAM_min( m, n ), index ] = min( BER_theta( 2, : ) ) ;       
        theta_16QAM( m , n ) = theta( index ) ; 
        
        ESNR_16QAM_min( m, n ) = ESNR_16QAM( m, n, index ) ;
        ESNR_dB_16QAM_min( m, n ) = ESNR_dB_16QAM( m, n, index ) ;
    end
end

save(  strcat( 'Effective_SNR_QPSK_16QAM', num2str(now), '.mat' ) ) ;


figure( 2000 )
semilogy( EbNo_dB0_QPSK, ber_QPSK, EbNo_dB0_QPSK, ber_QPSK_calc, '--' , ...
    EbNo_dB_QPSK , BER_QPSK(:, 1, 1) , '-.', EbNo_dB_QPSK, BER_QPSK(:, 2, 1) , '-.', EbNo_dB_QPSK, BER_QPSK(:, 3, 1) , '-.' , ...
    EbNo_dB0_16QAM, ber_16QAM, EbNo_dB0_16QAM, ber_QPSK_calc, '--' , ...
    EbNo_dB_16QAM, BER_16QAM(:, 1, 1) , '-.', EbNo_dB_16QAM, BER_16QAM(:, 2, 1) , '-.', EbNo_dB_16QAM, BER_16QAM(:, 3, 1) , '-.' ) ;
grid ;
ylim([10^(-7) 1])
xlabel( 'E_b/N_o' ) ;
ylabel( 'BER' ) ;

% figure( 2010 )
% semilogy( SNR0_dB, ber_QPSK.*2, SNR0_dB, ber_QPSK_calc.*2, '--' , SNR0_dB, SER_QPSK( :, 1, 1 ), '-o', SNR0_dB, SER_QPSK( :, 2, 1 ), '--') ;
% grid ;
% ylim([10^(-7) 1])
% xlabel( 'E_s/N_o' ) ;
% ylabel( 'SER' ) ;


figure( 2020 )
plot( SNR0_dB, SNR0_dB, SNR_dB , ESNR_dB_QPSK( :, 1, 1 ), '--' ) ;
hold on ;
figure( 2020 )
plot( SNR_dB , ESNR_dB_QPSK( :, 2, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 2), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 3, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 3), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 4, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 4), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 5, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 5), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 6, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 6), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 7, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 7), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 8, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 8), '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 9, 1 ), '--' , SNR_dB , ESNR_dB_QPSK_min(:, 9), '-.' ) ;
xlim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;


figure( 2025 )
plot( SNR_dB, SNR_dB, '-', SNR_dB , ESNR_dB_16QAM( :, 1, 1 ), '--' ) ;
hold on ;
figure( 2025 )
plot( SNR_dB , ESNR_dB_16QAM( :, 2, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 2), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 3, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 3), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 4, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 4), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 5, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 5), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 6, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 6), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 7, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 7), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 8, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 8), '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 9, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 9), '-.' ) ;
xlim([0 20]) ;
grid ;
xlabel( 'E_s/N_o' ) ;
% ylabel( 'Effective E_s/N_o' ) ;


figure( 2027 )
plot( SNR0_dB, SNR0_dB, SNR_dB , ESNR_dB_16QAM( :, 1, 1 ), '--' , SNR_dB , ESNR_dB_16QAM( :, 2, 1 ), '--' , SNR_dB , ESNR_dB_16QAM_min(:, 2), '-.' ) ;
xlim([0 20])
ylim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;


figure( 2030 )
plot(SNR_dB, 10.^( (ESNR_dB_QPSK( :, 1, 1 ) - SNR_dB')/10 )     , ...
    SNR_dB , 10.^( (ESNR_dB_QPSK( :, 2, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^( (ESNR_dB_QPSK_min(:, 2)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_QPSK( :, 3, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 3)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^( (ESNR_dB_QPSK( :, 4, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^( (ESNR_dB_QPSK_min(:, 4)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_QPSK( :, 5, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 5)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_QPSK( :, 6, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 6)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_QPSK( :, 7, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 7)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_QPSK( :, 8, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 8)  - (SNR_dB )')/10 )  , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_QPSK( :, 9, 1 ) - (SNR_dB )')/10 )  , '--' , SNR_dB , 10.^(( ESNR_dB_QPSK_min(:, 9)  - (SNR_dB )')/10 )  , '-.' ) ;
xlim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;

figure( 2040 )
plot( SNR_dB  , ESNR_dB_QPSK( :, 1, 1 ) - SNR_dB', ...
    SNR_dB , ESNR_dB_QPSK( :, 2, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 2) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 3, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 3) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 4, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 4) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 5, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 5) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 6, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 6) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 7, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 7) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 8, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 8) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_QPSK( :, 9, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_QPSK_min(:, 9) - (SNR_dB )', '-.' ) ;
xlim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Modulation Efficiency' ) ;


figure( 2050 )
plot(SNR_dB, 10.^( (ESNR_dB_16QAM( :, 1, 1 ) - SNR_dB')/10 )     , ...
    SNR_dB , 10.^( (ESNR_dB_16QAM( :, 2, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^( (ESNR_dB_16QAM_min(:, 2)  - (SNR_dB )')/10 )   , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 3, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 3)  - (SNR_dB )')/10 )   , '-.', ...
    SNR_dB , 10.^( (ESNR_dB_16QAM( :, 4, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^( (ESNR_dB_16QAM_min(:, 4)  - (SNR_dB )')/10 )   , '-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 5, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 5)  - (SNR_dB )')/10 )   , '-.', ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 6, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 6)  - (SNR_dB )')/10 )   ,'-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 7, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 7)  - (SNR_dB )')/10 )   ,'-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 8, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 8)  - (SNR_dB )')/10 )   ,'-.' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM( :, 9, 1 ) - (SNR_dB )')/10 )  , '--' , ...
    SNR_dB , 10.^(( ESNR_dB_16QAM_min(:, 9)  - (SNR_dB )')/10 )   ,'-.' ) ;
xlim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;


figure( 2055 )
plot(SNR_dB, 10.^( (ESNR_dB_16QAM( :, 1, 1 ) - SNR_dB')/10 ) ,  SNR_dB , 10.^( (ESNR_dB_16QAM( :, 2, 1 ) - (SNR_dB )')/10 )  , '--' ,  SNR_dB , 10.^( (ESNR_dB_16QAM_min(:, 2)  - (SNR_dB )')/10 )   , '-.' ,  SNR_dB , 10.^(( ESNR_dB_16QAM( :, 3, 1 ) - (SNR_dB )')/10 )  , '--' ) ;
xlim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Effective E_s/N_o' ) ;



figure( 2060 )
plot( SNR_dB  , ESNR_dB_16QAM( :, 1, 1 ) - SNR_dB', ...
    SNR_dB , ESNR_dB_16QAM( :, 2, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 2) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 3, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 3) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 4, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 4) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 5, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 5) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 6, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 6) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 7, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 7) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 8, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 8) - (SNR_dB )', '-.' , ...
    SNR_dB , ESNR_dB_16QAM( :, 9, 1 ) - (SNR_dB )', '--' , SNR_dB , ESNR_dB_16QAM_min(:, 9) - (SNR_dB )', '-.' ) ;
xlim([0 20])
grid ;
xlabel( 'E_s/N_o' ) ;
ylabel( 'Modulation Efficiency' ) ;
