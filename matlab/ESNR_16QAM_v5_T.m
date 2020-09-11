close all
clear all


N_1 = 16 ;  N_2 = 4 ;


b_16QAM = 0:(N_1-1) 
s_16QAM = qammod( b_16QAM, N_1, 0 ) ;
qamdemod( s_16QAM, N_1, 0 ) 

scatterplot( s_16QAM )
hold on; % Make sure the annotations go in the same figure.
for jj = 1:length( s_16QAM )
   text(  real(s_16QAM(jj)),imag(s_16QAM(jj)),[' ' num2str( b_16QAM(jj) ) ]  );
end
hold off;


b_QPSK = 0 : (N_2-1) ;
s_QPSK = qammod( b_QPSK, N_2, 0 ) 
qamdemod( s_QPSK, 4, 0 )


s_QPSK_QPSK = kron( s_QPSK*2, ones( size( s_QPSK ) ) )      + kron( ones( size( s_QPSK ) ), s_QPSK )
b_QPSK_QPSK = kron( b_QPSK*N_2, ones( size( b_QPSK ) ) )    + kron( ones( size( b_QPSK ) ), b_QPSK ) 
qamdemod( s_QPSK_QPSK, 16, 0 )

scatterplot( s_QPSK_QPSK )
hold on;    % Make sure the annotations go in the same figure.
for jj = 1:length( s_QPSK_QPSK )
   text(  real( s_QPSK_QPSK(jj) ), imag( s_QPSK_QPSK(jj) ), [ ' ' num2str( b_QPSK_QPSK(jj) )  ': '  num2str(jj-1) ]  );
end
hold off;

SNR0_dB = 0 : 1 : 20 ;
SNR0    = 10.^(SNR0_dB./10);

EbNo_dB0_QPSK  = SNR0_dB  - 10*log10( 2 )  ;
EbNo_dB0_16QAM = SNR0_dB  - 10*log10( 4 )  ;

ber_QPSK            = berawgn( EbNo_dB0_QPSK,  'qam' , 4  ) ;
ber_16QAM           = berawgn( EbNo_dB0_16QAM, 'qam' , 16 ) ;

MAX_RUNS = 1.0e5 ;

runs = ceil( ber_16QAM.^(-1) ).*1000 ;

SNR_dB  = 0:1:20 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0:0.05:1 ;
theta   = theta .* (pi/4) ;

S_QPSK = ones( MAX_RUNS, 1 )  * qammod( 0:3, 4, 0 ) ;

beta = 0 ;

for m = 1 : length( SNR_dB )
    
    snr_dB = SNR_dB( m ) ;
    RUNS = min( runs( m ), MAX_RUNS ) ;    
    
    for n = 1 : length( beta )
        
        for k = 1 : length( theta )
            
            zeta = exp( j*theta(k) ) * sqrt(0.4) ;
            
            b0 = randint( RUNS, 1, 4 ) ;    s0 = qammod( b0, 4, 0 ) ;
            b1 = randint( RUNS, 1, 4 ) ;    s1 = qammod( b1, 4, 0 ) * zeta ;
            
            S = kron( S_QPSK(1:RUNS, :), ones( 1, 4 ) ) + kron( ones(1, 4), S_QPSK(1:RUNS, :) ).*zeta ;
            
            r  = awgn( s0, snr_dB, 'measured' ) + s1 ;
            
            [ c, d ] = min( r*ones( 1, 16 ) - S , [], 2 ) ; 
            
            [ num, BER( k, n, m ) ]     = biterr( b0*4+b1, d-1 ) ; 
            [ num, BER_B( k, n, m ) ]   = biterr( b0, floor((d-1)/4) ) ; 
            [ num, BER_E( k, n, m ) ]   = biterr( b1, mod(d-1, 4) ) ;
        end
        
        BER_1( n, m ) = BER( 1, n, m ) ;
        BER_B_1( n ,m ) = BER_B( 1, n, m ) ;
        BER_E_1( n ,m ) = BER_E( 1, n, m ) ;
        
        [ BER_min( n, m ), index ]= min( BER( :, n, m )  ) ; 
        BER_B_min( n ,m ) = BER_B( index, n, m ) ;
        BER_E_min( n ,m ) = BER_E( index, n, m ) ;

        theta_opt( n, m ) = theta( index ) ;        
    end
end

ESNR_B_1        = qfuncinv( BER_B_1 ).^2 ;            
ESNR_dB_B_1     = log10( ESNR_B_1 )*10 ;
ESNR_E_1        = qfuncinv( BER_E_1 ).^2 ;
ESNR_dB_E_1     = log10( ESNR_E_1 ).*10 ;


ESNR_B_min      = qfuncinv( BER_B_min ).^2 ;            
ESNR_dB_B_min   = log10( ESNR_B_min )*10 ;
ESNR_E_min      = qfuncinv( BER_E_min ).^2 ;
ESNR_dB_E_min   = log10( ESNR_E_min ).*10 ;


figure( 1000 ) ;
semilogy( SNR0_dB, ber_QPSK, SNR0_dB, ber_16QAM, SNR_dB , BER_min( 1, : ) , SNR_dB , BER_B_min( 1, : ) ,'-s' , SNR_dB , BER_E_min( 1, : ), '-o' ) ;
grid ;
xlabel( 'SNR (dB)') ;
ylabel('BER') ;
legend( 'QPSK', '16QAM', 'Hierarchical', 'Base-Layer QPSK', 'Enhancement-Layer QPSK' ) ;


figure( 1010 ) ;
plot( SNR0_dB, ESNR_B_1( 1, : )./SNR0, SNR0_dB, ESNR_B_min( 1, : )./SNR0, '-o' , SNR_dB , ESNR_E_1( 1, : )./SNR0/zeta , SNR_dB , ESNR_E_min( 1, : )./SNR0/zeta ,'-s' , SNR_dB , ones(size(SNR_dB)) ) ;
grid ;
ylim([0 1.2] );
xlabel( 'SNR (dB)') ;
ylabel('BER') ;
legend( 'Base Layer', 'Enhanced Base Layer', 'Enhancement Layer', 'Enhanced Enhancement Layer', 'QPSK' ) ;

