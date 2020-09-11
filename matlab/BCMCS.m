clear all
close all

SNR_dB = -10:1:20 ;

SNR = 10.^(SNR_dB./10) ;

C(1,:) = log10( 1 + SNR ) ;

figure(100)
plot( SNR_dB, C(1,:) )
grid ;
xlabel('SNR (dB)')
ylabel('Spectral Efficiency')