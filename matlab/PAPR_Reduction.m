
clear all; close all ;

M = 4

L = 128 ;

bin_dB = -ceil( 10*log10(L) ) : 0.5 : ceil( 10*log10(L) ) ;
bin = 10.^(bin_dB./10) ;

trials = 10000 ;

theta   = 1.0/L : 1.0/L : 1 ;
theta   = theta .* (pi/2) ;


for k = 1 : trials

    b = unidrnd( M, L , 1 ) - 1 ;
    s0 = qammod( b, M ) ;
    % s0  = pskmod( b, M ) ;
    s1 = modnorm( s0, 'avpow' , 1 ) * s0 ;
    st1 = ifft( s1, L ).*sqrt(L) ;

    b = unidrnd( M, L , 1 ) - 1 ;

    for m = 1 : 5 
        
        s0 = qammod( b, M ) ;
        % s0  = pskmod( b, M ) ;
        s1 = modnorm( s0, 'avpow' , 1 ) * s0 ;
        st2 = ifft( s1.*m/10.0, L ).*sqrt(L) ;

        pt      = ( st1 + st2 ).* conj( st1 + st2 ) ;
        pt_m    = mean(pt) ;
        pt_max  = max( pt ) ;
        PAPR_max( k, m ) = 10*log10( pt_max / pt_m ) ;
        PAPR_max2( k, m ) = PAPR_max( k, m ) ;
        delay( k , m ) = 0 ;

        for l_ = 1 : L - 1

            st2 = circshift( st2, 1 ) ;

            pt      = ( st1 + st2 ).* conj( st1 + st2 ) ;
            pt_m    = mean(pt) ;
            pt_max  = max( pt ) ;
            temp = 10*log10( pt_max / pt_m ) ;

            if temp < PAPR_max2( k, m )

                PAPR_max2( k, m ) = temp ;
                delay( k , m ) = l_ ;

            end
        end
        
    end
    
    
 end


for m = 1 : 5
    pt_count = histc( PAPR_max( :, m ), bin_dB ) ;     
    pt_pr( :, m ) = pt_count ./ trials ;

    pt_count2 = histc( PAPR_max2( :, m ), bin_dB ) ;     
    pt_pr2( :, m ) = pt_count2 ./ trials ;
end



%         figure(100) ;
%         plot( 0:L-1, pt_dB ) ;
%         grid ;
%         xlim( [0 L] ) ;
%         xlabel( 'Subcarrier', 'FontSize', 24 ) ;
%         ylabel( 'Power Spectrum (dB)', 'FontSize', 24 ) ;
% 
%         figure(102) ;
%         bar( bin_dB, pt_pr ) 
%         xlabel('Peak-to-Average-Power Ratio (dB)', 'FontSize', 24) ;
%         ylabel('Probability Density Function', 'FontSize', 24) ;
%         grid
%         return
%         
%         %hold on
%         
%         figure(104) ;
%         plot( bin_dB, pt_pr, '<' ) ;
%         grid
% 
%         [ st_abs, I ] = sort( abs( st ) ) ; 
%         figure(106) ;
%         plot( real( st(I) ), imag( st(I) ) ) ;
%         xlabel(' Real Part of OFDM Signal', 'FontSize', 24)
%         ylabel(' Imaginary Part of OFDM Signal', 'FontSize', 24)
%         grid



len_bin = length( bin_dB ) ;
pt_CCDF = zeros( len_bin, 1 ) ;
pt_CCDF = pt_pr ;

pt_CCDF2 = pt_pr2 ;

% pt_CCDF4 = pt_pr4 ;
 
for k = len_bin -1 : -1 : 1
    for m = 1 : 5
        pt_CCDF( k, m ) = pt_CCDF( k + 1, m ) + pt_CCDF( k, m ) ;

        pt_CCDF2( k, m ) = pt_CCDF2( k + 1, m ) + pt_CCDF2( k, m ) ;
    end
end

save(  strcat( 'Enh_HM_PAPR', num2str(now), '.mat' ) ) ;
 
         figure( 200 ) ;
         semilogy( bin_dB , 1 - ( 1 - exp(-bin) ).^L ,  bin_dB , pt_CCDF(:,1), '-d',  bin_dB , pt_CCDF(:,3), '-d',  bin_dB , pt_CCDF(:,5), '-d'  , bin_dB , pt_CCDF2(:, 1) , '-o', bin_dB , pt_CCDF2(:, 2) , '-o', bin_dB , pt_CCDF2(:, 3) , '-o', bin_dB , pt_CCDF2(:, 4) , '-o', bin_dB , pt_CCDF2(:, 5) , '-o' ) ;
         xlabel('Peak-to-Average Power Ratio (dB)', 'FontSize', 24) ;
         ylabel('CCDF', 'FontSize', 24)
         xlim([0 12])
         ylim([10^(-4) 1 ])
         grid ;
         legend( 'Gaussian Approximation, L=128', 'Regular QPSK/QPSK, P_2/P_1=0.01, L=128', 'Regular QPSK/QPSK, P_2/P_1=0.09, L=128', 'Regular QPSK/QPSK, P_2/P_1=0.25, L=128', 'Enhanced QPSK/QPSK, P_2/P_1=0.01, L=128', 'Enhanced QPSK/QPSK, P_2/P_1=0.04, L=128', 'Enhanced QPSK/QPSK, P_2/P_1=0.09, L=128', 'Enhanced QPSK/QPSK, P_2/P_1=0.16, L=128', 'Enhanced QPSK/QPSK, P_2/P_1=0.25, L=128') ;
         
