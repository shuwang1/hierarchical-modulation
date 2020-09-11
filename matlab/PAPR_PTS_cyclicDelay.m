
clear all; close all ;

M = 4

L = 128 ;

bin_dB = -ceil( 10*log10(L) ) : 0.5 : ceil( 10*log10(L) ) ;
bin = 10.^(bin_dB./10) ;

trials = 100000 ;

for k = 1 : trials

    b = unidrnd( M, L , 1 ) - 1 ;
    s0 = qammod( b, M ) ;
    % s0  = pskmod( b, M ) ;
    s1 = modnorm( s0, 'avpow' , 1 ) * s0 ;
    % st      = ifft( s1, L ).*sqrt(L) ;
    
    st1 = ifft( [s1(1:L/2); zeros(L/2, 1) ], L ).*sqrt(L) ;
    st2 = ifft( [zeros(L/2, 1); s1((L/2+1):L);  ], L ).*sqrt(L) ;
    
    pt      = ( st1 + st2 ) .* conj( st1 + st2 ) ;
    pt_m    = mean(pt) ;
    pt_max  = max( pt ) ;
    PAPR_max(:, k) = 10*log10( pt_max / pt_m ) * ones( 3, 1) ;

    for l_ = 1 : L - 1
        
        st3 = circshift( st2, l_ ) ;
        pt      = ( st1 + st3 ).* conj( st1 + st3 ) ;
        pt_m    = mean(pt) ;
        pt_max  = max( pt ) ;
        temp    = 10*log10( pt_max / pt_m ) ;
        
        if temp < PAPR_max(2, k)            
            PAPR_max(2, k) = temp ;            
        end
        
        st3 = st2.*exp( -j*2*pi/L*l_ ) ;
        pt      = ( st1 + st3 ).* conj( st1 + st3 ) ;
        pt_m    = mean(pt) ;
        pt_max  = max( pt ) ;
        temp    = 10*log10( pt_max / pt_m ) ;
        
        if temp < PAPR_max(3, k)            
            PAPR_max(3, k) = temp ;            
        end
        
    end
end

pt_count = histc( PAPR_max(1,:), bin_dB ) ;     
pt_pr = pt_count ./ trials ;

pt_count2 = histc( PAPR_max(2,:), bin_dB ) ;     
pt_pr2 = pt_count2 ./ trials ;

pt_count3 = histc( PAPR_max(3,:), bin_dB ) ;     
pt_pr3 = pt_count3 ./ trials ;



len_bin = length( bin_dB ) ;
pt_CCDF = zeros( len_bin, 1 ) ;
pt_CCDF = pt_pr ;

pt_CCDF2 = pt_pr2 ;
pt_CCDF3 = pt_pr3 ;
 
for k = len_bin -1 : -1 : 1
    pt_CCDF( k ) = pt_CCDF( k + 1 ) + pt_CCDF( k ) ;
    pt_CCDF2( k ) = pt_CCDF2( k + 1 ) + pt_CCDF2( k ) ;
    pt_CCDF3( k ) = pt_CCDF3( k + 1 ) + pt_CCDF3( k ) ;
end
 
         figure( 200 ) ;
         semilogy( bin_dB , 1 - ( 1 - exp(-bin) ).^L , bin_dB , pt_CCDF , '-<' , bin_dB , pt_CCDF2 , '-o' , bin_dB , pt_CCDF3 , '-d' ) ;
         xlabel('Peak-to-Average Power Ratio (dB)')
         ylabel('CCDF');
         xlim([0 12])
         ylim([10^(-4) 1 ])
         grid ;
         legend('Guassian Approximation', ' Original', ' with Cyclic Delay, G=2 ', ' with PTS, G=2') ;