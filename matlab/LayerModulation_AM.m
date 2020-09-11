%% LayerModulation_AM

global ER N1 N2 Ns

s1 = ( 0:1:(N1-1) ) - (N1-1)./2
P1 = mean( s1.^2 )
A1 = sqrt(P1)

s2 = ( 0:1:(N2-1) ) - (N2-1)./2
P2 = mean( s2.^2 )
A2 = sqrt( P2 )

A_ER = sqrt(ER)
As = sqrt( 1 + ER ) ;

for m = 1:length(ER)
    for n = 1:1:N1
        s0( (1+(n-1)*N2):(n*N2) , m ) = s1(n)./A1.*A_ER + s2./A2;
    end
end

figure(100)
plot( s0, ones([1 Ns]), 'd' , s0, -ones([1 Ns]), 'o' )
grid

RUNS = 10000 ;
w0      = randn(Ns,RUNS) ;
w0_2    = abs(w0).^2 ;

SNR_dB  = 0 : 1 : 30 ;
SNR     = 10.^(SNR_dB./10);

ss = zeros([Ns 1]) ;

clock
tic

for n = 1 : length( ER )
    for k = 1 : N1
        s0( (1+(k-1)*N2):(k*N2), n ) = s1(k).*(A_ER/A1/As) + s2./(A2*As) ;
    end
    for r = 1 : length( SNR )
        ss = s0( :, n ) .* sqrt( SNR(r) ) ;
        for p = 1 : Ns
            E_Q(p) = mean( log2( sum( exp( - ( abs( ss(p,n) + w0 - ss(:,n)*ones(1,RUNS) ).^2 - w0_2 )./2 ), 1 ) ) );
        end
        Capacity( r, n) = log2( Ns ) - mean( E_Q(1:Ns) ) ;
    end
end
save(  strcat( 'layered_modualtion_AM_AM', num2str(now), '.mat' ) ) ;
clock
toc

figure(200)
plot( SNR_dB, log2(1+SNR)./2, '--' , SNR_dB, Capacity( :, 1), SNR_dB, log2(1+SNR), '-.' , SNR_dB, Capacity( :, 1).*2  )
ylim( [0  4.1] )
grid ;
