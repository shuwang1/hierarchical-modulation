%% LayerModulation_QAM

global ER N1 N2 Ns

SNR_dB  = 0 : 1 : 30 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0.0:0.001:1 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;

s1_ini_phase = 0 ;
s1 = qammod( [0:1:N1-1]' , N1, s1_ini_phase ) ;
P1 = mean( abs( s1 ).^2 ) ;
A1 = sqrt(P1) ;

s2_ini_phase = 0 ;
s2 = qammod( [0:1:N2-1]' , N2, s2_ini_phase ) ;
P2 = mean( abs( s2 ).^2 ) ;
A2 = sqrt( P2 ) ;

A_ER = sqrt(ER) ;
As   = sqrt( 1 + ER ) ;

s = zeros( Ns, L_theta ) ;
for m = 1:1:L_theta
    s2( :, m ) = qammod( [0:1:N2-1]' , N2, theta(m) ) ;
end
s = kron( s1, ones(N2,1) )*ones(1, L_theta) + kron( ones(N1,1), s2 )./2 ;

h = scatterplot( s(:,1), 1, 0, 'bs' ) ;
hold on
scatterplot( s1(:,1), 1, 0, 'ro', h ) ;
title('Constellation for Hierachical Modulation QPSK/QPSK');
grid ;
hold on
text( real( s(3,1) ),imag( s(3,1) ), ['   ' 'a'] , 'FontSize',20);
text( real( s(7,1) ),imag( s(7,1) ), ['   ' 'b'] , 'FontSize',20);
text( real( s(13,1) ),imag( s(13,1) ), ['   ' 'c'] , 'FontSize',20);

h = scatterplot( s(:,4), 1, 0, 'bs' ) ;
hold on
scatterplot( s1(:,1), 1, 0, 'ro', h ) ;
title('Constellation for Hierachical Modulation QPSK/QPSK');
grid ;
text( real( s(3,4) ),imag( s(3,4) ), ['   ' 'a'] , 'FontSize',20);
text( real( s(7,4) ),imag( s(7,4) ), ['   ' 'b'] , 'FontSize',20);
text( real( s(13,4) ),imag( s(13,4) ), ['   ' 'c'] , 'FontSize',20);

RUNS = 5000 ;
%w0      = ones(Ns,1)*( randn(1,RUNS) + j*randn(1,RUNS) ) ;
w0      = randn(Ns,RUNS) + j*randn(Ns,RUNS);
w0_2    = abs(w0).^2 ;
w0_2r   = real(w0).^2 ;
w0_2i   = imag(w0).^2 ;

for m = 1 : L_theta
    s2( :, m ) = qammod( [0:1:N2-1]' , N2, theta(m) ) ;
end

for m = 1 : length(ER)
    s3( :, :, m ) = kron( s1, ones(N2,1) )*ones(1,L_theta).*(A_ER(m)/As(m)/A1) + kron( ones(N1,1), s2./(As(m)*A2) ) ; 
end
    
clock
tic

for m = 1 : L_theta
    for n = 1 : length( ER )
        for r = 1 : length(SNR)
            s0 = s3( :, m , n ) .* sqrt( SNR( r )*2 )  ;
            for p = 1 : Ns
                w   = ones(Ns,1)*w0(p,:) ;
                w_2 = ones(Ns,1)*w0_2(p,:) ;
                E_Q(p) = mean( log2( sum( exp( - ( abs( s0(p) + w - s0*ones(1,RUNS) ).^2 - w_2 ) ./ 2 ), 1 ) ) );
            end
            Capacity( r, m, n) = log2( Ns ) - mean( E_Q(1:Ns) ) ;             
        end
    end
end


for m = 1 : L_theta
    for n = 1 : length( ER )
        for r = 1 : length(SNR)
            Cap( n, r, m ) = Capacity( r, m, n) ;             
        end
    end
end

for m = 1 : length(ER)
    for n = 1 : length(SNR)
       [ C_max( n, m ), index ] = max( Capacity( n, :, m)  ) ;
       theta_opt( n, m ) = theta( index ) ;
    end
end

save(  strcat( 'layered_modualtion_QPSK_QPSK', num2str(now), '.mat' ) ) ;
clock
toc

scatterplot( reshape( s3, prod(size(s3)),1 ) )
grid

scatterplot( reshape( w0, prod(size(w0)),1 ) )
grid

figure(110)
plot( SNR_dB, Capacity( :, 1, 1), SNR_dB, C_max( :, 1 ), SNR_dB, log2(1+SNR), '--', SNR_dB, C_max( :, 1 ) - Capacity( :, 1, 1) )
ylim( [0 4.1] ) ;
xlabel('SNR (dB)')
ylabel('Achievable Rates (bits/symbol)')
grid ;

figure(120)
plot( SNR_dB, Capacity( :, 1, 2), SNR_dB, C_max( :, 2 ), SNR_dB, log2(1+SNR), '--', SNR_dB, C_max( :, 2 ) - Capacity( :, 1, 2) )
ylim( [0 4.1] ) ;
xlabel('SNR (dB)')
ylabel('Achievable Rates (bits/symbol)')
grid ;

figure(130)
plot( SNR_dB, Capacity( :, 1, 3), SNR_dB, C_max( :, 3 ), SNR_dB, log2(1+SNR), '--', SNR_dB, C_max( :, 3 ) - Capacity( :, 1, 3) )
ylim( [0 4.1] ) ;
xlabel('SNR (dB)')
ylabel('Achievable Rates (bits/symbol)')
grid ;

figure(140)
plot( SNR_dB, Capacity( :, 1, 4), SNR_dB, C_max( :, 4 ), SNR_dB, log2(1+SNR), '--', SNR_dB, C_max( :, 4 ) - Capacity( :, 1, 4) )
ylim( [0 4.1] ) ;
xlabel('SNR (dB)')
ylabel('Achievable Rates (bits/symbol)')
grid ;

figure(200)
plot(SNR_dB, theta_opt(:,1)./pi*180, SNR_dB, theta_opt(:,2)./pi*180, SNR_dB, theta_opt(:,3)./pi*180, SNR_dB, theta_opt(:,4)./pi*180, SNR_dB, theta_opt(:,5)./pi*180, SNR_dB, theta_opt(:,6)./pi*180, SNR_dB, theta_opt(:,7)./pi*180, SNR_dB, theta_opt(:,8)./pi*180 )
xlabel('SNR (dB)')
ylabel('Optimal Rotation Angle (degree)')
grid ;

figure(300)
plot( SNR_dB, C_max( :, 1) - Capacity( :, 1, 1), SNR_dB, C_max( :, 2) - Capacity( :, 1, 2), SNR_dB, C_max( :, 3) - Capacity( :, 1, 3), SNR_dB, C_max( :, 4 ) - Capacity( :, 1, 4), SNR_dB, C_max( :, 5 ) - Capacity( :, 1, 5), SNR_dB, C_max( :, 6 ) - Capacity( :, 1, 6), SNR_dB, C_max( :, 7 ) - Capacity( :, 1, 7), SNR_dB, C_max( :, 8 ) - Capacity( :, 1, 8) )
xlabel('SNR (dB)')
ylabel('Achievable Capacity Gain (bits/symbol)')
grid ;

return

figure(120)
% plot( SNR_dB, Capacity( :, 1, 10), SNR_dB, C_max( :, 10 ) )
plot( SNR_dB, Capacity( :, 1, 3), SNR_dB, C_max( :, 3 ), SNR_dB, log2(1+SNR) )
ylim( [0 4.5] ) ;
grid ;