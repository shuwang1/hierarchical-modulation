clear all;
close all ;

alpha = 1 ;
gamma = 0.0000001:0.01:1 ;      %% The power splitting in terms of minimum Euclid distance
gamma_dB = -10:-2:-30 ;
% beta    = 10.^(gamma_dB./10)*alpha ;
beta = gamma*alpha
beta0   = 0.3 ;

SNR_dB  = 0:2:30 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0:0.01:1 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;


N_1 = 4 ;                   %% QPSK
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

h = scatterplot( s(:,1), 1, 0, 'bs' ) ;
hold on
scatterplot( s1(:,1), 1, 0, 'ro', h ) ;
%title('Constellation for Hierachical Modulation 16QAM/QPSK');
grid ;
hold on
text( real( s(3,1) ),imag( s(3,1) ), ['   ' 'a'] , 'FontSize',20);
text( real( s(7,1) ),imag( s(7,1) ), ['   ' 'b'] , 'FontSize',20);
text( real( s(13,1) ),imag( s(13,1) ), ['   ' 'c'] , 'FontSize',20);

h = scatterplot( s(:,4), 1, 0, 'bs' ) ;
hold on
scatterplot( s1(:,1), 1, 0, 'ro', h ) ;
%title('Constellation for Hierachical Modulation 16QAM/QPSK');
grid ;
text( real( s(3,4) ),imag( s(3,4) ), ['   ' 'a'] , 'FontSize',20);
text( real( s(7,4) ),imag( s(7,4) ), ['   ' 'b'] , 'FontSize',20);
text( real( s(13,4) ),imag( s(13,4) ), ['   ' 'c'] , 'FontSize',20);


% hold on ;
% for m = 1:1:N_1
%    text( real( s1(m) ),imag( s1(m) ), ['   ' num2str(m-1)] , 'FontSize',8);
% end

for m = 1:1:L_theta
    ss2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    for n = 1:1:length(beta)
       ss_A( m, n ) = s1(3) + ss2(1).*sqrt(P_1/P_2*beta(n)) ;
       ss_B( m, n ) = s1(3) + ss2(2).*sqrt(P_1/P_2*beta(n)) ;
       ss_C( m, n ) = s1(1) + ss2(4).*sqrt(P_1/P_2*beta(n)) ;
       minED( m, n) = min( norm(ss_A( m, n )-ss_C( m, n )), norm(ss_B( m, n )-ss_C( m, n )) ) ;
    end
end

h = scatterplot( ss_A(:,4), 1, 0, 'rd' ) ;
hold on ;
scatterplot( ss_B(:,4), 1, 0, 'bd', h ) ;

hold on ;
scatterplot( ss_C(:,4), 1, 0, 'yo', h ) ;
grid ;

figure(10) ;
mesh( theta./(pi/4)*45, beta,  minED') ;
xlabel('Rotation Angle \theta (degree)') ;
ylabel('Power Splitting Ratio \gamma') ;
zlabel('Minimum Euclid Distance \Delta = min\{|a-c|,|b-c|\}')

figure(20) ;
plot(  theta./(pi/4)*45, minED( :, 1), theta./(pi/4)*45, minED( :, 2), theta./(pi/4)*45, minED( :, 3), theta./(pi/4)*45, minED( :, 4),'-.' , theta./(pi/4)*45, minED( :, 5),'--', theta./(pi/4)*45, minED( :, 6),'--', theta./(pi/4)*45, minED( :, 7),'--'  , theta./(pi/4)*45, minED( :, 8),'--' )
xlabel('Rotation Angle \theta (degree)') ;
ylabel('Minimum Euclid Distance \Delta = min\{|a-c|,|b-c|\}') ;
legend('P_2/P_1=0.1','P_2/P_1=0.2', 'P_2/P_1=0.3','P_2/P_1=0.4','P_2/P_1=0.5','P_2/P_1=0.6','P_2/P_1=0.7','P_2/P_1=0.8')
grid ;

RUNS = 5000 ;
s0 = zeros( N_1*N_2, 1 ) ;
PP = zeros( N_1*N_2, 1 ) ;

for m = 1 : length(theta)
    s2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    for n = 1 : length(beta)
        for k = 1:1:N_1
            s0( (1+(k-1)*N_2):(k*N_2), 1 ) = s1(k,1) + s2.*sqrt(P_1/P_2*beta(n)) ;
        end
        
        for r = 1 : length(SNR)
            
            sigma = sqrt( P_1/SNR(r) );
            
            for p = 1 : N_s
                
                w = ( randn( N_s, RUNS ) + j*randn( N_s, RUNS ) ) * sigma ;
                
                PP(p) = mean( log2( sum( exp( ( ( abs(w) ).^2 - ( abs( s0(p) + w - s0*ones(1,RUNS) ) ).^2 ) ./ (2*sigma^2) ) ) ) ) ;
            end
            
            Capacity( r, m, n ) = 4 - mean( PP ) ;
            
        end
        
    end
end

for m = 1 : length(beta)
    for n = 1 : length(SNR)
        C_max( n , m ) = max( Capacity( n, :, m )  ) ;
    end
end

figure(300)
mesh( theta./(pi/4)*45, SNR_dB, Capacity( :, :, 7) ) ;

figure(400)
plot( SNR_dB, Capacity(:, 1, 2 ),'--', SNR_dB, C_max( :, 2 ),'--', SNR_dB, Capacity(:, 1, 4 ),'-.', SNR_dB, C_max( :, 4 ), '-.' , SNR_dB, Capacity(:, 1, 6 ), '-s', SNR_dB, C_max( :, 6 ), '-s' , SNR_dB, Capacity(:, 1, 7 ),'-o' , SNR_dB, C_max( :, 7 ),'-o' )
xlabel('SNR (dB)')
ylabel('Spetral Efficiency')
grid;


figure(401)
plot( SNR_dB, Capacity(:, 1, 2 ),'--', SNR_dB, C_max( :, 2 ),'--', SNR_dB, Capacity(:, 1, 3 ),'-.', SNR_dB, C_max( :, 3 ), '-.'  )
xlabel('SNR (dB)')
ylabel('Spetral Efficiency')
grid;
legend('P_2/P_1=0.1','P_2/P_1=0.2');

figure(402)
plot( SNR_dB, Capacity(:, 1, 4 ),'--', SNR_dB, C_max( :, 4 ),'-s', SNR_dB, Capacity(:, 1, 5 ),'-.', SNR_dB, C_max( :, 5 ), '-o'  )
xlabel('SNR (dB)')
ylabel('Bits/Symbol')
grid;
legend('regular, P_2/P_1=0.5','enhanced, P_2/P_1=0.5','regular, P_2/P_1=0.6','enhanced, P_2/P_1=0.6');


