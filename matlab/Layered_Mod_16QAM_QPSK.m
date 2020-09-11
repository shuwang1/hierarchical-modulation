clear all;
close all ;

alpha = 1 ;
gamma = 0.000001:0.001:0.1 ;      %% The power splitting in terms of minimum Euclid distance
gamma_dB = -10:-1:-30 ;
% beta    = 10.^(gamma_dB./10)*alpha ;
beta = gamma*alpha
beta0   = 0.025 ;

SNR_dB  = 0:3:30 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0:0.01:1 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;


N_1 = 16 ;                  %% 16QAM
L1  = zeros( N_1, 1 ) ;
N_2 = 4 ;                   %% QPSK
L2  = zeros( N_2, 1 ) ;
N_s = N_1 * N_2 ;           %% The total constellation size => 64

s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;

P_1 = mean( norm( s1, 2 )^2 ) ;
P_2 = mean( norm( qammod( [0:1:N_2-1]' , N_2, 0 ), 2 )^2 ) ;

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
text( real( s(33) ),imag( s(33) ), ['   ' 'a'] , 'FontSize',20);
text( real( s(49) ),imag( s(49) ), ['   ' 'b'] , 'FontSize',20);
text( real( s(5) ),imag( s(5) ), ['   ' 'c'] , 'FontSize',20);


h = scatterplot( s(:,40), 1, 0, 'bs' ) ;
hold on
scatterplot( s1(:,1), 1, 0, 'ro', h ) ;
%title('Constellation for Hierachical Modulation 16QAM/QPSK');
grid ;
text( real( s(33) ),imag( s(33) ), ['   ' 'a'] , 'FontSize',20);
text( real( s(49) ),imag( s(49) ), ['   ' 'b'] , 'FontSize',20);
text( real( s(5) ),imag( s(5) ), ['   ' 'c'] , 'FontSize',20);




% hold on ;
% for m = 1:1:N_1
%    text( real( s1(m) ),imag( s1(m) ), ['   ' num2str(m-1)] , 'FontSize',8);
% end

for m = 1:1:L_theta
    ss2 = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    for n = 1:1:length(beta)
       ss_A( m, n ) = s1(1) + ss2(3).*sqrt(P_1/P_2*beta(n)) ;
       ss_B( m, n ) = s1(1) + ss2(4).*sqrt(P_1/P_2*beta(n)) ;
       ss_C( m, n ) = s1(5) + ss2(1).*sqrt(P_1/P_2*beta(n)) ;
       minED( m, n) = min( norm(ss_A( m, n )-ss_C( m, n )), norm(ss_B( m, n )-ss_C( m, n )) ) ;
    end
end

h = scatterplot( ss_A(:,50), 1, 0, 'rd' ) ;
hold on ;
scatterplot( ss_B(:,50), 1, 0, 'bd', h ) ;

hold on ;
scatterplot( ss_C(:,50), 1, 0, 'yo', h ) ;
grid ;

figure(10) ;
mesh( theta./(pi/4)*45, beta,  minED') ;
xlabel('Rotation Angle \theta (degree)') ;
ylabel('Power Splitting Ratio \gamma') ;
zlabel('Minimum Euclid Distance \Delta = min\{|a-c|,|b-c|\}')

figure(20) ;
plot(  theta./(pi/4)*45, minED( :, 20), theta./(pi/4)*45, minED( :, 30), theta./(pi/4)*45, minED( :, 40), theta./(pi/4)*45, minED( :, 50),'-.' , theta./(pi/4)*45, minED( :, 60),'--', theta./(pi/4)*45, minED( :, 70),'--', theta./(pi/4)*45, minED( :, 80),'--'  , theta./(pi/4)*45, minED( :, 90),'--' )
xlabel('Rotation Angle \theta (degree)') ;
ylabel('Minimum Euclid Distance \Delta = min\{|a-c|,|b-c|\}') ;
legend('P_2/P_1=0.02','P_2/P_1=0.03', 'P_2/P_1=0.04','P_2/P_1=0.05','P_2/P_1=0.06','P_2/P_1=0.07','P_2/P_1=0.08','P_2/P_1=0.09')
grid ;
return


RUNS = 5000 ;

A2 = sqrt(beta0) ;
for n = 1 : length( theta )
    s2 = qammod( 0:1:N_2-1 , N_2, theta(n) ) ;
    
    for m = 1 : length( beta )
        
        
    end
    
end

    L20(1, :) = L1(1) * A2 * ones(1,L_theta) ; 
    L20(2, :) = L1(2) * A2 * ones(1,L_theta) ; 
    L20(3, :) = L1(3) * A2 * ones(1,L_theta) ;
    L20(4, :) = L1(4) * A2 * ones(1,L_theta) ; 


figure(1)
plot( x , y ) ;
grid;

figure(2)
plot( reshape( real(s), N_s*length(theta), 1 ), reshape( imag(s), N_s*length(theta), 1 ), '.' , real(L1), imag(L1),'o', real(s(:,1)), imag(s(:,1)), 'd'  );
grid;



for m = 1 : length(gamma)
        
    L22 = L2.*A1.*sqrt( gamma(m) ) ;
        
    D_Euclid( m, : ) = min( abs( L11(2).*A1 + L22(4,:) - ( L11(1) + L22( 2,:) ) ), abs( L11(2) + L22(4,:) - ( L11(1) + L22( 3,:) ) ) ) ;
    
end

figure(1000)
mesh( theta./(pi/4)*45, gamma, D_Euclid )

for m = 1 : length(SNR)
    SNR( m ) 
    for n = 1 : length(theta)
        % A2 = sqrt(beta0) ;
        sigma = A1/sqrt( SNR( m ) ) ;
 
        for p = 1 : N_s
            
            for t = 1 : RUNS
             
                w = ( randn(1) + j*randn(1) ) * sigma ;
                                
                for q = 1 : N_s
                    d2 = ( abs( s(p,n) + w - s(q,n) ) )^2 - ( abs(w) )^2 ;
                    QQ( t, q ) = exp( - d2 / ( 2*sigma^2 ) ) ;                  
                end
                
            end            
            
            E_Q(p) = mean( log2( sum( QQ, 2 ) ) ) ;                      
            
        end
        
               
        Capacity(m,n) = log2( N_s ) - mean( E_Q ) ; 
    
    end
end


figure(3)
mesh( theta./(pi/4)*45, SNR, Capacity ) ;

figure(4)
plot( SNR_dB, Capacity(:,1), SNR_dB, max( Capacity, [], 2  ), '-d', 1:10, log2(1+10.^((1:10)./10)*(1+beta0^2)),'--'  )
xlabel('SNR (dB)')
ylabel('Spetral Efficiency')
grid;