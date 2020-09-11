clear all ; close all ;

ER = [ 0 : 0.25 : 10 ]

alpha = 10 ;
gamma = 1./(1 + ER) ;      %% The power splitting in terms of minimum Euclid distance
beta = gamma*alpha ;
beta0   = 0.5*alpha ;

SNR_dB  = 0 : 1 : 30 ;
SNR     = 10.^(SNR_dB./10);

theta   = 0:0.02:1 ;
theta   = theta .* (pi/4) ;
L_theta = length( theta ) ;

N_1 = 4 ;                  %% QPSK
L1  = zeros( N_1, 1 ) ;
s1_ini_phase = 0 ;
s1  = qammod( [0:1:N_1-1]' , N_1, s1_ini_phase ) ;
P_1 = mean( abs( s1 ).^2 ) 

N_2 = 4 ;                   %% QPSK
L2  = zeros( N_2, 1 ) ;
N_s = N_1 * N_2 ;           %% The total constellation size => 16
P_2 = mean( abs( qammod( [0:1:N_2-1]' , N_2, 0 ) ).^2 ) 

s = zeros( N_1*N_2, L_theta ) ;
for m = 1:1:L_theta
    s2( :, m ) = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
    for n = 1:1:N_2
        s( (1+(n-1)*N_1):(n*N_1) , m ) = s1 + s2( n, m ).*sqrt(1.0/4) ;
    end
end

h = scatterplot( s(:,1), 1, 0, 'bs' ) ;
hold on
scatterplot( s1(:,1), 1, 0, 'ro', h ) ;
title('Constellation for Hierachical Modulation 16QAM/QPSK');
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

for n = 1 : length( beta )
    ER( n ) = ( (alpha-beta(n))/P_1 ) / ( beta(n)/P_2 ) ;   
end
ER 

RUNS    = 5000 ;
%w0      = ones(N_s,1) * ( randn(1,RUNS) + j*randn(1,RUNS) ) ;
w0      = randn(N_s,RUNS) + j*randn(N_s,RUNS) ;
w0_2    = w0.*conj(w0) ;

for m = 1 :L_theta
    s2(:,m) = qammod( [0:1:N_2-1]' , N_2, theta(m) ) ;
end
for m = 1 : length( ER )
    s3( :, :, m ) = kron( s1.*sqrt( ER(m)/P_1/(ER(m)+1) ), ones( N_2, 1 ) )*ones(1,L_theta) + kron( ones(N_1, 1), s2.*sqrt( 1/P_2/(ER(m)+1) ) ) ;
end

clock
tic

for m = 1 : length( theta )
    for n = 1 : length( ER )
        for r = 1 : length(SNR)
            ss = s3( :, m, n ) .* sqrt( SNR(r)*2 ) ;
            for p = 1 : N_s
                w = ones(N_s,1)*w0(p,:) ;
                w_2 = ones(N_s,1)*w0_2(p,:) ;
                E_Q(p) = mean( log2( sum( exp( ( -abs( ss(p) + w - ss*ones(1,RUNS) ).^2 + w_2 ) ./ 2 ), 1 ) ) );
            end
            Capacity( r, m, n) = log2( N_s ) - mean( E_Q(1:N_s) ) ;             
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

figure(110)
plot( SNR_dB, Capacity( :, 1, 1), SNR_dB, C_max( :, 1 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 1 ) - Capacity( :, 1, 1) )
ylim( [0 4.1] ) ;
grid ;

figure(120)
plot( SNR_dB, Capacity( :, 1, 4), SNR_dB, C_max( :, 4 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 4 ) - Capacity( :, 1, 4) )
ylim( [0 4.1] ) ;
grid ;

figure(130)
plot( SNR_dB, Capacity( :, 1, 7), SNR_dB, C_max( :, 7 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 7 ) - Capacity( :, 1, 7) )
ylim( [0 4.1] ) ;
grid ;

figure(140)
plot( SNR_dB, Capacity( :, 1, 10), SNR_dB, C_max( :, 10 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 10 ) - Capacity( :, 1, 10) )
ylim( [0 4.1] ) ;
grid ;

figure(150)
plot( SNR_dB, Capacity( :, 1, 13), SNR_dB, C_max( :, 13 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 13 ) - Capacity( :, 1, 13) )
ylim( [0 4.1] ) ;
grid ;

figure(160)
plot( SNR_dB, Capacity( :, 1, 16), SNR_dB, C_max( :, 16 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 16 ) - Capacity( :, 1, 16) )
ylim( [0 4.1] ) ;
grid ;

figure(170)
plot( SNR_dB, Capacity( :, 1, 19), SNR_dB, C_max( :, 19 ), SNR_dB, log2(1+SNR), SNR_dB, C_max( :, 19 ) - Capacity( :, 1, 19) )
ylim( [0 4.1] ) ;
grid ;

figure(100)
plot( SNR_dB, theta_opt(:,1)./pi*180 , SNR_dB, theta_opt(:,2)./pi*180 , SNR_dB, theta_opt(:,3)./pi*180 , SNR_dB, theta_opt(:,4)./pi*180 , SNR_dB, theta_opt(:,5)./pi*180 , SNR_dB, theta_opt(:,6)./pi*180 )
grid ;


figure(200)
plot( SNR_dB, log2(1+SNR),'--' ,SNR_dB, Capacity( :, 1, 1), SNR_dB, Capacity( :, 1, 6),'--' ,SNR_dB, Capacity( :, 1, 11),SNR_dB, Capacity( :, 1, 16),'--' , SNR_dB, Capacity( :, 1, 17), '-.' , SNR_dB, Capacity( :, 1, 21), SNR_dB, Capacity( :, 1, 26),'--' ,SNR_dB, Capacity( :, 1, 31), SNR_dB, Capacity( :, 1, 36),'--' ,SNR_dB, Capacity( :, 1, 41) )
ylim( [0 4.1] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Unconstraint Shannon Capacity','QPSK/QPSK, ER=0.0','QPSK/QPSK, ER=1.25', 'QPSK/QPSK, ER=2.5','QPSK/QPSK, ER=3.75', 'QPSK/QPSK, ER=4.00' , 'QPSK/QPSK, ER=5.0', 'QPSK/QPSK, ER=6.25', 'QPSK/QPSK, ER=7.5','QPSK/QPSK, ER=8.75','QPSK/QPSK, ER=10')

figure(300)
plot( SNR_dB, (log2(1+SNR))' -  Capacity( :, 1, 17), '--' , SNR_dB, Capacity( :, 1, 15) -  Capacity( :, 1, 17), SNR_dB, C_max( :, 15) -  Capacity( :, 1, 17),'--' , SNR_dB, Capacity( :, 1, 16)- Capacity( :, 1, 17), SNR_dB, C_max( :, 16)- Capacity( :, 1, 17), '--', SNR_dB, C_max( :, 17)- Capacity( :, 1, 17) , SNR_dB, Capacity( :, 1, 18)- Capacity( :, 1, 17),SNR_dB, C_max( :, 18)- Capacity( :, 1, 17),'--' , SNR_dB, Capacity( :, 1, 19)- Capacity( :, 1, 17),  SNR_dB, C_max( :, 19)- Capacity( :, 1, 17), '-.',  SNR_dB, Capacity( :, 1, 26)- Capacity( :, 1, 17) )
ylim( [-0.05 0.05] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('R_{Shannon}-R_{4.00}','R_{3.50}-R_{4.00}','R_{3.50}^{opt}-R_{4.00}','R_{3.75}-R_{4.00}','R_{3.75}^{opt}-R_{4.00}','R_{4.00}^{opt}-R_{4.00}' ,'R_{4.25}-R_{4.00}', 'R_{4.25}^{opt}-R_{4.00}','R_{4.50}-R_{4.00}','R_{4.50}^{opt}-R_{4.00}','R_{6.25}-R_{4.00}' )

figure(310)
plot( SNR_dB, (log2(1+SNR))' -  Capacity( :, 1, 17), '--' , SNR_dB, Capacity( :, 1, 15) -  Capacity( :, 1, 17), SNR_dB, C_max( :, 15) -  Capacity( :, 1, 17),'--' , SNR_dB, Capacity( :, 1, 16)- Capacity( :, 1, 17), SNR_dB, C_max( :, 16)- Capacity( :, 1, 17), '--', SNR_dB, C_max( :, 17)- Capacity( :, 1, 17) , SNR_dB, Capacity( :, 1, 18)- Capacity( :, 1, 17),SNR_dB, C_max( :, 18)- Capacity( :, 1, 17),'--' , SNR_dB, Capacity( :, 1, 19)- Capacity( :, 1, 17),  SNR_dB, C_max( :, 19)- Capacity( :, 1, 17), '-.' )
ylim( [-0.05 0.05] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('R_{Shannon}-R_{4.00}','R_{3.50}-R_{4.00}','R_{3.50}^{opt}-R_{4.00}','R_{3.75}-R_{4.00}','R_{3.75}^{opt}-R_{4.00}','R_{4.00}^{opt}-R_{4.00}' ,'R_{4.25}-R_{4.00}', 'R_{4.25}^{opt}-R_{4.00}','R_{4.50}-R_{4.00}','R_{4.50}^{opt}-R_{4.00}' )


SNR0_dB  = -10 : 0.1 : 30 ;
SNR0     = 10.^(SNR0_dB./10);

N_s = 4 
for r = 1 : length(SNR0)
    ss = s1 ./ sqrt( P_1 ) .* sqrt( SNR0(r)*2 ) ;
    for p = 1 : N_s
        w = ones(N_s,1)*w0(p,:) ;
        w_2 = ones(N_s,1)*w0_2(p,:) ;
        E_Q(p) = mean( log2( sum( exp( ( -abs( ss(p) + w - ss*ones(1,RUNS) ).^2 + w_2 ) ./ 2 ), 1 ) ) );
    end
    R_QPSK( r ) = log2( N_s ) - mean( E_Q(1:N_s) ) ;             
end
figure(1000)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK ) ;
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','QPSK')

SNR0_dB(36 : 10 : 336)
R_QPSK_E(:,1) = ( R_QPSK( 33 : 10 : 336 ) )'
R_QPSK_B(:,1) = Capacity( :, 1, 15)-R_QPSK_E(:,1) ;
R_QPSK_B_e(:,1) = C_max( :, 15)-R_QPSK_E(:,1) ;

SNR0_dB(33 : 10 : 333)
R_QPSK_E(:,2) = ( R_QPSK( 33 : 10 : 333 ) )'
R_QPSK_B(:,2) = Capacity( :, 1, 16)-R_QPSK_E(:,2) ;
R_QPSK_B_e(:,2) = C_max( :, 16)-R_QPSK_E(:,2) ;

SNR0_dB(31 : 10 : 331)
R_QPSK_E(:,3) = ( R_QPSK( 31 : 10 : 331 ) )'
R_QPSK_B(:,3) = Capacity( :, 1, 17)-R_QPSK_E(:,3) ;
R_QPSK_B_e(:,3) = C_max( :, 17)-R_QPSK_E(:,3) ;

SNR0_dB(29 : 10 : 329)
R_QPSK_E(:,4) = ( R_QPSK( 29 : 10 : 329 ) )' ;
R_QPSK_B(:,4) = Capacity( :, 1, 18)-R_QPSK_E(:,4) ;
R_QPSK_B_e(:,4) = C_max( :, 18)-R_QPSK_E(:,4) ;

SNR0_dB(27 : 10 : 327)
R_QPSK_E(:,5) = ( R_QPSK( 27 : 10 : 327 ) )'
R_QPSK_B(:,5) = Capacity( :, 1, 19)-R_QPSK_E(:,5) ;
R_QPSK_B_e(:,5) = C_max( :, 19)-R_QPSK_E(:,5) ;

SNR0_dB(50 : 10 : 350)
R_QPSK_E(:,6) = ( R_QPSK( 50 : 10 : 350 ) )'
R_QPSK_B(:,6) = Capacity( :, 1, 10)-R_QPSK_E(:,6) ;
R_QPSK_B_e(:,6) = C_max( :, 10)-R_QPSK_E(:,6) ;
R_QPSK_B_r(:,6) = ( R_QPSK( 85 : 10 : 385 ) )' ;

figure(1050)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17) , SNR_dB, Capacity( :, 1, 10) , SNR_dB, C_max( :, 10) , SNR_dB , R_QPSK_E(:,6), SNR_dB, R_QPSK_B(:,6), SNR_dB, R_QPSK_B_e(:,6), SNR_dB, R_QPSK_B(:,3) ) ;
xlim( [ 0 25 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}', 'R_{2.25}', 'R_{2.25}^{opt}' , 'R_{E;2.25}', 'R_{B;2.25}', 'R_{B;2.25}^{opt}','R_{B;4.0}')

figure(1100)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17) , SNR_dB, Capacity( :, 1, 15) , SNR_dB , R_QPSK_E(:,1), SNR_dB, R_QPSK_B(:,1) ) ;
xlim( [ 0 30 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}', 'R_{3.50}' , 'R_{E;3.50}', 'R_{B;3.50}')

figure(1200)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17), SNR_dB, Capacity( :, 1, 16) , SNR_dB , R_QPSK_E(:,2), SNR_dB, R_QPSK_B(:,2) ) ;
xlim( [ 0 30 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}' , 'R_{3.75}' ,'R_{E;3.75}', 'R_{B;3.75}')

figure(1300)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17) , SNR_dB , R_QPSK_E(:,3), SNR_dB, R_QPSK_B(:,3) ) ;
xlim( [ 0 30 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}' , 'R_{E;4.00}', 'R_{B;4.00}')

figure(1400)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17) , SNR_dB, Capacity( :, 1, 18) , SNR_dB , R_QPSK_E(:,4), SNR_dB, R_QPSK_B(:,4) ) ;
xlim( [ 0 30 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}', 'R_{4.25}' , 'R_{E;4.25}', 'R_{B;4.25}')

figure(1500)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17), SNR_dB, Capacity( :, 1, 19) , SNR_dB , R_QPSK_E(:,5), SNR_dB, R_QPSK_B(:,5) ) ;
xlim( [ 0 30 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}' , 'R_{4.50}', 'R_{E;4.25}', 'R_{B;4.25}')

figure(1600)
plot( SNR_dB, R_QPSK_B(:,1) -  R_QPSK_B(:,3), '--' , SNR_dB, R_QPSK_B_e(:,1) -  R_QPSK_B(:,3), SNR_dB, R_QPSK_B(:,2) -  R_QPSK_B(:,3),'--' , SNR_dB, R_QPSK_B_e(:,2) -  R_QPSK_B(:,3), SNR_dB, R_QPSK_B_e(:,3) -  R_QPSK_B(:,3), '--' , SNR_dB, R_QPSK_B(:,4) -  R_QPSK_B(:,3), SNR_dB, R_QPSK_B_e(:,4) -  R_QPSK_B(:,3),'--' , SNR_dB, R_QPSK_B(:,5) -  R_QPSK_B(:,3), SNR_dB, R_QPSK_B_e(:,5) -  R_QPSK_B(:,3) )
ylim( [-0.05 0.05] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('R_{B;3.50}-R_{B;4.00}','R_{B;3.50}^{opt}-R_{B;4.00}','R_{B;3.75}-R_{B;4.00}', 'R_{B;3.75}^{opt}-R_{B;4.00}', 'R_{B;4.00}^{opt}-R_{B;4.00}','R_{B;4.25}-R_{B;4.00}','R_{B;4.25}^{opt}-R_{B;4.00}','R_{B;4.50}-R_{B;4.00}','R_{B;4.50}^{opt}-R_{B;4.00}' )

figure(1610)
plot( SNR_dB, R_QPSK_E(:,1) -  R_QPSK_E(:,3), '-v' , SNR_dB, R_QPSK_E(:,2) -  R_QPSK_E(:,3),'-^' , SNR_dB, R_QPSK_E(:,4) -  R_QPSK_E(:,3),'--' , SNR_dB, R_QPSK_E(:,5) -  R_QPSK_E(:,3) )
ylim( [-0.05 0.05] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('R_{E;3.50}-R_{E;4.00}','R_{E;3.75}-R_{E;4.00}','R_{E;4.25}-R_{E;4.00}','R_{E;4.50}-R_{E;4.00}' )

figure(1620)
plot(SNR_dB, (log2(1+SNR)), SNR0_dB, R_QPSK, SNR_dB, R_QPSK_E(:,1), '-' , SNR_dB, R_QPSK_E(:,2), '-' , SNR_dB, R_QPSK_E(:,3), SNR_dB, R_QPSK_E(:,4),'--' , SNR_dB, R_QPSK_E(:,5) )
xlim( [0 25 ] )
ylim( [0.1 2.1] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Unconstraint','Constraint','R_{E;3.50}','R_{E;3.75}', 'R_{E;4.00}','R_{E;4.25}','R_{E;4.50}' )

figure(1630)
plot(SNR_dB, (log2(1+SNR)), SNR0_dB, R_QPSK, SNR_dB, R_QPSK_B(:,1), '-' ,SNR_dB, R_QPSK_B_e(:,1), '--' , SNR_dB, R_QPSK_B(:,2), '-' , SNR_dB, R_QPSK_B_e(:,2), '--' , SNR_dB, R_QPSK_B(:,3), SNR_dB, R_QPSK_B_e(:,3),'--' ,SNR_dB, R_QPSK_B(:,4), '-' ,SNR_dB, R_QPSK_B_e(:,4),'--' , SNR_dB, R_QPSK_B(:,5), SNR_dB, R_QPSK_B_e(:,5),'--' )
xlim( [0 25 ] )
ylim( [0.5 2.4] ) ;
grid ;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Unconstraint','Constraint','R_{B;3.50}', 'R_{B;3.50}^{opt}', 'R_{B;3.75}', 'R_{B;3.75}^{opt}', 'R_{B;4.00}','R_{B;4.00}^{opt}','R_{B;4.25}','R_{B;4.25}^{opt}','R_{B;4.50}','R_{B;4.50}^{opt}' )

save(  strcat( 'layered_modualtion_QPSK_QPSK', num2str(now), '.mat' ) ) ;
clock

figure(2000)
plot( SNR0_dB, (log2(1+SNR0))',  SNR0_dB, R_QPSK, SNR_dB, Capacity( :, 1, 17) , SNR_dB , R_QPSK_E(:,3), SNR_dB, R_QPSK_B(:,3) , SNR_dB - 10*log10(0.8), Capacity( :, 1, 17)./2 ,'--' ) ;
xlim( [ 0 30 ])
ylim( [ 0 4.1] )
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Spectral Efficiency (Bit/Symbol)');
legend('Shannon Capacity','R_{QPSK}', 'R_{16QAM}' , 'R_{E;4.00}', 'R_{B;4.00}','Enhancement Layer')

