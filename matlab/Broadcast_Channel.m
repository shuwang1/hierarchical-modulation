clear all ;
close all ;
P = 150 ;
beta = 0.0001:0.001: 1

g1 = 1 ;
g2 = 0.25;

sigma = 1 ;


for m = 1 : length( beta )
    
    P1 = P * beta( m ) ;
    P2 = P * ( 1 - beta( m ) ) ;
    
    
    R1(1, m) = log2( 1 + P1*g1/sigma ) ;
    R2(1, m) = log2( 1 + P2*g2/( sigma + P1*g2) ) ;        

    R1(2, m) = beta( m )*log2( 1 + P*g1/(sigma) ) ;
    R2(2, m) = ( 1 - beta( m ) )*log2( 1 + P*g2/sigma ) ;        
        
    R1(3, m) = beta(m)*log2( 1 + P1*g1/beta(m)/sigma ) ;
    R2(3, m) = (1-beta(m))*log2( 1 + P2*g2/(1-beta(m))/sigma ) ;     
    
    R1(4, m) = 0.1*log2( 1 + P1*g1/sigma/0.1 ) ;
    R2(4, m) = 0.9*log2( 1 + P2*g2/sigma/0.9 ) ;

    R1(5, m) = 0.2*log2( 1 + P1*g1/sigma/0.2 ) ;
    R2(5, m) = 0.8*log2( 1 + P2*g2/sigma/0.8 ) ;
    
    R1(6, m) = 0.4*log2( 1 + P1*g1/sigma/0.4 ) ;
    R2(6, m) = 0.6*log2( 1 + P2*g2/sigma/0.6 ) ;

    R1(7, m) = 0.6*log2( 1 + P1*g1/sigma/0.6 ) ;
    R2(7, m) = 0.4*log2( 1 + P2*g2/sigma/0.4 ) ;

    R1(8, m) = 0.8*log2( 1 + P1*g1/sigma/0.8 ) ;
    R2(8, m) = 0.2*log2( 1 + P2*g2/sigma/0.2 ) ;
    
    R1(9, m) = 0.9*log2( 1 + P1*g1/sigma/0.9 ) ;
    R2(9, m) = 0.1*log2( 1 + P2*g2/sigma/0.1 ) ;

    rE = P1*g1/sigma ;
    rB = P2*g1/sigma ;
    
    tempE(m, :) = roots( [rE 1 -1 ] ) ;
  
    muE( m )  = 0 ;
    for n = 1 : 2
        if ( tempE(m, n) > 0 ) && ( tempE(m, n) > muE( m ) ) && ( tempE(m, n) <= 1)           
            muE( m ) = tempE(m, n) ;
        end
    end
    R1(10, m) = log2( 1 + rE*muE( m ) ) ;

    rE = P1*g1/sigma ;
    tempE2(m, :) = roots( [ rE*rB  (rB*rE+rB+rE) 1 -1 ] ) ;
    muE( m )  = 0 ;
    for n = 1 : 3
        if ( tempE2(m, n) > 0 ) && ( tempE2(m, n) > muE( m ) ) && ( tempE2(m, n) <= 1)           
            muE( m ) = tempE2(m, n) ;
        end
    end
    R1(11, m) = log2( 1 + rE*muE( m ) ) ;
    
    %%-------------------------------------------%%    
    rE = P1*g1/sigma ;
    rB = P2*g1/sigma ;
    
    R1(12, m) = C_MMSE( rE, 5 ) ;
    R1(13, m) = C_OPT( rE, 5 ) ;
    
    rE = P1*g2/sigma ;
    rB = P2*g2/sigma ;
    %rB = P2*g2 / ( sigma + P1*g2 ) ;
    
    %tempB(m, :) = roots( [rB 1 -1 ] ) ;
    tempB(m, :) = roots( [ rE*rB  (rB*rE+rB+rE) 1 -1 ] ) ;
    
    muB( m )  = 0 ;
    for n = 1 : 3
        if ( tempB(m, n) > 0 ) && ( tempB(m, n) > muB( m ) ) && ( tempB(m, n) <= 1)           
            muB( m ) = tempB(m, n) ;
        end
    end
    
    R2(10, m) = log2( 1 + rB*muB( m ) ) ;   
    rB = P2*g2/(sigma+P1*g2) ;
    tempB2(m, :) = roots( [rB 1 -1 ] ) ;
    muB( m )  = 0 ;
    for n = 1 : 2
        if ( tempB2(m, n) > 0 ) && ( tempB2(m, n) > muB( m ) ) && ( tempB2(m, n) <= 1)           
            muB( m ) = tempB2(m, n) ;
        end
    end
    
    R2(11, m) = log2( 1 + rB*muB( m ) ) ;   
    
    rE = P1*g2/sigma ;
    rB = P2*g2/(1+rE) ;
    
    R2(12, m) = C_MMSE( rB, 5 ) ;
    R2(13, m) = C_OPT( rB, 5 ) ;

    
end


figure(100)
plot( R2(2,:), R1(2,:) )
grid;
xlabel('R_2: Achievable Rate at Coverage Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
legend('Time-Division Multiplexing')


figure(110)
plot( R2(2,:), R1(2,:), R2(6,:), R1(6,:),'-.' )
grid;
xlabel('R_2: Achievable Rate at Coverage Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
legend('Time-Division Multiplexing', 'Frequency-Division Multiplexing, B_1/B_2=2/3' )

%return


figure(115)
plot( R2(2,:), R1(2,:), R2(5,:), R1(5,:) , R2(6,:), R1(6,:),'-.', R2(7,:), R1(7,:),':', R2(8,:), R1(8,:),'--' )
grid;
xlabel('R_2: Achievable Rate at Coverage Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
legend('Time-Division Multiplexing', 'Frequency-Division Multiplexing, B_1/B_2=1/4', 'Frequency-Division Multiplexing, B_1/B_2=2/3' , 'Frequency-Division Multiplexing, B_1/B_2=3/2', 'Frequency-Division Multiplexing, B_1/B_2=4/1')



figure(120)
plot( R2(2,:), R1(2,:), '--', R2(4,:), R1(4,:) , R2(5,:), R1(5,:),'-.', R2(6,:), R1(6,:),':')
grid;
xlabel('R_2: Achievable Rate at Coverage Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
legend('Time-Division Multiplexing','Frequency-Division Multiplexing, B_1/B_2=1/4', 'Frequency-Division Multiplexing, B_1/B_2=2/3' , 'Frequency-Division Multiplexing, B_1/B_2=4')


figure(130)
plot( R2(1,:), R1(1,:), R2(2,:), R1(2,:), '--' )
grid;
xlabel('R_2: Achievable Rate at Cell Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend('Superposition Precoding','Time-Division Multiplexing')

figure(132)
plot( R2(1,:), R1(1,:), R2(2,:), R1(2,:), '--', R2(10,2:end ), R1(10, 2:end ), '-.', R2(12,2:end ), R1(12, 2:end ), '--' , R2(13,2:end ), R1(13, 2:end ), ':' )
grid;
xlabel('R_2: Achievable Rate at Cell Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend('Superposition Precoding','Time-Division Multiplexing', 'Quasi-Orthogonal CDM', 'SIC with MMSE Detection', 'SIC with Optimal Detection')

return

figure(134)
plot( R2(1,:), R2(1,:) + R1(1,:), R2(2,:), R2(2,:)+ R1(2,:), '--', R2(10, 2:end ), R2(10,2:end )+ R1(10, 2:end ), '-.' )
grid;
xlabel('R_2: Achievable Rate at Cell Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend('Superposition Precoding','Time-Division Multiplexing', 'Quasi-Orthogonal CDM')


figure(136)
plot( P.*(1-beta).*g2./sigma, R2(10,:), P.*(1-beta).*g2./sigma, R2(1,:), P.*(1-beta).*g2./sigma, R2(11,:) )
grid;
xlabel('R_2: Achievable Rate at Cell Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend( 'Quasi-Orthogonal CDM', 'Superposition Precoding')

figure(138)
plot( P.*beta.*g1./sigma, R1(10,:), P.*beta.*g1./sigma, R1(1,:), P.*beta.*g1./sigma, R1(11,:) )
grid;
xlabel('R_2: Achievable Rate at Cell Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend( 'Quasi-Orthogonal CDM', 'Superposition Precoding', 'CDM2')

return

figure(140)
plot( R2(1,:), R2(1,:), R2(1,:), R1(1,:), R2(2,:), R1(2,:), '--', R2(4,:), R1(4,:) , R2(5,:), R1(5,:),'-.', R2(6,:), R1(6,:),':' )
grid;
xlabel('R_2: Achievable Rate at Cell Edge at d_2 (bps/Hz) ')
ylabel('R_1-R_2: Achievable Rate Gain inside Cell at d_1 (bps/Hz) ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend('Superposition Precoding','Time-Division Multiplexing','Frequency-Division Multiplexing, B_1/B_2=1/4', 'Frequency-Division Multiplexing, B_1/B_2=2/3' , 'Frequency-Division Multiplexing, B_1/B_2=4')


figure(145)
plot( [R2(1,1), R2(1,:)], [ 0, R1(1,:) + R2(1,:)], [R2(2,1), R2(2,:)], [0, R1(2,:) + R2(2,:)], '--',[R2(4,1), R2(4,:) ], [0, R1(4,:) + R2(4,:)] , [R2(5,1), R2(5,:)], [0, R1(5,:) + R2(5,:)],'-.',[R2(6,1), R2(6,:)],[ 0, R1(6,:) + R2(6,:)],':' )
grid;
xlabel('R_2 (bps/Hz): Achievable Rate at the Cell Edge')
ylabel('R_1 (bps/Hz): Achievable Rate in the Inner Coverage ')
%legend('Superposition Precoding','Time/Freq.-Division Multiplexing','Time-Division Multiplexing')
legend('Superposition Precoding','Time-Division Multiplexing','Frequency-Division Multiplexing, B_1/B_2=1/4', 'Frequency-Division Multiplexing, B_1/B_2=2/3' , 'Frequency-Division Multiplexing, B_1/B_2=4')


figure(200)
%plot(1:length(muE), muE, 1:length(temp(:, 1)), 1.0./temp(:, 1), 1:length(temp(:, 2)), 1.0./temp(:, 2), 1:length(temp(:, 3)), 1.0./temp(:, 3))
plot(1:length(muE), muE, 1:length(muB), muB )
grid;
ylim( [0  1] )
return

tau = 0:0.01:1
R1 = 11.5
R2 = 4 
R = (1-tau).*R1 + (2*tau-1).*R2

figure(150)
plot(tau, R)

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log-Normal Shadowing, 3GPP2 EMD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d0 = 35

d = 35 : 1 : 1.25e3 ;

PL0_dB = 28.6 + 35*log10(d) ;
PL1_dB = PL0_dB  + 8.9*randn( size(d) ) ;

figure(200)
plot( d, -PL0_dB, '--', d, -PL1_dB , 'LineWidth', 6 ) ;
xlim([35 1.25e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( 'PL (dB), Log-Normal Shadowing' ) ;
grid on;
legend( 'Modified Hata Urban Propation Model at 1.9GH, COST231', 'Log-Normal Shadowing' )


d = 35 : 1 : 3e3 ;

SNR0 =  2.0e15 ;

SNR_dB = log10( SNR0 )*10 - 35*log10(d) - 28.6 ; % - 8.9*randn( size(d) ) ;

SNR = 10.^( SNR_dB./10 ) ;

mu = log2( 1 + SNR ) ;

figure(300)
loglog( d, mu) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu (bps/Hz), Achievable Spectral Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )

figure(310)
semilogx( d, mu.*d) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu*d (bps/Hz), Achievable Transport Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )

figure(320)
semilogx( d, mu.*(d.^2)) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu*d^2 (bps/Hz), Achievable Coverage Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )


figure(330)
loglog( d, mu, d, mu.*(d), d, mu.*(d.^2) ) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu, \mu*d, \mu*d^2 -- Achievable Capacity' ) ;
grid on;
legend( 'Channel Capacity', 'Transport Capacity', 'Coverage Capacity' )


beta = 0.1 ;
SNR1 = beta*SNR0 ;
SNR1_dB = SNR_dB + log10(beta)*10 ;
mu1 = log2( 1 + 10.^(SNR1_dB./10) )

SNR2 = (1-beta)*SNR0 ;
SNR2_dB = SNR_dB + log10(1-beta)*10 ;
mu2 = log2( 1 + 10.^(SNR2_dB./10) )


figure(400)
loglog( d, mu, d, mu*(1-beta), d, mu + mu1) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu (bps/Hz), Achievable Spectral Efficiency' ) ;
grid on;
legend( '', 'Path loss factor \alpha=3.5' )

return

figure(410)
semilogx( d, mu.*d, d, (mu+mu1)*d  ) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu*d (bps/Hz), Achievable Transport Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )

figure(320)
semilogx( d, mu.*(d.^2)) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu*d^2 (bps/Hz), Achievable Coverage Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )


figure(330)
loglog( d, mu, d, mu.*(d), d, mu.*(d.^2) ) ;
xlim([35 3e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu, \mu*d, \mu*d^2 -- Achievable Capacity' ) ;
grid on;
legend( 'Channel Capacity', 'Transport Capacity', 'Coverage Capacity' )



return

d = 35 : 1 : 4e3

SNR0 =  1.0e14 ;

SNR_dB = log10( SNR0 )*10 - 35*log10(d) - 28.6 ; %- 8.9*randn( size(d) ) ;

SNR = 10.^( ( SNR_dB )./10 ) ;

mu = log2( 1 + (sqrt(SNR) + sqrt(SNR(end:-1:1)).*6).^2 ) ;

figure(400)
loglog( d, mu) ;
xlim([35 1.25e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu (bps/Hz), Achievable Spectral Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )

figure(410)
loglog( d, mu.*(d)) ;
xlim([35 1.25e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu (bps/Hz), Achievable Rate-Transport Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )


figure(420)
loglog( d, mu.*(d.^2)) ;
xlim([35 1.25e3])
xlabel( 'd (m), the separate between BS and MS' ) ;
ylabel( '\mu (bps/Hz), Achievable Rate-Coverage Efficiency' ) ;
grid on;
legend( 'Path loss factor \alpha=3.5' )


SNR = 0:0.001:1

figure(500)
plot( SNR, 2*log2(1+SNR), SNR, log2(1+2*SNR)) ;

