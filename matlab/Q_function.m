x_dB = -10: 0.1 : 20

x = 10.^( x_dB./10 )


figure(100)
semilogy(x_dB, erfc(sqrt(x)./2)./2, x_dB, exp(-x./4)./sqrt(x.*pi), x_dB, exp(-x./4)./sqrt(x.*pi) - erfc(sqrt(x)./2)./2 )
grid;
xlabel('\gamma (dB)')
ylabel('QPSK SER')
legend('Q-Function','Approximation','Approximation Error')
ylim([0 5 ])

figure(200)
semilogy(x_dB, erfc(sqrt(x)./2)./2, x_dB, exp(-x./4)./sqrt(x.*pi))
grid;
xlabel('\gamma (dB)')
ylabel('QPSK SER')
legend('Q-Function','Approximation')
ylim([0 5 ])

