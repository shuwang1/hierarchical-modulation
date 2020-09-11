SNR = 1:2:20

C0 = log2(1+SNR)

C2 = 2*log2(1+SNR/2)

C4 = 4*log2(1+SNR/4)

C8 = 8*log2(1+SNR/8)

C16 = 16*log2(1+SNR/16)

figure(10)
plot(SNR, C0, SNR, C2, '-.',SNR, C4,'--',SNR, C8,'-o',SNR, C16,'-v')
grid
xlabel('SNR');
ylabel('Spectral Efficiency (bps/Hz)')
legend('N=1','N=2','N=4','N=8','N=16')