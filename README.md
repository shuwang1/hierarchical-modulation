# Enhance and Optimize Hierarchical Modulations

Hierarchical modulation, also called layered modulation, is one of the techniques for multiplexing and modulating multiple data streams into one single symbol stream, where the base-layer symbols and enhancement-layer symbols are
synchronously superimposed together before being transmitted. When hierarchical modulation is employed, users with good reception and advanced receiver can demodulate more than one layer of data streams. For a user with conventional
receiver or poor reception, it may be able to demodulate the data streams embedded in low layer(s), e.g, the base layer only. From an information-theoretical perspective, hierarchical modulation is one practical implementation of superposition
precoding, which can help achieve the maximum sum rate of Gaussian broadcast channel with employing interference cancellation by receivers. From a network operation perspective,
a network operator can seamlessly target users with different services or QoS’s and this technique. 

However, traditional hierarchical modulation suffers from inter-layer interference (ILI) so that the achievable rates by low-layer data streams, e.g.
the base-layer data stream, can be dented by the interference from high-layer signal(s). For example, for the hierarchically modulated two-layer symbols comprising a 16QAM base layer and a QPSK enhancement layer, the base-layer throughput
loss can be up to about 1.5bits/symbol with the total receive signal-to-noise ratio (SNR) of about 23dB. This means, due to ILI, there is about 1.5/4 = 37.5% loss of the base-layer
achievable throughput with 23dB SNR. The demodulation error rate of either the base-layer and enhancement-layer symbols increases too. 

From an implementation point-view, it is known that the severe amplitude and phase fluctuations of wireless channels can significantly degrade the receiver demodulation performance since the demodulator must scale the
received signal so that the result signals is within the dynamic range of the followed analog-to-digital convertor (ADC) or, more generally, the receiver processing region, mostly with
automatic gain control (AGC). Even though pilots may be available for assisting the receiver channel estimation and equalization, there are channel estimation errors, especially
when the channel coherent time is short. If the channel is estimated in errors, it can lead to improperly compensated signals and incorrect demodulation even in the absence of
noise. 

Furthermore, multicarrier transmission, e.g. orthogonal frequency-division multiplexing (OFDM), is widely used for BCMCS as well as next generation wireless systems, due to
its high diversity gain and high spectral efficiency with simple receiver design. However, the advantages of OFDM, specially when it is modulated by high-order signal constellations,
are counter-balanced by the high peak-to-average-power ratio (PAPR) issue. High PAPR of modulated signals can significantly reduce the average output power of the high-power amplifier (HPA) at the transmitter due to more back-offs. It
also increases the receiver demodulation and decoding errors and therefore limits the throughput of whole transceiver chain. Therefore it is important to understand and optimize regular
hierarchical modulations for the best achievable performance.
