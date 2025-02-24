# Programmable Signal Generator

## **Interpolator**
An **interpolator** increases the sampling rate of a signal by inserting new samples and applying a **low-pass filter** to reconstruct the missing information. This prevents aliasing and ensures a smooth transition between original and newly generated samples.

### **Multistage interpolation**
Instead of performing direct interpolation by a factor of $16$, the system uses a multistage interpolation process, with multiple **FIR filters** applied sequentially. The goal is to reduce computational complexity, as designing a single FIR filter for $\displaystyle 16\times$ interpolation would require an extremely high-order filter.

The process follows these steps:
- **Upsample the signal** (insert zeros between samples)
- **Filter with a FIR filter**
- **Repeat this process multiple times** until the desired sampling frequency is achieved

This staged approach reduces computational load because each filter operates at a lower interpolation factor, which significantly decreases the required filter order.

<p align="center">
    <img src="images/interpolator.png" alt="Interpolator Block Diagram" style="width: 85%;">
</p>

### **FIR Filters Design**
**Halfband filters** are computationally more efficient than traditional low-pass filters because every second coefficient is zero, reducing the number of required multiplications by nearly $50\%$. Additionally, their fixed stopband at $\displaystyle \frac{f_s}{4}$ makes them ideal for interpolation, enabling efficient multi-stage signal processing with minimal complexity.

#### **FIR Filter Requirements**
To design the FIR filters, we need to determine:
- The **passband frequency** $\displaystyle F_{pass}$:

    The passband frequency in a multi-stage interpolation system depends on the **maximum signal frequency** $\displaystyle f_{max}$, **the interpolation factor** $\displaystyle M$ and **the sampling frequency** $f_s$:
    $$F_{pass} = \frac{f_{max}}{Mf_s}$$
    This ensures that each interpolation stage preserves the desired signal bandwidth while preventing aliasing.
- The **stopband frequency** $\displaystyle F_{stop}$:

    For the halfband filters stopband frequency $\displaystyle F_{stop}$ is fixed:
$$F_{stop} = \frac{f_s}{4}$$
- The **filter order** $\displaystyle N$:

    For optimal FIR filter calculated using **Remez algorithm**, the required filter order $\displaystyle N$ is estimated using **Harris' formula**:
    $$N = \frac{2A_{dB}}{23(0.5 - 2F_{pass})}$$
    where:
    - $\displaystyle A_{dB}$ is a filter attenuation

    Since FIR filters require even orders, we ensure:
    $$N = 2 \times \left\lceil\frac{N}{2}\right\rceil$$

#### **Remez Algorithm**
The **Remez exchange algorithm** (also known as the **Parks-McClellan algorithm**) is an optimal FIR filter design technique based on the **Chebyshev approximation**. It finds the filter coefficients that minimize the maximum deviation from the ideal filter response, resulting in a **minimax optimal filter**.

Given a desired frequency response $\displaystyle H_d(\omega)$, the error function to minimize is:
$$E(\omega) = W(\omega)(H_d(\omega) - H(\omega))$$
where:
- $\displaystyle H(\omega)$ is the actual frequency response of the FIR filter
- $\displaystyle W(\omega)$ is a weighting function

The algorithm iteratively finds an optimal set of filter coefficients that distribute the error uniformly across frequency bands, leading to a flat error response in the passband and stopband. The designed filter must satisfy frequency constraints. The **passband ripple** is checked using:
$$\Delta_{pass} = 10^{-\frac{A_{dB}}{20}}$$
A filter satisfies the passband constraint if:
$$|H(\omega) - 1| \leq 2\Delta_{pass},\ \forall \omega \leq F_{pass}$$
If the filter does not meet these conditions, the order $\displaystyle N$ is increased, and the process repeats.

#### **Halfband Optimization**
A key property of halfband filters is that every second coefficient is zero except for the center tap, which is $\displaystyle 0.5$:

$$h[n] =
\begin{cases}
    0.5, & n = \frac{N}{2} \\
    h[n], & n\ is\ even \\
    0, & n\ is\ odd
\end{cases}
$$

#### Polyphase Decomposition
Using direct convolution for filtering would be computationally expensive. Instead, we use **polyphase decomposition**, which splits the filter into multiple parallel sub-filters that operate at a lower sampling rate.

A discrete-time **FIR filter** with impulse response $\displaystyle h[n]$ and upsampling factor $\displaystyle M$ can be decomposed as:
$$H(z) = \sum_{k=0}^{M - 1} z^{-1}H_k(z^M)$$
where:
- $\displaystyle H_k(z)$ are the **polyphase components**, which operate at the lower input sampling rate.

By restructuring the filter using polyphase decomposition, the system:
- **Processes the input at a lower sampling rate**, reducing computational cost.
- **Performs fewer multiplications per output sample**, making the system more efficient.

<p align="center">
    <img src="images/polyphase.png" alt="Polyphase Decomposition" style="width: 85%;">
</p>

## **Complex Mixer**
A **complex mixer** is a digital signal processing block that shifts the frequency of a complex signal $\displaystyle S[n] = I[n] + jQ[n]$ by multiplying it with a complex exponential:
$$e^{jn\Omega_0} = \cos(n\Omega_0) + j\sin(n\Omega_0)$$
This is implemented using an **NCO + CORDIC** system:
- **The NCO generates phase values** $\displaystyle \theta[n] = e^{jn\Omega_0}$
- **CORDIC computes** $\displaystyle \cos(\theta[n])$ and $j\displaystyle \sin(\theta[n])$
- The original signal is then multiplied by the computed sine and cosine:
$$\hat{I}[n] = I[n]\cos(\theta[n]) - Q[n]\sin(\theta[n])$$
$$\hat{Q}[n] = I[n]\sin(\theta[n]) + Q[n]\cos(\theta[n])$$

<p align="center">
    <img src="images/complex_mixer.png" alt="Polyphase Decomposition" style="width: 85%;">
</p>

### **Numerically Controlled Oscillator (NCO)**
**A Numerically Controlled Oscillator (NCO)** is a digital frequency synthesizer that generates sinusoidal waveforms with high frequency resolution. An NCO operates based on a phase accumulator which increments the phase value by a fixed step size $\displaystyle W$ each clock cycle. This accumulated phase is then wrapped around at $\displaystyle 2\pi$, forming a periodic waveform. Since the **CORDIC algorithm** which is used for computing sine/cosine takes multiple iterations, the NCO introduces a latency (delay).

**NCO** has two key parameter:
- **Phase Accumulator Word Length** ($\displaystyle L$) is the frequency resolution of the **NCO**. For the smallest frequency step the **NCO** must resolve (**frequency resolution** $\displaystyle f_{res}$), minimum phase word length $\displaystyle L$ is:
$$L = \left\lceil\log_2(\frac{f_{res}}{f_s})\right\rceil$$

- **Phase Increment** ($\displaystyle W$) is an update of the phase at each clock cycle and it is determined by the **output frequency** $f_{out}$:
$$\theta[n] = \theta[n-1] + W$$
$$W = \frac{f_{out}2^L}{f_s}$$

### **CORDIC algorithm**
The **CORDIC (COordinate Rotation DIgital Computer) algorithm** is an iterative algorithm used to efficiently compute sine and cosine values using only bit shifts and additions. It rotates a vector iteratively until it aligns with the desired phase $\displaystyle \theta[n]$ which comes from the **NCO** phase accumulator.

After fixed $\displaystyle N$ iterations, the algorithm converges.

#### Initialization
At the beginning of the **CORDIC** algorithm inputs are set to the initial values:
- $\displaystyle v_0 = x_0 + jy_0,\ (x_0,\ y_0) = (1,\ 0)$ - Initial complex unit vector
- $\displaystyle z_0 = \theta[n]$ - Initial phase error from the **NCO**
- $\displaystyle f_i = 1 + 2^{-i},\ i = \overline{1,\ N}$ - Precomputed rotation factors
- $\alpha_i = arg(f_i) = \arctan(2^{-i})$ - Precomputed rotation angles

#### Iterative Process
At each iteration $\displaystyle i$, the algorithm:
- **Rotate the vector forward or backward, based on the current residual phase $z_i$:**
$$v_{i+1} =
\begin{cases}
    v_i f_i, & 0 < z_i < \frac{W}{2} \rightarrow forward \\
    v_i f_{i}^{*}, & otherwise \rightarrow backward
\end{cases}
$$
- **Adjust the phase error $z_i$ to ensure that the angle wraps around within the valid range**:
$$z_{i+1} =
\begin{cases}
    mod((z_i - \alpha_i + W),\ W), & 0 < z_i < \frac{W}{2} \\
    mod((z_i + \alpha_i + W),\ W), & otherwise
\end{cases}
$$

#### Convergence
Each iteration reduces the **phase error** $\displaystyle z_i$, bringing it closer to $\displaystyle 0$. After $\displaystyle N$ iterations the input vector has rotated by desired phase $\displaystyle \theta[n]$:
$$v_N = \cos(\theta[n]) + j\sin(\theta[n]) = e^{jn\Omega_0}$$

## **$\displaystyle sinc(x)$ Compensation Filter**
In a **digital-to-analog converter (DAC)**, when the output is held constant between sample points (**Zero-Order Hold**), the frequency response of the system is affected by a **sinc function** distortion:
$$H(f) = \frac{\sin(\pi fT_s)}{\pi fT_s}$$
This low-pass filtering effect causes amplitude attenuation at higher frequencies, especially near the Nyquist frequency. Sinc compensation is used to correct this distortion by pre-emphasizing the higher frequencies before conversion. To counteract this effect, we apply a pre-filter that boosts the high frequencies by the **inverse sinc function** before sending the signal to the **DAC**:
$$H_{compensation}(f) = \frac{\pi fT_s}{\sin(\pi fT_s)}$$
However, because the inverse sinc function is not perfectly realizable, we need to approximate it using a finite-order digital filter. The iterative design process starts with a low-order filter, evaluates its frequency response against the ideal inverse sinc function, and incrementally increases the filter order until the error falls below an error margin.

## **Digital-to-Analog Converter**
A **Digital-to-Analog Converter (DAC)** converts discrete digital samples into a continuous analog signal by holding each sample value for a time $\displaystyle T_s$ (sampling period). The continuous output can be represented as:
$$x(t) = \sum_{n} x[n]h(t - nT_s)$$
where:
- $\displaystyle x[n]$ are the digital samples
- $\displaystyle h(t)$ is the **DAC's** impulse response, which determines how each sample is translated into an analog signal

The output spectrum is affected by the **zero-order hold** effect, introducing a sinc-shaped frequency response:
$$H_{DAC}(f) = \frac{\sin(\pi fT_s)}{\pi fT_s}$$
which attenuates high-frequency components.

### **Non-Return-to-Zero (NRZ) Mode**
In **NRZ** mode, each digital sample is held constant for the full sampling period $\displaystyle T_s$, meaning the **DAC** output is a staircase-like waveform. The impulse response of the **DAC** in **NRZ** mode can be modeled as:
$$h_{NRZ}(t) =
\begin{cases}
    1, & 0 \leq t < T_s \\
    0, & otherwise
\end{cases}
$$
while the frequency response of the **DAC** in **NRZ** mode can be modeled as:
$$H_{NRZ}(f) = T_s sinc(fT_s)$$
This means:
- **The first zero occurs at** $\displaystyle f_s$
- **The main lobe extends from DC to** $\displaystyle \frac{f_s}{2}$
- **Attenuation is the lowest, hence signal is usable in the first Nyquist zone** $\displaystyle f\in [0,\ \frac{f_s}{2}]$

### **Bipolar Non-Return-to-Zero (RF) Mode**
**RF** Mode is a **DAC** output mode where each sample is held at $\displaystyle +V$ for $\displaystyle \frac{T_s}{2}$, and at $\displaystyle -V$ for the second half. This results in a balanced, alternating waveform that improves spectral properties by reducing even-order harmonics. The impulse response of the **DAC** in **RF** mode can be modeled as:
$$h_{RF}(t) =
\begin{cases}
    1, & 0 \leq t < \frac{T_s}{2} \\
    -1, & \frac{T_s}{2} \leq t < T_s \\
    0, & otherwise
\end{cases}
$$
while the frequency response of the **DAC** in **RF** mode can be modeled as:
$$H_{RF}(f) = T_s sinc(\frac{fT_s}{2})\sin(\frac{fT_s}{2})e^{-j(\pi f T_s - \frac{\pi}{2})}$$
This means:
- **The first zero occurs at** $\displaystyle 0$
- **The main lobe extends from DC to** $\displaystyle f_s$
- **Attenuation is the lowest, hence signal is usable in the second Nyquist zone** $\displaystyle f\in [\frac{f_s}{2},\ f_s]$

## Conclusion
By understanding the multi-stage **interpolation** approach, we see how **polyphase filtering** reduces computational complexity while achieving a smooth, high-resolution signal. The **complex mixer**, implemented using an **NCO** and **CORDIC**, enables precise frequency shifts while maintaining signal integrity. The effects of **DAC** operation in **NRZ** and **RF** modes illustrate how spectral shaping influences signal reconstruction and how **sinc compensation** corrects amplitude roll-off.

For details on testing and validation methods used in this system, see [Testing & Validation](TESTING.md).

For an overview of the system’s structure, see the [System Architecture](ARCHITECTURE.md).
