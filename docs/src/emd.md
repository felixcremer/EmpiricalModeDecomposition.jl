# Empirical Mode Decomposition
The empirical mode decomposition was developed by Huang in [Huang 1988].
It decomposes the signal into several Intrinsic Mode Functions (IMF).
But unlike the Fourier or wavelet transformation the IMFs are computed adaptively from the data.
The algorithms works as follows:

1. Detect the local extrema of the original dataset $x(t)$
2. Compute a spline $s_{min},s_{max}$ through the local minima and maxima respectively
3. Compute the mean of $s_{min}$ and $s_{max}$ and subtract it from $x(t)$.  

    ```math
    h_{1,k}(t)=x(t)-\frac{s_{min}-s_{max}}{2}
    ```  
4. Repeat the previous steps with $h_{1,k}$ as $x(t)$ until $h_{1,k}$ suffice a stopping criteria. $h_{1,k}$ is the first IMF $c_1(t)$ of $x(t)$
5. Repeat the steps 1 till 4 with the difference between $x(t)$ and $c_1(t)$ until there is no IMF to extract

This procedure decomposes the dataset $x(t)$ into N IMFs and one residual:
```math
x(t)=\sum_{i=1}^N c_i(t)+r(t)
```
