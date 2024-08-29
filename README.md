# Digital-Modulation-Channel-Coding-Simulation-with-MATLAB

I went through three steps during the development process to make the structure more obvious. 

The first step is to create a signal with a fixed Eb/No value, do not deploy the filter, and debug the vectors and signals in it.

The second step is to add a function to draw BER curve for the first step of the code, at this time the Eb/No value is an array, which faces a lot of problems, such as mod_bits variable 1 into a multidimensional matrix, which is very troublesome when calculating.

The third step is to add filters and compare BER curves with and without filters

Finally, we get two file to compare our different BER Curve.
