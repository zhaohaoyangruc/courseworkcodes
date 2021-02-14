% simulation and coefficients
b10 = sampleregress(10,10000);
b25 = sampleregress(25,10000);
b50 = sampleregress(50,10000);
b100 = sampleregress(100,10000);

% Plot the density
normplot(b10)
normplot(b25)
normplot(b50)
normplot(b100)