rng(100)

n=1e6; 

a = 0;
b = 1;

m = 1/2;
P = 1/12; 

x = m + sqrt(P)*randn(n,1);
y = a + (b - a)*rand(n,1);

gaussianMean = mean(x); 
gaussianStd = std(x);
uniformMean = mean(y);
uniformStd = std(y);

fprintf("Gaussian Mean %.4e \n" + ...
        "Gaussian Std  %.4e \n" + ...
        "Uniform Mean  %.4e \n" + ...
        "Uniform Std   %.4e\n", [gaussianMean, gaussianStd, uniformMean, uniformStd])

%%
figure(name="Uniform vs Gaussian Distributions")

tiledlayout
nexttile
histogram(x, 100)
title("Gaussian Distribution")
nexttile
histogram(y,100)
title("Uniform Distribution")
