X = random("Uniform", -0.5, 0.5, [1, 20]);
convolved = conv(X, X);
convolved2 = conv(convolved, X);
convolved3 = conv(convolved2, X);
convolved4 = conv(convolved3, X);
convolved5 = conv(convolved4, X);

mu = 0;
n = length(convolved);
sigma = sqrt(n) * std(X);
x1 = linspace(-4*sigma, 4*sigma, 40);
g1 = normpdf(x1, mu, sigma);

n = length(convolved2);
sigma = sqrt(n) * std(X);
x2 = linspace(-4*sigma, 4*sigma, 40);
g2 = normpdf(x2, mu, sigma);

n = length(convolved3);
sigma = sqrt(n) * std(X);
x3 = linspace(-4*sigma, 4*sigma, 40);
g3 = normpdf(x3, mu, sigma);

n = length(convolved4);
sigma = sqrt(n) * std(X);
x4 = linspace(-4*sigma, 4*sigma, 40);
g4 = normpdf(x4, mu, sigma);

n = length(convolved5);
sigma = sqrt(n) * std(X);
x5 = linspace(-4*sigma, 4*sigma, 40);
g5 = normpdf(x5, mu, sigma);


goldenrod = [0.85 0.65 0.13];
[ah, ha, H] = krijetem(subnum(3, 2));


axes(ah(1))
hold on
histogram(X, 10, 'FaceColor', goldenrod,'EdgeColor','none', "Normalization","pdf")
title("Uniform Distribution")


axes(ah(2))
hold on
plot(x1, g1)
histogram(convolved,10, 'FaceColor', goldenrod,'EdgeColor','none', "Normalization","pdf")
hold off

axes(ah(3))
hold on
plot(x2, g2)
histogram(convolved2,10, 'FaceColor', goldenrod,'EdgeColor','none', "Normalization","pdf")
hold off


axes(ah(4))
hold on
plot(x3, g3)
histogram(convolved3, 10,'FaceColor', goldenrod,'EdgeColor','none', "Normalization","pdf")
hold off

axes(ah(5))
hold on
plot(x4, g4)
histogram(convolved4,10, 'FaceColor', goldenrod,'EdgeColor','none', "Normalization","pdf")
hold off

axes(ah(6))
hold on
plot(x5, g5)
histogram(convolved5,10, 'FaceColor', goldenrod,'EdgeColor','none', "Normalization","pdf")
hold off