library(pcapack)

set.seed(1234)

n <- 10
p <- 3
x <- matrix(rnorm(n*p), n, p)



truth <- 
structure(c(-2.68584033968119, -2.88415245190642, -2.91587343172424, 
-2.32829867631603, -1.41937447005234, -1.64398662890595, -3.1349175063739, 
-2.29805740265251, -2.94920342990902, 0, -1.51462777283308, -1.05660447643105, 
-0.346417436049099, -1.84359220411886, 0.137199864132763, -0.757113114142974, 
-0.919181604623999, -1.78997626102861, -1.34742335871751, 0, 
-0.327506007657822, -1.94509617914423, -2.38979056971792, 0.880701961829639, 
-1.37344315191914, -2.30583120696983, -0.524840458425345, -1.64610751216655, 
-1.01123725567031, 0), .Dim = c(10L, 3L))


truth

cma(x)
