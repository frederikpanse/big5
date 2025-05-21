library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# read the data I downloaded from the Open-Source Psychometrics Project
dat <- read.csv("data.csv", header = TRUE, sep = "\t")
# remove all missing values
dat[dat == 0] <- NA
dat$gender[dat$gender == 3] <- NA
dat <- na.omit(dat)
# correct the score for reverse items
rev_items <- c( "E2", "E4", "E6", "E8", "E10", "N2", "N4", "A1", "A3", "A5", 
	       "A7", "C2", "C4", "C6", "C8", "O2", "O4", "O6")
dat[rev_items] <- 6 - dat[rev_items]

# get the score for each personality dimension for each subject
dat$E <- dat[ ,  8:17] |> rowMeans()  # Extraversion
dat$N <- dat[ , 18:27] |> rowMeans()  # Neuroticism
dat$A <- dat[ , 28:37] |> rowMeans()  # Agreeableness
dat$C <- dat[ , 38:47] |> rowMeans()  # Conscientiousness
dat$O <- dat[ , 48:57] |> rowMeans()  # Openess (Intellect/Imagination)

# compare densities in Extraversion for male (1) and female (2) participants
# Extraversion
dat[dat$gender == 1,]$E |> density() |> plot(col = "blue", 
					     xlab = "Extraversion score", 
					     main = "")
dat[dat$gender == 2,]$E |> density() |> lines(col = "pink")
# Neuroticism
dat[dat$gender == 2,]$N |> density() |> plot(col = "pink", 
					     xlab = "Neuroticism score", 
					     main = "")
dat[dat$gender == 1,]$N |> density() |> lines(col = "blue")

set.seed(1909)  # for replicability

##### Building the model #####
## Model 1
# start with one question on Extraversion, no predictor, and only 100 subjects
dat1 <- dat[sample(nrow(dat), 100), c("gender", "E1")]
table(dat1$E1) |> barplot()
# translate frequencies to log-cumulative-odds
pr_k <- table(dat1$E1) / nrow(dat1)
cum_pr_k <- cumsum(pr_k)
plot(1:5 , cum_pr_k, type="b", xlab="response", ylab="cumulative proportion", 
     ylim=c(0,1))
round(lco <- qlogis(cum_pr_k) , 2)  # apply logit
plot(1:5 , lco, type="b", xlab="response", ylab="log-cumulative-odds")

datlist <- list(N = nrow(dat1), 
		K = 5, 
		y = dat1$E1)

s <- "
data {
  int<lower=2>                   K;
  int<lower=0>                   N;
  array[N] int<lower=1, upper=K> y;
}
parameters {
  ordered[K - 1] c;
}
model {
  for (n in 1:N) {
    y[n] ~ ordered_probit(0, c);  // no x variable
  }
}
"
m1 <- stan(model_code = s, data = datlist)
print(m1, probs = c(.025, .975))

cum_pr_k_hat <- extract(m1)$c |> colMeans() |> plogis()

# function to turn any cumulative probabilities into normal probabilities
cum_to_pr <- function(cum_prob) {
  pr_k_hat <- matrix(nrow = nrow(cum_prob), ncol = 5)
  pr_k_hat[, 1] <- cum_prob[, 1]
  pr_k_hat[, 2] <- cum_prob[, 2] - cum_prob[, 1]
  pr_k_hat[, 3] <- cum_prob[, 3] - cum_prob[, 2]
  pr_k_hat[, 4] <- cum_prob[, 4] - cum_prob[, 3]
  pr_k_hat[, 5] <- 1 - cum_prob[, 4]
  pr_k_hat
}

pr_k |> plot(ylim = c(0, .5))
points(pr_k_hat)

# get some uncertainty
sample_cum_prob <- extract(m1)$c |> plogis() 
sample_prob <- apply(sample_cum_prob, 1, cum_to_pr)
pr_k |> plot(ylim = c(0, .5))
for (i in seq_along(idx)) {
  points(sample_prob[ , idx[i]], col = "lightgrey")
}

## Model 2
# all 10 Extraversion questions, no predictor, complete pooling
select <- c("gender", paste0("E", 1:10))
dat2_wide <- dat[sample(nrow(dat), 100), select]
dat2 <- reshape(dat2, varying = paste0("E", 1:10), v.names = "score", 
        timevar = "q", idvar = "id", direction = "long")

datlist <- list(N = nrow(dat2), 
		K = 5, 
		y = dat2$score)

s <- "
data {
  int<lower=2>                   K;
  int<lower=0>                   N;
  array[N] int<lower=1, upper=K> y;
}
parameters {
  ordered[K - 1] c;
}
model {
  for (n in 1:N) {
    y[n] ~ ordered_logistic(0, c);  // no x variable
}
"
m2 <- stan(model_code = s, data = datlist)
print(m2, probs = c(.025, .975))

## Model 3
# all 10 extraversion questions, no predictor, no pooling with regard to scale
datlist <- list(N = nrow(dat2), 
		K = 5, 
		nq = length(unique(dat2$q)),
		q = dat2$q,
		y = dat2$score)

s <- "
data {
  int<lower=2>                         K;
  int<lower=0>                         N;
  int<lower=0>                         nq;
  array[N] int<lower=1, upper=nq>      q;
  array[N] int<lower=1, upper=K>       y;
}
parameters {
  ordered[K - 1] c_s[nq];
}
model {
  for (n in 1:N) {
    y[n] ~ ordered_logistic(0, c_s[q[n]]);
  }
}
"
m3 <- stan(model_code = s, data = datlist)
print(m3, probs = c(.025, .975))

## Model 4
# all 10 extraversion questions, no predictor, partial pooling regarding scale
datlist <- list(N = nrow(dat2), 
		K = 5, 
		nq = length(unique(dat2$q)),
		q = dat2$q,
		y = dat2$score)

s <- "
data {
  int<lower=2>                         K;
  int<lower=0>                         N;
  int<lower=0>                         nq;
  array[N] int<lower=1, upper=nq>      q;
  array[N] int<lower=1, upper=K>       y;
}
parameters {
  ordered[K - 1] c; // c over all questions
  ordered[K - 1] c_s[nq];
}
model {
  for (s in 1:nq) {
    c_s[s] ~ normal(c, 2);
  }
  for (n in 1:N) {
    y[n] ~ ordered_logistic(0, c_s[q[n]]);
  }
}
"
m4 <- stan(model_code = s, data = datlist)
print(m4, probs = c(.025, .975))

cum_pr_k_hat <- extract(m4)$c |> colMeans() |> plogis()
cum_pr_k_hat |> cum_to_pr() |> as.table() |> plot(ylim = c(0, 0.5))
lines(seq(0,5,.1), 
      dnorm(seq(0,5,.1), mean = mean(dat2$score), sd = sd(dat2$score)))

## Model 5
# add weakly informative prior for c and now use 250 people
# -> decision: just run the one model and change the input scale -> we assume 
# that the scales are independent
dat3_wide <- dat[sample(nrow(dat), 250), c(4, 8:57)]
dat3 <- reshape(dat3_wide, varying = list(2:11, 12:21, 22:31, 32:41, 42:51), 
		v.names = c("E", "N", "A", "C", "O"), 
		timevar = "q", idvar = "id", direction = "long")

datlist_E <- list(N = nrow(dat3), 
		K = 5, 
		nq = length(unique(dat3$q)),
		q = dat3$q,
		y = dat3$E)

# add a prior distribution for c
# What would be an appropriate prior distribution?
pr = c(.15, .20, .30, .20, .15)   # more extreme answers are less likely
rdirichlet <- function(n, alpha) {
  k <- length(alpha)
  x <- matrix(rgamma(n*k, alpha), ncol=k, byrow=TRUE)
  sm <- rowSums(x)
  return(x/sm)
}
# visualize prior predictive check
prior <- rdirichlet(1000, c(4, 5, 6, 5, 4)) |> colMeans() |> as.table()
pdf("prior.pdf", height=3.375, width=3.375, pointsize=10)
par(mai = c(0.5, 0.5, 0.1, 0.1), mgp = c(2, 0.7, 0))
plot(x = 1:5, y = prior, type = "h", ylim = c(0, 0.25), lwd = 10, lend = 2,
     axes = FALSE, ylab = "Probability", 
     xlab = "Response Category")
axis(1, at = c(0, 6))
axis(1, at = 1:5, labels = 1:5)
axis(2, at = c(-1, 0.3))
axis(2, at = seq(0, 0.25, by = 0.05))
dev.off()

s <- "
data {
  int<lower=2>                         K;
  int<lower=0>                         N;
  int<lower=0>                         nq;
  array[N] int<lower=1, upper=nq>      q;
  array[N] int<lower=1, upper=K>       y;
}
parameters {
  ordered[K - 1] c_s[nq];
  simplex[K] theta;  // probability vector
}
transformed parameters {
  ordered[K - 1] c;  // c over all questions
  
  for (k in 1:(K - 1)) {
    c[k] = logit(sum(theta[1:k]));
  }
}
model {
  theta ~ dirichlet([4, 5, 6, 5, 4]');  // Dirichlet prior
  for (s in 1:nq) {
    c_s[s] ~ normal(c, 2);
  }
  for (n in 1:N) {
    y[n] ~ ordered_logistic(0, c_s[q[n]]);
  }
}
"
m5_E <- stan(model_code = s, data = datlist_E)
print(m5_E, pars = "c", probs = c(.025, .975))

cum_pr_k_hat <- extract(m5_E)$c |> colMeans() |> plogis()
cum_pr_k_hat |> cum_to_pr() |> as.table() |> plot(ylim = c(0, 0.5))
lines(seq(0,5,.1), dnorm(seq(0,5,.1), mean = mean(dat3$E), sd = sd(dat3$E)))


## Model 6
# add gender as a predictor 
datlist_E <- list(N = nrow(dat3), 
		 K = 5, 
		nq = length(unique(dat3$q)),
		 q = dat3$q,
		 y = dat3$E,
                 x = dat3$gender)

s <- "
data {
  int<lower=2>                    K;
  int<lower=0>                    N;
  int<lower=0>                    nq;
  array[N] int<lower=1, upper=nq> q;
  array[N] int<lower=1, upper=K>  y;
  array[N] int<lower=1, upper=2>  x;
}
parameters {
  ordered[K - 1] c_s[nq];
  simplex[K] theta;  // probability vector
  real beta;
}
transformed parameters {
  ordered[K - 1] c;  // c over all questions
  
  for (k in 1:(K - 1)) {
    c[k] = logit(sum(theta[1:k]));
  }
}
model {
  theta ~ dirichlet([4, 5, 6, 5, 4]');  // Dirichlet prior
  for (s in 1:nq) {
    c_s[s] ~ normal(c, 2);
  }
  for (n in 1:N) {
    y[n] ~ ordered_logistic(x[n] * beta, c_s[q[n]]);
  }
}
"
m6_E <- stan(model_code = s, data = datlist_E)

# change datlist to analyse each trait
datlist_N <- datlist_A <- datlist_C <- datlist_O <- datlist_E
datlist_N$y = dat3$N
datlist_A$y = dat3$A
datlist_C$y = dat3$C
datlist_O$y = dat3$O
m6_N <- stan(model_code = s, data = datlist_N)
m6_A <- stan(model_code = s, data = datlist_A)
m6_C <- stan(model_code = s, data = datlist_C)
m6_O <- stan(model_code = s, data = datlist_O)
print(m6_E, pars = c("c", "beta"), probs = c(.025, .975))
print(m6_N, pars = c("c", "beta"), probs = c(.025, .975))
print(m6_A, pars = c("c", "beta"), probs = c(.025, .975))
print(m6_C, pars = c("c", "beta"), probs = c(.025, .975))
print(m6_O, pars = c("c", "beta"), probs = c(.025, .975))

# get the response probabilities
get_resp_prob <- function(model) {
  c_samples <- extract(model)$c
  beta_samples <- extract(model)$beta

  cum_pr_k_hat_m <- plogis(c_samples)
  cum_pr_k_hat_f <- plogis(sweep(c_samples, 1, beta_samples, "-"))

  resp_prob_m <- cum_to_pr(cum_pr_k_hat_m)
  resp_prob_f <- cum_to_pr(cum_pr_k_hat_f)

  return(list(resp_prob_m = resp_prob_m, resp_prob_f = resp_prob_f))
}

# visualize the posterior distribution
pdf("main.pdf", height=6, width=4, pointsize=10)
par(mfrow = c(5, 1), mai = c(.1, .4, 0, .1), mgp = c(2, .3, 0),
    omi = c(.3, .2, 0, 0))
for (i in 1:5) {
  scale <- c("E", "N", "A", "C", "O")
  trait_names <- c("Extraversion", "Neuroticism", "Agreeableness", 
		   "Conscientiousness", "Openness")
  model <- paste0("m6_", scale)
  results <- get_resp_prob(get(model[i]))

  # metric for comparison
  mean_m <-  mean(dat3[[scale[i]]][dat3$gender == 1])
  sd_m   <-  sd(dat3[[scale[i]]][dat3$gender == 1])
  mean_f <-  mean(dat3[[scale[i]]][dat3$gender == 2])
  sd_f   <-  sd(dat3[[scale[i]]][dat3$gender == 2])

  plot(x = .90:4.90, y = as.table(colMeans(results$resp_prob_m)), 
       type = "h", ylim = c(0, dnorm(mean_f, mean_f, sd_f)*1.6), lwd = 10,
       axes = FALSE, ylab = "", xlim = c(.95, 5.05), 
       col = "deepskyblue2", lend = 2)
  axis(1, labels = FALSE, lwd.ticks = 0, at = c(0, 6))
  mtext(trait_names[i], side = 2, line = 1, cex = 0.8)
  points(x = 1.10:5.10, y = as.table(colMeans(results$resp_prob_f)), 
	 type = "h", col = rgb(249/255, 32/255, 68/255), lwd = 10, lend = 2)
  lines(seq(0, 6, .1), col = "deepskyblue2",
	dnorm(seq(0, 6, .1), mean = mean_m, sd = sd_m))
  lines(seq(0, 6, .1), col = rgb(249/255, 32/255, 68/255), 
	dnorm(seq(0, 6, .1), mean = mean_f, sd = sd_f))

  # get uncertainty
  sample_idx <- sample(4000, 30)
  for (j in seq_along(sample_idx)) {
    points(x = 1.10:5.10, y = as.table(results$resp_prob_f[sample_idx[j], ]), 
	   type = "p", col = rgb(255/255, 192/255, 203/255, alpha = 0.7))
    points(x = 0.90:4.90, y = as.table(results$resp_prob_m[sample_idx[j], ]), 
	   type = "p", col = rgb(156/255, 216/255, 237/255, alpha = 0.7))
  }
}
mtext("Response Category", 1, 1, outer = TRUE, adj = 0.57)
mtext("Posterior Frequency of Responses", 2, 0, outer = TRUE)
axis(1, labels = 1:5, lwd.ticks = 0, at = 1:5, cex = 3)
legend("topleft", 
       legend = c("Women", "Men"),
       col = c(rgb(249/255, 32/255, 68/255), "deepskyblue2"),
       pch = 15,  # Filled square
       pt.cex = 1.5,  # Size of the square
       cex = 1,  # Text size
       bty = "n",  # No box around the legend
       inset = c(0.05, 0.05))  # Adjust position if needed
dev.off()

## Model 7
# Sensitivity analysis for extraversion scale; flat prior c
s <- "
data {
  int<lower=2>                    K;
  int<lower=0>                    N;
  int<lower=0>                    nq;
  array[N] int<lower=1, upper=nq> q;
  array[N] int<lower=1, upper=K>  y;
  array[N] int<lower=1, upper=2>  x;
}
parameters {
  ordered[K - 1] c_s[nq];
  simplex[K] theta;  // probability vector
  real beta;
}
transformed parameters {
  ordered[K - 1] c;  // c over all questions
  
  for (k in 1:(K - 1)) {
    c[k] = logit(sum(theta[1:k]));
  }
}
model {
  theta ~ dirichlet([1, 1, 1, 1, 1]');  // Dirichlet prior
  for (s in 1:nq) {
    c_s[s] ~ normal(c, 2);
  }
  for (n in 1:N) {
    y[n] ~ ordered_logistic(x[n] * beta, c_s[q[n]]);
  }
}
"
m7_E_sensitivity <- stan(model_code = s, data = datlist_E)
print(m7_E_sensitivity, pars = c("c", "beta"), probs = c(.025, .975))


########## For the analysis with the metric model

# calculating Cohen´s d
cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  
  mean_diff <- mean(x) - mean(y)
  
  s1 <- var(x)
  s2 <- var(y)
  
  pooled_sd <- sqrt(((nx - 1) * s1 + (ny - 1) * s2) / (nx + ny - 2))
  
  d <- mean_diff / pooled_sd
  return(d)
}

aggregate(cbind(E, N, A, C, O) ~ gender, dat3, mean)
aggregate(cbind(E, N, A, C, O) ~ gender, dat3, sd)

cohens_d(dat3$E[dat3$gender == 2], dat3$E[dat3$gender == 1])
cohens_d(dat3$N[dat3$gender == 2], dat3$N[dat3$gender == 1])
cohens_d(dat3$A[dat3$gender == 2], dat3$A[dat3$gender == 1])
cohens_d(dat3$C[dat3$gender == 2], dat3$C[dat3$gender == 1])
cohens_d(dat3$O[dat3$gender == 2], dat3$O[dat3$gender == 1])

t.test(dat3$E[dat3$gender == 1], dat3$E[dat3$gender == 2])
t.test(dat3$N[dat3$gender == 1], dat3$N[dat3$gender == 2])
t.test(dat3$A[dat3$gender == 1], dat3$A[dat3$gender == 2])
t.test(dat3$C[dat3$gender == 1], dat3$C[dat3$gender == 2])
t.test(dat3$O[dat3$gender == 1], dat3$O[dat3$gender == 2])
