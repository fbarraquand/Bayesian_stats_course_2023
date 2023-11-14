# load packages
library(tidyverse)
theme_set(theme_light(base_size = 16))
library(gganimate)
library(magick)

# deer data, 19 "success" out of 57 "attempts"
survived <- 19
released <- 57

# log-likelihood function
loglikelihood <- function(x, p){
  dbinom(x = x, size = released, prob = p, log = TRUE)
}

# prior density
logprior <- function(p){
  dunif(x = p, min = 0, max = 1, log = TRUE)
}

# posterior density function (log scale)
posterior <- function(x, p){
  loglikelihood(x, p) + logprior(p) # - log(Pr(data))
}

# propose candidate value
move <- function(x, away = .2){ 
  logitx <- log(x / (1 - x))
  logit_candidate <- logitx + rnorm(1, 0, away)
  candidate <- plogis(logit_candidate)
  return(candidate)
}

metropolis <- function(steps = 100, inits = 0.5){
  
  # pre-alloc memory
  theta.post <- rep(NA, steps)
  
  # start
  theta.post[1] <- inits
  
  for (t in 2:steps){
    
    # propose candidate value for prob of success
    theta_star <- move(theta.post[t-1])
    
    # calculate ratio R
    pstar <- posterior(survived, p = theta_star)  
    pprev <- posterior(survived, p = theta.post[t-1])
    logR <- pstar - pprev
    R <- exp(logR)
    
    # decide to accept candidate value or to keep current value
    accept <- rbinom(1, 1, prob = min(R, 1))
    theta.post[t] <- ifelse(accept == 1, theta_star, theta.post[t-1])
  }
  theta.post
}

#---------- apply Metropolis

steps <- 1000
chain1 <- metropolis(steps = steps, inits = 0.2)
chain2 <- metropolis(steps = steps, inits = 0.5)
chain3 <- metropolis(steps = steps, inits = 0.7)

df <- data.frame(iter = rep(1:steps, 3), 
                 value = c(chain1, chain2, chain3),
                 chain = c(rep("chain1", steps), 
                           rep("chain2", steps), 
                           rep("chain3", steps)))

#---------- time series
static_tsplot <- df %>%
  mutate(posterior_mean = mean(value)) %>%
  ggplot(aes(x = iter, y = value, group = chain, color = chain)) +
  geom_line(size = 1, alpha = 0.7) + 
  geom_hline(aes(yintercept = posterior_mean, linetype = "posterior mean")) + 
  scale_linetype_manual(name = "", values = c(2,2)) + 
  labs(color = "", x = "iteration", y = "probability of success")
static_tsplot  

# animate
animated_tsplot <- static_tsplot +
  transition_reveal(along = iter, 
                    range = as.integer(c(1, max(df$iter) + 50))) # trick to pause
animated_tsplot  

# save
a_gif <- animate(animated_tsplot,
                 width = 600, 
                 height = 300)


#---------- histogram

# histogram
static_hist <- df %>% 
  mutate(posterior_mean = mean(value)) %>%
  split(.$iter) %>% 
  accumulate(~ bind_rows(.x, .y)) %>% 
  bind_rows(.id = "frame") %>% 
  mutate(frame = as.integer(frame)) %>%
  ggplot(aes(x = value, fill = chain)) +
  geom_histogram(color = "white", bins = 15, alpha = 0.7, position = "identity") + 
  labs(x = "probability of success", fill = "") +
  geom_vline(aes(xintercept = posterior_mean), lty = 2)
static_hist

# animate
anim_hist <- static_hist + 
  transition_manual(frame) +
  ease_aes("linear") +
  enter_fade() +
  exit_fade()
anim_hist

# save
b_gif <- animate(anim_hist,
                 width = 600, 
                 height = 300)


#---------- put side-by-side
new_gif <- image_append(c(a_gif[1], b_gif[1]))
for(i in 2:100){
  combined <- image_append(c(a_gif[i], b_gif[i]))
  new_gif <- c(new_gif, combined)
}

new_gif

