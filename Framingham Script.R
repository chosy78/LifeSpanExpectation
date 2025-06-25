library(tidyverse)
library(readr)
library(caret)
library(kableExtra)
library(car)
library(FSA)

set.seed(5)

# Read in the data (set working directory to data location on your computer)
raw_data <- read.csv("fhs_03MAR23.csv", header = TRUE)
raw_data[raw_data == 'M'] <- NA

# Baseline covariates of interest: age, sex 
# Additional covariates of interest: bmi, sbp, dbp, cholesterol

################################################################################
##################         PART 1: Variable Selection         ##################
################################################################################


# Create flags for extreme valued covariates
framingham <- raw_data %>% mutate(surv_time = as.numeric(DATEDTH), age = as.numeric(age_or_ex1), sex = female, 
                                    sbp = as.numeric(SBP_or_ex1_1exam), dbp = as.numeric(DBP_or_ex1_1exam), 
                                    bmi = as.numeric(bmi_or_ex1), chol = as.numeric(Chol_or_ex1)) %>%
  select(surv_time, age, sex, sbp, dbp, bmi, chol) %>%
  filter(complete.cases(.)) %>% 
  mutate(sex = as.factor(ifelse(sex == 1, "Female", ifelse(sex == 0, "Male", "Missing"))),
         surv_time = as.numeric(surv_time)/365.25, bmi_flag = ifelse(bmi >= 30, 1, 0), sbp_flag = ifelse(sbp >= 140, 1, 0),
         dbp_flag = ifelse(dbp >= 90, 1, 0), chol_flag = ifelse(chol >= 200, 1, 0))
 
# Check correlations 
framingham %>% mutate(sex = ifelse(sex == "Female", 1, 0)) %>% 
  select(surv_time, age, sex, sbp, dbp, bmi, chol) %>%
  cor()
  
# Create 10 folds for cross validation
rand <- sample(nrow(framingham))
samples <- list()
for (s in 1:10){
  samples[[s]] <- rand[rand %% 10 + 1 == s]
}

# Create empty dataframes for MSE and bias
MSE <- data.frame(matrix(ncol=3, nrow=10))
colnames(MSE) <- c('Model 1', 'Model 2', 'Model 3')
bias <- data.frame(matrix(ncol=3, nrow=10))
colnames(bias) <- c('Model 1', 'Model 2', 'Model 3')

# Fit all three models
for (s in 1:10){
  train_set <- framingham[-samples[[s]],]
  test_set <- framingham[samples[[s]],] %>% filter(bmi_flag == 1 | sbp_flag == 1 | dbp_flag == 1 | chol_flag == 1)
  
  test_fv1 <- vector()
  for (i in 1:nrow(test_set)){
    bmi_ind <- test_set[i,]$bmi_flag
    bp_ind <- test_set[i,]$sbp_flag+test_set[i,]$dbp_flag
    chol_ind <- test_set[i,]$chol_flag
    if (bmi_ind == 1 & bp_ind >= 1 & chol_ind == 1){
      model1 <- lm(data = train_set, surv_time ~ age + sex + bmi + sbp + dbp + chol)
    } else if (bmi_ind == 1 & bp_ind >= 1 & chol_ind == 0){
      model1 <- lm(data = train_set, surv_time ~ age + sex + bmi + sbp + dbp)
    } else if (bmi_ind == 1 & bp_ind == 0 & chol_ind == 1){
      model1 <- lm(data = train_set, surv_time ~ age + sex + bmi + chol)
    } else if (bmi_ind == 1 & bp_ind == 0 & chol_ind == 0){
      model1 <- lm(data = train_set, surv_time ~ age + sex + bmi)
    } else if (bmi_ind == 0 & bp_ind >= 1 & chol_ind == 1){
      model1 <- lm(data = train_set, surv_time ~ age + sex + sbp + dbp + chol)
    } else if (bmi_ind == 0 & bp_ind >= 1 & chol_ind == 0){
      model1 <- lm(data = train_set, surv_time ~ age + sex + sbp + dbp)
    } else if (bmi_ind == 0 & bp_ind == 0 & chol_ind == 1){
      model1 <- lm(data = train_set, surv_time ~ age + sex + chol)
    } else if (bmi_ind == 0 & bp_ind == 0 & chol_ind == 0){
      model1 <- lm(data = train_set, surv_time ~ age + sex)
    } else {
      model1 <- NA
    }
    test_fv1[i] <- predict(model1, newdata = test_set[i,])
  }
  
  model2 <- lm(data = train_set, surv_time ~ age + sex)
  # summary(model2)
  
  model3 <- lm(data = train_set, surv_time ~ age + sex + bmi + sbp + dbp + chol)
  # summary(model3)
  
  test_fv2 <- predict(model2, newdata = test_set)
  test_fv3 <- predict(model3, newdata = test_set)
  
  MSE[s,1] <- mean((test_fv1 - test_set$surv_time)^2)
  MSE[s,2] <- mean((test_fv2 - test_set$surv_time)^2)
  MSE[s,3] <- mean((test_fv3 - test_set$surv_time)^2)
  
  bias[s,1] <- mean(test_fv1 - test_set$surv_time)
  bias[s,2] <- mean(test_fv2 - test_set$surv_time)
  bias[s,3] <- mean(test_fv3 - test_set$surv_time)
}

# F-test model 2 vs model 3
anova(model2, model3, test="F")

# Format output table
MSE_bias_table <- round(matrix(c(mean(MSE[,1]), mean(MSE[,2]), mean(MSE[,3]), 
                                 mean(bias[,1]), mean(bias[,2]), mean(bias[,3])), ncol = 2),2)
colnames(MSE_bias_table) <- c('MSE','Bias')
rownames(MSE_bias_table) <- c('Model 1', 'Model 2', 'Model 3')
kbl(MSE_bias_table) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  add_header_above(c("Variable Selection" = 3))


################################################################################
##################          PART 2: Measurement Error         ##################
################################################################################

################################################################################
#################          PART 2a: Adding Random Noise        #################
################################################################################

# Find mean value of covariate over first few visits
# This will be used to empirically model mean and sd for noise

sbp_sd <- raw_data %>% select(SBP_or_ex1_adm, SBP_or_ex1_1exam, SBP_or_ex1_2exam,
                              SBP_or_ex2_adm, SBP_or_ex2_1exam, SBP_or_ex2_2exam,
                              SBP_or_ex3_adm, SBP_or_ex3_1exam, SBP_or_ex3_2exam) %>% 
          filter(complete.cases(.)) %>% mutate_all(as.numeric) %>% as.matrix() %>% apply(., 1, sd) %>% mean()

dbp_sd <- raw_data %>% select(DBP_or_ex1_adm, DBP_or_ex1_1exam, DBP_or_ex1_2exam,
                              DBP_or_ex2_adm, DBP_or_ex2_1exam, DBP_or_ex2_2exam,
                              DBP_or_ex3_adm, DBP_or_ex3_1exam, DBP_or_ex3_2exam) %>% 
          filter(complete.cases(.)) %>% mutate_all(as.numeric) %>% as.matrix() %>% apply(., 1, sd) %>% mean()

chol_sd <- raw_data %>% select(Chol_or_ex1, Chol_or_ex2, Chol_or_ex3, Chol_or_ex4, Chol_or_ex5) %>%
           filter(complete.cases(.)) %>% mutate_all(as.numeric) %>% as.matrix() %>% apply(., 1, sd) %>% mean()

bmi_sd <- raw_data %>% select(bmi_or_ex1, bmi_or_ex4, bmi_or_ex5) %>%
          filter(complete.cases(.)) %>% mutate_all(as.numeric) %>% as.matrix() %>% apply(., 1, sd) %>% mean()


# Create empty data frames for MSE and bias
mse_df <- data.frame(matrix(ncol=10, nrow=48))
bias_df <- data.frame(matrix(ncol=10, nrow=48))
f_test <- list()

# Conduct variable selection process but now with random noise
counter <- 0
for (s in 1:10){
  train_set <- framingham[-samples[[s]],] %>% 
    mutate(sbp.me = sbp, dbp.me = dbp, chol.me = chol)
  
  mse_vec <- vector()
  bias_vec <- vector()
  mse_list <- list()
  bias_list <- list()
  j <- 1
  
  # Assess different levels of random noise
  for (m in c(0, 0.5, 1, 2)){
    for (v in c(0, 0.5, 1, 2)){
      test_fv_temp1 <- vector()
      test_fv_temp2 <- vector()
      test_fv_temp3 <- vector()
      framingham_me <- framingham %>% mutate(sbp_temp = sbp + rnorm(length(framingham$sbp), m*sbp_sd, v*sbp_sd))
      names(framingham_me)[ncol(framingham_me)] <- paste0("sbp",j)
      pred.sbp <- paste0("sbp",j)
      framingham_me <- framingham_me %>% mutate(dbp_temp = dbp + rnorm(length(framingham$dbp), m*dbp_sd, v*dbp_sd))
      names(framingham_me)[ncol(framingham_me)] <- paste0("dbp",j)
      pred.dbp <- paste0("dbp",j)
      framingham_me <- framingham_me %>% mutate(chol_temp = chol + rnorm(length(framingham$chol), m*chol_sd, v*chol_sd))
      names(framingham_me)[ncol(framingham_me)] <- paste0("chol",j)
      pred.chol <- paste0("chol",j)
      test_set_me <- framingham_me[samples[[s]],] %>% filter(bmi_flag == 1 | sbp_flag == 1 | dbp_flag == 1 | chol_flag == 1)
      names(train_set)[(ncol(train_set)-2):ncol(train_set)] <- c(paste0("sbp",j), paste0("dbp",j), paste0("chol",j))
      for (i in 1:nrow(test_set_me)){
        bmi_ind <- test_set_me[i,]$bmi_flag
        bp_ind <- ifelse((test_set_me[i,ncol(test_set_me)-2] >= 140 | test_set_me[i,ncol(test_set_me)-1] >= 90) , 1, 0)
        chol_ind <- ifelse(test_set_me[i,ncol(test_set_me)] >= 200 , 1, 0)
        if (bmi_ind == 1 & bp_ind == 1 & chol_ind == 1){
          model1 <- lm(data = train_set, surv_time ~ age + sex + bmi + get(pred.sbp) + get(pred.dbp) + get(pred.chol))
        } else if (bmi_ind == 1 & bp_ind == 1 & chol_ind == 0){
          model1 <- lm(data = train_set, surv_time ~ age + sex + bmi + get(pred.sbp) + get(pred.dbp))
        } else if (bmi_ind == 1 & bp_ind == 0 & chol_ind == 1){
          model1 <- lm(data = train_set, surv_time ~ age + sex + bmi + get(pred.chol))
        } else if (bmi_ind == 1 & bp_ind == 0 & chol_ind == 0){
          model1 <- lm(data = train_set, surv_time ~ age + sex + bmi)
        } else if (bmi_ind == 0 & bp_ind == 1 & chol_ind == 1){
          model1 <- lm(data = train_set, surv_time ~ age + sex + get(pred.sbp) + get(pred.dbp) + get(pred.chol))
        } else if (bmi_ind == 0 & bp_ind == 1 & chol_ind == 0){
          model1 <- lm(data = train_set, surv_time ~ age + sex + get(pred.sbp) + get(pred.dbp))
        } else if (bmi_ind == 0 & bp_ind == 0 & chol_ind == 1){
          model1 <- lm(data = train_set, surv_time ~ age + sex + get(pred.chol))
        } else if (bmi_ind == 0 & bp_ind == 0 & chol_ind == 0){
          model1 <- lm(data = train_set, surv_time ~ age + sex)
        } else {
          model1 <- NA
        }
        test_fv_temp1[i] <- predict(model1, newdata = test_set_me[i,])
      }
      
      model2 <- lm(data = train_set, surv_time ~ age + sex)
      # summary(model2)
      
      model3 <- lm(data = train_set, surv_time ~ age + sex + bmi + get(pred.sbp) + get(pred.dbp) + get(pred.chol))
      # summary(model3)
      
      
      test_fv_temp2 <- predict(model2, newdata = test_set_me)
      test_fv_temp3 <- predict(model3, newdata = test_set_me)
      
      mse_vec[3*j-2] <- mean((test_fv_temp1 - test_set_me$surv_time)^2)
      mse_vec[3*j-1] <- mean((test_fv_temp2 - test_set_me$surv_time)^2)
      mse_vec[3*j] <- mean((test_fv_temp3 - test_set_me$surv_time)^2)
      
      bias_vec[3*j-2] <- mean(test_fv_temp1 - test_set_me$surv_time)
      bias_vec[3*j-1] <- mean(test_fv_temp2 - test_set_me$surv_time)
      bias_vec[3*j] <- mean(test_fv_temp3 - test_set_me$surv_time)
      
      mse_list[[j]] <- mse_vec
      bias_list[[j]] <- bias_vec
      j <- j+1
    }
  }
  mse_df[,s] <- c(mse_vec)
  bias_df[,s] <- c(bias_vec)
  f_test[[s]] <- anova(model2, model3, test="F")
  counter <- counter + 1
  print(counter)
  
}

# Format output tables
mse_table <- round(matrix(rowMeans(mse_df), nrow = 16, ncol = 3, byrow=TRUE),2)
bias_table <- round(matrix(rowMeans(bias_df), nrow = 16, ncol = 3, byrow=TRUE),2)

rownames(mse_table) <- c('No Measurement Error', 'Mean: 0 | SD: 0.1', 'Mean: 0 | SD: 0.5', 'Mean: 0 | SD: 1.0',
                         'Mean: 0.1 | SD: 0', 'Mean: 0.1 | SD: 0.1', 'Mean: 0.1 | SD: 0.5', 'Mean: 0.1 | SD: 1.0',
                         'Mean: 0.2 | SD: 0', 'Mean: 0.2 | SD: 0.1', 'Mean: 0.2 | SD: 0.5', 'Mean: 0.2 | SD: 1.0',
                         'Mean: 0.5 | SD: 0', 'Mean: 0.5 | SD: 0.1', 'Mean: 0.5 | SD: 0.5', 'Mean: 0.5 | SD: 1.0' )
colnames(mse_table) <- c('Model 1', 'Model 2', 'Model 3')
rownames(bias_table) <- c('No Measurement Error', 'Mean: 0 | SD: 0.1', 'Mean: 0 | SD: 0.5', 'Mean: 0 | SD: 1.0',
                          'Mean: 0.1 | SD: 0', 'Mean: 0.1 | SD: 0.1', 'Mean: 0.1 | SD: 0.5', 'Mean: 0.1 | SD: 1.0',
                          'Mean: 0.2 | SD: 0', 'Mean: 0.2 | SD: 0.1', 'Mean: 0.2 | SD: 0.5', 'Mean: 0.2 | SD: 1.0',
                          'Mean: 0.5 | SD: 0', 'Mean: 0.5 | SD: 0.1', 'Mean: 0.5 | SD: 0.5', 'Mean: 0.5 | SD: 1.0' )
colnames(bias_table) <- c('Model 1', 'Model 2', 'Model 3')

# kbl(mse_table) %>%
#   kable_styling(bootstrap_options = "striped", full_width = F)
# kbl(bias_table) %>%
#   kable_styling(bootstrap_options = "striped", full_width = F)

combined_table <- cbind(mse_table, bias_table)

kbl(combined_table) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  add_header_above(c(" " = 1, "MSE" = 3, "Bias" = 3))



# Create plots for measurement error results
plot_df <- as.data.frame(combined_table)
colnames(plot_df) <- c('m1', 'm2', 'm3', 'b1', 'b2', 'b3') 

mse_plot_df <- plot_df[,1:3] %>% as.data.frame() %>%
  mutate(sd = c(0, 0.5, 1.0, 2.0, 0, 0.5, 1.0, 2.0, 0, 0.5, 1.0, 2.0, 0, 0.5, 1.0, 2.0)) %>%
  mutate(mu = c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0)) %>%
  pivot_longer(cols = c('m1', 'm2', 'm3'), names_to = 'model', values_to = 'mse') %>%
  mutate(model = as.numeric(gsub('m', '', model))) %>%
  mutate(mu = as.factor(mu))

bias_plot_df <- plot_df[,4:6] %>% as.data.frame() %>%
  mutate(sd = c(0, 0.5, 1.0, 2.0, 0, 0.5, 1.0, 2.0, 0, 0.5, 1.0, 2.0, 0, 0.5, 1.0, 2.0)) %>%
  mutate(mu = c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0)) %>%
  pivot_longer(cols = c('b1', 'b2', 'b3'), names_to = 'model', values_to = 'bias') %>%
  mutate(model = as.numeric(gsub('b', '', model))) %>%
  mutate(mu = as.factor(mu))



ggplot(mse_plot_df, aes(x = sd, y = mse, color=as.factor(model))) +
  geom_line(lwd = 1) +
  geom_point() +
  facet_wrap( ~ mu, ncol = 4, 
              labeller = labeller(mu = c("0" = "Mean: 0",
                              "0.5" = "Mean: 0.5",
                              "1" = "Mean: 1.0",
                              "2" = "Mean: 2.0"))) +
  #ggtitle("Measurement Error Effect of MSE") +
  xlab("Measurement Error Variance") +
  ylab("Mean Squared Error") +
  scale_color_discrete(name="Model", labels=c("Variable Selection Model", "Baseline Model", "Saturated Model")) +
  theme_minimal() +
  theme(legend.position="right", plot.title = element_text(hjust = 0.5))

ggplot(bias_plot_df, aes(x = sd, y = bias, color=as.factor(model))) +
  geom_line(lwd = 1) +
  geom_point() +
  geom_hline(yintercept=0, linetype='dotted', lwd = 1) +
  facet_wrap( ~ mu, ncol = 4, 
              labeller = labeller(mu = c("0" = "Mean: 0",
                                         "0.5" = "Mean: 0.5",
                                         "1" = "Mean: 1.0",
                                         "2" = "Mean: 2.0"))) +
  #ggtitle("Measurement Error Effect of Bias") +
  xlab("Measurement Error Variance") +
  ylab("Bias") +
  scale_color_discrete(name="Model", labels=c("Variable Selection Model", "Baseline Model", "Saturated Model")) +
  theme_minimal() +
  theme(legend.position="right", plot.title = element_text(hjust = 0.5))






################################################################################
###############      PART 2b: Taking the Maximum Measurement    ################
################################################################################


# Find maximum measurement for all covariates across different exams
framingham_max <- raw_data %>% mutate(surv_time = as.numeric(DATEDTH), age = as.numeric(age_or_ex1), sex = female, 
                                  sbp = as.numeric(SBP_or_ex1_1exam), dbp = as.numeric(DBP_or_ex1_1exam), 
                                  bmi = as.numeric(bmi_or_ex1), chol = as.numeric(Chol_or_ex1)) %>%
  mutate(sbp_max = pmax(SBP_or_ex1_adm, SBP_or_ex1_1exam, SBP_or_ex1_2exam,
                        SBP_or_ex2_adm, SBP_or_ex2_1exam, SBP_or_ex2_2exam,
                        SBP_or_ex3_adm, SBP_or_ex3_1exam, SBP_or_ex3_2exam, na.rm = TRUE)) %>%
  mutate(dbp_max = pmax(DBP_or_ex1_adm, DBP_or_ex1_1exam, DBP_or_ex1_2exam,
                        DBP_or_ex2_adm, DBP_or_ex2_1exam, DBP_or_ex2_2exam,
                        DBP_or_ex3_adm, DBP_or_ex3_1exam, DBP_or_ex3_2exam, na.rm = TRUE)) %>%
  mutate(chol_max = pmax(Chol_or_ex1, Chol_or_ex2, Chol_or_ex3, Chol_or_ex4, Chol_or_ex5, na.rm = TRUE)) %>%
  mutate(bmi_max = pmax(bmi_or_ex1, bmi_or_ex4, bmi_or_ex3, bmi_or_ex4, bmi_or_ex5, na.rm = TRUE)) %>%
  mutate(sbp_max = as.numeric(sbp_max), dbp_max = as.numeric(dbp_max), 
         bmi_max = as.numeric(bmi_max), chol_max = as.numeric(chol_max)) %>%
  select(surv_time, age, sex, sbp, dbp, bmi, chol, sbp_max, dbp_max, chol_max, bmi_max) %>%
  filter(complete.cases(.)) %>% 
  mutate(sex = as.factor(ifelse(sex == 1, "Female", ifelse(sex == 0, "Male", "Missing"))),
         surv_time = as.numeric(surv_time)/365.25, bmi_flag = ifelse(bmi >= 30, 1, 0), 
         sbp_flag = ifelse(sbp >= 140, 1, 0), dbp_flag = ifelse(dbp >= 90, 1, 0), 
         chol_flag = ifelse(chol >= 200, 1, 0), bmi_max_flag = ifelse(bmi_max >= 30, 1, 0), 
         sbp_max_flag = ifelse(sbp_max >= 140, 1, 0), dbp_max_flag = ifelse(dbp_max >= 90, 1, 0), 
         chol_max_flag = ifelse(chol_max >= 200, 1, 0))


# Define empty dataframes for MSE and bias
MSE_max <- data.frame(matrix(ncol=3, nrow=10))
colnames(MSE_max) <- c('Model 1', 'Model 2', 'Model 3')
bias_max <- data.frame(matrix(ncol=3, nrow=10))
colnames(bias_max) <- c('Model 1', 'Model 2', 'Model 3')

# fit all three models
for (s in 1:10){
  train_set_max <- framingham_max[-samples[[s]],] %>%
    mutate(sbp_max = sbp, dbp_max = dbp, chol_max = chol, bmi_max = bmi)
  test_set_max <- framingham_max[samples[[s]],] %>% filter(bmi_flag == 1 | sbp_flag == 1 | dbp_flag == 1 | chol_flag == 1)
  
  
  test_fv1_max <- vector()
  for (i in 1:nrow(test_set_max)){
    bmi_max_ind <- test_set_max[i,]$bmi_max_flag
    bp_max_ind <- test_set_max[i,]$sbp_max_flag+test_set_max[i,]$dbp_max_flag
    chol_max_ind <- test_set_max[i,]$chol_max_flag
    if (bmi_max_ind == 1 & bp_max_ind >= 1 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + bmi_max + sbp_max + dbp_max + chol_max)
    } else if (bmi_max_ind == 1 & bp_max_ind >= 1 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + bmi_max + sbp_max + dbp_max)
    } else if (bmi_max_ind == 1 & bp_max_ind == 0 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + bmi_max + chol_max)
    } else if (bmi_max_ind == 1 & bp_max_ind == 0 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + bmi_max)
    } else if (bmi_max_ind == 0 & bp_max_ind >= 1 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + sbp_max + dbp_max + chol_max)
    } else if (bmi_max_ind == 0 & bp_max_ind >= 1 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + sbp_max + dbp_max)
    } else if (bmi_max_ind == 0 & bp_max_ind == 0 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex + chol_max)
    } else if (bmi_max_ind == 0 & bp_max_ind == 0 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max, surv_time ~ age + sex)
    } else {
      model1_max <- NA
    }
    test_fv1_max[i] <- predict(model1_max, newdata = test_set_max[i,])
  }
  
  model2_max <- lm(data = train_set_max, surv_time ~ age + sex)
  # summary(model2)
  
  model3_max <- lm(data = train_set_max, surv_time ~ age + sex + bmi_max + sbp_max + dbp_max + chol_max)
  # summary(model3)
  
  
  test_fv2_max <- predict(model2_max, newdata = test_set_max)
  test_fv3_max <- predict(model3_max, newdata = test_set_max)
  
  
  MSE_max[s,1] <- mean((test_fv1_max - test_set_max$surv_time)^2)
  MSE_max[s,2] <- mean((test_fv2_max - test_set_max$surv_time)^2)
  MSE_max[s,3] <- mean((test_fv3_max - test_set_max$surv_time)^2)
  
  bias_max[s,1] <- mean(test_fv1_max - test_set_max$surv_time)
  bias_max[s,2] <- mean(test_fv2_max - test_set_max$surv_time)
  bias_max[s,3] <- mean(test_fv3_max - test_set_max$surv_time)
}

# Format output tables 
MSE_bias_table_max <- round(matrix(c(mean(MSE_max[,1]), mean(MSE_max[,2]), mean(MSE_max[,3]), 
                                 mean(bias_max[,1]), mean(bias_max[,2]), mean(bias_max[,3])), ncol = 2),2)
colnames(MSE_bias_table_max) <- c('MSE','Bias')
rownames(MSE_bias_table_max) <- c('Model 1', 'Model 2', 'Model 3')
kbl(MSE_bias_table_max) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  add_header_above(c("Variable Selection" = 3))









################################################################################
###########      PART 2c: Combining Measurement Error Paradigms      ###########
################################################################################


# Add random noise and then find maximum measurement for all covariates across different exams
framingham_max_me <- raw_data %>% mutate(surv_time = as.numeric(DATEDTH), age = as.numeric(age_or_ex1), sex = female, 
                                      sbp = as.numeric(SBP_or_ex1_1exam), dbp = as.numeric(DBP_or_ex1_1exam), 
                                      bmi = as.numeric(bmi_or_ex1), chol = as.numeric(Chol_or_ex1)) %>%
  mutate(sbp_max = pmax(as.numeric(SBP_or_ex1_adm) + rnorm(1, 0, sbp_sd), 
                        as.numeric(SBP_or_ex1_1exam) + rnorm(1, 0, sbp_sd), 
                        as.numeric(SBP_or_ex1_2exam) + rnorm(1, 0, sbp_sd),
                        as.numeric(SBP_or_ex2_adm) + rnorm(1, 0, sbp_sd), 
                        as.numeric(SBP_or_ex2_1exam) + rnorm(1, 0, sbp_sd), 
                        as.numeric(SBP_or_ex2_2exam) + rnorm(1, 0, sbp_sd),
                        as.numeric(SBP_or_ex3_adm) + rnorm(1, 0, sbp_sd), 
                        as.numeric(SBP_or_ex3_1exam) + rnorm(1, 0, sbp_sd), 
                        as.numeric(SBP_or_ex3_2exam) + rnorm(1, 0, sbp_sd), na.rm = TRUE)) %>%
  mutate(dbp_max = pmax(as.numeric(DBP_or_ex1_adm) + rnorm(1, 0, dbp_sd), 
                        as.numeric(DBP_or_ex1_1exam) + rnorm(1, 0, dbp_sd), 
                        as.numeric(DBP_or_ex1_2exam) + rnorm(1, 0, dbp_sd),
                        as.numeric(DBP_or_ex2_adm) + rnorm(1, 0, dbp_sd), 
                        as.numeric(DBP_or_ex2_1exam) + rnorm(1, 0, dbp_sd), 
                        as.numeric(DBP_or_ex2_2exam) + rnorm(1, 0, dbp_sd),
                        as.numeric(DBP_or_ex3_adm) + rnorm(1, 0, dbp_sd), 
                        as.numeric(DBP_or_ex3_1exam) + rnorm(1, 0, dbp_sd), 
                        as.numeric(DBP_or_ex3_2exam) + rnorm(1, 0, dbp_sd), na.rm = TRUE)) %>%
  mutate(chol_max = pmax(as.numeric(Chol_or_ex1) + rnorm(1, 0, chol_sd), 
                         as.numeric(Chol_or_ex2) + rnorm(1, 0, chol_sd), 
                         as.numeric(Chol_or_ex3) + rnorm(1, 0, chol_sd), 
                         as.numeric(Chol_or_ex4) + rnorm(1, 0, chol_sd), 
                         as.numeric(Chol_or_ex5) + rnorm(1, 0, chol_sd), na.rm = TRUE)) %>%
  mutate(bmi_max = pmax(as.numeric(bmi_or_ex1), 
                        as.numeric(bmi_or_ex4), 
                        as.numeric(bmi_or_ex3), 
                        as.numeric(bmi_or_ex4), 
                        as.numeric(bmi_or_ex5), na.rm = TRUE)) %>%
  mutate(sbp_max = as.numeric(sbp_max), dbp_max = as.numeric(dbp_max), 
         bmi_max = as.numeric(bmi_max), chol_max = as.numeric(chol_max)) %>%
  select(surv_time, age, sex, sbp, dbp, bmi, chol, sbp_max, dbp_max, chol_max, bmi_max) %>%
  filter(complete.cases(.)) %>% 
  mutate(sex = as.factor(ifelse(sex == 1, "Female", ifelse(sex == 0, "Male", "Missing"))),
         surv_time = as.numeric(surv_time)/365.25, bmi_flag = ifelse(bmi >= 30, 1, 0), 
         sbp_flag = ifelse(sbp >= 140, 1, 0), dbp_flag = ifelse(dbp >= 90, 1, 0), 
         chol_flag = ifelse(chol >= 200, 1, 0), bmi_max_flag = ifelse(bmi_max >= 30, 1, 0), 
         sbp_max_flag = ifelse(sbp_max >= 140, 1, 0), dbp_max_flag = ifelse(dbp_max >= 90, 1, 0), 
         chol_max_flag = ifelse(chol_max >= 200, 1, 0))


# Define empty dataframes for MSE and bias
MSE_max_me <- data.frame(matrix(ncol=3, nrow=10))
colnames(MSE_max_me) <- c('Model 1', 'Model 2', 'Model 3')
bias_max_me <- data.frame(matrix(ncol=3, nrow=10))
colnames(bias_max_me) <- c('Model 1', 'Model 2', 'Model 3')

# fit all three models
for (s in 1:10){
  train_set_max_me <- framingham_max_me[-samples[[s]],] %>%
    mutate(sbp_max = sbp, dbp_max = dbp, chol_max = chol, bmi_max = bmi)
  test_set_max_me <- framingham_max_me[samples[[s]],] %>% filter(bmi_flag == 1 | sbp_flag == 1 | dbp_flag == 1 | chol_flag == 1)
  
  
  test_fv1_max_me <- vector()
  for (i in 1:nrow(test_set_max_me)){
    bmi_max_ind <- test_set_max_me[i,]$bmi_max_flag
    bp_max_ind <- test_set_max_me[i,]$sbp_max_flag + test_set_max_me[i,]$dbp_max_flag
    chol_max_ind <- test_set_max_me[i,]$chol_max_flag
    if (bmi_max_ind == 1 & bp_max_ind >= 1 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + bmi_max + sbp_max + dbp_max + chol_max)
    } else if (bmi_max_ind == 1 & bp_max_ind >= 1 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + bmi_max + sbp_max + dbp_max)
    } else if (bmi_max_ind == 1 & bp_max_ind == 0 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + bmi_max + chol_max)
    } else if (bmi_max_ind == 1 & bp_max_ind == 0 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + bmi_max)
    } else if (bmi_max_ind == 0 & bp_max_ind >= 1 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + sbp_max + dbp_max + chol_max)
    } else if (bmi_max_ind == 0 & bp_max_ind >= 1 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + sbp_max + dbp_max)
    } else if (bmi_max_ind == 0 & bp_max_ind == 0 & chol_max_ind == 1){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex + chol_max)
    } else if (bmi_max_ind == 0 & bp_max_ind == 0 & chol_max_ind == 0){
      model1_max <- lm(data = train_set_max_me, surv_time ~ age + sex)
    } else {
      model1_max <- NA
    }
    test_fv1_max_me[i] <- predict(model1_max, newdata = test_set_max_me[i,])
  }
  
  model2_max_me <- lm(data = train_set_max_me, surv_time ~ age + sex)
  
  model3_max_me <- lm(data = train_set_max_me, surv_time ~ age + sex + bmi_max + sbp_max + dbp_max + chol_max)
  
  
  test_fv2_max_me <- predict(model2_max, newdata = test_set_max_me)
  test_fv3_max_me <- predict(model3_max, newdata = test_set_max_me)
  
  
  MSE_max_me[s,1] <- mean((test_fv1_max_me - test_set_max_me$surv_time)^2)
  MSE_max_me[s,2] <- mean((test_fv2_max_me - test_set_max_me$surv_time)^2)
  MSE_max_me[s,3] <- mean((test_fv3_max_me - test_set_max_me$surv_time)^2)
  
  bias_max_me[s,1] <- mean(test_fv1_max_me - test_set_max_me$surv_time)
  bias_max_me[s,2] <- mean(test_fv2_max_me - test_set_max_me$surv_time)
  bias_max_me[s,3] <- mean(test_fv3_max_me - test_set_max_me$surv_time)
}

# Format output tables 
MSE_bias_table_max_me <- round(matrix(c(mean(MSE_max_me[,1]), mean(MSE_max_me[,2]), mean(MSE_max_me[,3]), 
                                     mean(bias_max_me[,1]), mean(bias_max_me[,2]), mean(bias_max_me[,3])), ncol = 2),2)
colnames(MSE_bias_table_max_me) <- c('MSE','Bias')
rownames(MSE_bias_table_max_me) <- c('Model 1', 'Model 2', 'Model 3')
kbl(MSE_bias_table_max_me) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  add_header_above(c("Variable Selection" = 3))




