## ----library, include=FALSE---------------------------------------------
# Download core libraries
if (!require(tidyverse)) {
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")
}
if (!require(caret)) {
  install.packages("caret", repos = "http://cran.us.r-project.org")
}

library(tidyverse)
library(caret)
options(digits = 3)


# Dataset ****************************************************************


## ----download-data, include=FALSE---------------------------------------
# Download age tsv
dl <- tempfile()
download.file("https://genomics.senescence.info/species/dataset.zip", dl)

age <- read_tsv(unzip(dl, "anage_data.txt"))

# The indexes below give us numerical columns, for each col,
# do a 5-number summary and count num of NAs. Then save it for the appendix.
summary_age <-
  sapply(age[, c(10:21, 26:30)], function(x)
    summary(x))
summary_age <- as_tibble(t(summary_age), rownames = "Var")
save(summary_age, file = "summary_age.rda")

rm(dl, summary_age)
unlink("anage_data.txt")


## ----sample-row---------------------------------------------------------
# Because we haven't cleaned and format dataset yet, I'm temporarily
# showing Mammal instead of Mammalia from original dataset
# just to simulate the final format
age %>%
  filter(`Common name` == "Bowhead whale") %>%
  select(
    Class,
    `Common name`,
    `Maximum longevity (yrs)`,
    `Adult weight (g)`,
    `Female maturity (days)`,
    `Gestation/Incubation (days)`,
    `Litter/Clutch size`
  ) %>%
  mutate(Class = "Mammal")


# Data cleaning **********************************************************


## ----clean-row, include=FALSE-------------------------------------------
age <- age %>%
  filter(
    # keep birds and mammals, remove all other classes
    Class %in% c("Aves", "Mammalia") &
      # keep captivity and wild, remove unknown origin
      `Specimen origin` %in% c("captivity", "wild") &
      # keep small, medium and huge, remove tiny
      `Sample size` != "tiny" &
      # keep acceptable and high, remove low and questionable
      `Data quality` %in% c("acceptable", "high")
  ) %>%
  mutate(
    # now the renaming happens
    Class = factor(
      Class,
      levels = c("Aves", "Mammalia"),
      labels = c("Bird", "Mammal")
    ))


## ----clean-col-step1, include=FALSE-------------------------------------
# keep information: class, name
# keep outcome: max-longevity
# keep 5 candidate predictors
age <- age %>%
  # At this stage don't require rows to have no NAs for the candidates
  # Since if a candidate is dropped at the end,
  # we would have required more than we need, which leads to fewer rows/data
  filter(
    !is.na(`Maximum longevity (yrs)`) # outcome definitely can't be NA
  ) %>%
  select(
    Class, # Bird or Mammal
    `Common name`, # species name, e.g. Golden eagle
    `Maximum longevity (yrs)`, # highest num of years a species can live
    `Adult weight (g)`, # adult weight in grams
    `Female maturity (days)`, # female age of sexual maturity
    `Male maturity (days)`, # male age of sexual maturity
    `Gestation/Incubation (days)`, # pregnancy time
    `Litter/Clutch size` # number of young
  )


## ----clean-col-step2-nzv------------------------------------------------
# Eliminate predictors with near-zero variance, otherwise models may crash
# Luckily all 5 predictors pass this test.
# Remove index 1:3 because predictors start at col 4
nearZeroVar(age[, -(1:3)], saveMetrics = TRUE)


## ----clean-col-step2-cor------------------------------------------------
# Remove highly correlated predictors to improve model performance
if (!require(corrplot)) {
  install.packages("corrplot", repos = "http://cran.us.r-project.org")
}
library(corrplot)

# The correlation matrix
cor_matrix <- cor(
  age[, -(1:3)],
  use = "complete.obs" # Ignores cells with NAs for now
)

# Temp shorten variable names for the plot
dimnames(cor_matrix) <- list(
  c("wt", "fmat", "mmat", "ges", "lit"),
  c("wt", "fmat", "mmat", "ges", "lit")
)

# Correlation plot has cor coefficients, and do some pretty styling
corrplot(
  cor_matrix, method = "shade", shade.col = NA,
  cl.pos = "n", tl.col = "black", tl.srt = 45, diag = FALSE,
  addCoef.col = "white", col = COL2(diverging = "PiYG", n = 2))

dimnames(cor_matrix) <- list( # rename variables back
  names(age[,-(1:3)]),
  names(age[,-(1:3)])
)

# Which predictor is highly correlated, define cutoff to be 0.75? Male maturity
# Can also see same information in plot above
highly_cor_vars <- findCorrelation(cor_matrix, cutoff = .75, names = TRUE, verbose = TRUE)
highly_cor_vars


## ----clean-col-step3, include=FALSE-------------------------------------
# Now that we remove male maturity, and have finalized 4 predictors
# we require each row to not have any NA for the 4 final predictors
age <- age %>%
  select(-all_of(highly_cor_vars)) %>%
  filter(
    !is.na(`Adult weight (g)`),
    !is.na(`Female maturity (days)`),
    !is.na(`Gestation/Incubation (days)`),
    !is.na(`Litter/Clutch size`)
  ) %>%
  mutate( # Make these 2 in correct integer type, from double
    `Female maturity (days)` =
      as.integer(`Female maturity (days)`),
    `Gestation/Incubation (days)` =
      as.integer(`Gestation/Incubation (days)`),
  )

rm(cor_matrix, highly_cor_vars)


## ----class-summary------------------------------------------------------
# How many birds and mammals do we have?
age %>%
  group_by(Class) %>%
  count()


# Data splitting *********************************************************


## ----split-data, include=FALSE------------------------------------------
# Want 3 datasets: training, test, validation. Do 2 80%-20% splits.
# training set: develop models; test set: evaluate performance;
# validation set: ONLY used in very end, evaluate final rmses
set.seed(1)
training_index <- # 80% available to us, 20% validation
  createDataPartition(
    age$`Maximum longevity (yrs)`,
    times = 1,
    p = 0.8,
    list = FALSE
  )

available_set <- age[training_index, ]
validation_set <- age[-training_index, ]

# train: 80%*80%=64% original, test: 20%*80%=16%
set.seed(1)
training_index <-
  createDataPartition(
    available_set$`Maximum longevity (yrs)`,
    times = 1,
    p = 0.8,
    list = FALSE
  )

train_set <- available_set[training_index, ]
test_set <- available_set[-training_index, ]

rm(training_index, age, available_set)


# Data visualization *****************************************************


## ----summary-longevity--------------------------------------------------
# 5-number summary of max-longevity
quantile(train_set$`Maximum longevity (yrs)`, c(0, 0.25, 0.5, 0.75, 1))


## ----plot-longevity-----------------------------------------------------
# Boxplot of max-longevity, by class
train_set %>%
  ggplot(aes(Class, `Maximum longevity (yrs)`, fill = Class)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") + # no fill legend
  ggtitle("Max-longevity")


## ----min-max-longevity--------------------------------------------------
# Longest lifespan animals in each class, slice max
train_set %>%
  group_by(Class) %>%
  slice_max(order_by = `Maximum longevity (yrs)`, n = 1, with_ties = FALSE)

# Shortest lifespan animals in each class, slice min
train_set %>%
  group_by(Class) %>%
  slice_min(order_by = `Maximum longevity (yrs)`, n = 1, with_ties = FALSE)


## ----plot-wt------------------------------------------------------------
# Density plot of adult weight, by class
train_set %>%
  ggplot(aes(`Adult weight (g)`, fill = Class)) +
  geom_density(position = "identity", alpha = 0.5) +
  geom_vline( # dotted lines to highlight where density peaks are
    xintercept = c(25, 500, 7000), alpha = 0.5, linetype = 2
  ) +
  scale_x_log10() +
  annotation_logticks(sides = "b") + # a log ruler on bottom
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Adult weight")


## ----plot-mat-ges-lit---------------------------------------------------
# Scatter plot of pregnancy time vs. maturity and num of young, by class
train_set %>%
  ggplot(aes(`Gestation/Incubation (days)`, `Female maturity (days)`)) +
  geom_point(aes(color = Class, size = `Litter/Clutch size`)) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() +
  scale_color_brewer(palette = "Set1") +
  labs(size = paste( # breaks name from 1 line to 2 lines
    str_wrap("Litter/Clutch size", width = 5),
    collapse = "\n")) +
  ggtitle("Female maturity vs. pregnancy time and num of young")


## ----plot-longevity-ges-lit---------------------------------------------
# Scatter plot of max-longevity vs. pregnancy time and num of young, by class
# Added density contours showing point dense areas, and lm line
train_set %>%
  ggplot(aes(`Gestation/Incubation (days)`, `Maximum longevity (yrs)`)) +
  geom_point(aes(color = Class, size = `Litter/Clutch size`)) +
  geom_smooth(method = "lm", color = "gold") +
  geom_density_2d(color = "lavender", alpha = 0.5, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(color = "white") +
  scale_color_brewer(palette = "Pastel1") +
  # Can't see dark dots on dark theme, change to light color
  guides(size = guide_legend(override.aes = list(color = "lavender"))) +
  labs(size = paste(
    str_wrap("Litter/Clutch size", width = 5),
    collapse = "\n")) +
  theme_dark() +
  ggtitle("Max-longevity vs. pregnancy time and num of young")


## ----plot-longevity-mat-wt----------------------------------------------
# Scatter plot of max-longevity vs. adult weight and female maturity
# Added density contours showing point dense areas, and lm line
train_set %>%
  # log to distinguish colors apart, can see easily
  ggplot(aes(`Adult weight (g)`, `Maximum longevity (yrs)`,
             color = log10(`Female maturity (days)`))) +
  geom_point() +
  geom_smooth(method = "lm", color = "gold") +
  geom_density_2d(color = "lemonchiffon", alpha = 0.5, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(color = "white") +
  scale_color_distiller(palette = "YlGnBu") +
  labs(
    color = paste(
      str_wrap("log of Female maturity (days)", width = 5),
      collapse = "\n")
  ) +
  theme_dark() +
  ggtitle("Max-longevity vs. adult weight and female maturity")


## ----caret-format, include=FALSE----------------------------------------
# Stick with caret conventions, put predictor columns into a matrix named x
convert_matrix <- function(data) {
  as.matrix(
    data,
    dimnames = list(
      NULL,
      c(`Adult weight (g)`, `Female maturity (days)`,
        `Gestation/Incubation (days)`, `Litter/Clutch size`)))
}

# Outcome renamed to y
train_set <- train_set %>%
  rename(y = `Maximum longevity (yrs)`) %>%
  mutate(x = convert_matrix(train_set[, -(1:3)])) %>%
  select(x, y, Class, `Common name`)

test_set <- test_set %>%
  rename(y = `Maximum longevity (yrs)`) %>%
  mutate(x = convert_matrix(test_set[, -(1:3)])) %>%
  select(x, y, Class, `Common name`)

validation_set <- validation_set %>%
  rename(y = `Maximum longevity (yrs)`) %>%
  mutate(x = convert_matrix(validation_set[, -(1:3)])) %>%
  select(x, y, Class, `Common name`)

rm(convert_matrix)


## ----parallel, include=FALSE--------------------------------------------
# Parallel speeds up model development
if (!require(doParallel)) {
  install.packages("doParallel", repos = "http://cran.us.r-project.org")
}
library(doParallel)

cl <- makePSOCKcluster(4)
registerDoParallel(cl)


# Linear methods *********************************************************
# linear methods generate linear equations to predict animals y_hat


## ----download-linear, include=FALSE-------------------------------------
# Download linear methods libraries
if (!require(Cubist)) {
  install.packages("Cubist", repos = "http://cran.us.r-project.org")
}
library(Cubist)


## ----cubist-train, echo=TRUE--------------------------------------------
# Select the best tuning combo which has the min training rmse
# The fit is created based on best combo
# Cubist generate many linear models, each model has many linear equations
set.seed(1)
fit <- train(
  x = train_set$x,
  y = train_set$y,
  method = "cubist",
  tuneGrid = expand.grid( # expand grid does cross join on tuning params
    committees = seq(1, 9, 2), # number of models
    neighbors = seq(5, 9, 1) # used to adjust equations
  )
)


## ----cubist-plot--------------------------------------------------------
# Save the best tuning combo for later use
fit$bestTune
best_tunes <- list(cubist = fit$bestTune)

# Plot the full tuning result, training rmses of all combos
ggplot(fit, highlight = TRUE) +
  scale_x_continuous(breaks = seq(1, 9, 2)) +
  ggtitle("Tuning result cubist")


## ----cubist-predict-----------------------------------------------------
# Predict test set outcomes, save the list of y_hat
set.seed(1)
linear_models <- tibble(
  cubist = predict(fit, test_set$x)
)

# Calculate rmse, save it to compare with other methods later
rmses_linear <- tibble(
  method = "cubist",
  rmse = RMSE(linear_models$cubist, test_set$y)
)
rmses_linear

rm(fit)


## ----lm-----------------------------------------------------------------
# lm makes one linear equation to fit all data
# No tuning parameters for lm
set.seed(1)
fit <- train(
  x = train_set$x,
  y = train_set$y,
  method = "lm"
)

set.seed(1)
linear_models <- linear_models %>% mutate(
  lm = predict(fit, test_set$x)
)

# The coefficients for predictors, which show what the linear equation is
fit$finalModel$coefficients

rmses_linear <- rmses_linear %>% add_row(
  method = "lm",
  rmse = RMSE(linear_models$lm, test_set$y)
)
rmses_linear

rm(fit)


## ----linear-ensemble----------------------------------------------------
# Create ensemble by taking the mean of cubist and lm
linear_models <- linear_models %>%
  rowwise() %>% # take one row at a time
  mutate(
    ensemble = mean(c(cubist, lm)) # take the mean
  ) %>%
  ungroup() # put rows back together

rmses_linear <- rmses_linear %>% add_row(
  method = "ensemble",
  rmse = RMSE(linear_models$ensemble, test_set$y)
)

# Rank methods based on performance
rmses_linear <- rmses_linear %>% arrange(rmse)
rmses_linear


## ----linear-plot--------------------------------------------------------
# Residuals plot compare method performance for cubist, lm, and ensemble
# residual = true outcome - predicted outcome
# r = 0 perfect prediction, r > 0 underestimate, r < 0 overestimate
# From left to right, methods subplots arranged from low to high rmse
linear_models %>%
  transmute( # transmute is mutate + only select the specified columns
    true = test_set$y,
    r_cubist = true - cubist,
    r_lm = true - lm,
    r_ensemble = true - ensemble
  ) %>%
  pivot_longer( # create method and residual column, and spreads the data
    cols = -true,
    names_to = "method",
    names_pattern = "r_(.*)", # remove the r_ prefix from r_cubist, etc.
    values_to = "residual"
  ) %>%
  # Neat trick is to use rmses_linear$method which we ranked already
  mutate(method = factor(
    method,
    levels = rmses_linear$method,
    ordered = TRUE
  )) %>%
  ggplot(aes(true, residual)) +
  # ndensity to give each subplot a total density of 1, which is correct
  # otherwise the entire plot has a total density of 1, which is wrong
  geom_density_2d_filled(contour_var = "ndensity") +
  geom_hline(yintercept = 0, size = 1, alpha = 0.75, color = "white") +
  geom_point(color = "white", alpha = 0.25) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  theme_dark() +
  ggtitle("Linear models residuals") +
  facet_wrap(method~.)


## ----linear-mistakes----------------------------------------------------
# Top 3 residual mistakes made by the ensemble
test_set %>%
  mutate(
    residual = y - linear_models$ensemble
  ) %>%
  select(Class, `Common name`, y, residual) %>%
  slice_max(order_by = abs(residual), n = 3)

rm(linear_models, rmses_linear)


# Non-linear methods *****************************************************
# non-linear methods group animals, and animals in same group have same y_hat


## ----download-nonlinear, include=FALSE----------------------------------
# Download non-linear methods libraries
if (!require(randomForest)) {
  install.packages("randomForest", repos = "http://cran.us.r-project.org")
}
if (!require(rpart)) {
  install.packages("rpart", repos = "http://cran.us.r-project.org")
}
if (!require(ggdendro)) {
  install.packages("ggdendro", repos = "http://cran.us.r-project.org")
}
library(randomForest)
library(rpart)
library(ggdendro) # used to plot a decision tree


## ----rpart-train--------------------------------------------------------
# rpart partitions to make a decision tree
# cp: if rpart can't improve by cp amount by doing the next partition, it stops
set.seed(1)
fit <- train(
  x = train_set$x,
  y = train_set$y,
  method = "rpart",
  tuneGrid = data.frame(cp = seq(0.01, 0.1, 0.01))
)

fit$bestTune

best_tunes <- c(best_tunes, list(rpart = fit$bestTune))

ggplot(fit, highlight = TRUE) +
  scale_x_continuous(breaks = seq(0.02, 0.1, 0.02)) +
  labs(title = "Tuning result rpart")

set.seed(1)
nonlinear_models <- tibble(
  rpart = predict(fit, test_set$x)
)

rmses_nonlinear <- tibble(
  method = "rpart",
  rmse = RMSE(nonlinear_models$rpart, test_set$y)
)
rmses_nonlinear


## ----rpart-plot---------------------------------------------------------
# rpart generates a decision tree, plot it
data <- dendro_data(fit$finalModel)
# shorten variable names, otherwise plot won't fit
data$labels$label <- str_replace_all(data$labels$label, "Fem.*\\)", "mat ")
data$labels$label <- str_replace_all(data$labels$label, "Adult.*\\)", "wt ")
data$labels$label <- str_replace_all(data$labels$label, "Ges.*\\)", "ges ")
data$labels$label <- str_replace_all(data$labels$label, "Lit.*\\)", "lit ")
data$leaf_labels$label <- as.numeric(data$leaf_labels$label) # correct the type

ggplot() +
  geom_segment( # draw branches
    data = segment(data),
    aes(x = x, y = y, xend = xend, yend = yend), alpha = 0.5
  ) +
  geom_label( # draw nodes: splits on predictor, e.g. weight < 200
    data = label(data),
    aes(x = x, y = y, label = label), fill = "lightskyblue"
  ) +
  geom_label( # draw leafs: predicted outcome, e.g. 10 years
    data = leaf_label(data),
    aes(x = x, y = y, label = label,), fill = "yellowgreen"
  ) +
  theme_dendro() +

  ggtitle("rpart tree (training set)")

rm(fit)


## ----rf-----------------------------------------------------------------
# rf randomly choose predictors to partition, mtry is number of predictors
# rf suggests that mtry = sqrt(num of predictors) = sqrt(4) = 2,
# so take mtry = 2 and won't tune other values
set.seed(1)
fit <- train(
  x = train_set$x,
  y = train_set$y,
  method = "rf",
  tuneGrid = data.frame(mtry = 2)
)

set.seed(1)
nonlinear_models <- nonlinear_models %>% mutate(
  rf = predict(fit, test_set$x)
)

rmses_nonlinear <- rmses_nonlinear %>% add_row(
  method = "rf",
  rmse = RMSE(nonlinear_models$rf, test_set$y)
)
rmses_nonlinear

rm(fit)


## ----nonlinear-ensemble-------------------------------------------------
# Create non-linear ensemble by taking the mean of rf and rpart
nonlinear_models <- nonlinear_models %>%
  rowwise() %>%
  mutate(ensemble = mean(c(rf, rpart))) %>%
  ungroup()

rmses_nonlinear <- rmses_nonlinear %>% add_row(
  method = "ensemble",
  rmse = RMSE(nonlinear_models$ensemble, test_set$y)
)

# arrange models from low to high rmse
rmses_nonlinear <- rmses_nonlinear %>% arrange(rmse)
rmses_nonlinear


## ----nonlinear-plot-----------------------------------------------------
# Residual plot for non-linear methods, this is just like the residual
# plot done before for linear methods, so will skip comments here
nonlinear_models %>%
  transmute(
    true = test_set$y,
    r_rf = true - rf,
    r_rpart = true - rpart,
    r_ensemble = true - ensemble
  ) %>%
  pivot_longer(
    cols = -true,
    names_to = "method",
    names_pattern = "r_(.*)",
    values_to = "residual"
  ) %>%
  mutate(method = factor(
    method,
    levels = rmses_nonlinear$method,
    ordered = TRUE
  )) %>%
  ggplot(aes(true, residual)) +
  geom_density_2d_filled(contour_var = "ndensity") +
  geom_hline(yintercept = 0, size = 1, alpha = 0.75, color = "white") +
  geom_point(color = "white", alpha = 0.25) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  theme_dark() +
  ggtitle("Non-linear models residuals") +
  facet_wrap(method~.)


## ----nonlinear-mistakes-------------------------------------------------
# Top 3 residual mistakes of non-linear ensemble
test_set %>%
  mutate(
    residual = y - nonlinear_models$ensemble
  ) %>%
  select(Class, `Common name`, y, residual) %>%
  slice_max(order_by = abs(residual), n = 3)

rm(nonlinear_models, rmses_nonlinear)


# Final 2 models *********************************************************
# Create final 2 ensembles using all data we are allowed to use,
# available set = train + test set, which is all data except validation set


## ----final-linear, include=FALSE----------------------------------------
# form the available set by putting train + test together
available_set <- rbind(train_set, test_set)
rm(train_set, test_set)

set.seed(1)
fit <- train(
  x = available_set$x,
  y = available_set$y,
  method = "cubist",
  tuneGrid = best_tunes$cubist # use the stored best tuning combo
)

set.seed(1)
final_models <- tibble(cubist = predict(fit, validation_set$x))

cubist <- fit$finalModel
save(cubist, file = "cubist.rda") # save linear models details for appendix
rm(fit, cubist)

set.seed(1)
fit <- train(
  x = available_set$x,
  y = available_set$y,
  method = "lm"
)

set.seed(1)
final_models <- final_models %>% mutate(
  lm = predict(fit, validation_set$x)
)

lm <- fit$finalModel$coefficients
save(lm, file = "lm.rda") # save coefficients of linear eq for appendix
rm(fit, lm)

# form the linear ensemble
final_models <- final_models %>%
  rowwise() %>%
  mutate(linear = mean(c(lm, cubist))) %>%
  ungroup() %>%
  select(-c(lm, cubist))


## ----final-nonlinear, include=FALSE-------------------------------------
set.seed(1)
fit <- train(
  x = available_set$x,
  y = available_set$y,
  method = "rpart",
  tuneGrid = best_tunes$rpart
)

set.seed(1)
final_models <- final_models %>% mutate(
  rpart = predict(fit, validation_set$x)
)

rpart <- fit$finalModel
save(rpart, file = "rpart.rda") # save to plot decision tree for appendix
rm(fit, rpart)

set.seed(1)
fit <- train(
  x = available_set$x,
  y = available_set$y,
  method = "rf",
  tuneGrid = data.frame(mtry = 2)
)

set.seed(1)
final_models <- final_models %>% mutate(
  rf = predict(fit, validation_set$x)
)

rm(fit)

# form the non-linear ensemble
final_models <- final_models %>%
  rowwise() %>%
  mutate(nonlinear = mean(c(rpart, rf))) %>%
  ungroup() %>%
  select(-c(rpart, rf))

rm(best_tunes, available_set)

# we won't need parallel computing anymore now that we're done
stopCluster(cl)
rm(cl)


# Results ****************************************************************


## ----final-rmse---------------------------------------------------------
# the final rmses
rmses_final <- tibble(
  method = c("linear", "nonlinear"),
  rmse = c(
    RMSE(final_models$linear, validation_set$y),
    RMSE(final_models$nonlinear, validation_set$y)
  ))

# rank methods based on rmse from low to high
rmses_final <- rmses_final %>% arrange(rmse)
rmses_final


## ----final-boxplot------------------------------------------------------
# 5-number summaries for true outcomes y, linear y_hat, and non-linear y_hat
summary_outcome <- sapply(
  c(validation_set %>% transmute(true = y), final_models), # rename y to true
  function(o) summary(o))
as_tibble(t(summary_outcome), rownames = "Source") # formatting
rm(summary_outcome)

# Boxplot of final ensembles models pred vs. true max-longevity
# From left to right, want methods subplots arranged from low to high rmse
# true outcomes are not predictions, but still put it on the left since 0 error
final_models %>%
  mutate(
    true = validation_set$y # add true outcomes to our dataframe
  ) %>%
  pivot_longer( # create 2 columns, source and max-longevity, spread the data
    cols = everything(),
    names_to = "source",
    values_to = "max_longevity"
  ) %>%
  mutate(source = factor(
    source, # since true is not inside rmses_final$method, add it to the factor
    levels = c("true", rmses_final$method),
    ordered = TRUE
  )) %>%
  ggplot(aes(source, max_longevity, fill = source)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  ggtitle("Final ensemble models pred vs. true max-longevity")


## ----final-residual-----------------------------------------------------
# Residuals plot facet by class and colored by method
# Want to compare how models did for birds vs. mammals
# Want to compare linear vs. non-linear ensemble
validation_set %>%
  mutate(
    true = y,
    r_linear = true - final_models$linear,
    r_nonlinear = true - final_models$nonlinear
  ) %>%
  select(Class, true, r_linear, r_nonlinear) %>%
  pivot_longer(
    cols = c(r_linear, r_nonlinear),
    names_to = "method",
    names_pattern = "r_(.*)",
    values_to = "residual"
  ) %>%
  mutate(method = factor(
    method,
    levels = rmses_final$method,
    ordered = TRUE
  )) %>%
  ggplot(aes(true, residual, color = method)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = method), se = FALSE) + # a loess line
  geom_hline(yintercept = 0, size = 1, alpha = 0.75) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  scale_color_brewer(palette = "Pastel1") +
  theme_dark() +
  ggtitle("Final ensemble models residuals by class") +
  facet_wrap(Class~.)


## ----final-mistakes-----------------------------------------------------
# Top 5 mistakes made by linear ensemble
validation_set %>%
  mutate(
    linear = final_models$linear,
    residual = y - linear
  ) %>%
  select(Class, `Common name`, y, linear, residual) %>%
  slice_max(order_by = abs(residual), n = 5)

# Top 5 mistakes made by non-linear ensemble
validation_set %>%
  mutate(
    nonlinear = final_models$nonlinear,
    residual = y - nonlinear
  ) %>%
  select(Class, `Common name`, y, nonlinear, residual) %>%
  slice_max(order_by = abs(residual), n = 5)


## ----cleanup, include=FALSE---------------------------------------------
rm(final_models, rmses_final, validation_set)


# Appendix ***************************************************************


## ----summary_age, echo=FALSE--------------------------------------------
# 5-number summary and number of NAs for original age dataset numerical cols
load("summary_age.rda")
summary_age

rm(summary_age)
unlink("summary_age.rda")


## ----summary-lm---------------------------------------------------------
# final lm coefficients
load("lm.rda")
lm

rm(lm)
unlink("lm.rda")


## ----summary-rpart------------------------------------------------------
# final rpart decision tree
load("rpart.rda")
data <- dendro_data(rpart)
data$labels$label <- str_replace_all(data$labels$label, "Fem.*\\)", "mat ")
data$labels$label <- str_replace_all(data$labels$label, "Adult.*\\)", "wt ")
data$labels$label <- str_replace_all(data$labels$label, "Ges.*\\)", "ges ")
data$labels$label <- str_replace_all(data$labels$label, "Lit.*\\)", "lit ")
data$leaf_labels$label <- as.numeric(data$leaf_labels$label)

ggplot() +
  geom_segment(
    data = segment(data),
    aes(x = x, y = y, xend = xend, yend = yend), alpha = 0.5
  ) +
  geom_label(
    data = label(data),
    aes(x = x, y = y, label = label), fill = "lightskyblue"
  ) +
  geom_label(
    data = leaf_label(data),
    aes(x = x, y = y, label = label,), fill = "yellowgreen"
  ) +
  theme_dendro() +
  ggtitle("rpart tree (final)")

rm(rpart, data)
unlink("rpart.rda")


## ----summary-cubist-----------------------------------------------------
# final cubist linear models equations
load("cubist.rda")
summary(cubist)

rm(cubist)
unlink("cubist.rda")

