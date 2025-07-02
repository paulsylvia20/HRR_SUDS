SUDs Overdose LMM Main
================
Paul Sylvia
2025-06-30

The data preparation procedures in Step 1 of the analysis yielded a
fully disaggregated dataset that captures both county- and cluster-level
measures related to drug overdose, socioeconomics, and healthcare
utilization. In line with the exploratory nature of the study, the
present analysis generates preliminary step-forward regression scripts
that enable the researcher to step in predictors at L1 or L2 and as
either fixed or random effects as specified. The resulting model is then
tested and interpreted in terms of it’s predictive power for the target
outcome variable, overdose deaths. The analysis aims to examine how
socioeconomics (at level 1) and hospital attributes (at level 2) are
related and my cause overdose deaths at level 1. While we attempt to
explain as much variabilty in the outcome as possible, we will conclude
by interpreting the primary effects as found in the Final controlled
model.

``` r
library.list <- c("dplyr", # For piping and data manipulation
                  "lme4", # For building HLM models
                  "lmerTest", # For what?
                  "ggplot2",
                  "car" #For assessing multicolinearity
                  )

#Here we load packages from libary.list
for (i in 1:length(library.list)) {
  if (!library.list[i] %in% rownames(installed.packages())) {
    install.packages(library.list[i])
  }
  library(library.list[i], character.only = TRUE)}

Final <- read.csv("Final.csv", )[-1]


Final$population <- as.numeric(scale(Final$population))
Final$population_grand <- as.numeric(scale(Final$population_grand))
```

### Examination of L1 relationships

In any analysis, it is best practice to investigate the relational
structure of predictors and avoid including variables which introduce
multicolinearity or limit desirable suppression effects. Below, we see
strong correlations between socioeconomic predictors like median income,
poverty rate, and social capital. However, we considered this
acceptable, especially as we progress using a step-up procedure. While
strong relationships among predictors are not desirable, they are not a
major concern and we can see that some of these even have opposite
relationships to the outcome variable (see the simple regression
figures). For this reason we do not recommend the exclusion of any
variables.

``` r
library(ggcorrplot)
ggcorrplot(cor(Final[3:14]),
           type = "upper",
           lab = TRUE,
           lab_size = 2,
           colors = c("blue", "white", "red"),
           title = "Level 1 Corr Matrix",
           tl.cex = 6) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

ggcorrplot(cor(Final[27:42]),
           type = "upper",
           lab = TRUE,
           lab_size = 2,
           colors = c("blue", "white", "red"),
           title = "Level 2 Correlation Matrix",
           tl.cex = 6) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))
```

<img src="HLM_analysis_files/figure-gfm/initial L1 examination 1-1.png" width="50%" /><img src="HLM_analysis_files/figure-gfm/initial L1 examination 1-2.png" width="50%" />

``` r
L1_pred <- Final[-4][3:14]
L1_out <- log(Final[4])

for (i in 1:11){
  graph_dat <- cbind(out = L1_out, L1_pred)
  graph <- ggplot(data = graph_dat, aes(y=.data[[names(graph_dat)[1]]], x=.data[[names(graph_dat)[i+1]]])) +
    geom_point() +
    geom_smooth(method = lm)
    
  assign(paste("graph_", i, sep=""), graph)
  }

library("ggpubr")
print(ggarrange(graph_1, graph_2, graph_3, graph_4, graph_5, graph_6, 
          graph_7, graph_8, graph_9, graph_10, graph_11, ncol = 3))
```

<img src="HLM_analysis_files/figure-gfm/initial L1 examination 2-1.png" style="display: block; margin: auto;" /><img src="HLM_analysis_files/figure-gfm/initial L1 examination 2-2.png" style="display: block; margin: auto;" /><img src="HLM_analysis_files/figure-gfm/initial L1 examination 2-3.png" style="display: block; margin: auto;" /><img src="HLM_analysis_files/figure-gfm/initial L1 examination 2-4.png" style="display: block; margin: auto;" />

### Outcome transformation

Another challenge faced in this study is that the chosen outcome
variable measures a somewhat rare occurrence and exhibits a strong floor
effect. To remedy this, we will subject it to log transformation. In the
figures below we see how the residuals change from being oddly shaped
with extreme outliers (left) to classically and desirably shaped like a
normal data cloud (right).

<img src="HLM_analysis_files/figure-gfm/Transform Check-1.png" style="display: block; margin: auto;" />

As a result, we find a few important values to be carried into the
remaining analysis. First, we find an unconditional grand mean
(gamma_00) of -6.7 which indicates that the average county across all
clusters is predicted to have ~2 overdose deaths per 10 million
individuals in given year. This value is statistically significantly
greater than zero. Moreover, we have found an unconditional interclass
correlation coefficient (ICC) of 3.3% which indicates the proportion of
total variability that exists between clusters. This is a rather low ICC
which, at the outset, implies that HRR’s do not have a strong aggregate
effect on surrounding county’s overdose deaths, although the aggregate
effects do appear to be non-zero. This number tells us how much
variability is available for our analysis to explain.

    ## Gamma_00 : -6.694953  p-value =  1.77387e-170

    ## Total L1 variance remaining:  3.322391

    ## Total L2 variance remaining:  0.1133314

    ## ICC:  0.03298619

### Building the stepwise function

The next step will define a function for building and extracting the
results of HLM models in a flexible, stepwise manner. We will apply this
function to conduct a stepwise regression procedure based on deviance
and reliability. We will build a model that has the best fit to the
data, while balancing against model reliability. In HLM, reliability is
an indication of a models capacity to capture the typical individual. A
model which is reliable does not tend to seriously mispredict the
observed data. This is distinct from model fit because a model with
strong fit may do so at the cost of failing for a smaller subset of the
data.

Please note that I have utilized ChatGPT to add comments and improve
readability of the script. But all coding (and writing) has been
conducted without the assistance of LLM’s.

``` r
#' Compare function to evaluate the effect of adding new variables to a multilevel model
#' Supports both fixed effects and random slopes
#' Returns a summary matrix with model comparison metrics

compare <- function(prev_mod, vars, fit.random.slopes = FALSE, random.intercept = TRUE) {
  
  # --- 1. Identify variables in the previous model ---
  prev_vars_fixed <- row.names(summary(prev_mod)$coefficients)[-1]  # Fixed effects (excluding intercept)
  prev_vars_random <- row.names(summary(prev_mod)$varcor$HRR)[-1]   # Random slope terms
  
  # --- 2. Determine candidate variables to step in (for fixed effects) ---
  if (length(prev_vars_fixed) == 0) {
    new_vars_fixed <- names(vars)  # If no previous fixed effects, use all candidate vars
  } else {
    if (sum(is.na(match(prev_vars_fixed, names(vars)))) == length(prev_vars_fixed)) {
      new_vars_fixed <- names(vars)
    } else {
      new_vars_fixed <- names(vars)[-na.omit(match(prev_vars_fixed, names(vars)))]  # Exclude already used vars
    }
  }

  # --- 3. Determine candidate variables for random slopes ---
  if (length(prev_vars_random) == 0) {
    new_vars_random <- prev_vars_fixed  # If no prior random slopes, use fixed terms
  } else {
    new_vars_random <- prev_vars_fixed[-match(prev_vars_random, prev_vars_fixed)]  # Exclude already random
  }

  # --- 4. Run model comparisons (Fixed effects only case) ---
  if (fit.random.slopes == FALSE) {
    mod_change <- matrix(ncol = 7, nrow = length(new_vars_fixed))
    colnames(mod_change) <- c("Chi_sq", "p_val", "reliability", "L1_expl", "L2_expl", "Singular?", "Warnings")
    rownames(mod_change) <- new_vars_fixed

    for (i in new_vars_fixed) {
      # Add one new predictor at a time to the model
      f <- as.formula(paste(". ~ . +", i))
      mod <- update(prev_mod, formula = f)

      # Extract variance components for reliability & explained variance
      tau_00 <- as.data.frame(VarCorr(mod))[1, 4]     # L2 variance (intercept)
      sigma2 <- as.data.frame(VarCorr(mod))[2, 4]     # L1 residual variance
      rel <- tau_00 / (tau_00 + sigma2 / 11.7)         # Estimate of reliability

      # Model comparison with likelihood ratio test
      test <- suppressMessages(anova(mod, prev_mod))

      # Proportion of explained variance
      expl_L1 <- (sigma2_null - sigma2) / sigma2_null
      expl_L2 <- (tau_00_null - tau_00) / tau_00_null

      # Store results in matrix
      mod_change[i, ] <- c(
        signif(test$Chisq[2], 5),
        signif(test$`Pr(>Chisq)`[2], 5),
        signif(rel, 5),
        signif(expl_L1, 5),
        signif(expl_L2, 5),
        as.character(isSingular(mod)),
        "None"  # Placeholder, replaced below if needed
      )

      # Capture and clean warnings
      temp <- mod@optinfo$conv$lme4$messages
      if (length(temp > 0)) {
        ind <- as.numeric(lapply(temp, function(x) x != "boundary (singular) fit: see help('isSingular')"))
        mod_change[i, "Warnings"] <- if (sum(ind) > 0) paste(temp[ind], collapse = ", ") else "None"
      }
    }
  }

  # --- 5. Run model comparisons (Random slopes version) ---
  else if (fit.random.slopes == TRUE) {
    mod_change <- matrix(ncol = 7, nrow = length(new_vars_random))
    colnames(mod_change) <- c("Chi_sq", "p_val", "reliability", "L1_expl", "L2_expl", "Singular?", "Warnings")
    rownames(mod_change) <- new_vars_random

    for (i in new_vars_random) {
      # Build model formula with new random slope
      fixed_f <- paste(prev_vars_fixed, collapse = " + ")
      random_f <- paste(c(prev_vars_random, i), collapse = " + ")
      f <- paste("overdose_t ~ 1 +", fixed_f, "+ (1 +", random_f, "| HRR)")
      dat <- prev_mod@call$data

      # Fit model quietly
      suppressWarnings(
        suppressMessages(
          mod <- eval(bquote(lmer(as.formula(.(f)), data = .(dat))))
        )
      )

      # Extract variances and reliability
      VarCorr <- as.data.frame(VarCorr(mod))
      tau_00 <- VarCorr[1, 4]
      sigma2 <- VarCorr[nrow(VarCorr), 4]
      rel <- tau_00 / (tau_00 + sigma2 / 11.7)

      # Likelihood ratio test
      test <- suppressMessages(anova(mod, prev_mod))

      # Calculate explained variances
      expl_L1 <- (sigma2_null - sigma2) / sigma2_null
      expl_L2 <- (tau_00_null - tau_00) / tau_00_null

      # Store results
      mod_change[i, ] <- c(
        signif(test$Chisq[2], 5),
        signif(test$`Pr(>Chisq)`[2], 5),
        signif(rel, 5),
        signif(expl_L1, 5),
        signif(expl_L2, 5),
        as.character(isSingular(mod)),
        "None"
      )

      # Clean and store warnings
      temp <- mod@optinfo$conv$lme4$messages
      if (length(temp > 0)) {
        ind <- as.numeric(lapply(temp, function(x) x != "boundary (singular) fit: see help('isSingular')"))
        mod_change[i, "Warnings"] <- if (sum(ind) > 0) paste(temp[ind], collapse = ", ") else "None"
      }
    }
  }

  # --- 6. Return the summary matrix ---
  return(mod_change)
}
```

*Example Output*

    ##                      Chi_sq    p_val        reliability L1_expl      
    ## population_grand     "251.73"  "1.091e-56"  "0.17528"   "0.067311"   
    ## hispanic_grand       "6.259"   "0.012356"   "0.34584"   "0.0079815"  
    ## older_grand          "145.88"  "1.3758e-33" "0.31567"   "0.048266"   
    ## younger_grand        "30.439"  "3.4462e-08" "0.32682"   "0.013618"   
    ## black_grand          "9.3053"  "0.0022849"  "0.26978"   "0.0011064"  
    ## rural_grand          "265.33"  "1.1811e-59" "0.25585"   "0.078394"   
    ## gini_grand           "0.63122" "0.42691"    "0.28624"   "-1.9417e-05"
    ## med.inc_grand        "104.86"  "1.31e-24"   "0.28068"   "0.032233"   
    ## employment_grand     "2.461"   "0.11671"    "0.28475"   "0.00041498" 
    ## social_capital_grand "12.493"  "0.00040841" "0.30465"   "0.0056378"  
    ## poverty_rate_grand   "20.134"  "7.2211e-06" "0.29504"   "0.0070869"  
    ##                      L2_expl      Singular? Warnings
    ## population_grand     "0.50332"    "FALSE"   "None"  
    ## hispanic_grand       "-0.31407"   "FALSE"   "None"  
    ## older_grand          "-0.099999"  "FALSE"   "None"  
    ## younger_grand        "-0.1999"    "FALSE"   "None"  
    ## black_grand          "0.075305"   "FALSE"   "None"  
    ## rural_grand          "0.20606"    "FALSE"   "None"  
    ## gini_grand           "-0.0048327" "FALSE"   "None"  
    ## med.inc_grand        "0.053799"   "FALSE"   "None"  
    ## employment_grand     "0.0029035"  "FALSE"   "None"  
    ## social_capital_grand "-0.091594"  "FALSE"   "None"  
    ## poverty_rate_grand   "-0.041199"  "FALSE"   "None"

## Model Building Stage 1: The Level 1 model with Fixed Effects

The step-up procedure begins by developing and applying a custom
function for comparing fixed effect options. This function, compare(),
will allow us to sequentially introduce fixed effects which minimize
deviance and therefore produce the greatest improvements to model fit at
each stage. However, occasionally I have prioritized fixed effects which
are of theoretical interest (socio-economic variables) or I have avoided
introducing fixed effects that are detrimental to overall reliability.
In the first step, as an example, I have first introduced rural_grand.
Looking at the table just above, we find it has the largest chi^2 which
corresponds to the greatest reduction in deviance in the associated
likelihood ratio test. However, in the second step I have stepped in
med.inc_grand even though older_grand had a lower deviance. This is
because I prioritized it as a variable of interest, meanwhile their
deviances were essentially identical. This procedure is somewhat
subjective, but is important for generating useful results.

``` r
compare(null, vars = predictors)
L1_step1 <- lmer(data = Final, overdose_t ~ 1 + rural_grand #stepped in rural
                   + (1 |HRR))

compare(L1_step1, vars = predictors)
L1_step2 <- lmer(data = Final, overdose_t ~ 1 + rural_grand #stepped in median income
                   + med.inc_grand + (1 |HRR))

compare(L1_step2, vars = predictors)
L1_step3 <- lmer(data = Final, overdose_t ~ 1 + rural_grand #stepped in gini coefficient
                   + med.inc_grand + gini_grand + (1 |HRR))

compare(L1_step3, vars = predictors)
L1_step4 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand #stepped in %older
                   + med.inc_grand + gini_grand + (1 |HRR))

compare(L1_step4, vars = predictors)
L1_step5 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand #stepped in employment
                   + med.inc_grand + gini_grand + employment_grand + social_capital_grand + (1 |HRR))

compare(L1_step5, vars = predictors)
L1_step6 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand #stepped in social capital
                   + med.inc_grand + gini_grand + employment_grand + social_capital_grand + (1 |HRR))

compare(L1_step6, vars = predictors)
L1_step7 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand + population_grand #stepped in population
                   + med.inc_grand + gini_grand + employment_grand + social_capital_grand + (1 |HRR))
```

*Intermediate results for Level 1 model with fixed effects*

    ## Total L1 variance remaining:  2.938011

    ## Total L2 variance remaining:  0.05414736

    ## Conditional ICC:  0.01809642

    ## % L1 variance explained:  0.1156936

    ## % L2 variance explained:  0.5222209

In continuation of the previous step, we are interested in whether any
aggregated effects predict or control the outcome. This is something
that is anticipated given the considerable L2 variance explained so far.
This generally occurs when a level 1 variable means something else in
its aggregate form. For example, socioeconomic status of a county is
distinct conceptually from socioeconomic status of an entire region.
I.e., we might expect greater overdose deaths in low socioeconomic
status counties, but not necessarily in low socioeconomic status regions
(or HRRs). We were able to step in two aggregated effects, %Hispanic and
median income, leading to an additional 30% of the L2 variance being
explained by the Level 1 model with fixed effects. At this stage, we
have found compelling evidence by the ~80% L2 variance explained that
different hospital referral networks are not associated with
considerable changes in overdose deaths which is a surprising finding
given the role of medical providers in preventing overdose on the
ground.

``` r
agg_predictors <- Final[c(30,32:41)]

compare(L1_step7, vars=agg_predictors, fit.random.slopes = FALSE)
L1_step8 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand + population_grand
                + med.inc_grand + gini_grand + employment_grand + social_capital_grand 
                + hispanic_grand.mean
                + (1 |HRR))

compare(L1_step8, vars=agg_predictors, fit.random.slopes = FALSE)
L1_step9 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand + population_grand
                + med.inc_grand + gini_grand + employment_grand + social_capital_grand 
                + hispanic_grand.mean + med.inc_grand.mean
                + (1 |HRR))
```

*Final result for Level 1 model with fixed effects*

    ## Total L1 variance remaining:  2.93618

    ## Total L2 variance remaining:  0.01933297

    ## Conditional ICC:  0.006541325

    ## % L1 variance explained:  0.1162446

    ## % L2 variance explained:  0.829412

## Model Building Stage 2: Level 1 Model with Random Effects

The following procedure investigates two questions. First, to what
extent do hospital systems alleviate the consequence of substance use
disorder and which features are associated with overdose deaths in
surrounding counties? This question will be addressed by examining fixed
effects at L2 and by testing cross-level moderation effects on
gini_coefficient (we may have tested all moderating effects in grid
search fashion, but decided instead to focus this stage of the study on
variables with theoretical importance). Second, previous research has
identified social inequality as a potential driver of substance use
related deaths. However, these findings have been inconsistent across
studies. We are interested in whether differences between hospital
system correspond to a change in relationship between inequality and
overdose deaths. In particular, do any hospital system attributes
moderate the relationship of income inequality and overdose deaths? We
apply the same grid search approach and identify only 1 random effect.

``` r
compare(L1_step9, vars=predictors, fit.random.slopes = TRUE)
L1r_step1 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand + population_grand #stepped in gini_grand as a random effect
                + med.inc_grand + gini_grand + employment_grand + social_capital_grand 
                + hispanic_grand.mean + med.inc_grand.mean
                + (1 + gini_grand |HRR))
```

*Results for the final Level 1 model with Random Effects*

    ## Total L1 variance remaining:  2.891508

    ## Total L2 variance remaining:  0.01997934

    ## Conditional ICC:  0.006862246

    ## % L1 variance explained:  0.1296907

    ## % L2 variance explained:  0.8237086

    ## Psuedo R^2 0.1712722 
    ## 

    ## Total variance in the gini coeficient slope to be explained:  6.345003

### Building the stepwise function

The final model building procedure attempts to do two things: predict
the outcome and predict the coefficient for the gini index (income
inequality) using hospital system variables. But before building the
model further we will write a new compare function that works similarly
to the previous one. The decision criteria are very similar, however,
this time we are more concerned with statistical significance of model
coefficients, removing variables that are not significant at the p \<
.05 level.

``` r
#' Compare2 function to evaluate the effect of adding predictors or interactions to a model
#' Supports main effects or moderator (interaction with gini_grand) testing
#' Returns a summary matrix with model comparison metrics and coefficient details

compare2 <- function(mod, predictors, moderators = FALSE) {

  # --- 1. Initialize storage matrix ---
  mod_change <- matrix(ncol = 4, nrow = length(predictors))
  rownames(mod_change) <- predictors
  colnames(mod_change) <- c("Chi_sq", "p_val", "coef", "coef_pval")

  # --- 2. Loop over predictors ---
  for (i in seq_along(predictors)) {

    if (moderators == FALSE) {
      # --- 2a. Main effect only ---
      
      # Construct updated formula with added predictor
      f <- as.formula(paste(". ~ . +", predictors[[i]]))
    
    } else {
      # --- 2b. Interaction term with gini_grand ---
      
      # Build interaction string
      interaction_term <- paste(predictors[[i]], "*gini_grand", sep = "")
      
      # Construct formula including both main effect and interaction
      f <- as.formula(paste(". ~ . +", predictors[[i]], "+", interaction_term))
    }

    # --- 3. Fit updated model and compare to base model ---
    new_mod <- suppressWarnings(update(mod, formula = f))
    test <- suppressWarnings(suppressMessages(anova(new_mod, mod)))

    # --- 4. Extract coefficient table and results ---
    coef_table <- summary(new_mod)$coefficients
    last_row <- nrow(coef_table)  # Assumes new term is last row

    mod_change[i, "Chi_sq"]   <- test$Chisq[2]              # Likelihood ratio test statistic
    mod_change[i, "p_val"]    <- test$`Pr(>Chisq)`[2]        # p-value for model improvement
    mod_change[i, "coef"]     <- coef_table[last_row, 1]     # Estimate of added term
    mod_change[i, "coef_pval"]<- coef_table[last_row, 5]     # p-value of added term
  }

  # --- 5. Return results matrix ---
  return(mod_change)
}
```

*Example Output*

    ##                   Chi_sq        p_val         coef    coef_pval
    ## sub_beds        5.998143 0.0498333198 2.239186e-02 0.0698357537
    ## emergency_hosp 17.583574 0.0001519762 1.058837e-01 0.0003000428
    ## medicaid_dc    16.319528 0.0002859299 5.456023e-05 0.0008491446

## Model Building Stage 3: the Level 2 Model

Now we step in any L2 variables that function either as main effects or
moderators. Initially we find no significant main effects. However, we
find one moderator with a significant coefficient and significant
likelihood ratio test:

``` r
#Rescale the L2 predictors to facilitate model fitting
Final$medicaid_dc_std <- as.numeric(scale(Final$medicaid_dc/Final$cluster_population))
Final$sub_beds_std <- as.numeric(scale(Final$sub_beds/Final$cluster_population))
Final$emergency_hosp_std <- as.numeric(scale(Final$emergency_hosp/Final$cluster_population))

L2_predictors <- Final[27:29]
```

``` r
compare2(L1r_step1, names(L2_predictors), moderators=FALSE) #no significant main effects at level 2

compare2(L1r_step1, names(L2_predictors), moderators=TRUE)
L2_step1 <- lmer(data = Final, overdose_t ~ 1 + rural_grand + older_grand + population_grand #stepped in gini_grand as a random effect
                + med.inc_grand + gini_grand + employment_grand + social_capital_grand 
                + hispanic_grand.mean + med.inc_grand.mean
                + emergency_hosp_std
                + emergency_hosp_std*gini_grand
                + (1 + gini_grand |HRR))
```

# Findings

The finding of the study is twofold. First, we were able to identify
employment as a statistically significant predictor of overdose deaths.
Specifically, we found that a movement in employment from the least
employed county to the most employed county in the sample led to a
predicted 1.34 decrease in the *order of magnitude* of overdose deaths
for the county at the grand mean for all of the controlling variables
(ex. population, rurality, median income, etc.). For the grand mean
county, this meant about a 74% reduction in the number of overdose
deaths, however, this percent change varies depending on the initial
order of magnitude. We also detect a moderating effect of the total
number of emergency rooms in a cluster on the relationship between
income inequality and overdose deaths, controlling for relevant factors
like population, urbanicity, age and ethnicity. Below we describe this
finding more closely. For now we can note that about 27% of variability
in the relationship between income inequality and overdose deaths is
attributable to the availability of emergency rooms.

    ## Total L1 variance remaining:  2.892758

    ## Total L2 variance remaining:  0.02333052

    ## Conditional ICC:  0.008000621 
    ## 

    ## % L1 variance explained:  0.1293143

    ## % L2 variance explained:  0.7941389 
    ## 

    ## Psuedo R^2 0.1698166 
    ## 

    ## Gini Slope variance remaining 4.678746

    ## %Gini Slope variance explained 0.2626093

First, we identify the relevant effect sizes by rerunning the model on a
standardized dataset and we find that the main effect for employment and
the interaction effect for income inequality and emergency room count
have practically significant effect sizes for the social sciences where
an effect of .1 is considered small.

``` r
std_dat <- as.data.frame(cbind(HRR = Final$HRR, sapply(cbind(
                  overdose_t = Final$overdose_t, Final[row.names(summary(L2_step1)$coefficients)[-c(1,12)]]), 
                  function(x) as.numeric(scale(x)))))

std_mod <- lmer(formula(L2_step1), data = std_dat)
```

    ## Main effect for employment:

    ##     Estimate   Std. Error           df      t value     Pr(>|t|) 
    ##  -0.04878225   0.02113476 782.84550793  -2.30815265   0.02125046

    ## 
    ##  Main effect for income innequality (not statistically significantly > 0):

    ##     Estimate   Std. Error           df      t value     Pr(>|t|) 
    ##  -0.01590534   0.02115272 120.95075521  -0.75192870   0.45355469

    ## 
    ##  Main effect for emergency room count (not statistically significantly > 0):

    ##     Estimate   Std. Error           df      t value     Pr(>|t|) 
    ##   0.01449643   0.02049607 134.92941489   0.70727870   0.48061308

    ## 
    ##  Interaction effect for emergency room count and income innequality (not statistically significantly > 0):

    ##     Estimate   Std. Error           df      t value     Pr(>|t|) 
    ##  0.053623543  0.019652105 74.595758795  2.728641228  0.007927009

Lastly, we can interpret the direction of the interaction. Recall that
we are interested as to why previous research has failed to consistently
show a relationship between substance use deaths and income inequality.
Here we suggest that differences across hospital systems are to blame.
Our finding shows that hospital systems with many emergency rooms show
an *increased* effect of inequality on substance use deaths (made
apparent in the figure below). While more research is required to
explain why this might be the case, my guess is that hospital systems
with more emergency rooms see greater utilization by more affluent
individuals. If true, we would expect more socially unequal communities
to see greater overdoses. On the other hand, in hospital systems with
fewer emergency rooms, their might be less disparity in utilization.

    ## `geom_smooth()` using formula = 'y ~ x'

<img src="HLM_analysis_files/figure-gfm/mod graph-1.png" style="display: block; margin: auto;" />
