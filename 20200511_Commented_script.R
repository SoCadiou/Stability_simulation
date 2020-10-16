################################################################################
##this code allows to perform a simulation to assess performance in terms of
##stability, sensitivity and specificity of prespecified statistical methods
##used to find the true predictors of an health among the exposome it contains 5
##parts:

##1. defining the functions allowing to generate a realistic dataset of exposome
#and an outcome linearly related to some variables of this exposome. 

#It needs a real exposome dataset as input, as well as parameters allowing to
#define the #association between the exposome and the outcome (number of
#predictors, #variability explained, correlation..)

##2. defining the methods assessed

##3. defining some functions used to assess methods performance

##4. defining the simulation function, which, for a given scenario, generates
##the datasets, applies the methods and assess their performance. This function
##allows to parallelize the simulation.

##5. runnning the simulation itself with parallelization, repeating X times the
##function defined in 4. for each scenario and saving the results.
################################################################################

library(mvtnorm)
library(boot)
library(parallel)
library(reshape)
library(glmnet)
library(DSA)
library(OmicsMarkeR)
library(Rcpp)
library(stringr)

###########################################################
######### 1. Generating datasets - functions ##############
###########################################################

##Define functions used to generate datasets function which generates a
##dataset of exposures with same number of variables and individuals and similar
##correlation structure than a real exposome matrix provided and an outcome
##linearly generated

simulator_2layers <- function(E_true,
                              #real exposome data
                              R2_tot = 0.1 ,
                              #total variability explained all predictors
                              n_Ey = 5,
                              #number of predictors
                              BetaEy = 0.01,
                              #Beta coefficient for each predictors. Can be a
                              #vector of values or a unique value
                              test_and_training = TRUE,
                              #generate a dataset of the same size of E_true (if
                              #FALSE) or double the number of individuals (if
                              #TRUE)
                              pos_and_neg = FALSE,
                              #if FALSE, all effects are positive; if TRUE, half
                              #are negative
                              corr = F,
                              #if TRUE, the correlation between the predictors
                              #is controled
                              range_corr = c(0, 1)) {
  #if corr = TRUE, range of correlation for true
  #predictors
  
  ##creating a new exposome dataset by bootstraping
  data.X <- as.data.frame(E_true)
  names_row <- rownames(data.X)
  data.X <- data.X[sample(1:nrow(data.X), 2 * nrow(data.X), replace = TRUE),]
  rownames(data.X) <- c(names_row, sprintf('boot%s', names_row))
  dataExp <- data.X
  
  remove(data.X)
  
  
  ##creating a linearly generated outcome
  ##defining the vector of Beta coefficients
  if (length(BetaEy) == 1) {
    if (pos_and_neg == FALSE) {
      Betapred_yE <- rep(BetaEy, n_Ey)
    } else{
      Betapred_yE <- rep(BetaEy, round(n_Ey / 2))
      Betapred_yE <- c(Betapred_yE, rep(-BetaEy, n_Ey - length(Betapred_yE)))
    }
  } else{
    if (length(BetaEy) != n_Ey) {
      stop(
        "error: Betas for M explaining Y not explained by E are not 
        consistent with the number of predictors"
      )
    }
    Betapred_yE <- BetaEy
  }
  ##creating the linear combination
  yE <- simResponseSimple(
    met = dataExp,
    Nmet = n_Ey,
    beta = Betapred_yE,
    corr = corr,
    range_corr = range_corr
  )
  
  ##adding a gaussian to Y to reach the wanted level of variability explained by
  ##the linear combination of predictors
  Y <- yE$resp
  if (!is.na(R2_tot)) {
    if (((R2_tot) != 0)) {
      sigma <- var(Y) * (1 / R2_tot - 1)
    } else{
      R2 = 0.00000001
      sigma <- var(Y) * (1 / R2_tot - 1)
    }
    Y <-
      as.matrix(Y + rnorm(length(Y), mean(Y), sqrt(sigma)), ncol = 1)
    Y <- as.data.frame(Y)
  }
  
  ##estimating the R2
  R2 <- estimatedR2(dataExp, yE$predictors, Y)$r.squared
  ##results to return
  resultats <-  list(
    Y_train = Y[1:(nrow(dataExp) / 2), , drop = FALSE],
    ##train part of the generated Y vector
    E_train = dataExp[1:(nrow(dataExp) / 2), , drop = FALSE],
    ##train part of the generated exposome dataset
    Y_test = Y[(nrow(dataExp) / 2):nrow(dataExp), , drop =
                 FALSE],
    ##test part of the generated Y vector
    E_test = dataExp[(nrow(dataExp) / 2):nrow(dataExp), , drop = FALSE],
    ##test part of the generated exposome dataset
    yE = yE,
    ##yE (output of the simResponseSimple) object containing the list of
    ##predictors and the Betas
    R2 = R2,
    ##estimated R2
    list_predictor = as.character(yE$predictors) ##vectors of predictors
  )
  return(resultats)
}


####function to generate a linear response####
simResponseSimple <-function(met,
                             ##dataframe of potential predictors
                             Nmet = NA,
                             ##number of predictors
                             beta = NULL,
                             ##Betas coefficient for predictors. Can be a vector of values or a
                             ##unique value
                             cpg = NULL,
                             ##name of forced predictors if necessary
                             corr = FALSE,
                             #if TRUE, the correlation between the predictors is controled
                             range_corr = c(0, 1)) {
  #if corr = TRUE, range of correlation for true predictors
  if (all(c(is.na(Nmet), is.null(cpg))) == TRUE) {
    ##case with no link between the response and the dataset of potential
    ##predictors
    return (list(
      resp = as.matrix(rep(0, nrow(met)), ncol = 1),
      beta = NA,
      predictors = NA
    ))
  }
  if (corr == FALSE | Nmet == 1) {
    ##case of only 1 predictors and no correlation cotnrol
    temp <- Nmet - length(cpg)
    if (temp != 0) {
      if (length(cpg) == 0) {
        wh <- sample((1:ncol(met)), temp) ##drawing predictors
      } else{
        ##drawing predictors while conserving those specified if some were
        ##specified as input
        wh <- sample((1:ncol(met)[-cpg]), temp) 
        wh <- c(cpg, wh)
      }
    } else{
      wh <- cpg
    }
  } else{
    if (length(cpg) != 0) {
      stop("set corr to true is only possible when names of predictors
             are not provided")
    }
    ##if a specified correlation between predictors is set by the users,
    ##selecting the predictors to be in this specified range
    wh <- submatFindSimpl(Mat <- as.matrix(met), range = range_corr, Nvar = Nmet)
    wh <- which(colnames(met) %in% wh)
  }
  ##defining a matrix of predictors
  CovMat <- as.matrix(met[, wh])
  colnames(CovMat) <- colnames(met)[wh]
  # computing the response
  mean <- CovMat %*% matrix(beta, ncol = 1)
  rownames(mean) <- rownames(met)
  names(beta) <- colnames(CovMat)
  return (list(
    resp = mean,
    ##response vector
    beta = beta,
    ##Betas coefficient vector
    predictors = colnames(CovMat) ##vector of predictors
  ))
}

####function to choose a set of predictors among a set of variables with a
####constraint on the correlation range between predictors####
submatFindSimpl <- function(Mat, range = c(0, 1), Nvar) {
  # verifying formats and values of inputs
  if (Nvar > ncol(Mat))
    stop("No matrix of the correct size meeting the range criterion")
  if (Nvar < 2)
    stop("Nvar must be at least 2")
  # computing the correlation matrix
  Mat <- abs(cor(Mat))
  diag(Mat) <- NA
  # removing rows with no correlation value in the given range
  wh <- which(apply(Mat, 1, min, na.rm = T) > range[2] |
                apply(Mat, 1, max, na.rm = T) < range[1])
  if (length(wh) > 0)
    Mat <- Mat[-wh, -wh]
  # iteratively selecting and testing samples for the correct correlation
  Res <- NA
  samp1 <- sample(1:ncol(Mat), size = ncol(Mat))
  t1 <- 1
  while (t1 <= ncol(Mat) & all(is.na(Res))) {
    var <- array(NA, 0)
    t2 <- samp1[t1]
    while (all(is.na(Res)) & length(t2) > 0) {
      if (length(t2) == 1) var <- c(var, t2)
      if (length(t2) > 1) var <- c(var, sample(t2, size = 1))
      t2 <- (1:ncol(Mat))[-var]
      if (length(t2) > 1)
        t2 <- t2[which(apply(as.matrix(Mat[t2, var] <= range[2] &
                                         Mat[t2, var] >= range[1] &
                                         Mat[t2, var] < 1), 1, min) == 1)]
      if (length(t2) == 1){
        if (min(c(Mat[t2, var] <= range[2], Mat[t2, var] >= range[1], 
                  Mat[t2, var] < 1)) != 1)
          t2 <- NA
      }
      
      if (length(var) == Nvar)
        Res <- var
    }
    t1 <- t1 + 1
  }
  if (all(is.na(Res)))
    stop("Not enough variable with the given correlation range")
  return(colnames(Mat)[Res])
}

####function which estimates R2 from a dataset of potential predictors, the list
####of true predictors and a vector of outcome####
estimatedR2 <- function(X, truepred, Y) {
  if ("y" %in% truepred) {
    stop("error: one of the true predictors is named y")
  }
  if (ncol(Y) != 1) {
    stop("error:Y is multidimensionnal")
  }
  if (nrow(X) != nrow(Y)) {
    stop("error: not the same number of rows")
  }
  if (isTRUE(all.equal(rownames(X), rownames(Y))) == FALSE) {
    stop("error: individuals are not ordered similarly in X and Y")
  }
  if (all(truepred %in% colnames(X))) {
    data <- X[, colnames(X) %in% truepred, drop = FALSE]
    data <- cbind(Y, data)
    colnames(data)[1] <- "y"
    mod <- lm(y ~ ., as.data.frame(data))
    toselect.x <- summary(mod)$coeff[-1, 4]
    r <- list(summary(mod)$r.squared,
              summary(mod)$adj.r.squared,
              names(toselect.x)[toselect.x == TRUE])
    names(r) <- c("r.squared", "adj.r.squared", "pred")
    return(r)
  } else{
    stop("error: X does not countain all true predictors")
    
  }
}

################################################################################
##an other function simulator_2layers can be used if one want to control the
##correlation of the overall dataset in this case, the simulated exposome matrix
##is no longer obtained by bootstrapping the real exposome but from a
##correlation matrix the correlation matrix must be provided by the user or
##specified as nulll uncomment the section below and adapt the input of
##functions f0 (section 4.) and clusterApply (section 5.) to use it

# simulator_2layers <-
#   function(names_rows_true,
#            cormat,
#            R2_tot = 0.1 ,
#            n_Ey = 5,
#            BetaEy = 0.01,
#            test_and_training = TRUE,
#            pos_and_neg = FALSE,
#            corr_all = "real",
#            corr_pred = TRUE,
#            range_corr = c(0, 1)) {
#     ##generation of exposome dataset DIFFERENT FROM THE OTHER FUNCTION 
#       simulator_2layers
#     if (corr_all == "real") {
#       data.X <- data.frame(rmvnorm(1173 * 2, rep(0, ncol(cormat)),
#                                    cormat))
#     } else{
#       if (corr_all == "null") {
#         cormat2 <- diag(ncol(cormat))
#         data.X <- data.frame(rmvnorm(1173 * 2, rep(0, ncol(cormat)),
#                                      cormat2))
#       } else{
#         stop("Correlation for the whole exposome (null or real) must 
#               be specifed")
#       }
#     }
#     colnames(data.X) <- colnames(cormat)
#     rownames(data.X) <-
#       c(names_rows_true, sprintf('boot%s', (names_rows_true)))
#     dataExp <- data.X
#     remove(data.X)
# 
#     ##FROM HERE, THE FUNCTION IS IDENTICAL TO THE OTHER simulator_2layers
#     function ##creating a linearly generated outcome ##defining the vector of
#     Beta coefficients
#
#     if (length(BetaEy) == 1) {
#       if (pos_and_neg == FALSE) {
#         Betapred_yE <- rep(BetaEy, n_Ey)
#       } else{
#         Betapred_yE <- rep(BetaEy, round(n_Ey / 2))
#         Betapred_yE <-
#           c(Betapred_yE, rep(-BetaEy, n_Ey - length(Betapred_yE)))
#       }
#     } else{
#       if (length(BetaEy) != n_Ey) {
#         stop(
#           "error: Betas for M explaining Y not explained by E are not 
#             consistent with the number of predictors"
#         )
#       }
#       Betapred_yE <- BetaEy
#     }
#     ##creating the linear combination
#     yE <-
#       simResponseSimple(
#         met = dataExp,
#         Nmet = n_Ey,
#         beta = Betapred_yE,
#         corr = corr,
#         range_corr = range_corr
#       )
#     
#     ##adding a gaussian to Y to reach the wanted level of variability 
#     explained by the linear combination of predictors
#     Y <- yE$resp
#     if (!is.na(R2_tot)) {
#       if (((R2_tot) != 0)) {
#         sigma <- var(Y) * (1 / R2_tot - 1)
#       } else{
#         R2 = 0.00000001
#         sigma <- var(Y) * (1 / R2_tot - 1)
#       }
#       Y <-
#         as.matrix(Y + rnorm(length(Y), mean(Y), sqrt(sigma)), ncol = 1)
#       Y <- as.data.frame(Y)
#     }
#     
#     ##estimating the R2
#     R2 <- estimatedR2(dataExp, yE$predictors, Y)$r.squared
#     ##results to return
#     resultats <-
#       list(
#         Y_train = Y[1:(nrow(dataExp) / 2), , drop = FALSE],
#         ##train part of the generated Y vector
#         E_train = dataExp[1:(nrow(dataExp) / 2), , drop = FALSE],
#         ##train part of the generated exposome dataset
#         Y_test = Y[(nrow(dataExp) / 2):nrow(dataExp), , drop =
#                      FALSE],
#         ##test part of the generated Y vector
#         E_test = dataExp[(nrow(dataExp) / 2):nrow(dataExp), , drop = FALSE],
#         ##test part of the generated exposome dataset
#         yE = yE,
#         ##yE (output of the simResponseSimple) object containing the list of
#         ##predictors and the Betas
#         R2 = R2,
#         ##estimated R2
#         list_predictor = as.character(yE$predictors) ##vectors of predictors
#       )
#     return(resultats)
#   }

###########################################################
########### 2. Methods to be tested - functions ###########
###########################################################


####a function used to compute residuals from a linear model if covariates are
####part of the inputs of any the function of the methods tested####
getresiduals_2df <-
  function(data_Y_in, data_covar_in, name_Y, covar) {
    data_covar <-
      data_covar_in[, colnames(data_covar_in) %in% covar, drop = FALSE]
    data_Y <-
      data_Y_in[rownames(data_Y_in) %in% rownames(data_covar), 
                colnames(data_Y_in) ==
                  name_Y, drop = FALSE]
    data_covar <-
      data_covar[rownames(data_covar) %in% rownames(data_Y), , drop = FALSE]
    data_covar <- data_covar[rownames(data_Y), , drop = FALSE]
    data_output <- data_Y
    data <- cbind(data_Y, data_covar)
    mod <- lm(data = data)
    data_output[, 1] <- as.data.frame(residuals(mod))
    return(data_output)
  }

####ExWAS####
ewas <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0),
           ##if necessary, vector of covariates name
           corr = "BY",
           ##name of multiple testing correction to be applied ("BH" or "Bon" or
           ##"BY" or "None")
           ntest = NULL) {
    ##if ntest is a numeric, correction of multiple testing will be a Bonferroni
    ##correction considering ntest as the number of tests performed
    require(parallel)
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar<-data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                                  rownames(data_covar_in) %in% rownames(data_Xs_in), 
                                colnames(data_covar_in) %in% covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in),
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y), , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), , drop = FALSE]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), 
                   , drop = FALSE]
      data_Xs <- as.data.frame(data_Xs[rownames(data_Y), ])
      colnames(data_Xs) <- colnames(data_Xs_in)
    }
    ##checking consistency of the datasets
    if (is.null(data_Y) == TRUE |
        is.null(data_Xs) == TRUE |
        !(name_Y %in% colnames(data_Y))) {
      stop("Données incohérentes entre elles")
    }
    
    ##applying univariate regression for each exposure outcome association
    p.values <- mclapply(1:ncol(data_Xs), function(x, data_Xs) {
      c(colnames(data_Xs)[x], summary(lm(Y ~ .,
                                         data = data.frame(
                                           cbind(var1 = data_Xs[, x],
                                                 Y = data_Y[, 1])
                                         )))$coefficients[2, ])
    },
    data_Xs)
    if (length(p.values) == 1) {
      p.values <-
        as.matrix(as.vector(unlist(p.values[[1]])), ncol = 5, byrow = TRUE)[-4, ]
      p.values <- t(as.data.frame(p.values))
    } else{
      p.values <-
        cbind(matrix(unlist(p.values), ncol = 5, byrow = TRUE)[, -4])
    }
    p.values <- as.data.frame(p.values)
    colnames(p.values) <- c("var", "Est", "Sd", "pVal")
    p.values <- p.values[p.values$var != "Intercept", ]
    p.values$pVal <- as.numeric(as.character(p.values$pVal))
    p.values.adj <- p.values
    pVal <- as.numeric(as.character(p.values$pVal))
    ##applying correction for multiple testing
    if (corr == "None") {
      wh <- which(pVal <= 0.05)
      p.values.adj$pVal_adj <- pVal
    }
    if (corr == "Bon") {
      wh <- which(pVal <= 0.05 / nrow(p.values))
      p.values.adj$pVal_adj <- pVal * nrow(p.values)
    }
    if (corr == "BH") {
      wh <- which(p.adjust(pVal, "BH") <= 0.05)
      p.values.adj$pVal_adj <- p.adjust(pVal, "BH")
    }
    if (corr == "BY") {
      wh <- which(p.adjust(pVal, "BY") <= 0.05)
      p.values.adj$pVal_adj <- p.adjust(pVal, "BY")
    }
    if (!corr %in% c("Bon", "BH", "BY", "", "None"))
      stop("Please specify a known correction method for
                                                 multiple testing")
    if (!is.null(ntest)) {
      p.values.adj$pVal_adj <- pVal * ntest
    }
    wh_num <- wh
    wh <- p.values$var[wh]
    a <- list(wh, wh_num, p.values.adj)
    ##returning selected exposures and pvalues
    names(a) <- c("selected", "indices_selected", "pval")
    return(a)
  }


####lasso - basic implementation#### this includes a basic 10-fold cross
##validation process as implemented in the CVglmnet package
lasso <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0)) {
    ##if necessary, vector of covariates name
    
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), 
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in),
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), 
                   , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) ==  name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    ##applying lasso (as implemented in glmnet package, a path of penalization
    ##parameter lambda  according to the MSE computed is computed by 10-fold
    ##cross-validation)
    cvfit <- cv.glmnet(data_Xs, data_Y, family = "gaussian",
                       alpha = 1)
    
    ##Compute predicted Y "Y_predit"
    Y_predit <-
      predict(cvfit, newx = data_Xs, s = "lambda.min")
    ##selecting the model with the penalization parameter minimizing MSE
    Y_predit <- Y_predit[rownames(Y_predit), ]
    
    ##dataframe of predictors selected
    tmp_coeffs <-
      coef(cvfit, s = "lambda.min") 
    ##selecting the model with the penalization parameter minimizing MSE
    cg_select <-
      data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                 coefficient = tmp_coeffs@x)
    
    cg_select <- cg_select$name[cg_select$name != "(Intercept)"]
    a <- list()
    if (length(cg_select) != 0) {
      a <- list("selected" = cg_select,
                "prediction" = Y_predit,
                "null" = "nul")
    } else{
      a <-
        list(
          "selected" = character(),
          "prediction" = "pas_de_prediction",
          "null" = "nul"
        )
    }
    return(a)
    
  }

####lasso - basic implementation but using lambda.1se instead of lambda.min
#### this includes a basic 10-fold cross
##validation process as implemented in the CVglmnet package
lasso_1SE <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0)) {
    ##if necessary, vector of covariates name
    
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), 
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in),
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), 
                   , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) ==  name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    ##applying lasso (as implemented in glmnet package, a path of penalization
    ##parameter lambda  according to the MSE computed is computed by 10-fold
    ##cross-validation)
    cvfit <- cv.glmnet(data_Xs, data_Y, family = "gaussian",
                       alpha = 1)
    
    ##Compute predicted Y "Y_predit"
    Y_predit <-
      predict(cvfit, newx = data_Xs, s = "lambda.1se")
    ##selecting the model within 1MSE of the model using
    ##the penalization parameter minimizing MSE
    Y_predit <- Y_predit[rownames(Y_predit), ]
    
    ##dataframe of predictors selected
    tmp_coeffs <-
      coef(cvfit, s = "lambda.min") 
    ##selecting the model with the penalization parameter minimizing MSE
    cg_select <-
      data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                 coefficient = tmp_coeffs@x)
    
    cg_select <- cg_select$name[cg_select$name != "(Intercept)"]
    a <- list()
    if (length(cg_select) != 0) {
      a <- list("selected" = cg_select,
                "prediction" = Y_predit,
                "null" = "nul")
    } else{
      a <-
        list(
          "selected" = character(),
          "prediction" = "pas_de_prediction",
          "null" = "nul"
        )
    }
    return(a)
    
  }


####lasso_stab : implementation of Meinshausen 2010
#### repeating lasso on
##subsamples and performing the selection according to empirical probability of
##selection computed over the repeated runs using a threshold specified by user
lasso_stab_Meinshausen <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0),
           ##if necessary, vector of covariates name
           prop = 0.85) {
    #minimal threshold for an empirical probability to make the corresponding
    #variable selected
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in),
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in),
                   , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    if (is.null(data_Y)) {
      stop("Données incohérentes entre elles _ Y")
    }
    if (is.null(data_Xs)) {
      stop("Données incohérentes entre elles _ Xs")
    }
    if (!(name_Y %in% colnames(data_Y))) {
      stop("Données incohérentes entre elles _ nom Y")
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    list_iter <- list()
    
    length_lambdas <- numeric()
    n_iter_stab <-  100 
    ##setting the number of repetitions from which the empirical probabilities
    ##will be computed
    for (k2 in (1:n_iter_stab)) {
      #randomly selecting a subsample containing 50% of individuals
      selec <-
        sample(rownames(data_Xs), round(1 * nrow(data_Xs) / 2)) #
      i_df <- data.matrix(data_Xs[rownames(data_Xs) %in% selec, ])
      block_pheno <-
        data_Y[rownames(data_Y) %in% rownames(i_df), , drop = FALSE]
      block_pheno <- block_pheno[rownames(i_df), , drop = FALSE]
      ##applying lasso to the subsample
      model.lasso <-
        glmnet(
          x = i_df,
          y = data.matrix(block_pheno),
          family = "gaussian",
          alpha = 1
        )
      #saving the selection for each lambda in a dataframe
      selection_pour_une_iter <-
        as.data.frame(cbind(rownames(model.lasso$beta), as.numeric(tabulate(
          model.lasso$beta@i + 1
        ))))
      colnames(selection_pour_une_iter) <-
        c("variables", paste("iter", k2))
      ##adding this dataframe to the list of dataframes computed for all
      ##precedent subsamples
      list_iter <- c(list_iter, list(selection_pour_une_iter))
      ##saving the length of penalization path
      length_lambdas <-
        c(length_lambdas, length(model.lasso$lambda))
    }
    ##merging all dataframes of selection by variable
    M <-
      as.data.frame(Reduce(function(x, y)
        merge(x, y, by = "variables"), list_iter))
    M[, 2:ncol(M)] <-
      lapply(M[, 2:ncol(M)], function(x)
        as.numeric(as.character(x)))
    ##summing the frequencies of selection by variables
    M$sum <- rowSums(M[, 2:ncol(M)])
    M1 <- M[, c(1, ncol(M))]
    ##computing the empirical probabilities from the sum of frequencies of
    ##selection by variables and the number of possible selections
    M1$proba_estimee <- M1$sum / (sum(length_lambdas))
    ##selecting variables from the empirical probabilities and the threshold
    cg_select <- M1$variables[M1$proba_estimee >= prop]
    ##returning the selection
    a <- list()
    if (length(cg_select) != 0) {
      a <-
        list(
          "selected" = cg_select,
          "selection_iter" = list(M, M1),
          "iteration" = n_iter_stab
        )
    } else{
      a <-
        list(
          "selected" = character(0),
          "selection_iter" = list(M, M1),
          "iteration" = n_iter_stab
        )
    }
    return(a)
    
  }

####lasso_stab : second implementation of Meinshausen 2010
#### selection of a
##set of lambda parameters on a subsample then repeating lasso on subsamples
##with those lambdas and performing the selection according to empirical
##probability of selection computed over the repeated runs using a threshold
##specified by user
lasso_stab_Meinshausen2 <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0),
           ##if necessary, vector of covariates names
           prop = 0.95) {
    #minimal threshold for an empirical probability to make the corresponding
    #variable selected
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), 
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in),
                  colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    if (is.null(data_Y)) {
      stop("Données incohérentes entre elles _ Y")
    }
    if (is.null(data_Xs)) {
      stop("Données incohérentes entre elles _ Xs")
    }
    if (!(name_Y %in% colnames(data_Y))) {
      stop("Données incohérentes entre elles _ nom Y")
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    
    
    list_iter <- list()
    ###defining a vector of penalized parameters lambda on a subsample
    selec <- sample(rownames(data_Xs), round(1 * nrow(data_Xs) / 2))
    i_df <- data.matrix(data_Xs[rownames(data_Xs) %in% selec, ])
    block_pheno <-
      data_Y[rownames(data_Y) %in% rownames(i_df), , drop = FALSE]
    block_pheno <- block_pheno[rownames(i_df), , drop = FALSE]
    temp <-
      glmnet(
        x = i_df,
        y = block_pheno,
        family = "gaussian",
        alpha = 1
      ) 
    ##lasso is applied; the path of lambdas will be used for all other
    ##subsamples
    lambdas <- temp$lambda
    ###repeating lasso fo all subsamples
    n_iter_stab <- 100 
    ##setting the number of repetitions from which the empirical probabilities
    ##will be computed
    for (k2 in (1:n_iter_stab)) {
      selec <-
        sample(rownames(data_Xs), round(1 * nrow(data_Xs) / 2))
      #randomly selecting a subsample containing 50% of individuals
      i_df <- data.matrix(data_Xs[rownames(data_Xs) %in% selec, ])
      block_pheno <-
        data_Y[rownames(data_Y) %in% rownames(i_df), , drop = FALSE]
      block_pheno <- block_pheno[rownames(i_df), , drop = FALSE]
      ##applying lasso to the subsample
      model.lasso <-
        glmnet(
          x = i_df,
          y = data.matrix(block_pheno),
          family = "gaussian",
          alpha = 1,
          lambda = as.numeric(lambdas)
        )
      #saving the selection for each lambda in a dataframe
      selection_pour_une_iter <-
        as.data.frame(cbind(rownames(model.lasso$beta), as.numeric(tabulate(
          model.lasso$beta@i + 1
        ))))
      colnames(selection_pour_une_iter) <-
        c("variables", paste("iter ", k2))
      ##adding this dataframe to the list of dataframes computed for all
      ##precedent subsamples
      list_iter <- c(list_iter, list(selection_pour_une_iter))
    }
    ##merging all dataframe of selection by variable
    M <-
      as.data.frame(Reduce(function(x, y)
        merge(x, y, by = "variables"), list_iter))
    M[, 2:ncol(M)] <-
      lapply(M[, 2:ncol(M)], function(x)
        as.numeric(as.character(x)))
    ##summing the frequencies of selection by variables
    M$sum <- rowSums(M[, 2:ncol(M)])
    M1 <- M[, c(1, ncol(M))]
    ##computing the empirical probabilities from the sum of frequencies of
    ##selection by variables and the number of possible selections
    M1$proba_estimee <- M1$sum / (length(lambdas) * n_iter_stab)
    M1$proba_estimee <- M1$sum / (sum(length_lambdas))
    ##selecting variables from the empirical probabilities and the threshold
    cg_select <- M1$variables[M1$proba_estimee >= prop]
    ##returning the selection
    a <- list()
    if (length(cg_select) != 0) {
      a <-
        list(
          "selected" = cg_select,
          "selection_iter" = list(M, M1),
          "iteration" = n_iter_stab
        )
    } else{
      a <-
        list(
          "selected" = character(0),
          "selection_iter" = list(M, M1),
          "iteration" = n_iter_stab
        )
    }
    return(a)
    
  }

####LASSO_CV2####

### loop applying 100 times lasso with the 10-fold validations procedures on the
### whole dataset the average of the mean error curves (MSE as a function of the
### MSE) gives the penalization parameter which minimizes this averaged MSE and
### which will be used in the final model
lasso_moy_MSE <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0)) {
    ##if necessary, vector of covariates name
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), 
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in),
                  colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), 
                   , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    lambdas <- numeric(0)
    MSEs <- data.frame(matrix(NA, nrow = 100, ncol = 100))
    ##repeating 100 times the cross-validation process
    for (i in 1:100) {
      cv <- cv.glmnet(x = data_Xs, y = data_Y, alpha = 1) ##applying lasso
      MSEs[1:length(cv$lambda), i] <-
        cv$cvm 
      ##saving the cross-validation path (ie MSE values for each lambdas) as a
      ##colum of MSEs dataframe
      if (length(cv$lambda) > length(lambdas)) {
        lambdas <-
          cv$lambda ##saving the vector of lambdas of maximum length
      }
      print(cv$lambda[1:10])
      
    }
    ##restricting MSEs to non empty rows
    MSEs <- MSEs[1:length(lambdas), ]
    rownames(MSEs) <- lambdas
    ##choosing the lambda minimimizing the mean of MSE computed across all
    ##lambdas
    lambda.min <-
      as.numeric(names(which.min(rowMeans(MSEs, na.rm = TRUE))))
    
    ##applying lasso with this lambda as forced penalization parameter
    model.enet <- glmnet(
      data_Xs,
      data_Y,
      family = "gaussian",
      alpha = 1,
      lambda = lambda.min
    )
    
    ##selecting variables
    tmp_coeffs <- coef(model.enet)
    cg_select <-
      data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1],
                 coefficient = tmp_coeffs@x)
    cg_select <- cg_select$name[cg_select$name != "(Intercept)"]
    ##returning selected variables
    a <- list()
    if (length(cg_select) != 0) {
      a <-
        list("selected" = cg_select,
             "prediction" = "prediction",
             "null" = "nul")
    } else{
      a <-
        list(
          "selected" = character(),
          "prediction" = "pas_de_prediction",
          "null" = "nul"
        )
    }
    return(a)
    
  }


####LASSO_CV1####
### loop applying 100 times lasso with the 10-fold validations procedures on the
### whole dataset the average of the penalization parameters minimizing the MSE
### for each run gives the penalization parameter which will be used in the
### final model

lasso_moy_lambda <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0)) {
    ##if necessary, vector of covariates name
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), 
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in),
                   , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    
    ##initialization
    lambdas <- NULL
    ##repeating 100 times the cross-validation process
    for (i in 1:100) {
      cv <- cv.glmnet(x = data_Xs, y = data_Y, alpha = 1)
      lambdas <- c(lambdas, cv$lambda.min)
    }
    ##averaging the optimal parameters across repetitions
    lambda_min <- mean(lambdas)
    ##applying lasso with this lambda as forced penalization parameter
    model.enet <- glmnet(
      data_Xs,
      data_Y,
      family = "gaussian",
      alpha = 1,
      lambda = lambda_min
    )
    
    ##selecting variables
    tmp_coeffs <- coef(model.enet)
    cg_select <-
      data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1],
                 coefficient = tmp_coeffs@x)
    
    cg_select <- cg_select$name[cg_select$name != "(Intercept)"]
    ##returning selected variables
    a <- list()
    if (length(cg_select) != 0) {
      a <-
        list("selected" = cg_select,
             "prediction" = "prediction",
             "null" = "nul")
    } else{
      a <-
        list(
          "selected" = character(),
          "prediction" = "pas_de_prediction",
          "null" = "nul"
        )
    }
    return(a)
  }




####Mix Method####
##Computing empirical probabilities derived from selection frequencies when
##running cross-validated lasso on subsamples
lasso_moy_Meinshausen <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0),
           ##if necessary, vector of covariates name
           prop = 0.5) {
    ##threshold to select variables from their empirical probabilities
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in),
                      colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in),
                   , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), 
                  colnames(data_Y_in) ==  name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    if (is.null(data_Y)) {
      stop("Données incohérentes entre elles _ Y")
    }
    if (is.null(data_Xs)) {
      stop("Données incohérentes entre elles _ Xs")
    }
    if (!(name_Y %in% colnames(data_Y))) {
      stop("Données incohérentes entre elles _ nom Y")
    }
    ##initialization
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    vector_selected <- character()
    n_iter_stab <- 100  
    #setting the number of repetitions from which the empirical probabilities
    #will be computed
    for (k2 in (1:n_iter_stab)) {
      #randomly selecting a subsample containing 50% of individuals
      selec <-
        sample(rownames(data_Xs), round(1 * nrow(data_Xs) / 2))
      i_df <- data.matrix(data_Xs[rownames(data_Xs) %in% selec, ])
      block_pheno <-
        data_Y[rownames(data_Y) %in% rownames(i_df), , drop = FALSE]
      block_pheno <- block_pheno[rownames(i_df), , drop = FALSE]
      ##applying cross_validated lasso to the subsample
      model.lasso <-
        cv.glmnet(
          x = i_df,
          y = data.matrix(block_pheno),
          family = "gaussian",
          alpha = 1
        )
      ##selecting the model minimizing the MSE
      tmp_coeffs <- coef(model.lasso, s = "lambda.min")
      cg_select <-
        data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
      cg_select <- cg_select$name[cg_select$name != "(Intercept)"]
      #saving the variables selected for this iteration
      vector_selected <- c(vector_selected, as.character(cg_select))
    }
    ##computing frequencies of selection
    tab <- as.data.frame(table(vector_selected))
    colnames(tab) <- c("exp", "proba")
    ##computing empirical probabilities
    tab$proba <- tab$proba / n_iter_stab
    ###selecting variables from the empirical probabilities and the threshold
    cg_select_tot <- tab$exp[tab$proba >= prop]
    ##returning the selection
    a <- list()
    if (length(cg_select_tot) != 0) {
      a <-
        list("selected" = cg_select_tot,
             "freq" = tab,
             "iteration" = n_iter_stab)
    } else{
      a <- list(
        "selected" = character(0),
        "freq" = tab,
        "iteration" = n_iter_stab
      )
    }
    return(a)
    
  }


####Elastic-Net####
Enet <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0)) {
    ##if necessary, vector of covariates name
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    ##preparation of CV folds (which must be the same for alpha and lambda CV)
    nfolds=10
    
    n <- nrow(data_Xs)
    folds <- rep(1:nfolds,length.out=n)[sample(n,n)] 
    #step 1: do all crossvalidations for each alpha in range 0.1 - 1.0 
    alphasOfInterest <- seq(1,0.1,-0.1)
    cvs <- lapply(alphasOfInterest, 
                  function(curAlpha){cv.glmnet(data_Xs,data_Y,alpha=curAlpha,
                                               family="gaussian",nfolds=nfolds,
                                               foldid=folds,standardize=FALSE)})
    #step 2: collect the optimum lambda for each alpha 
    optimumPerAlpha <- sapply(seq_along(alphasOfInterest), function(curi){
      curcvs <- cvs[[curi]] 
      curAlpha <- alphasOfInterest[curi]
      indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
      return(c(lam=curcvs$lambda.min,alph=curAlpha,cvm=curcvs$cvm[indOfMin],
               cvup=curcvs$cvup[indOfMin])) })
    #step 3: find the overall optimum 
    posOfOptimum <- which.min(optimumPerAlpha["cvm",]) 
    overall.alpha.min <- optimumPerAlpha["alph",posOfOptimum]
    overall.lambda.min <- optimumPerAlpha["lam",posOfOptimum]
    overall.criterionthreshold <- optimumPerAlpha["cvup",posOfOptimum]
    
    
    ##final model
    model <-
      glmnet(x=data_Xs,y=data_Y,alpha=overall.alpha.min,
             lambda=overall.lambda.min,family="gaussian",standardize=FALSE)
    ##selecting variables in the model minimizig the MSE
    tmp_coeffs <- coef(model, s =as.numeric(overall.lambda.min))
    cg_select<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], 
                          coefficient = tmp_coeffs@x)
    cg_select<-cg_select$name[cg_select$name!="(Intercept)"&!is.na(cg_select$name)&cg_select$name!="<NA>"]
    ##returning selection
    a<-list()
    if (length(cg_select)!=0){
      a<-list("selected"=cg_select,"alpha_final"=overall.alpha.min,
              "lambda_final"=overall.lambda.min)  
    }else{
      a<-list("selected"=character(0),"alpha_final"=overall.alpha.min,
              "lambda_final"=overall.lambda.min)
    }
    return(a)
    
  }


####Elastic-Net stabilized by repeating the CV process####
Enet_CV <-
  function(data_Xs_in,
           ##dataset of explanatory variables ("exposures")
           data_Y_in,
           ##dataset of univariate variable of interest ("outcome")
           name_Y,
           ##variable of interest name
           data_covar_in = NULL,
           ##if neccessary, dataset of covariates ("confounders")
           covar = character(0)) {
    ##if necessary, vector of covariates name
    if (length(covar) > 0) {
      ##if necessary, computing residuals of the linear model explaining the
      ##variable of interest by the covariates
      data_covar <-
        data_covar_in[rownames(data_covar_in) %in% rownames(data_Y_in) &
                        rownames(data_covar_in) %in% rownames(data_Xs_in), colnames(data_covar_in) %in%
                        covar, drop = FALSE]
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_covar) &
                    rownames(data_Y_in) %in% rownames(data_Xs_in), colnames(data_Y_in) == name_Y, drop =
                    FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_covar) &
                     rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_covar <- data_covar[rownames(data_Y), ]
      data_Xs <- data_Xs[rownames(data_Y), ]
      data_Y <- getresiduals_2df(data_Y, data_covar, name_Y, covar)
    } else{
      data_Y <-
        data_Y_in[rownames(data_Y_in) %in% rownames(data_Xs_in), colnames(data_Y_in) ==
                    name_Y, drop = FALSE]
      data_Xs <-
        data_Xs_in[rownames(data_Xs_in) %in% rownames(data_Y_in), , drop = FALSE]
      data_Xs <- data_Xs[rownames(data_Y), ]
    }
    data_Y <- data.matrix(data_Y)
    data_Xs <- data.matrix(data_Xs)
    
    
    nfolds=10
    
    ##the cross-validation process is repeated 100 times for alpha
    all_alpha_min<-numeric(0)
    n <- nrow(data_Xs)
    optimumAlpha<-NULL
    for (i in 1:100){
      print(i)
      ##preparation of CV folds (which must be the same for alpha and lambda CV)
      folds <- rep(1:nfolds,length.out=n)[sample(n,n)] 
      #step 1: do all crossvalidations for each alpha in range 0.1 - 1.0 
      alphasOfInterest <- seq(1,0.1,-0.1)
      cvs <- lapply(alphasOfInterest, 
                    function(curAlpha){cv.glmnet(data_Xs,data_Y,alpha=curAlpha,
                                                 family="gaussian",nfolds=nfolds,
                                                 foldid=folds,standardize=FALSE)})
      #step 2: collect the optimum lambda for each alpha 
      optimumPerAlpha <- sapply(seq_along(alphasOfInterest), function(curi){
        curcvs <- cvs[[curi]] 
        curAlpha <- alphasOfInterest[curi]
        indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
        return(c(lam=curcvs$lambda.min,alph=curAlpha,cvm=curcvs$cvm[indOfMin],
                 cvup=curcvs$cvup[indOfMin])) })
      posOfOptimum <- which.min(optimumPerAlpha["cvm",]) 
      overall.alpha.min <- optimumPerAlpha["alph",posOfOptimum]
      #list_optimumPerAlpha<-cbind(list_optimumPerAlpha,list(optimumPerAlpha))
      optimumAlpha<-c(optimumAlpha,overall.alpha.min)
    }
    ##the final alpha value is computed
    alpha_final<-mean(optimumAlpha)
    ##comment : it is not possible to obtain simultaneously the average values of
    ##lambda and alpha, as lambda must be computed according to alpha
    
    ##the cross-validation process is repeated 100 times for lambda
    all_lambda_min<-numeric(0)
    for (i in 1:100){
      overall.lambda.min<-cv.glmnet(x=data_Xs,y=data_Y,alpha=alpha_final,family="gaussian")
      plot(overall.lambda.min)
      
      all_lambda_min<-c(all_lambda_min,overall.lambda.min)
    }
    lambda_final<-mean(all_lambda_min)
    
    ##final model
    model <-
      glmnet(x=data_Xs,y=data_Y,alpha=alpha_final,lambda=lambda_final,family="gaussian",standardize=FALSE)
    ##selecting variables in the model minimizig the MSE
    tmp_coeffs <- coef(model, s =as.numeric(overall.lambda.min))
    cg_select<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
    cg_select<-cg_select$name[cg_select$name!="(Intercept)"&!is.na(cg_select$name)&cg_select$name!="<NA>"]
    ##returning selection
    a<-list()
    if (length(cg_select)!=0){
      a<-list("selected"=cg_select,"alpha_final"=overall.alpha.min,"lambda_final"=overall.lambda.min)  
    }else{
      a<-list("selected"=character(0),"alpha_final"=overall.alpha.min,"lambda_final"=overall.lambda.min)
    }
    return(a)
    
  }


####DSA####
##DSA is an iterative linear regression model search algorithm (Sinisi and van
##der Laan 2004) following three constraints: maximum order of interaction
##amongst predictors, the maximum power for a given predictor, and the maximum
##model size.
DSAreg <-
  function(Exp,
           ####dataset of explanatory variables ("exposures")
           resp,
           ##dataset of univariate variable of interest ("outcome")
           family = gaussian,
           ##family of the outcome
           maxsize = 15,
           ##maximum size of the model
           maxsumofpow = 2,
           ##maximum power for a given predictor
           maxorderint = 2) {
    ##maximum order of interaction
    Exp <-
      data.frame(cbind(resp = resp, data.frame(Exp))) 
    ##merging exp and resp in a unique dataframe
    colnames(Exp)[1] <- "resp"
    ##applying DSA function with 5 fold split and 1 cross-validation process
    res <-
      DSA(
        resp ~ 1,
        data = Exp,
        family = family,
        maxsize = maxsize,
        maxsumofpow
        = maxsumofpow,
        maxorderint = maxorderint ,
        nsplits = 1,
        usersplits = NULL
      )
    ##extracting the selected variables in case there are power or interaction
    form <- gsub("I[(]", "", colnames(coefficients(res)))
    form <-
      gsub("[*]", ":", gsub("[)]", "", gsub("[:^:]1", "", form)))
    if (length(grep(":", form)) > 0) {
      nam <- strsplit(form[grep("[:]", form)], ":")
      for (j in 1:length(nam)) {
        nam[[j]] <- gsub("[[:space:]]", "", nam[[j]])
        name <- nam[[j]][1]
        for (k in 2:length(nam[[j]]))
          name <- paste(name, ":", nam[[j]][k], sep = "")
        Exp <- cbind(Exp, name = apply(Exp[, nam[[j]]], 1, prod))
      }
    }
    form2 <- "resp~1"
    if (length(form) > 1)
      for (i in 2:length(form))
        form2 <- paste(form2, "+", form[i])
    ##putting the selected variablesin a linear model
    res2 <- lm(form2, data = data.frame(Exp))
    #pred <- predict(res2,Exp)
    ##obtaining beta coefficients
    coef <- summary(res2)$coefficients
    coef <-
      as.character(rownames(coef)[rownames(coef) != "Intercept"])
    ##returning selecting
    return(list(
      selected = coef[coef != "(Intercept)"],
      pred = "prediction",
      null = "null"
    ))
  }
##multivariate regression
multi <- function(Exp, resp) {
  ##Exp is the dataset of potential predictors ("the exposome")
  ##resp is the variable to explain ("the outcome")
  var <- colnames(Exp)
  Exp <- data.frame(cbind(resp = resp, data.frame(Exp)))
  colnames(Exp)[1] <- "resp"
  formula <-
    as.formula(paste("resp ~", paste(var, collapse = " + ")))
  model <- lm(formula = formula, data = Exp)
  selection <- as.data.frame(summary(model)$coeff[-1, 4])
  colnames(selection) <- "pVal"
  selection$pVal_adj <- p.adjust(selection$pVal, "BH")
  selected <-
    as.character(unique(rownames(selection)[selection$pVal_adj <= 0.05]))
  return(list(
    selected = selected,
    pred = "prediction",
    null = "null"
  ))
}

#########################################################
########### 3. Methods assessment - functions ###########
#########################################################


####computing sensitivity####
sensitivity <- function(truepred, predfound) {
  return(length(truepred[truepred %in% predfound]) / length(truepred))
}
####computing false discovery proportion####
fdp <- function(truepred, predfound) {
  if (length(predfound) == 0) {
    return(0)
  } else{
    return(length(predfound[!predfound %in% truepred]) / length(predfound))
  }
}
####computing specificity####
specificity <- function(truepred, predfound, n_base) {
  return((n_base - length(truepred) - length(predfound[!predfound %in% truepred])) /
           (n_base - length(truepred)))
}

################################################################################
###### 4. Function to parallelize - performing simulation and assessment #######
################################################################################

##take only the iteration number as input and use it as a seed
f0 <- function(x) {
  ##setting seed
  set.seed(x)
  ##generating datasets
  simu <-
    simulator_2layers(
      E_true = dataExp_true,
      R2_tot = R2_fixed,
      n_Ey = n_Ey,
      BetaEy = 0.1,
      test_and_training = TRUE,
      pos_and_neg = pos_and_neg,
      corr = corr,
      range_corr = corr_range
    )
  simu$Y_train <- scale(simu$Y_train)
  simu$Y_test <- scale(simu$Y_test)
  list_predBMI_E <- list()
  ##creating a list where to save methods results
  predBMI_E <-
    lapply(1:n_method, function(i)
      lapply(1:7, function(x)
        list()))
  ###repeating methods on the datasets
  for (j_stab in (1:n_iter_stab)) {
    ##setting seed
    set.seed(x + j_stab)
    ##measuring computation time
    start_time_meth <- Sys.time()
    ##ExWAS
    pred_iter <-
      list(ewas_BH = ewas(
        as.data.frame(simu$E_train),
        as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train)),
        corr = "BH"
      ))
    end_time_meth <- Sys.time()
    pred_iter$ewas_BH[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    print("ewas")
    ###LASSO
    start_time_meth <- Sys.time()
    predlasso <-
      lasso(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train))
      )
    end_time_meth <- Sys.time()
    predlasso[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    print("lasso")
    ###LASSO Meinshausen 1
    start_time_meth <- Sys.time()
    predlasso_stab_Meinshausen <-
      lasso_stab_Meinshausen(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train)),
        prop = 0.85
      )
    end_time_meth <- Sys.time()
    ###LASSO Mix
    predlasso_stab_Meinshausen[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    print("lasso_stab")
    start_time_meth <- Sys.time()
    predlasso_moy_Meinshausen <-
      lasso_moy_Meinshausen(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train)),
        prop = 0.5
      )
    end_time_meth <- Sys.time()
    predlasso_moy_Meinshausen[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    
    start_time_meth <- Sys.time()
    ###LASSO Meinshausen 2
    predlasso_stab_Meinshausen2 <-
      lasso_stab_Meinshausen2(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train)),
        prop = 0.95
      )
    end_time_meth <- Sys.time()
    predlasso_stab_Meinshausen2[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    
    start_time_meth <- Sys.time()
    ###LASSO CV2
    predlasso_moy_MSE <-
      lasso_moy_MSE(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train))
      )
    end_time_meth <- Sys.time()
    predlasso_moy_MSE[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    
    ###LASSO CV1
    start_time_meth <- Sys.time()
    predlasso_moy_lambda <-
      lasso_moy_lambda(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train))
      )
    end_time_meth <- Sys.time()
    predlasso_moy_lambda[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    
    
    ##ElasticNet
    start_time_meth <- Sys.time()
    predEnet <-
      Enet(
        data_Xs_in = as.data.frame(simu$E_train),
        data_Y_in = as.data.frame(simu$Y_train),
        colnames(as.data.frame(simu$Y_train))
      )
    end_time_meth <- Sys.time()
    predEnet[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    print("enet")
    
    start_time_meth <- Sys.time()
    
    predDSA <-
      DSAreg(
        Exp = simu$E_train,
        resp = simu$Y_train,
        maxsize = floor(ncol(simu$E_train) / 10),
        maxsumofpow = 1,
        maxorderint = 1
      )
    end_time_meth <- Sys.time()
    predDSA[[2]] <-
      as.numeric(difftime(end_time_meth, start_time_meth,  units = c("secs")))
    print("DSA effectué")
    ##combining results of the different methods for this iteration
    pred_iter <-
      c(
        pred_iter,
        lasso = list(predlasso),
        lasso_stab_Meinshausen = list(predlasso_stab_Meinshausen),
        lasso_stab_Meinshausen2 = list(predlasso_stab_Meinshausen2),
        lasso_moy_Meinshausen = list(predlasso_moy_Meinshausen),
        lasso_moy_MSE = list(predlasso_moy_MSE),
        lasso_moy_lambda = list(predlasso_moy_lambda),
        Enet = list(predEnet,
                    DSA = list(predDSA))
      )
    
    #pred_iter<-c(pred_iter,lasso=list(predlasso),
    #lasso_stab=list(predlasso_stab),
    #lasso_stab_plus_spec=list(predlasso_stab_ps),
    #Enet=list(predEnet),DSA=list(predDSA))
    #pred_iter<-c(pred_iter,lasso=list(predlasso)) #assessing performance of
    #each method
    truepred <- simu$yE$predictors
    for (k1 in (1:length(pred_iter))) {
      predfound <- pred_iter[[k1]]$selected
      if (exists("predfound") & exists("truepred")) {
        if (length(predfound) == 0) {
          print("no predictors found")
        }
        a <- sensitivity(truepred, predfound)
        b <- specificity(truepred, predfound, ncol(simu$E_train))
        c <- fdp(truepred, predfound)
        d <-
          estimatedR2(simu$E_test, predfound, simu$Y_test)$r.squared
        # print(a)
        # print(b)
        # print(c)
        # print(d)
        pred_iter[[k1]] <-
          c(
            pred_iter[[k1]],
            sens = a,
            spec = b,
            fdp = c,
            R2_test = d
          )
        remove(a)
        remove(b)
        remove
        remove(d)
        remove(predfound)
      }
      k1 <- k1 + 1
    }
    ##reshaping results saves
    for (k2 in (1:length(pred_iter))) {
      for (k3 in (1:length(pred_iter[[k2]]))) {
        predBMI_E[[k2]][[k3]] <-
          c(predBMI_E[[k2]][[k3]], list(pred_iter[[k2]][[k3]]))
        names(predBMI_E[[k2]]) <- names(pred_iter[[k2]])
      }
      names(predBMI_E) <- names(pred_iter)
    }
    print(j_stab)
    j_stab <- j_stab + 1
    
  }
  list_predBMI_E <- c(list_predBMI_E, list(predBMI_E))
  performance <- list()
  ##Assessing stability
  ###Assessing the frequency of selection of each hits
  for (k4 in (1:length(predBMI_E))) {
    tab_freq <-
      as.data.frame(table(table(as.character(
        unlist(predBMI_E[[k4]][[1]])
      ))))
    tab_freq$Var1 <-
      as.numeric(as.character(tab_freq$Var1)) / n_iter_stab
    colnames(tab_freq) <- c("freq", "nb_exp")
    tab_freq <- tab_freq[order(tab_freq$freq), ]
    ###Computing average Sorensen index as a measure of stability
    c_sor <- numeric(0)
    for (i1 in 1:(length(predBMI_E[[k4]][[1]]) - 1)) {
      for (i2 in (i1 + 1):length(predBMI_E[[k4]][[1]])) {
        if ((length(predBMI_E[[k4]][[1]][[i1]]) == 0) &
            length(predBMI_E[[k4]][[1]][[i2]]) == 0) {
          sor <- 1
        } else{
          sor <-
            sorensen(as.character(predBMI_E[[k4]][[1]][[i1]]), 
                     as.character((predBMI_E[[k4]][[1]][[i2]])))
        }
        c_sor <- c(c_sor, sor)
      }
    }
    
    sor_moy <- mean(c_sor)
    ##Saving performance measurement
    performance[[k4]] <-
      list(
        mean_nb_selec = mean(unlist(
          lapply(predBMI_E[[k4]][[1]], function(X)
            length(X))
        )),
        mean_sens = mean(unlist(predBMI_E[[k4]][[4]]), na.rm =
                           TRUE),
        mean_spec = mean(unlist(predBMI_E[[k4]][[5]]), na.rm =
                           TRUE),
        mean_fdp = mean(unlist(predBMI_E[[k4]][[6]]), na.rm =
                          TRUE),
        mean_R2_test = mean(unlist(predBMI_E[[k4]][[7]]), na.rm =
                              TRUE),
        nb_sup_20percent = sum(tab_freq$nb_exp[(which(tab_freq$freq >=
                                                        0.2))]),
        nb_sup_60percent = sum(tab_freq$nb_exp[(which(tab_freq$freq >=
                                                        0.6))]),
        sorensen = sor_moy,
        sd_nb_selec = sd(unlist(
          lapply(predBMI_E[[k4]][[1]], function(X)
            length(X))
        ), na.rm = TRUE),
        sd_sens = sd(unlist(predBMI_E[[k4]][[4]]), na.rm =
                       TRUE),
        sd_spec =   sd(unlist(predBMI_E[[k4]][[5]]), na.rm =
                         TRUE),
        sd_fdp =   sd(unlist(predBMI_E[[k4]][[6]]), na.rm =
                        TRUE),
        sd_R2_test =   sd(unlist(predBMI_E[[k4]][[7]]), na.rm =
                            TRUE),
        run_time <-
          mean(unlist(predBMI_E[[k4]][[2]]), na.rm = TRUE)
      )
    
    
  }
  remove(sor_moy)
  remove(tab_freq)
  names(performance) = names(predBMI_E)
  ##returning an object containing the simulated datasets with their
  ##characteristics, the results of each method and the performance measurements
  A <-
    list(simu = simu,
         performance = performance,
         list_predBMI_E = list_predBMI_E)
  remove(simu)
  remove(performance)
  remove(list_predBMI_E)
  gc()
  return(A)
  
  
}

########################################################
################# 5. SIMULATIONS #######################
########################################################

####loading input needed to generate datasets#####
dataExp_true <- readRDS("20190205 Exposome simu borne.rds")
##initialization
list_list_list_predBMI_E <- list()
list_list_performance <- list()
n_iter <- 30 ##number of iterations for each scenarios
n_iter_stab <- 15 ##number of iterations for each dataset
n_method <- 8 ##number of methods tested
c_n_R2_fixed <-
  c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8)
##values tested for R2 (variability of the outcome explained)
c_n_Ey <-
  c(1, 3, 10, 25, 50, 100) ##values tested for the number of predictors

#list_corr<-list(all=c(0,1),low=c(0,0.2),middle=c(0.2,0.4),high=c(0.5,1))
list_corr <- list(all = c(0, 1))
neg <- FALSE
#neg<-c(FALSE,TRUE)

seed = 22
##initalizing the dataset in which the performance of each method will be
##summarized by scenarios
comp_stability_method <- data.frame(
  Methods = vector(),
  Nb_true_predictors_of_BMI_in_E = numeric(0),
  Total_variability_of_BMI_explained_by_E =
    numeric(0),
  Mean_measured_variability_of_BMI_explained_by_E =
    numeric(0),
  Number_iterations = numeric(0),
  Number_iterations_stab = numeric(0),
  Mean_mean_number_predictors_found = numeric(0),
  Mean_mean_sensitivity = numeric(0),
  Mean_mean_specificity = numeric(0),
  Mean_mean_fdp = numeric(0),
  Mean_mean_R2_test = numeric(0),
  Mean_nb_sup_20percent = numeric(0),
  Mean_nb_sup_60percent = numeric(0),
  Mean_Sorensen = numeric(0),
  Mean_sd_number_predictors_found = numeric(0),
  Mean_sd_sensitivity = numeric(0),
  Mean_sd_specificity = numeric(0),
  Mean_sd_fdp = numeric(0),
  Mean_sd_R2_test = numeric(0),
  SD_mean_number_predictors_found = numeric(0),
  SD_mean_sensitivity = numeric(0),
  SD_mean_specificity = numeric(0),
  SD_mean_fdp = numeric(0),
  SD_mean_R2_test = numeric(0),
  SD_nb_sup_20percent = numeric(0),
  SD_nb_sup_60percent = numeric(0),
  SD_sd_number_predictors_found = numeric(0),
  SD_sd_sensitivity = numeric(0),
  SD_sd_specificity = numeric(0),
  SD_sd_fdp = numeric(0),
  SD_sd_R2_test = numeric(0),
  SD_Sorensen = numeric(0),
  which_scenario = numeric(0),
  Correlation = vector(),
  Negative_coefficient = vector(),
  Run_time = numeric(0)
)

##############################
####Run of the simulation#####
##############################
n = 1
##looping on each scenario ie looping on each list of parameters which define
##scenarios
for (i4 in 1:length(neg)) {
  pos_and_neg <- neg[i4]
  for (i3 in 1:length(list_corr)) {
    corr_range <- list_corr[[i3]]
    corr = TRUE
    if (corr_range[1] == 0 & corr_range[2] == 1) {
      corr = FALSE
    }
    for (i1 in 1:length(c_n_Ey)) {
      n_Ey <- c_n_Ey[i1]
      for (i2 in 1:length(c_n_R2_fixed)) {
        R2_fixed <- c_n_R2_fixed[i2]
        
        ##inside the loops a scenario is set
        
        test_correl <- tryCatch({
          simu <-
            simulator_2layers(
              E_true = dataExp_true,
              R2_tot = R2_fixed,
              n_Ey = n_Ey,
              BetaEy = 0.1,
              test_and_training = TRUE,
              pos_and_neg = pos_and_neg,
              corr = corr,
              range_corr = corr_range
            )
        }, silent = FALSE)
        if (class(test_correl) == "try-error") {
          print(n)
          
        } else{
          n_row = nrow(comp_stability_method)
          list_list_predBMI_E <- list()
          list_performance <- list()
          print("cluster")
          
          ##parallelization
          start_time <- Sys.time()
          cl <-
            makeCluster(getOption("cl.cores", round(detectCores())))
          clusterExport(
            cl,
            list(
              "dataExp_true",
              "simulator_2layers",
              "simResponseSimple",
              "estimatedR2",
              "getresiduals_2df",
              "ewas",
              "lasso",
              "lasso_stab_Meinshausen2",
              "lasso_moy_lambda",
              "lasso_moy_MSE",
              "lasso_moy_Meinshausen",
              "lasso_stab_Meinshausen",
              "Enet",
              "wqs",
              "sensitivity",
              "fdp",
              "specificity",
              "f0",
              "R2_fixed",
              "n_Ey",
              "n_iter_stab",
              "submatFindSimpl",
              "DSAreg",
              "corr",
              "corr_range",
              "pos_and_neg",
              "n_method"
            )
          )
          clusterEvalQ(
            cl,
            list(
              library("boot"),
              library("reshape"),
              library("glmnet"),
              library("DSA"),
              library("mvtnorm"),
              library("gWQS"),
              library("OmicsMarkeR"),
              library("Rcpp")
            )
          )
          ##applying f0 in parallel
          results_1_jeu <- clusterApply(cl, 1:n_iter, f0)
          stopCluster(cl)
          ##getting results for this scenario
          simulated_data <- lapply(results_1_jeu, function(x)
            x$simu)
          ##formatting results for this scenario
          list_list_predBMI_E <-
            lapply(results_1_jeu, function(x)
              x$list_predBMI_E)
          list_performance <- lapply(results_1_jeu, function(x)
            x$performance)
          remove(results_1_jeu)
          
          ##saving simulation parameters for this scenario
          param_simu <-
            data.frame(
              Parameters = vector(),
              Fixed_or_measured = vector(),
              Value = numeric(0)
            )
          param_simu[1, ] <-
            c("Nb_true predictors of BMI in E",
              "Fixed",
              mean(unlist(
                lapply(simulated_data, function(X)
                  length(X$yE$beta))
              )))
          param_simu[2, ] <-
            c("Total variability of BMI explained by E",
              "Fixed",
              R2_fixed)
          param_simu[3, ] <-
            c("Mean variability of BMI explained by E",
              "Measured",
              mean(unlist(
                lapply(simulated_data, function(X)
                  (X$R2))
              )))
          param_simu[4, ] <- c("Number_iterations", "Fixed", n_iter)
          param_simu[5, ] <-
            c("Number_iterations_on_a_same_dataset",
              "Fixed",
              n_iter_stab)
          param_simu[6, ] <-
            c("Correlation_range",
              "Fixed",
              str_c(
                as.character(corr_range[1]),
                as.character(corr_range[2])
              ))
          param_simu[7, ] <-
            c("Negative_coefficients",
              "Fixed",
              ifelse(pos_and_neg == TRUE, "yes", "no"))
          
          ##saving methods performance for this scenario
          for (i1 in (1:length(list_performance[[1]]))) {
            comp_stability_method[i1 + n_row, ] <-
              c(
                as.character(names(list_performance[[1]][i1])),
                param_simu[1, 3],
                R2_fixed,
                param_simu[3, 3],
                n_iter,
                n_iter_stab,
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[1]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[2]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[3]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[4]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[5]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[6]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[7]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[8]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[9]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[10]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[11]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[12]]))
                ), na.rm = TRUE),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[13]]))
                ), na.rm = TRUE),
                
                
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[1]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[2]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[3]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[4]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[5]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[6]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[7]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[8]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[9]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[10]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[11]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[12]]))
                ), na.rm = TRUE),
                sd(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[13]]))
                ), na.rm = TRUE),
                n,
                param_simu[6, 3],
                ifelse(pos_and_neg == TRUE, "yes", "no"),
                mean(unlist(
                  lapply(list_performance, function(X)
                    (X[[i1]][[14]]))
                ), na.rm = TRUE)
                
                
              )
            i1 <- i1 + 1
          }
          
          print(n)
          
          list_list_list_predBMI_E <-
            c(list_list_list_predBMI_E,
              list(list_list_predBMI_E))
          end_time <- Sys.time()
          end_time - start_time
          ##saving results externaly
          saveRDS(comp_stability_method,
                  "comp_stability_method_n30_15.Rds")
          saveRDS(
            simulated_data,
            file = paste(n,
                         '_simulated_data_scenario_iteration_stability_n30_15.Rds')
          )
          saveRDS(
            list_list_list_predBMI_E,
            "list_list_list_predBMI_E_stability_n30_15.Rds"
          )
          saveRDS(list_list_performance,
                  "list_list_performance_stability_n30_15.Rds")
          remove(simulated_data)
          remove(list_list_predBMI_E)
          remove(list_performance)
          
        }
        n = n + 1
        
      }
    }
  }
}


