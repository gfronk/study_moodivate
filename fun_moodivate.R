# fun_moodivate.R

# this script contains functions that are called 
# in the moodivate analysis workflow (.qmd files)

make_splits <- function(d, cv_resample_type, 
                        cv_resample = NULL, 
                        cv_outer_resample = NULL, 
                        cv_inner_resample = NULL, 
                        the_seed = NULL,
                        cv_group = NULL) {
  
  # d: (training) dataset to be resampled 
  # cv_resample_type: can be boot, kfold, or nested
  # resample: specifies for repeats and folds for CV (1_x_10; 10_x_10) or num splits for bootstrapping (100)
  # inner_resample: specifies repeats/folds or num bootstrap splits for nested cv inner loop - same format as above
  # outer_resample: specifies repeats/folds for outer nested cv loop - cannot use bootstrapping here

  if(is.null(the_seed)) {
    error("make_splits() requires a seed")
  } else set.seed(the_seed)
  
  # bootstrap splits
  if (cv_resample_type == "boot") {
    splits <- d %>% 
      bootstraps(times = cv_resample)
  }
  
  # kfold - includes repeated kfold
  if (cv_resample_type == "kfold") {
    # get number of repeats and folds
    n_repeats <- as.numeric(str_remove(cv_resample, "_x_\\d{1,2}"))
    n_folds <- as.numeric(str_remove(cv_resample, "\\d{1,3}_x_"))
    
      splits <- d %>% 
        vfold_cv(v = n_folds, repeats = n_repeats) 
  }
  
  
  # nested   
  if (cv_resample_type == "nested") {
    # get number of repeats and folds for outer cv loop
    outer_n_repeats <- as.numeric(str_remove(cv_outer_resample, "_x_\\d{1,2}"))
    outer_n_folds <- as.numeric(str_remove(cv_outer_resample, "\\d{1,3}_x_"))
    
    # get repeats/folds or bootstrap splits for inner loop
    if (str_detect(cv_inner_resample, "_x_")) {
      inner_n_repeats <- as.numeric(str_remove(cv_inner_resample, "_x_\\d{1,2}"))
      inner_n_folds <- as.numeric(str_remove(cv_inner_resample, "\\d{1,3}_x_"))
    } else {
      inner_boot_splits <- cv_inner_resample
    }
    
    # create splits for ungrouped nested cv with kfold inner
    if (is.null(cv_group) & str_detect(cv_inner_resample, "_x_")) {
      
      splits <- d %>% 
        nested_cv(outside = vfold_cv(v = outer_n_folds, repeats = outer_n_repeats) , 
                  inside = vfold_cv(v = inner_n_folds, repeats = inner_n_repeats))
    }

  }
  
  return(splits)
}

make_rset <- function(splits, cv_resample_type, split_num = NULL, 
                      inner_split_num = NULL, outer_split_num = NULL) {
  # used to make an rset object that contains a single split for use in 
  # tuning glmnet on CHTC  
  
  if (cv_resample_type == "nested") {
    split <- splits$inner_resamples[[outer_split_num]] %>% 
      slice(inner_split_num) 
  }
  
  if (cv_resample_type == "kfold") {
    split <- splits %>% 
      slice(split_num)
  }
  
  if (cv_resample_type == "boot") {
    stop("Make rset does not work for bootstrap resamples")
  }
  
  rset <- manual_rset(split$splits, split$id)
  return(rset)
}

tune_model <- function(config, rec, splits, ml_mode, cv_resample_type, hp2_glmnet_min = NULL,
                       hp2_glmnet_max = NULL, hp2_glmnet_out = NULL, y_level_pos = NULL) {
  # config: single-row config-specific tibble from jobs
  # splits: rset object that contains all resamples
  # rec: recipe (created manually or via build_recipe() function)
  
  # set metrics for regression or classification
  if (ml_mode == "regression") {
    mode_metrics <- metric_set(mae, rmse, rsq)
  }
  
  if (ml_mode == "classification") {
    mode_metrics <- metric_set(roc_auc, accuracy,
                               sens, yardstick::spec, ppv, npv)
  }
  
  if (algorithm == "glmnet") {
    grid_penalty <- expand_grid(penalty = exp(seq(hp2_glmnet_min, hp2_glmnet_max, 
                                                  length.out = hp2_glmnet_out)))
    
    # make rset for single held-in/held_out split
    # does not work for bootstrapping
    split <- make_rset(splits, cv_resample_type = cv_resample_type, 
                       inner_split_num = config$inner_split_num, 
                       outer_split_num = config$outer_split_num)
    
    # backward compatible for tune controls that didnt set family 
    if (!exists("glm_family")) glm_family <- if_else(ml_mode == "regression", "gaussian", "binomial")
    
    if (ml_mode == "classification") {
      models <- logistic_reg(penalty = tune(),
                             mixture = config$hp1) %>%
        set_engine("glmnet", glm_family = glm_family) %>%
        set_mode("classification") %>%
        tune_grid(preprocessor = rec,
                  resamples = split,
                  grid = grid_penalty,
                  # metrics assume that positive event it first level
                  # make sure this is true in recipe
                  metrics = mode_metrics)
    } else {
      
      models <- linear_reg(penalty = tune(),
                           mixture = config$hp1) %>%
        set_engine("glmnet", family = glm_family) %>%
        set_mode("regression") %>%
        tune_grid(preprocessor = rec,
                  resamples = split,
                  grid = grid_penalty,
                  # metrics assume that positive event is first level
                  # make sure this is true in recipe
                  metrics = mode_metrics)
    }
    # create tibble of penalty and metrics returned 
    results <- collect_metrics(models, summarize = FALSE) %>% 
      rename(hp2 = penalty) %>% 
      select(hp2, .metric, .estimate) %>% 
      pivot_wider(., names_from = ".metric",
                  values_from = ".estimate") %>%  
      bind_cols(config %>% select(-hp2), .) %>% 
      relocate(hp2, .after = hp1) 
    
    return(results)
  }
  
}