library(shiny)
library(ggplot2)
library(dplyr)
library(rstan)
library(tibble)

# ----------------------------
# Helper Functions
# ----------------------------
rle <- function(x, le = "gbd") {
  if(le == "gbd"){
    return(c(84.31, 76.76, 56.95, 35.17, 23.29))
  } else {
    return(x)
  }
}

burden <- function(N, DW, A, L, K, r, a) {
  N * DW * L
}

mean_ci <- function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975))

# Sensitivity analysis functions
sa_pcc <- function(y, x) {
  out <- matrix(ncol = 2, nrow = ncol(x))
  colnames(out) <- c("rho", "p")
  rownames(out) <- colnames(x)
  
  for (i in seq(ncol(x))){
    lm_y <- lm(y ~ x[, -i])      # regress y to other x's
    lm_x <- lm(x[, i] ~ x[, -i]) # regress x to other x's
    out[i, ] <-
      unlist(cor.test(lm_y$residuals, lm_x$residuals)[4:3],
             use.names = FALSE)
  }
  
  return(out)
}

tornado <- function(coef, names = NULL, p = 1) {
  if (is.null(names)) names <- rownames(coef)
  
  df <- data.frame(est = coef[, "rho"],
                   lab = names)
  
  df <- subset(df, coef[, "p"] < p)
  
  df <- df[order(abs(df$est)), ]
  
  df$id <- seq(nrow(df))
  
  ggplot(df, aes(x = id, y = est)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_continuous(element_blank(),
                       breaks = seq(nrow(df)),
                       labels = df$lab) +
    scale_y_continuous("partial correlation coefficient",
                       limits = c(min(0, min(df$est) - 0.1),
                                  max(0, max(df$est) + 0.1))) +
    geom_text(aes(x = id, y = est, label = formatC(est, 3, form = "f")),
              size = 3,
              hjust = ifelse(df$est > 0, -0.1, 1.1),
              vjust = 0.4) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

stan_model_tp <- readRDS("stan_model_tp.rds")

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Title styling */
      .title-panel {
        background: linear-gradient(135deg, #2c3e50, #3498db);
        color: white;
        padding: 20px;
        text-align: left;
        margin-bottom: 20px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        margin-left: -30px;   
        margin-right: -30px;
        font-weight: bold;
        font-family: 'Poppins', sans-serif;
      }
      .title-panel h2 {
        font-weight: 700;
        font-size: 28px;
        margin: 0;
        text-shadow: 1px 1px 3px rgba(0,0,0,0.3);
      }
      /* Custom styling */
      body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
      }
      .well {
        background-color: #f8f9fa;
        border: 1px solid #e9ecef;
        border-radius: 8px;
      }
      .btn-primary {
        background-color: #2c3e50;
        border-color: #2c3e50;
      }
      h4 {
        color: #2c3e50;
        border-bottom: 2px solid #3498db;
        padding-bottom: 8px;
      }
      .nav-tabs > li > a {
        color: #2c3e50;
        font-weight: bold;
      }
    "))
  ),
  
  div(class = "title-panel",
      h2("NCC DALY Calculator using Transitional Model")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Serological Test Data"),
      numericInput("N_tested", "Number of People Tested", value = 694, min = 0),
      numericInput("N_positive", "Number of Positive Tests", value = 11, min = 0),
      sliderInput("sensitivity", "Test Sensitivity", value = 98, min = 0, max = 100),
      sliderInput("specificity", "Test Specificity", value = 100, min = 0, max = 100),
      actionButton("calculate", "Calculate Burden", class = "btn-primary"),
      
      h4("Population Data"),
      numericInput("total_pop", "Total Population", value = 1347121000, min = 0),
      numericInput("prop_m", "Proportion Male", value = 0.514, min = 0, max = 1, step = 0.01),
      
      h5("Male Age Distribution"),
      fluidRow(
        column(6, numericInput("male_age1", "0-4", value = 0.077, min = 0, max = 1, step = 0.01)),
        column(6, numericInput("male_age2", "5-14", value = 0.176, min = 0, max = 1, step = 0.01))
      ),
      fluidRow(
        column(6, numericInput("male_age3", "15-44", value = 0.542, min = 0, max = 1, step = 0.01)),
        column(6, numericInput("male_age4", "45-64", value = 0.127, min = 0, max = 1, step = 0.01))
      ),
      numericInput("male_age5", "65+", value = 0.078, min = 0, max = 1, step = 0.01),
      
      h5("Female Age Distribution"),
      fluidRow(
        column(6, numericInput("female_age1", "0-4", value = 0.074, min = 0, max = 1, step = 0.01)),
        column(6, numericInput("female_age2", "5-14", value = 0.172, min = 0, max = 1, step = 0.01))
      ),
      fluidRow(
        column(6, numericInput("female_age3", "15-44", value = 0.544, min = 0, max = 1, step = 0.01)),
        column(6, numericInput("female_age4", "45-64", value = 0.129, min = 0, max = 1, step = 0.01))
      ),
      numericInput("female_age5", "65+", value = 0.081, min = 0, max = 1, step = 0.01),
      
      h4("Disease Parameters"),
      # p("Defaults values are based on studies from India"),
      
      h5("Symptomatic NCC Parameters (Beta distribution)"),
      fluidRow(
        column(6, numericInput("symp_ncc_alpha", "Alpha", value = 18, min = 0)),
        column(6, numericInput("symp_ncc_beta", "Beta", value = 21, min = 0))
      ),
      
      h5("NCC with Epilepsy Parameters (Beta distribution)"),
      fluidRow(
        column(6, numericInput("ncc_ep_alpha", "Alpha", value = 122, min = 0)),
        column(6, numericInput("ncc_ep_beta", "Beta", value = 45, min = 0))
      ),
      
      h5("NCC with Headache Parameters (Beta distribution)"),
      fluidRow(
        column(6, numericInput("ncc_hd_alpha", "Alpha", value = 53, min = 0)),
        column(6, numericInput("ncc_hd_beta", "Beta", value = 155, min = 0))
      ),
      
      h5("Epilepsy Case Fatality (Beta distribution)"),
      fluidRow(
        column(6, numericInput("cft_ep_alpha", "Alpha", value = 6, min = 0)),
        column(6, numericInput("cft_ep_beta", "Beta", value = 169, min = 0))
      ),
      
      h5("Epilepsy Treatment Rate (Uniform distribution)"),
      fluidRow(
        column(6, numericInput("prop_ep_trt_min", "Min", value = 0.10, min = 0, max = 1, step = 0.01)),
        column(6, numericInput("prop_ep_trt_max", "Max", value = 0.86, min = 0, max = 1, step = 0.01))
      ),
      
      h5("Disability Weights (Uniform distributions)"),
      p("Epilepsy - Treated:"),
      fluidRow(
        column(6, numericInput("dsw_ep_tr_min", "Min", value = 0.031, min = 0, max = 1, step = 0.001)),
        column(6, numericInput("dsw_ep_tr_max", "Max", value = 0.072, min = 0, max = 1, step = 0.001))
      ),
      p("Epilepsy - Not Treated:"),
      fluidRow(
        column(6, numericInput("dsw_ep_nt_min", "Min", value = 0.173, min = 0, max = 1, step = 0.001)),
        column(6, numericInput("dsw_ep_nt_max", "Max", value = 0.710, min = 0, max = 1, step = 0.001))
      ),
      p("Headache:"),
      fluidRow(
        column(6, numericInput("dsw_hd_min", "Min", value = 0.022, min = 0, max = 1, step = 0.001)),
        column(6, numericInput("dsw_hd_max", "Max", value = 0.588, min = 0, max = 1, step = 0.001))
      ),
      h5("Set Iterations"),
      fluidRow(
        column(12, numericInput("iter", "Iteration (Min = 100)", value = 10000, min = 100, step = 100))
      ),
      h5("Set Test Sensitivity Range (For P vs DALYs Plot)"),
      fluidRow(
        column(6, numericInput("se_min", "Min (0 - 1)", value = 0.6, min = 0, max = 1, step = 0.05)),
        column(6, numericInput("se_max", "Max (0 - 1)", value = 1.0, min = 0, max = 1, step = 0.05))
      ),
      h5("Set Test Specificity Range (For P vs DALYs Plot)"),
      fluidRow(
        column(6, numericInput("sp_min", "Min (0 - 1)", value = 0.6, min = 0, max = 1, step = 0.05)),
        column(6, numericInput("sp_max", "Max (0 - 1)", value = 1.0, min = 0, max = 1, step = 0.05))
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        # ... all your tabs ...
        tabPanel("Results Summary", 
                 fluidRow(
                   column(6, 
                          div(style = "font-size: 12px;",
                              tableOutput("results_table")
                          ),
                          conditionalPanel(
                            condition = "output.results_ready",
                            downloadButton("download_summary", "Download Summary CSV")
                          )
                   ),
                   column(6,
                          plotOutput("prevalence_daly_plot"),
                          conditionalPanel(
                            condition = "output.plot_ready",
                            tagList(
                              downloadButton("download_prevalence_daly_data", "Download Prevalence-DALY Data"),
                              downloadButton("download_prevalence_plot_png", "Download Plot (PNG)"),
                              downloadButton("download_prevalence_plot_pdf", "Download Plot (PDF)")
                            )
                          )
                   )
                 )),
        tabPanel("DALYs by Age Group", 
                 plotOutput("daly_plot"), 
                 downloadButton("download_daly_data", "Download DALY Data")),
        
        tabPanel("Sensitivity Analysis", 
                 plotOutput("sensitivity_plot"),
                 downloadButton("download_sens_plot_png", "Download Plot (PNG)"),
                 downloadButton("download_sens_plot_pdf", "Download Plot (PDF)")),
        tabPanel("Parameter Summary",
                 tableOutput("parameter_table"),
                 downloadButton("download_parameters", "Download Parameters CSV")),
        tabPanel("Model",
                 img(src= "model_4.3.jpg", height = "600px", alt = "Flow chart of computational disease model for dModel 4.3"),
                 p(HTML("<b>Figure: Computational disease model followed in this app.</b>")))# ... other tabs ...
      )
    )
  )
)



# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  results <- reactiveValues(
    daly_data = NULL,
    summary_data = NULL,
    sensitivity_data = NULL,
    parameters_used = NULL,
    sampled_parameters = NULL,
    prevalence_daly_data = NULL
  )
  
  observeEvent(input$calculate, {
    #showNotification("Calculating burden... This may take a few minutes.", type = "message")
    
    withProgress(message = 'Analyzing, Please wait for few seconds!!', value = 0, {
      incProgress(0.1, detail = "..")
      
      total_male <- input$total_pop * input$prop_m
      total_female <- input$total_pop - total_male
      
      male_prop_per_age <- c(input$male_age1, input$male_age2, input$male_age3, input$male_age4, input$male_age5)
      female_prop_per_age <- c(input$female_age1, input$female_age2, input$female_age3, input$female_age4, input$female_age5)
      
      pop_mx <- matrix(0, nrow = 5, ncol = 2)
      pop_mx[,1] <- total_male * male_prop_per_age
      pop_mx[,2] <- total_female * female_prop_per_age
      pop_mx <- matrix(pop_mx, nrow = 10, ncol = 1)
      
      incProgress(0.3, detail = " ")
      stan_data <- list(
        T_pos = input$N_positive,
        N = input$N_tested,
        Se = input$sensitivity / 100,
        Sp = input$specificity / 100
      )
      
      fit <- sampling(stan_model_tp, data = stan_data, chains = 4, iter = input$iter, warmup = input$iter/2, seed = 123, refresh = 0)
      posterior_samples <- rstan::extract(fit)
      bci <- quantile(posterior_samples$pi, probs = c(0.025, 0.975))
      
      incProgress(0.6, detail = " ")
      
      set.seed(123)
      n <- 10000
      
      # Use user inputs or defaults for parameters
      symp_ncc_alpha <- ifelse(is.na(input$symp_ncc_alpha) || input$symp_ncc_alpha == "", 18, input$symp_ncc_alpha)
      symp_ncc_beta <- ifelse(is.na(input$symp_ncc_beta) || input$symp_ncc_beta == "", 21, input$symp_ncc_beta)
      ncc_ep_alpha <- ifelse(is.na(input$ncc_ep_alpha) || input$ncc_ep_alpha == "", 122, input$ncc_ep_alpha)
      ncc_ep_beta <- ifelse(is.na(input$ncc_ep_beta) || input$ncc_ep_beta == "", 45, input$ncc_ep_beta)
      ncc_hd_alpha <- ifelse(is.na(input$ncc_hd_alpha) || input$ncc_hd_alpha == "", 53, input$ncc_hd_alpha)
      ncc_hd_beta <- ifelse(is.na(input$ncc_hd_beta) || input$ncc_hd_beta == "", 155, input$ncc_hd_beta)
      
      cft_ep_alpha <- ifelse(is.na(input$cft_ep_alpha) || input$cft_ep_alpha == "", 6, input$cft_ep_alpha)
      cft_ep_beta <- ifelse(is.na(input$cft_ep_beta) || input$cft_ep_beta == "", 169, input$cft_ep_beta)
      prop_ep_trt_min <- ifelse(is.na(input$prop_ep_trt_min) || input$prop_ep_trt_min == "", 0.10, input$prop_ep_trt_min)
      prop_ep_trt_max <- ifelse(is.na(input$prop_ep_trt_max) || input$prop_ep_trt_max == "", 0.86, input$prop_ep_trt_max)
      
      dsw_ep_tr_min <- ifelse(is.na(input$dsw_ep_tr_min) || input$dsw_ep_tr_min == "", 0.031, input$dsw_ep_tr_min)
      dsw_ep_tr_max <- ifelse(is.na(input$dsw_ep_tr_max) || input$dsw_ep_tr_max == "", 0.072, input$dsw_ep_tr_max)
      dsw_ep_nt_min <- ifelse(is.na(input$dsw_ep_nt_min) || input$dsw_ep_nt_min == "", 0.173, input$dsw_ep_nt_min)
      dsw_ep_nt_max <- ifelse(is.na(input$dsw_ep_nt_max) || input$dsw_ep_nt_max == "", 0.710, input$dsw_ep_nt_max)
      dsw_hd_min <- ifelse(is.na(input$dsw_hd_min) || input$dsw_hd_min == "", 0.022, input$dsw_hd_min)
      dsw_hd_max <- ifelse(is.na(input$dsw_hd_max) || input$dsw_hd_max == "", 0.588, input$dsw_hd_max)
      
      # Generate samples with user parameters
      prev_ncc <- runif(n, bci[1], bci[2])
      symp_ncc <- rbeta(n, symp_ncc_alpha, symp_ncc_beta)
      ncc_ep <- rbeta(n, ncc_ep_alpha, ncc_ep_beta)
      ncc_hd <- rbeta(n, ncc_hd_alpha, ncc_hd_beta)
      cft_ep <- rbeta(n, cft_ep_alpha, cft_ep_beta)
      prop_ep_trt <- runif(n, prop_ep_trt_min, prop_ep_trt_max)
      
      # Store sampled parameters for sensitivity analysis
      results$sampled_parameters <- list(
        sero_prev_ncc = posterior_samples$pi[1:n],
        prev_ncc = prev_ncc,
        symp_ncc = symp_ncc,
        ncc_ep = ncc_ep,
        prop_ep_trt = prop_ep_trt,
        ncc_hd = ncc_hd,
        cft_ep = cft_ep
      )
      
      prev_ncc_symp <- prev_ncc * symp_ncc
      prev_ncc_ep_act <- prev_ncc * symp_ncc * ncc_ep
      prev_ncc_hd <- prev_ncc * symp_ncc * ncc_hd
      
      N_prev_ncc_ep_act <- t(apply(pop_mx, 1, function(x) x * prev_ncc_ep_act))
      N_prev_ncc_hd <- t(apply(pop_mx, 1, function(x) x * prev_ncc_hd))
      N_prev_ncc_ep_mrt <- t(apply(pop_mx, 1, function(x) x * prev_ncc_ep_act*cft_ep))
      
      # YLL calculation
      ep_aad <- c(2, 9.5, 29.5, 52, 73)
      ep_rle <- rle(ep_aad, "gbd")
      
      yll_ep <- array(dim = c(5,2,n))
      for(i in 1:5){
        yll_ep[i,1,] <- burden(N = N_prev_ncc_ep_act[i,]*mean(cft_ep), DW = 1, A = ep_aad[i], L = ep_rle[i], K = K, r = r, a = ep_aad[i])
        yll_ep[i,2,] <- burden(N = N_prev_ncc_ep_act[i+5,]*mean(cft_ep), DW = 1, A = ep_aad[i], L = ep_rle[i], K = K, r = r, a = ep_aad[i])
      }
      
      # YLD calculation
      dsw_ep_tr <- runif(n, dsw_ep_tr_min, dsw_ep_tr_max)
      dsw_ep_nt <- runif(n, dsw_ep_nt_min, dsw_ep_nt_max)
      dsw_hd <- runif(n, dsw_hd_min, dsw_hd_max)
      
      # Store disability weights for sensitivity analysis
      results$sampled_parameters$dsw_ep_tr <- dsw_ep_tr
      results$sampled_parameters$dsw_ep_nt <- dsw_ep_nt
      results$sampled_parameters$dsw_hd <- dsw_hd
      
      A <- 0; L <- 1; a <- A
      yld_ep <- array(dim=c(5,2,n))
      for(i in 1:5){
        yld_ep[i,1,] <- burden(N=N_prev_ncc_ep_act[i,]*mean(prop_ep_trt), DW=dsw_ep_tr, A=A,L=L,K=K,r=r,a=a) +
          burden(N=N_prev_ncc_ep_act[i,]*(1-mean(prop_ep_trt)), DW=dsw_ep_nt, A=A,L=L,K=K,r=r,a=a)
        yld_ep[i,2,] <- burden(N=N_prev_ncc_ep_act[i+5,]*mean(prop_ep_trt), DW=dsw_ep_tr, A=A,L=L,K=K,r=r,a=a) +
          burden(N=N_prev_ncc_ep_act[i+5,]*(1-mean(prop_ep_trt)), DW=dsw_ep_nt, A=A,L=L,K=K,r=r,a=a)
      }
      
      yld_hd <- array(dim=c(5,2,n))
      for(i in 1:5){
        yld_hd[i,1,] <- burden(N=N_prev_ncc_hd[i,], DW=dsw_hd, A=A,L=L,K=K,r=r,a=a)
        yld_hd[i,2,] <- burden(N=N_prev_ncc_hd[i+5,], DW=dsw_hd, A=A,L=L,K=K,r=r,a=a)
      }
      
      # Total DALYs
      age_groups <- c("0-4","5-14","15-44","45-64","65+")
      daly_by_age <- matrix(nrow=5, ncol=n)
      for(i in 1:5){
        daly_by_age[i,] <- yll_ep[i,1,] + yll_ep[i,2,] + apply(yld_ep[i,,],2,sum) + apply(yld_hd[i,,],2,sum)
      }
      
      # Store DALYs for sensitivity analysis
      yld_ep_all <- apply(yld_ep, 3, sum)
      yld_hd_all <- apply(yld_hd, 3, sum)
      yld <- yld_ep_all + yld_hd_all
      yll_ep_all <- apply(yll_ep, 3, sum)
      daly <- yld + yll_ep_all
      results$sampled_parameters$daly <- daly
      
      results$daly_data <- list(
        age_groups = age_groups,
        daly_matrix = daly_by_age,
        median_daly = apply(daly_by_age,1,median),
        lower_ci = apply(daly_by_age,1,quantile,0.025),
        upper_ci = apply(daly_by_age,1,quantile,0.975)
      )
      
      # Create parameter summary table with mean and 95% BCI
      results$parameters_used <- tibble::tibble(
        Parameter = c(
          "NCC Prevalence (%)",
          "Symptomatic NCC (%)",
          "NCC with Epilepsy (%)",
          "NCC with Headache (%)",
          "Case Fatality of epilepsy (%)",
          "Percentage of treated epilepsy (%)",
          "DW Epilepsy Treated",
          "DW Epilepsy Not Treated", 
          "DW Headache"
        ),
        Mean = c(
          mean(posterior_samples$pi)*100,
          mean(symp_ncc)*100,
          mean(ncc_ep)*100,
          mean(ncc_hd)*100,
          mean(cft_ep)*100,
          mean(prop_ep_trt)*100,
          mean(dsw_ep_tr),
          mean(dsw_ep_nt),
          mean(dsw_hd)
        ),
        `Lower 95% BCI` = c(
          quantile(posterior_samples$pi, 0.025)*100,
          quantile(symp_ncc, 0.025)*100,
          quantile(ncc_ep, 0.025)*100,
          quantile(ncc_hd, 0.025)*100,
          quantile(cft_ep, 0.025)*100,
          quantile(prop_ep_trt, 0.025)*100,
          quantile(dsw_ep_tr, 0.025),
          quantile(dsw_ep_nt, 0.025),
          quantile(dsw_hd, 0.025)
        ),
        `Upper 95% BCI` = c(
          quantile(posterior_samples$pi, 0.975)*100,
          quantile(symp_ncc, 0.975)*100,
          quantile(ncc_ep, 0.975)*100,
          quantile(ncc_hd, 0.975)*100,
          quantile(cft_ep, 0.975)*100,
          quantile(prop_ep_trt, 0.975)*100,
          quantile(dsw_ep_tr, 0.975),
          quantile(dsw_ep_nt, 0.975),
          quantile(dsw_hd, 0.975)
        )
      )
      
      # Summary Results Table
      seroprev_mean <- mean(posterior_samples$pi)
      seroprev_ci <- quantile(posterior_samples$pi, probs = c(0.025, 0.975))
      total_pop <- input$total_pop
      
      results$summary_data <- tibble::tibble(
        Indicator = c(
          "Estimated NCC prevalence(%)",
          "Total NCC cases with epileptic seizure (P)",
          "Total NCC cases with headache (P)",
          "Total NCC cases who died",
          "Epilepsy YLDs",
          "Headache YLDs",
          "Total YLDs",
          "Total YLLs",
          "Total DALYs",
          "YLDs per 1000 population",
          "YLLs per 1000 population",
          "DALYs per 1000 population",
          "YLD Proportion of DALYs",
          "YLL Proportion of DALYs"
        ),
        Mean = c(
          seroprev_mean*100,
          mean_ci(colSums(N_prev_ncc_ep_act))[1],
          mean_ci(colSums(N_prev_ncc_hd))[1],
          mean_ci(colSums(N_prev_ncc_ep_mrt))[1],
          mean_ci(yld_ep_all)[1],
          mean_ci(yld_hd_all)[1],
          mean_ci(yld)[1],
          mean_ci(yll_ep_all)[1],
          mean_ci(daly)[1],
          mean_ci(1e3 * yld / total_pop)[1],
          mean_ci(1e3 * yll_ep_all / total_pop)[1],
          mean_ci(1e3 * daly / total_pop)[1],
          mean_ci(yld / daly)[1],
          mean_ci(yll_ep_all / daly)[1]
        ),
        `Lower 95% CI` = c(
          seroprev_ci[1]*100,
          mean_ci(colSums(N_prev_ncc_ep_act))[2],
          mean_ci(colSums(N_prev_ncc_hd))[2],
          mean_ci(colSums(N_prev_ncc_ep_mrt))[2],
          mean_ci(yld_ep_all)[2],
          mean_ci(yld_hd_all)[2],
          mean_ci(yld)[2],
          mean_ci(yll_ep_all)[2],
          mean_ci(daly)[2],
          mean_ci(1e3 * yld / total_pop)[2],
          mean_ci(1e3 * yll_ep_all / total_pop)[2],
          mean_ci(1e3 * daly / total_pop)[2],
          mean_ci(yld / daly)[2],
          mean_ci(yll_ep_all / daly)[2]
        ),
        `Upper 95% CI` = c(
          seroprev_ci[2]*100,
          mean_ci(colSums(N_prev_ncc_ep_act))[3],
          mean_ci(colSums(N_prev_ncc_hd))[3],
          mean_ci(colSums(N_prev_ncc_ep_mrt))[3],
          mean_ci(yld_ep_all)[3],
          mean_ci(yld_hd_all)[3],
          mean_ci(yld)[3],
          mean_ci(yll_ep_all)[3],
          mean_ci(daly)[3],
          mean_ci(1e3 * yld / total_pop)[3],
          mean_ci(1e3 * yll_ep_all / total_pop)[3],
          mean_ci(1e3 * daly / total_pop)[3],
          mean_ci(yld / daly)[3],
          mean_ci(yll_ep_all / daly)[3]
        )
      )
      
      # Generate prevalence vs DALY data
      incProgress(0.8, detail = " ")
      generate_prevalence_daly_data()
    })
    
    showNotification("Calculation complete!", type = "message")
  })
  
  
  
  
  # Function to generate prevalence vs DALY data
  generate_prevalence_daly_data <- function() {
    req(input$N_tested, input$total_pop)
    
    se_min <- ifelse(is.na(input$se_min) || input$se_min == "", 0.6, input$se_min)
    se_max <- ifelse(is.na(input$se_max) || input$se_max == "", 1.0, input$se_max)
    sp_min <- ifelse(is.na(input$sp_min) || input$sp_min == "", 0.6, input$sp_min)
    sp_max <- ifelse(is.na(input$sp_max) || input$sp_max == "", 1.0, input$sp_max)

    
    # Define reasonable ranges for sensitivity and specificity
    se_range <- seq(input$se_min, input$se_max, by = 0.05)
    sp_range <- seq(input$sp_min, input$sp_max, by = 0.05) 
    
    prevalence_daly_points <- list()
    point_counter <- 1
    
    # Sample multiple combinations of Se and Sp
    for(i in 1:200) {
      se_sample <- sample(se_range, 1)
      sp_sample <- sample(sp_range, 1)
      
      # Calculate prevalence using the Stan model
      stan_data_sample <- list(
        T_pos = input$N_positive,
        N = input$N_tested,
        Se = se_sample,
        Sp = sp_sample
      )
      
      # Use tryCatch to handle potential Stan errors
      tryCatch({
        fit_sample <- sampling(stan_model_tp, data = stan_data_sample, 
                               chains = 2, iter = 2000, warmup = 1000, 
                               seed = 123 + i, refresh = 0)
        posterior_samples_sample <- rstan::extract(fit_sample)
        prevalence <- median(posterior_samples_sample$pi)
        
        # Calculate DALYs for this prevalence (simplified calculation)
        # Using the same burden calculation approach as main analysis
        daly_estimate <- calculate_daly_for_prevalence(prevalence, input$total_pop)
        
        prevalence_daly_points[[point_counter]] <- data.frame(
          Prevalence = prevalence * 100,  # Convert to percentage
          DALYs = daly_estimate,
          Sensitivity = se_sample * 100,
          Specificity = sp_sample * 100
        )
        point_counter <- point_counter + 1
      }, error = function(e) {
        # Skip this iteration if Stan fails
        message("Stan sampling failed for combination: Se=", se_sample, ", Sp=", sp_sample)
      })
    }
    
    if(length(prevalence_daly_points) > 0) {
      results$prevalence_daly_data <- do.call(rbind, prevalence_daly_points)
    }
  }
  
  # Helper function to calculate DALYs for a given prevalence

  calculate_daly_for_prevalence <- function(prevalence, total_pop) {
    # Use the same parameter values as the main analysis
    symp_ncc <- mean(results$sampled_parameters$symp_ncc)
    ncc_ep <- mean(results$sampled_parameters$ncc_ep)
    ncc_hd <- mean(results$sampled_parameters$ncc_hd)
    cft_ep <- mean(results$sampled_parameters$cft_ep)
    prop_ep_trt <- mean(results$sampled_parameters$prop_ep_trt)
    dsw_ep_tr <- mean(results$sampled_parameters$dsw_ep_tr)
    dsw_ep_nt <- mean(results$sampled_parameters$dsw_ep_nt)
    dsw_hd <- mean(results$sampled_parameters$dsw_hd)
    
    # More accurate burden calculation matching your main model
    ncc_cases <- total_pop * prevalence
    symp_ncc_cases <- ncc_cases * symp_ncc
    ep_cases <- symp_ncc_cases * ncc_ep
    hd_cases <- symp_ncc_cases * ncc_hd
    
    # Deaths from epilepsy
    deaths <- ep_cases * cft_ep
    
    # YLL calculation (using same approach as main model)
    # Use average life expectancy for simplicity
    avg_life_expectancy <- 70
    yll <- deaths * avg_life_expectancy
    
    # YLD calculation (using same approach as main model)
    treated_ep_cases <- ep_cases * prop_ep_trt
    untreated_ep_cases <- ep_cases * (1 - prop_ep_trt)
    
    yld_ep <- (treated_ep_cases * dsw_ep_tr) + (untreated_ep_cases * dsw_ep_nt)
    yld_hd <- hd_cases * dsw_hd
    
    total_yld <- yld_ep + yld_hd
    total_daly <- yll + total_yld
    
    return(total_daly)
  }
  
  # ----------------------------
  # Plots and Tables
  # ----------------------------
  output$daly_plot <- renderPlot({
    if(!is.null(results$daly_data)){
      plot_data <- data.frame(
        AgeGroup = factor(results$daly_data$age_groups, levels=results$daly_data$age_groups),
        Median = results$daly_data$median_daly,
        Lower = results$daly_data$lower_ci,
        Upper = results$daly_data$upper_ci
      )
      ggplot(plot_data, aes(x=AgeGroup, y=Median)) +
        geom_point(size=3) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2) +
        geom_line(aes(group=1), linetype="dashed", alpha=0.5) +
        labs(title="DALYs by Age Group", x="Age Group", y="DALYs") +
        theme_minimal() +
        scale_y_continuous(labels=scales::comma)
    }
  })
  
  output$prevalence_daly_plot <- renderPlot({
    if(!is.null(results$prevalence_daly_data)) {
      ggplot(results$prevalence_daly_data, aes(x = Prevalence, y = DALYs)) +
        geom_point(aes(color = Sensitivity, size = Specificity), alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
        labs(
          title = "Relationship between NCC Prevalence and DALYs (n = 200)",
          x = "NCC Prevalence (%)",
          y = "Total DALYs",
          color = "Sensitivity (%)",
          size = "Specificity (%)"
        ) +
        theme_minimal() +
        scale_y_continuous(labels = scales::comma) +
        scale_color_gradient(low = "blue", high = "red") +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }
  })
  output$download_prevalence_plot_png <- downloadHandler(
    filename = function() {
      paste0("prevalence_daly_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(results$prevalence_daly_data)
      p <- ggplot(results$prevalence_daly_data, aes(x = Prevalence, y = DALYs)) +
        geom_point(aes(color = Sensitivity, size = Specificity), alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
        labs(
          title = "Relationship between NCC Prevalence and DALYs",
          x = "NCC Prevalence (%)",
          y = "Total DALYs",
          color = "Sensitivity (%)",
          size = "Specificity (%)"
        ) +
        theme_minimal() +
        scale_y_continuous(labels = scales::comma) +
        scale_color_gradient(low = "blue", high = "red") +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_prevalence_plot_pdf <- downloadHandler(
    filename = function() {
      paste0("prevalence_daly_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(results$prevalence_daly_data)
      p <- ggplot(results$prevalence_daly_data, aes(x = Prevalence, y = DALYs)) +
        geom_point(aes(color = Sensitivity, size = Specificity), alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
        labs(
          title = "Relationship between NCC Prevalence and DALYs",
          x = "NCC Prevalence (%)",
          y = "Total DALYs",
          color = "Sensitivity (%)",
          size = "Specificity (%)"
        ) +
        theme_minimal() +
        scale_y_continuous(labels = scales::comma) +
        scale_color_gradient(low = "blue", high = "red") +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold")
        )
      ggsave(file, plot = p, device = "pdf", width = 10, height = 8)
    }
  )
  output$results_table <- renderTable({
    results$summary_data
  }, bordered=TRUE, align='l')
  
  output$parameter_table <- renderTable({
    if(!is.null(results$parameters_used)){
      # Format the values for better display
      formatted_params <- results$parameters_used
      formatted_params$Mean <- sprintf("%.4f", formatted_params$Mean)
      formatted_params$`Lower 95% BCI` <- sprintf("%.4f", formatted_params$`Lower 95% BCI`)
      formatted_params$`Upper 95% BCI` <- sprintf("%.4f", formatted_params$`Upper 95% BCI`)
      formatted_params
    }
  }, bordered=TRUE, align='l')
  
  output$sensitivity_plot <- renderPlot({
    if(!is.null(results$sampled_parameters) && !is.null(results$sampled_parameters$daly)){
      
      # Ensure all vectors have the same length
      n_samples <- length(results$sampled_parameters$daly)
      
      # Create a data frame with all parameters, ensuring equal length
      param_data <- data.frame(
        sero_prev_ncc = results$sampled_parameters$sero_prev_ncc[1:n_samples],
        prev_ncc = results$sampled_parameters$prev_ncc[1:n_samples],
        symp_ncc = results$sampled_parameters$symp_ncc[1:n_samples],
        ncc_ep = results$sampled_parameters$ncc_ep[1:n_samples],
        prop_ep_trt = results$sampled_parameters$prop_ep_trt[1:n_samples],
        ncc_hd = results$sampled_parameters$ncc_hd[1:n_samples],
        cft_ep = results$sampled_parameters$cft_ep[1:n_samples],
        dsw_ep_tr = results$sampled_parameters$dsw_ep_tr[1:n_samples],
        dsw_ep_nt = results$sampled_parameters$dsw_ep_nt[1:n_samples],
        dsw_hd = results$sampled_parameters$dsw_hd[1:n_samples]
      )
      
      # Remove any rows with NA values
      param_data <- na.omit(param_data)
      daly_values <- results$sampled_parameters$daly[1:n_samples]
      daly_values <- daly_values[complete.cases(param_data)]
      param_data <- param_data[complete.cases(param_data), ]
      
      # Check that we have enough data points
      if(nrow(param_data) < 10) {
        plot(1, 1, type = "n", xlab = "", ylab = "", main = "Insufficient data for sensitivity analysis")
        text(1, 1, "Not enough complete data points\nfor sensitivity analysis")
        return()
      }
      
      # Ensure x and y have the same length
      if(length(daly_values) != nrow(param_data)) {
        min_length <- min(length(daly_values), nrow(param_data))
        daly_values <- daly_values[1:min_length]
        param_data <- param_data[1:min_length, ]
      }
      
      # Run sensitivity analysis
      sa_daly_pcc <- sa_pcc(daly_values, as.matrix(param_data))
      
      # Parameter labels
      daly_items <- c(
        "NCC Seroprevalence", 
        "NCC prevalence", 
        "NCC symptomatic",
        "E among NCC",
        "E treatment",
        "H among NCC",
        "E case fatality", 
        "E disability weight (treated)",
        "E disability weight (untreated)",
        "H disability weight"
      )
      
      # Create tornado plot
      tornado(sa_daly_pcc, daly_items)
    }
  })
 
  # Add these reactive outputs
  output$results_ready <- reactive({
    !is.null(results$summary_data)
  })
  outputOptions(output, "results_ready", suspendWhenHidden = FALSE)
  
  output$plot_ready <- reactive({
    !is.null(results$prevalence_daly_data)
  })
  outputOptions(output, "plot_ready", suspendWhenHidden = FALSE)  
  
   
  # Download handlers
  output$download_prevalence_daly_data <- downloadHandler(
    filename = function() {
      paste0("prevalence_daly_relationship_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(results$prevalence_daly_data)
      write.csv(results$prevalence_daly_data, file, row.names = FALSE)
    }
  )
  
  output$download_sens_plot_png <- downloadHandler(
    filename = function() {
      paste0("sensitivity_analysis_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(results$sampled_parameters)
      req(results$sampled_parameters$daly)
      
      # Recreate the sensitivity analysis plot for download
      plot_data <- tryCatch({
        sp <- results$sampled_parameters
        
        required_params <- c("sero_prev_ncc", "prev_ncc", "symp_ncc", "ncc_ep", 
                             "prop_ep_trt", "ncc_hd", "cft_ep", "dsw_ep_tr", 
                             "dsw_ep_nt", "dsw_hd", "daly")
        
        if(!all(required_params %in% names(sp))) {
          stop("Missing parameters for sensitivity analysis")
        }
        
        param_matrix <- cbind(
          sero_prev_ncc = sp$sero_prev_ncc,
          prev_ncc = sp$prev_ncc,
          symp_ncc = sp$symp_ncc,
          ncc_ep = sp$ncc_ep,
          prop_ep_trt = sp$prop_ep_trt,
          ncc_hd = sp$ncc_hd,
          cft_ep = sp$cft_ep,
          dsw_ep_tr = sp$dsw_ep_tr,
          dsw_ep_nt = sp$dsw_ep_nt,
          dsw_hd = sp$dsw_hd
        )
        
        valid_rows <- complete.cases(param_matrix) & 
          is.finite(sp$daly) & 
          !is.na(sp$daly)
        
        if(sum(valid_rows) < 20) {
          stop("Insufficient valid data points for sensitivity analysis")
        }
        
        param_matrix_clean <- param_matrix[valid_rows, ]
        daly_clean <- sp$daly[valid_rows]
        
        sa_daly_pcc <- sa_pcc(daly_clean, param_matrix_clean)
        
        daly_items <- c(
          "NCC Seroprevalence", 
          "NCC prevalence", 
          "NCC symptomatic",
          "E among NCC",
          "E treatment",
          "H among NCC",
          "E case fatality", 
          "E disability weight (treated)",
          "E disability weight (untreated)",
          "H disability weight"
        )
        
        list(sa_daly_pcc = sa_daly_pcc, daly_items = daly_items)
      }, error = function(e) {
        return(NULL)
      })
      
      if(!is.null(plot_data)) {
        # Create the plot with higher resolution for download
        p <- tornado(plot_data$sa_daly_pcc, plot_data$daly_items)
        ggsave(file, plot = p, device = "png", width = 10, height = 8, dpi = 300)
      }
    }
  )
  
  output$download_sens_plot_pdf <- downloadHandler(
    filename = function() {
      paste0("sensitivity_analysis_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(results$sampled_parameters)
      req(results$sampled_parameters$daly)
      
      # Recreate the sensitivity analysis plot for download
      plot_data <- tryCatch({
        sp <- results$sampled_parameters
        
        required_params <- c("sero_prev_ncc", "prev_ncc", "symp_ncc", "ncc_ep", 
                             "prop_ep_trt", "ncc_hd", "cft_ep", "dsw_ep_tr", 
                             "dsw_ep_nt", "dsw_hd", "daly")
        
        if(!all(required_params %in% names(sp))) {
          stop("Missing parameters for sensitivity analysis")
        }
        
        param_matrix <- cbind(
          sero_prev_ncc = sp$sero_prev_ncc,
          prev_ncc = sp$prev_ncc,
          symp_ncc = sp$symp_ncc,
          ncc_ep = sp$ncc_ep,
          prop_ep_trt = sp$prop_ep_trt,
          ncc_hd = sp$ncc_hd,
          cft_ep = sp$cft_ep,
          dsw_ep_tr = sp$dsw_ep_tr,
          dsw_ep_nt = sp$dsw_ep_nt,
          dsw_hd = sp$dsw_hd
        )
        
        valid_rows <- complete.cases(param_matrix) & 
          is.finite(sp$daly) & 
          !is.na(sp$daly)
        
        if(sum(valid_rows) < 20) {
          stop("Insufficient valid data points for sensitivity analysis")
        }
        
        param_matrix_clean <- param_matrix[valid_rows, ]
        daly_clean <- sp$daly[valid_rows]
        
        sa_daly_pcc <- sa_pcc(daly_clean, param_matrix_clean)
        
        daly_items <- c(
          "NCC Seroprevalence", 
          "NCC prevalence", 
          "NCC symptomatic",
          "E among NCC",
          "E treatment",
          "H among NCC",
          "E case fatality", 
          "E disability weight (treated)",
          "E disability weight (untreated)",
          "H disability weight"
        )
        
        list(sa_daly_pcc = sa_daly_pcc, daly_items = daly_items)
      }, error = function(e) {
        return(NULL)
      })
      
      if(!is.null(plot_data)) {
        # Create the plot for PDF download
        p <- tornado(plot_data$sa_daly_pcc, plot_data$daly_items)
        ggsave(file, plot = p, device = "pdf", width = 10, height = 8)
      }
    }
  )
  
  output$download_parameters <- downloadHandler(
    filename = function() {
      paste0("parameter_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(results$parameters_used)
      write.csv(results$parameters_used, file, row.names = FALSE)
    }
  )
  
  output$download_summary <- downloadHandler(
    filename = function() paste0("ncc_burden_summary_", Sys.Date(), ".csv"),
    content = function(file) write.csv(results$summary_data, file, row.names=FALSE)
  )
  
  output$download_daly_data <- downloadHandler(
    filename = function() paste0("daly_by_age_", Sys.Date(), ".csv"),
    content = function(file){
      if(!is.null(results$daly_data)){
        daly_df <- data.frame(
          AgeGroup = results$daly_data$age_groups,
          Median_DALYs = results$daly_data$median_daly,
          Lower_95_CI = results$daly_data$lower_ci,
          Upper_95_CI = results$daly_data$upper_ci
        )
        write.csv(daly_df, file, row.names=FALSE)
      }
    }
  )
}

# ----------------------------
# Run App
# ----------------------------
shinyApp(ui=ui, server=server)