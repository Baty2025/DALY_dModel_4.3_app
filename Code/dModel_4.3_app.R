library(shiny)
library(ggplot2)
library(dplyr)
library(rstan)
library(tibble)

# ----------------------------
# Helper Functions
# ----------------------------
rle <- function(x, le = "gbd") {
  # Example life expectancy function
  # Replace with actual function if different
  if(le == "gbd"){
    return(c(86.02, 85.21, 81.25, 74.25, 64.75))
  } else {
    return(x)
  }
}

burden <- function(N, DW, A, L, K, r, a) {
  # Simplified DALY calculation (replace with your full function)
  N * DW * L
}

mean_ci <- function(x) c(mean(x), quantile(x, 0.025), quantile(x, 0.975))

# ----------------------------
# Precompile Stan Model
# ----------------------------
stan_model_tp <- stan_model(model_code = "
data {
  int<lower=0> T_pos;
  int<lower=0> N;
  real<lower=0, upper=1> Se;
  real<lower=0, upper=1> Sp;
}
parameters {
  real<lower=0, upper=1> pi;
}
model {
  T_pos ~ binomial(N, pi * Se + (1 - pi) * (1 - Sp));
  pi ~ beta(1, 1);
}
")

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  titlePanel("NCC Burden of Disease Calculator"),
  sidebarLayout(
    sidebarPanel(
      h4("Serological Test Data"),
      numericInput("N_tested", "Number of People Tested", value = 1063, min = 0),
      numericInput("N_positive", "Number of Positive Tests", value = 169, min = 0),
      sliderInput("sensitivity", "Test Sensitivity", value = 69.3, min = 0, max = 100),
      sliderInput("specificity", "Test Specificity", value = 97, min = 0, max = 100),
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
      numericInput("female_age5", "65+", value = 0.081, min = 0, max = 1, step = 0.01)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("DALYs by Age Group",
                 plotOutput("daly_plot"),
                 downloadButton("download_daly_data", "Download DALY Data")),
        tabPanel("Results Summary",
                 tableOutput("results_table"),
                 downloadButton("download_summary", "Download Summary CSV")),
        tabPanel("Sensitivity Analysis",
                 plotOutput("sensitivity_plot"))
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
    sensitivity_data = NULL
  )
  
  observeEvent(input$calculate, {
    showNotification("Calculating burden... This may take a few minutes.", type = "message")
    
    withProgress(message = 'Running analysis', value = 0, {
      
      incProgress(0.1, detail = "Processing population data")
      total_male <- input$total_pop * input$prop_m
      total_female <- input$total_pop - total_male
      male_prop_per_age <- c(input$male_age1, input$male_age2, input$male_age3, input$male_age4, input$male_age5)
      female_prop_per_age <- c(input$female_age1, input$female_age2, input$female_age3, input$female_age4, input$female_age5)
      pop_mx <- matrix(0, nrow = 5, ncol = 2)
      pop_mx[,1] <- total_male * male_prop_per_age 
      pop_mx[,2] <- total_female * female_prop_per_age 
      pop_mx <- matrix(pop_mx, nrow = 10, ncol = 1)
      
      incProgress(0.3, detail = "Sampling seroprevalence with Stan")
      stan_data <- list(
        T_pos = input$N_positive,
        N = input$N_tested,
        Se = input$sensitivity / 100,
        Sp = input$specificity / 100
      )
      fit <- sampling(stan_model_tp, data = stan_data, chains = 1, iter = 200, warmup = 100, seed = 123, refresh = 0)
      posterior_samples <- rstan::extract(fit)
      bci <- quantile(posterior_samples$pi, probs = c(0.025, 0.975))
      
      incProgress(0.6, detail = "Calculating burden of disease")
      set.seed(123)
      n <- 500  # reduced for Shiny speed
      prev_ncc <- runif(n, bci[1], bci[2])
      symp_ncc <- rbeta(n, 11, 21)
      prev_ncc_symp <- prev_ncc * symp_ncc
      ncc_ep <- rbeta(n, 122, 45)
      prev_ncc_ep_act <- prev_ncc * symp_ncc * ncc_ep
      ncc_hd <- rbeta(n, 53, 155)
      prev_ncc_hd <- prev_ncc * symp_ncc * ncc_hd
      N_prev_ncc_ep_act <- t(apply(pop_mx, 1, function(x) x * prev_ncc_ep_act))
      N_prev_ncc_hd <- t(apply(pop_mx, 1, function(x) x * prev_ncc_hd))
      
      cft_ep <- runif(n, 0.033, 0.316)
      
      N_prev_ncc_ep_mrt <- t(apply(pop_mx, 1, function(x) x * prev_ncc_ep_act*cft_ep))
      
      prop_ep_trt <- runif(n, 0.25, 0.90)
      
      # YLL calculation
      ep_aad <- c(2, 9.5, 29.5, 52, 65)
      ep_rle <- rle(ep_aad, "gbd")
      K <- 0; r <- 0
      yll_ep <- array(dim = c(5,2,n))
      for(i in 1:5){
        yll_ep[i,1,] <- burden(N = N_prev_ncc_ep_act[i,]*mean(cft_ep), DW = 1, A = ep_aad[i], L = ep_rle[i], K = K, r = r, a = ep_aad[i])
        yll_ep[i,2,] <- burden(N = N_prev_ncc_ep_act[i+5,]*mean(cft_ep), DW = 1, A = ep_aad[i], L = ep_rle[i], K = K, r = r, a = ep_aad[i])
      }
      
      # YLD calculation
      dsw_ep_tr <- runif(n, 0.031, 0.072)
      dsw_ep_nt <- runif(n, 0.173, 0.710)
      dsw_hd <- runif(n, 0.022, 0.588)
      A <- 0; L <- 1; K <- 0; r <- 0; a <- A
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
      results$daly_data <- list(
        age_groups = age_groups,
        daly_matrix = daly_by_age,
        median_daly = apply(daly_by_age,1,median),
        lower_ci = apply(daly_by_age,1,quantile,0.025),
        upper_ci = apply(daly_by_age,1,quantile,0.975)
      )
      
      # ----------------------------
      # Summary Results Table
      # ----------------------------
      # Estimated seroprevalence
      seroprev_mean <- mean(posterior_samples$pi)
      seroprev_ci <- quantile(posterior_samples$pi, probs = c(0.025, 0.975))
      
      
      yld_ep_all <- apply(yld_ep, 3, sum)
      yld_hd_all <- apply(yld_hd, 3, sum)
      yld <- yld_ep_all + yld_hd_all
      yll_ep_all <- apply(yll_ep, 3, sum)
      daly <- yld + yll_ep_all
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
          mean_ci(colSums(N_prev_ncc_ep_mrt))[1],  # placeholder for deaths
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
      
    })
    
    showNotification("Calculation complete!", type = "message")
  })
  
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
  
  output$results_table <- renderTable({
    results$summary_data
  }, bordered=TRUE, align='l')
  
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
  
  output$sensitivity_plot <- renderPlot({
    if(!is.null(results$daly_data)){
      sens_data <- data.frame(
        Parameter = c("Sensitivity","Specificity","Population"),
        Value = c(input$sensitivity,input$specificity,log10(input$total_pop)),
        Impact = c(0.8,0.6,0.9)  # placeholder
      )
      ggplot(sens_data, aes(x=Parameter, y=Impact)) +
        geom_col(fill="steelblue") +
        labs(title="Parameter Impact on DALY Estimates", y="Relative Impact") +
        theme_minimal()
    }
  })
}

# ----------------------------
# Run App
# ----------------------------
shinyApp(ui=ui, server=server)
