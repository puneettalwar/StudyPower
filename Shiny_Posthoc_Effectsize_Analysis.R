# app.R
library(shiny)
library(ggplot2)
library(grid)
library(pwr)

# ---------------- Power and CI Functions ----------------

compute_power <- function(f2, u, v, F_crit) {
  lambda <- f2 * (v + 1)
  1 - pf(F_crit, df1 = u, df2 = v, ncp = lambda)
}

find_f2 <- function(target_power, u, v, F_crit) {
  f2_vals <- seq(0.001, 5, by = 0.0001)
  for (f2 in f2_vals) {
    pwr <- compute_power(f2, u, v, F_crit)
    if (pwr >= target_power) {
      return(list(f2 = f2, power = pwr))
    }
  }
  return(NULL)
}

ci_r2_noncentral <- function(f2, N, k, conf.level = 0.99) {
  df1 <- k
  df2 <- N - k - 1
  alpha <- 1 - conf.level
  lambda <- f2 * (df2 + 1)
  
  if (df2 <= 0) stop("Degrees of freedom (df2) must be > 0")
  
  F_lower <- tryCatch(qf(alpha / 2, df1, df2, ncp = lambda), error = function(e) NA)
  F_upper <- tryCatch(qf(1 - alpha / 2, df1, df2, ncp = lambda), error = function(e) NA)
  
  if (is.na(F_lower) || is.na(F_upper)) {
    return(c(lower = NA, upper = NA))
  }
  
  R2_lower <- (df1 * F_lower) / (df1 * F_lower + df2)
  R2_upper <- (df1 * F_upper) / (df1 * F_upper + df2)
  
  R2_lower <- max(0, min(1, R2_lower))
  R2_upper <- max(0, min(1, R2_upper))
  
  return(c(lower = R2_lower, upper = R2_upper))
}

ci_f2_noncentral <- function(R2_ci) {
  f2_lower <- R2_ci[1] / (1 - R2_ci[1])
  f2_upper <- R2_ci[2] / (1 - R2_ci[2])
  return(c(lower = f2_lower, upper = f2_upper))
}

ci_r_from_r2 <- function(R2_ci) {
  c(lower = sqrt(max(0, R2_ci[1])), upper = sqrt(min(1, R2_ci[2])))
}

generate_nomogram_plot <- function(alpha_input) {
  f2_vals <- seq(0.01, 5, by = 0.001)
  u <- 1
  N_vals <- c(30, 50, 80, 100, 150, 200, 300, 500)
  
  df <- expand.grid(f2 = f2_vals, N = N_vals)
  df$v <- df$N - u - 1
  df$F_crit <- qf(1 - alpha_input, df1 = u, df2 = df$v)
  df$lambda <- df$f2 * (df$v + 1)
  
  df$power <- mapply(function(f2, v, F_crit) {
    lambda <- f2 * (v + 1)
    1 - pf(F_crit, df1 = u, df2 = v, ncp = lambda)
  }, df$f2, df$v, df$F_crit)
  
  df <- df[df$power >= 0 & df$power <= 1, ]
  
  ggplot(df, aes(x = f2, y = power, color = factor(N), group = N)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    scale_color_viridis_d(name = "Sample Size (N)", option = "B") +
    labs(
      title = "Power vs. Effect Size (f²)",
      subtitle = paste("α =", alpha_input, "| 1 predictor tested"),
      x = "Effect Size (f²)",
      y = "Power"
    ) +
    theme_minimal(base_size = 14)
}

# ------------------------- UI --------------------------

ui <- fluidPage(
  titlePanel("Specify Input Parameters"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("N", "Sample size (N):", value = 500, min = 10),
      numericInput("k", "Total predictors (k):", value = 6, min = 1),
      numericInput("u", "Tested predictors (u):", value = 6, min = 1),
      numericInput("alpha", "Significance level (alpha):", value = 0.01, min = 0.0001, max = 0.2, step = 0.001),
      numericInput("conf_level", "Confidence level for CI:", value = 0.99, min = 0.80, max = 0.999, step = 0.01),
      
      checkboxInput("use_observed", "Use observed effect size (post-hoc)", value = FALSE),
      
      conditionalPanel(
        condition = "input.use_observed == false",
        numericInput("power_target", "Target power:", value = 0.80, min = 0.01, max = 0.99, step = 0.01)
      ),
      conditionalPanel(
        condition = "input.use_observed == true",
        radioButtons("obs_type", "Observed effect size type:",
                     choices = c("R²" = "r2", "f²" = "f2"), inline = TRUE),
        conditionalPanel(
          condition = "input.obs_type == 'r2'",
          numericInput("observed_r2", "Observed R²:", value = 0.05, min = 0, max = 1, step = 0.001)
        ),
        conditionalPanel(
          condition = "input.obs_type == 'f2'",
          numericInput("observed_f2", "Observed f²:", value = 0.05, min = 0, max = 10, step = 0.001)
        )
      ),
      
      actionButton("compute", "Run Analysis")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions",
                 fluidRow(
                   column(12,
                          div(style = "font-size:16px; line-height:1.6;",
                              h2("Linear Multiple Regression: Power & Sensitivity Analysis"),
                              p("This app replicates and extends G*Power for post-hoc power and sensitivity analysis in linear multiple regression."),
                              
                              h3("Inputs"),
                              tags$ul(
                                tags$li(strong("Sample size (N): "), "Number of participants (must be > k + 1)."),
                                tags$li(strong("Total predictors (k): "), "Total predictors in the model."),
                                tags$li(strong("Tested predictors (u): "), "Predictors of interest being tested (k for the full model))."),
                                tags$li(strong("Significance level (alpha): "), "Type I error probability."),
                                tags$li(strong("Confidence level: "), "Confidence interval for effect sizes."),
                                tags$li(strong("Analysis type: "), "Choose between Target Power Analysis or Observed R² Analysis."),
                                tags$li(strong("Run Analysis: "), "Press the button to compute results.")
                              ),
                              
                              h3("Outputs"),
                              tags$ul(
                                tags$li(strong("Numerical Results: "), "Displays effect sizes (f, f², R², r, d) with confidence intervals."),
                                #tags$li(strong("Numerical Results: "), "Displays effect sizes (f, fÂ², RÂ², r, d) with confidence intervals."),
                                tags$li(strong("Power Curve: "), "Shows power vs effect size (f²)."),
                                tags$li(strong("Nomogram: "), "Shows power curves across different sample sizes.")
                              ),
                              
                              h3("How to Use"),
                              tags$ol(
                                tags$li("For sensitivity analysis: Enter N, k, u, alpha, and target power."),
                                tags$li("For post-hoc analysis: Check 'Use observed R² or f², enter R² or f²."),
                                tags$li("Click 'Run Analysis'."),
                                tags$li("Interpret results and plots for your study context.")
                              ),
                              
                              h3("Notes"),
                              tags$ul(
                                tags$li("Ensure N > k + 1."),
                                tags$li("Nomogram assumes 1 tested predictor (u = 1)."),
                                tags$li("Confidence intervals are based on non-central F distributions.")
                              ),
                              
                              hr(),
                              p("Use this app for replication of G*Power analyses in R, with added confidence interval reporting."),
                          )                  
                   )
                 )
        ),
        tabPanel("Results", verbatimTextOutput("results")),
        tabPanel("Power Curve", plotOutput("powerPlot")),
        tabPanel("Nomogram", plotOutput("nomogram"))
        
        
      ),
      #verbatimTextOutput("results")
    )
    
  )
)

# ----------------------- Server -------------------------

server <- function(input, output) {
  observeEvent(input$compute, {
    req(input$N > input$k + 1)
    v <- input$N - input$k - 1
    F_crit <- qf(1 - input$alpha, df1 = input$u, df2 = v)
    
    posthoc_power <- NA_real_
    
    if (!input$use_observed) {
      # --- MDE path ---
      result <- find_f2(input$power_target, input$u, v, F_crit)
      if (is.null(result)) {
        output$results <- renderText("Power target not achievable with given parameters.")
        return()
      }
      f2 <- result$f2
      R2 <- f2 / (1 + f2)
      f  <- sqrt(f2)
      r  <- sqrt(R2)
      d  <- sqrt(f2) * 2
      
      ci_r2 <- ci_r2_noncentral(f2, input$N, input$k, conf.level = input$conf_level)
      
    } else {
      # --- Observed effect size path ---
      if (input$obs_type == "r2") {
        R2 <- input$observed_r2
        f2 <- R2 / (1 - R2)
      } else {
        f2 <- input$observed_f2
        R2 <- f2 / (1 + f2)
      }
      f <- sqrt(f2)
      r <- sqrt(R2)
      d <- sqrt(f2) * 2
      
      ci_r2 <- ci_r2_noncentral(f2, input$N, input$k, conf.level = input$conf_level)
      posthoc_power <- compute_power(f2, input$u, v, F_crit)
    }
    
    ci_r <- ci_r_from_r2(ci_r2)
    ci_f2 <- ci_f2_noncentral(ci_r2)
    
    output$results <- renderText({
      header <- ifelse(input$use_observed, "Observed effect size analysis\n", "Minimum detectable effect size:\n")
      power_line <- if (input$use_observed) {
        paste0("Post-hoc power (α = ", input$alpha, "): ", round(posthoc_power, 3), "\n")
      } else {
        paste0("Target power: ", round(input$power_target, 3), "\n")
      }
      
      paste0(
        header,
        power_line,
        "f = ", round(f, 3), "\n",
        "f² = ", round(f2, 3), "\n",
        "R² = ", round(R2, 3), "\n",
        "r = ", round(r, 3), "\n",
        "d = ", round(d, 3), "\n\n",
        input$conf_level * 100, "% CI for f²: [", round(ci_f2[1], 3), ", ", round(ci_f2[2], 3), "]\n",
        input$conf_level * 100, "% CI for R²: [", round(ci_r2[1], 3), ", ", round(ci_r2[2], 3), "]\n",
        input$conf_level * 100, "% CI for r:  [", round(ci_r[1], 3), ", ", round(ci_r[2], 3), "]\n\n\n",
        
        "Effect size thresholds:\n\n",
        "f² = 0.02 (small effect), 0.15 (medium effect), and 0.35 (large effect)\n",
        "r² = 0.01 (small effect), 0.09 (medium effect), and 0.25 (large effect)\n",
        "f = 0.10 (small effect), 0.25 (medium effect), and 0.40 (large effect)\n",
        "r = 0.10 (small effect), 0.30 (medium effect), and 0.50 (large effect)\n",
        "d = 0.2 (small effect), 0.5 (medium effect), and 0.8 (large effect)\n"
      )
    })
    
    output$nomogram <- renderPlot({
      generate_nomogram_plot(input$alpha)
    })
    
    output$powerPlot <- renderPlot({
      f2_vals <- seq(0.001, 5, by = 0.001)
      power_vals <- sapply(f2_vals, function(f2_val) compute_power(f2_val, input$u, v, F_crit))
      
      plot(f2_vals, power_vals, type = "l", lwd = 2,col = "blue",
           xlab = "Cohen's f²", ylab = "Power", main = "Power vs f²",
           ylim = c(0, 1))
      
      if (!input$use_observed) {
        abline(h = input$power_target, lty = 2)
        abline(v = f2, col = "red", lty = 2)
        legend("bottomright",
               legend = c("Power Curve", "Target Power", "Min Detectable f²"),
               lty = c(1, 2, 2), col = c("blue", "black", "red"), bty = "n")
      } else {
        abline(h = posthoc_power, lty = 2)
        abline(v = f2, col = "red", lty = 2)
        legend("bottomright",
               legend = c("Power Curve", "Post-hoc Power", "Observed f²"),
               lty = c(1, 2, 2), col = c("blue", "black", "red"), bty = "n")
      }
    })
  })
}

shinyApp(ui, server)
