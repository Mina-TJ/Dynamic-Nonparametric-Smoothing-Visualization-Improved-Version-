library(shiny)
library(ggplot2)
library(MASS)
library(shinythemes)

# UI definition
ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  tags$head(
    tags$style(HTML("
      .well { background-color: #ffffff; }
      .nav-tabs { margin-bottom: 20px; }
      .control-label { font-weight: 500; }
      .shiny-output-error { color: #ff0000; }
      .btn { margin: 5px; }
      .progress-container { margin-top: 15px; }
      .metric-box { 
        background: #f8f9fa; 
        padding: 10px; 
        border-radius: 5px; 
        margin: 5px 0;
      }
    "))
  ),
  titlePanel("Dynamic Nonparametric Smoothing Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tabsetPanel(
        tabPanel("Data",
                 br(),
                 selectInput("data_type", "Select Data Type",
                             choices = c("Sine Wave" = "sine",
                                         "Linear Trend" = "linear",
                                         "Quadratic" = "quadratic",
                                         "Volatile Series" = "volatile",
                                         "Step Function" = "step",
                                         "Custom Upload" = "custom"),
                             selected = "sine"),
                 
                 conditionalPanel(
                   condition = "input.data_type == 'custom'",
                   fileInput("data_file", "Upload CSV",
                             accept = c(".csv"))
                 ),
                 
                 numericInput("n_points", "Number of Points",
                              value = 300, min = 50, max = 1000),
                 
                 sliderInput("noise_level", "Noise Level",
                             min = 0, max = 2, value = 0.25, step = 0.05)
        ),
        
        tabPanel("Smoothing",
                 br(),
                 selectInput("smooth_method", "Smoothing Method",
                             choices = c("Mean" = "mean",
                                         "Median" = "median",
                                         "Trimmed Mean" = "trimmed",
                                         "Local Regression" = "loess",
                                         "Kernel" = "kernel"),
                             selected = "kernel"),
                 
                 sliderInput("window_size", "Window Size",
                             min = 3, max = 101, value = 21, step = 2),
                 
                 conditionalPanel(
                   condition = "input.smooth_method == 'trimmed'",
                   sliderInput("trim_percent", "Trim Percentage",
                               min = 0, max = 0.4, value = 0.1, step = 0.05)
                 ),
                 
                 conditionalPanel(
                   condition = "input.smooth_method == 'loess'",
                   sliderInput("degree", "Polynomial Degree",
                               min = 0, max = 2, value = 1, step = 1)
                 ),
                 
                 conditionalPanel(
                   condition = "input.smooth_method == 'kernel'",
                   selectInput("kernel_type", "Kernel Type",
                               choices = c("Gaussian" = "gaussian",
                                           "Epanechnikov" = "epanechnikov",
                                           "Triangular" = "triangular",
                                           "Uniform" = "uniform"),
                               selected = "gaussian"),
                   sliderInput("bandwidth", "Bandwidth Multiplier",
                               min = 0.2, max = 3, value = 1, step = 0.1)
                 ),
                 
                 checkboxInput("show_confidence", "Show Confidence Band", FALSE),
                 checkboxInput("compare_methods", "Compare Methods", FALSE)
        ),
        
        tabPanel("Animation",
                 br(),
                 sliderInput("animation_speed", "Speed (ms)",
                             min = 10, max = 500, value = 50, step = 10),
                 
                 numericInput("jump_to", "Jump to Position",
                              value = 1, min = 1),
                 
                 actionButton("jump", "Jump", class = "btn-info btn-sm"),
                 
                 hr(),
                 
                 actionButton("start", "Start/Stop", 
                              class = "btn-success btn-block"),
                 actionButton("reset", "Reset", 
                              class = "btn-warning btn-block"),
                 actionButton("complete", "Complete", 
                              class = "btn-primary btn-block"),
                 
                 hr(),
                 
                 div(class = "progress-container",
                     htmlOutput("progress_bar")
                 )
        )
      ),
      
      hr(),
      
      verbatimTextOutput("stats"),
      
      conditionalPanel(
        condition = "input.compare_methods == true",
        div(class = "metric-box",
            htmlOutput("comparison_metrics")
        )
      ),
      
      tags$div(
        style = "margin-top: 20px; padding-top: 15px; border-top: 1px solid #e5e5e5;",
        tags$p(
          style = "color: #666; font-size: 0.9em; margin-bottom: 5px;",
          "Dr. Aydede - SMU",
          tags$br(),
          "Prepared for MBAN Students"
        ),
        tags$p(
          style = "color: #666; font-size: 0.8em; font-style: italic;",
          paste("©️ Copyright", format(Sys.Date(), "%Y"), "- All Rights Reserved")
        )
      )
    ),
    
    mainPanel(
      width = 9,
      plotOutput("smoothing_plot", height = "550px"),
      
      conditionalPanel(
        condition = "input.compare_methods == true",
        plotOutput("comparison_plot", height = "300px")
      ),
      
      wellPanel(
        uiOutput("method_explanation")
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Define kernel functions
  kernel_functions <- list(
    gaussian = function(x) dnorm(x, sd = 1),
    epanechnikov = function(x) ifelse(abs(x) <= 1, 3/4 * (1 - x^2), 0),
    triangular = function(x) ifelse(abs(x) <= 1, 1 - abs(x), 0),
    uniform = function(x) ifelse(abs(x) <= 1, 0.5, 0)
  )
  
  # Initialize reactive values
  rv <- reactiveValues(
    current_pos = 1,
    is_running = FALSE,
    true_signal = NULL
  )
  
  # Generate data based on inputs
  data <- reactive({
    if (input$data_type == "custom" && !is.null(input$data_file)) {
      df <- read.csv(input$data_file$datapath)
      if (ncol(df) >= 2) {
        df <- df[, 1:2]
        names(df) <- c("time", "y")
        rv$true_signal <- NULL
        return(df)
      }
    }
    
    n <- input$n_points
    set.seed(42)
    x <- sort(runif(n) * 2 * pi)
    
    y_true <- switch(input$data_type,
                     "sine" = sin(x),
                     "linear" = x/pi - 1,
                     "quadratic" = (x/pi - 1)^2 - 1,
                     "step" = ifelse(x < pi, -0.5, 
                                     ifelse(x < 1.5*pi, 0.5, -0.5)),
                     "volatile" = {
                       base <- cumsum(rnorm(n, 0, 0.3))
                       spikes <- numeric(n)
                       spike_positions <- seq(20, n, by = round(n/10))
                       spikes[spike_positions] <- rnorm(length(spike_positions), 0, 2)
                       shifts <- cumsum(sample(c(0, 1, -1), n, 
                                               prob = c(0.98, 0.01, 0.01), 
                                               replace = TRUE))
                       base + spikes + shifts
                     })
    
    rv$true_signal <- data.frame(time = 1:n, y = y_true)
    y <- y_true + rnorm(n) * input$noise_level
    data.frame(time = 1:n, y = y)
  })
  
  # Store smoothed values
  smoothed_data <- reactiveVal(data.frame(time = numeric(0), y = numeric(0)))
  comparison_data <- reactiveVal(list())
  
  # Calculate smoothed value with confidence interval
  calculate_smooth <- function(window_data, return_variance = FALSE) {
    values <- window_data$y
    times <- window_data$time
    
    result <- switch(input$smooth_method,
                     "mean" = list(fitted = mean(values), 
                                   variance = var(values) / length(values)),
                     "median" = list(fitted = median(values),
                                     variance = 1.57 * var(values) / length(values)),
                     "trimmed" = {
                       trimmed_mean <- mean(values, trim = input$trim_percent)
                       list(fitted = trimmed_mean,
                            variance = var(values) / length(values))
                     },
                     "kernel" = {
                       kernel_fn <- kernel_functions[[input$kernel_type]]
                       center_time <- mean(times)
                       h <- input$bandwidth * sd(times)
                       x_scaled <- (times - center_time) / h
                       weights <- kernel_fn(x_scaled)
                       weights <- weights / sum(weights)
                       fitted_val <- sum(values * weights)
                       variance_val <- sum(weights^2 * (values - fitted_val)^2)
                       list(fitted = fitted_val, variance = variance_val)
                     },
                     "loess" = {
                       if (length(values) >= 4) {
                         if (input$degree == 0) {
                           kernel_fn <- kernel_functions[["gaussian"]]
                           center_time <- mean(times)
                           h <- input$window_size/4
                           x_scaled <- (times - center_time) / h
                           weights <- kernel_fn(x_scaled)
                           weights <- weights / sum(weights)
                           fitted_val <- sum(values * weights)
                           list(fitted = fitted_val,
                                variance = var(values) / length(values))
                         } else {
                           fit <- if (input$degree == 1) {
                             tryCatch(MASS::rlm(values ~ times), 
                                      error = function(e) lm(values ~ times))
                           } else {
                             lm(values ~ poly(times, 2))
                           }
                           fitted_val <- predict(fit, newdata = data.frame(times = mean(times)))
                           residual_var <- var(residuals(fit))
                           list(fitted = fitted_val,
                                variance = residual_var / length(values))
                         }
                       } else {
                         list(fitted = mean(values),
                              variance = var(values) / length(values))
                       }
                     })
    
    if (return_variance) {
      return(result)
    } else {
      return(result$fitted)
    }
  }
  
  # Calculate all comparison methods
  calculate_all_methods <- function(window_data) {
    methods <- c("mean", "median", "kernel")
    results <- list()
    
    for (method in methods) {
      temp_input <- input
      temp_input$smooth_method <- method
      
      values <- window_data$y
      times <- window_data$time
      
      results[[method]] <- switch(method,
                                  "mean" = mean(values),
                                  "median" = median(values),
                                  "kernel" = {
                                    kernel_fn <- kernel_functions[["gaussian"]]
                                    center_time <- mean(times)
                                    h <- sd(times)
                                    x_scaled <- (times - center_time) / h
                                    weights <- kernel_fn(x_scaled)
                                    weights <- weights / sum(weights)
                                    sum(values * weights)
                                  }
      )
    }
    results
  }
  
  # Handle animation timing
  observe({
    if (rv$is_running) {
      invalidateLater(input$animation_speed)
      
      isolate({
        current_data <- data()
        window <- input$window_size
        max_pos <- nrow(current_data) - window + 1
        
        if (rv$current_pos <= max_pos) {
          window_data <- current_data[rv$current_pos:(rv$current_pos + window - 1), ]
          
          smooth_result <- calculate_smooth(window_data, return_variance = TRUE)
          
          new_point <- data.frame(
            time = rv$current_pos + (window/2 - 1),
            y = smooth_result$fitted,
            se = sqrt(smooth_result$variance)
          )
          
          current_smoothed <- smoothed_data()
          smoothed_data(rbind(current_smoothed, new_point))
          
          # Calculate comparison if enabled
          if (input$compare_methods) {
            comp_results <- calculate_all_methods(window_data)
            current_comp <- comparison_data()
            for (method in names(comp_results)) {
              if (is.null(current_comp[[method]])) {
                current_comp[[method]] <- data.frame(time = numeric(0), y = numeric(0))
              }
              current_comp[[method]] <- rbind(
                current_comp[[method]],
                data.frame(time = new_point$time, y = comp_results[[method]])
              )
            }
            comparison_data(current_comp)
          }
          
          rv$current_pos <- rv$current_pos + 1
        } else {
          rv$is_running <- FALSE
        }
      })
    }
  })
  
  # Handle buttons
  observeEvent(input$start, {
    rv$is_running <- !rv$is_running
  })
  
  observeEvent(input$reset, {
    rv$current_pos <- 1
    rv$is_running <- FALSE
    smoothed_data(data.frame(time = numeric(0), y = numeric(0)))
    comparison_data(list())
  })
  
  observeEvent(input$complete, {
    rv$is_running <- FALSE
    current_data <- data()
    window <- input$window_size
    max_pos <- nrow(current_data) - window + 1
    
    all_smoothed <- data.frame(time = numeric(0), y = numeric(0), se = numeric(0))
    all_comp <- list()
    
    withProgress(message = 'Computing smoothed values...', value = 0, {
      for (pos in 1:max_pos) {
        window_data <- current_data[pos:(pos + window - 1), ]
        smooth_result <- calculate_smooth(window_data, return_variance = TRUE)
        
        all_smoothed <- rbind(all_smoothed, data.frame(
          time = pos + (window/2 - 1),
          y = smooth_result$fitted,
          se = sqrt(smooth_result$variance)
        ))
        
        if (input$compare_methods) {
          comp_results <- calculate_all_methods(window_data)
          for (method in names(comp_results)) {
            if (is.null(all_comp[[method]])) {
              all_comp[[method]] <- data.frame(time = numeric(0), y = numeric(0))
            }
            all_comp[[method]] <- rbind(
              all_comp[[method]],
              data.frame(time = pos + (window/2 - 1), y = comp_results[[method]])
            )
          }
        }
        
        incProgress(1/max_pos)
      }
    })
    
    smoothed_data(all_smoothed)
    comparison_data(all_comp)
    rv$current_pos <- max_pos + 1
  })
  
  observeEvent(input$jump, {
    max_pos <- nrow(data()) - input$window_size + 1
    rv$current_pos <- min(max(1, input$jump_to), max_pos)
  })
  
  # Update jump_to max value
  observe({
    max_pos <- nrow(data()) - input$window_size + 1
    updateNumericInput(session, "jump_to", max = max_pos)
  })
  
  # Render the main plot
  output$smoothing_plot <- renderPlot({
    current_data <- data()
    current_smoothed <- smoothed_data()
    
    p <- ggplot(current_data, aes(x = time, y = y)) +
      geom_point(alpha = 0.4, size = 1.5) +
      theme_minimal(base_size = 13) +
      labs(title = "Dynamic Nonparametric Smoothing",
           x = "Time", y = "Value") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # Add true signal if available
    if (!is.null(rv$true_signal)) {
      p <- p + geom_line(data = rv$true_signal,
                         aes(x = time, y = y),
                         color = "darkgreen", linetype = "dashed",
                         linewidth = 0.8, alpha = 0.6)
    }
    
    max_pos <- nrow(current_data) - input$window_size + 1
    if (rv$current_pos <= max_pos) {
      window_start <- rv$current_pos
      window_end <- rv$current_pos + input$window_size - 1
      window_data <- current_data[window_start:window_end, ]
      
      # Highlight window region
      p <- p + 
        annotate("rect", 
                 xmin = window_start, xmax = window_end,
                 ymin = -Inf, ymax = Inf,
                 fill = "lightblue", alpha = 0.1) +
        geom_vline(xintercept = c(window_start, window_end), 
                   color = "blue", alpha = 0.4, linetype = "dashed")
      
      # Method-specific visualizations
      if (input$smooth_method == "kernel") {
        kernel_fn <- kernel_functions[[input$kernel_type]]
        mid_point <- (window_start + window_end) / 2
        h <- input$bandwidth * sd(window_data$time)
        
        x_scaled_data <- (window_data$time - mid_point) / h
        weights <- kernel_fn(x_scaled_data)
        weights <- weights / sum(weights)
        weighted_mean <- sum(window_data$y * weights)
        
        p <- p + 
          geom_hline(yintercept = weighted_mean,
                     color = "blue", linetype = "dashed", linewidth = 0.8) +
          geom_point(data = data.frame(
            time = window_data$time,
            y = window_data$y,
            weight = weights
          ),
          aes(x = time, y = y, size = weight),
          color = "blue", alpha = 0.7) +
          scale_size_continuous(range = c(1, 6), guide = "none")
      } else if (input$smooth_method == "loess" && input$degree > 0) {
        if (input$degree == 1) {
          fit <- tryCatch(MASS::rlm(y ~ time, data = window_data),
                          error = function(e) lm(y ~ time, data = window_data))
        } else {
          fit <- lm(y ~ poly(time, 2), data = window_data)
        }
        
        pred_times <- seq(window_start, window_end, length.out = 50)
        fitted_df <- data.frame(
          time = pred_times,
          y = predict(fit, newdata = data.frame(time = pred_times))
        )
        p <- p + geom_line(data = fitted_df, aes(x = time, y = y),
                           color = "blue", linewidth = 1.2)
      } else {
        smooth_val <- calculate_smooth(window_data)
        p <- p + geom_segment(
          aes(x = window_start, xend = window_end,
              y = smooth_val, yend = smooth_val),
          color = "blue", linewidth = 1.2)
      }
    }
    
    # Add smoothed line with confidence band
    if (nrow(current_smoothed) > 0) {
      if (input$show_confidence && "se" %in% names(current_smoothed)) {
        p <- p + geom_ribbon(
          data = current_smoothed,
          aes(x = time, ymin = y - 1.96*se, ymax = y + 1.96*se),
          fill = "red", alpha = 0.2
        )
      }
      
      p <- p + geom_line(data = current_smoothed, 
                         aes(x = time, y = y),
                         color = "red", linewidth = 1.2)
    }
    
    p
  })
  
  # Comparison plot
  output$comparison_plot <- renderPlot({
    if (!input$compare_methods) return(NULL)
    
    comp_data <- comparison_data()
    if (length(comp_data) == 0) return(NULL)
    
    combined <- do.call(rbind, lapply(names(comp_data), function(method) {
      df <- comp_data[[method]]
      if (nrow(df) > 0) {
        df$method <- method
        df
      } else {
        NULL
      }
    }))
    
    if (is.null(combined) || nrow(combined) == 0) return(NULL)
    
    ggplot(combined, aes(x = time, y = y, color = method)) +
      geom_line(linewidth = 1) +
      scale_color_brewer(palette = "Set1", name = "Method") +
      theme_minimal(base_size = 12) +
      labs(title = "Method Comparison",
           x = "Time", y = "Smoothed Value") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom")
  })
  
  # Progress bar
  output$progress_bar <- renderUI({
    current_data <- data()
    max_pos <- nrow(current_data) - input$window_size + 1
    progress_pct <- min(100, round((rv$current_pos / max_pos) * 100))
    
    HTML(sprintf('
      <div class="progress">
        <div class="progress-bar progress-bar-striped" 
             role="progressbar" 
             style="width: %d%%"
             aria-valuenow="%d" 
             aria-valuemin="0" 
             aria-valuemax="100">
          %d%%
        </div>
      </div>
    ', progress_pct, progress_pct, progress_pct))
  })
  
  # Comparison metrics
  output$comparison_metrics <- renderUI({
    if (!input$compare_methods || is.null(rv$true_signal)) return(NULL)
    
    current_smoothed <- smoothed_data()
    comp_data <- comparison_data()
    
    if (nrow(current_smoothed) < 10) return(HTML("<p>Computing...</p>"))
    
    # Calculate MSE for current method
    merged <- merge(current_smoothed, rv$true_signal, by = "time")
    mse_current <- mean((merged$y.x - merged$y.y)^2)
    
    # Calculate MSE for comparison methods
    mse_comp <- sapply(names(comp_data), function(method) {
      if (nrow(comp_data[[method]]) > 0) {
        merged <- merge(comp_data[[method]], rv$true_signal, by = "time")
        if (nrow(merged) > 0) {
          mean((merged$y.x - merged$y.y)^2)
        } else NA
      } else NA
    })
    
    HTML(sprintf('
      <strong>MSE vs True Signal:</strong><br/>
      Current (%s): %.4f<br/>
      %s
    ', input$smooth_method, mse_current,
                 paste(sprintf("%s: %.4f", names(mse_comp), mse_comp), collapse = "<br/>")))
  })
  
  # Statistics output
  output$stats <- renderText({
    method_info <- if(input$smooth_method == "kernel") {
      sprintf("%s kernel (BW: %.1f)", 
              tools::toTitleCase(input$kernel_type), 
              input$bandwidth)
    } else if(input$smooth_method == "loess") {
      sprintf("LOESS (degree %d)", input$degree)
    } else {
      tools::toTitleCase(input$smooth_method)
    }
    
    max_pos <- nrow(data()) - input$window_size + 1
    
    sprintf("Position: %d / %d\nWindow: %d points\nMethod: %s\nPoints smoothed: %d",
            rv$current_pos,
            max_pos,
            input$window_size,
            method_info,
            nrow(smoothed_data()))
  })
  
  # Method explanations
  output$method_explanation <- renderUI({
    explanation <- switch(input$smooth_method,
                          "mean" = HTML("
        <h4>Moving Average</h4>
        <p>Calculates the arithmetic mean of all points within the window. 
        Each point receives equal weight (1/n). Simple and intuitive, but 
        sensitive to outliers. The blue line shows the current window's average.</p>
      "),
                          "median" = HTML("
        <h4>Moving Median</h4>
        <p>Takes the middle value within the window. Highly robust to outliers
        and extreme values. May produce less smooth results than the mean, but
        better preserves sharp transitions and edges in the data.</p>
      "),
                          "trimmed" = HTML(sprintf("
        <h4>Trimmed Mean (%.0f%% trim)</h4>
        <p>Removes %.0f%% of extreme values from each end before averaging.
        Combines smoothness of the mean with robustness of the median. Good
        balance for data with occasional outliers.</p>
      ", input$trim_percent * 100, input$trim_percent * 100)),
                          "loess" = HTML(sprintf("
        <h4>Local Regression (Polynomial Degree %d)</h4>
        <p><strong>Degree 0:</strong> Weighted constant (like kernel smoothing)<br/>
        <strong>Degree 1:</strong> Fits local linear trend (good for varying slopes)<br/>
        <strong>Degree 2:</strong> Fits local quadratic (captures curvature)</p>
        <p>The blue %s shows how the method adapts to local patterns.</p>
      ", input$degree, ifelse(input$degree == 0, "line", 
                              ifelse(input$degree == 1, "line", "curve")))),
                          "kernel" = HTML(sprintf("
        <h4>Kernel Smoothing (%s Kernel)</h4>
        <p>Computes a weighted average where points near the window center
        receive higher weights. The point size shows each observation's weight.
        The dashed line shows the weighted average result.</p>
        <p><strong>Bandwidth:</strong> Controls weight decay. Higher = smoother
        but may miss features. Lower = tracks data closely but more variable.</p>
        <p><strong>Kernel shapes:</strong></p>
        <ul style='margin-left: 20px;'>
          <li><strong>Gaussian:</strong> Smooth bell curve, never zero</li>
          <li><strong>Epanechnikov:</strong> Parabolic, optimal efficiency</li>
          <li><strong>Triangular:</strong> Linear decay to zero</li>
          <li><strong>Uniform:</strong> Equal weights within bandwidth</li>
        </ul>
      ", tools::toTitleCase(input$kernel_type)))
    )
    
    additional_info <- ""
    if (!is.null(rv$true_signal)) {
      additional_info <- "<p><em>Green dashed line shows the true underlying signal (no noise).</em></p>"
    }
    
    tagList(
      h3("Method Explanation"),
      explanation,
      HTML(additional_info)
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)