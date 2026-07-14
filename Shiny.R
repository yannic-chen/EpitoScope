# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Script Name:        Shiny.R
# Purpose:            This script creates a shiny app for easy immunopeptidomics analysis
# Author:             Yannic Chen
# Date Created:       2025-11-19
# Last Modified:      2026-03-20
# Version:            0.3
# R Version:          4.5.2
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Details:
# This script contains the Shiny app skeleton.
# Helper and plotting functions are mainly in the global.R
# The data loading function in the app is disabled, since we are working with named dataframes.
#
# IMPORTANT:
# The Shiny App consists of 4 mandatory scripts and 1 optional script. For simplicity sake, keep them in the same folder.
# Mandatory:
# Shiny.R - This contain the server functions and the executable shiny app.
# ui.R - This contains the script for the UI
# global.R - This contains the functions to draw plots and various helper functions.
# report.RmD - This script generated a pdf report. Technically optional if you never generate reports.
# Optional:
# data_loading.R - contains Example on how to load the data for use in shiny app.

source("global.R") #global.R must be in the same folder. Otherwise change this path.
source("ui.R") #ui.R must be in the same folder. Otherwise change this path.

options(shiny.maxRequestSize = 5*1024^3) #Increase upload limit (in bytes) if needed. 1024^3 = 1 GB
options(width=10000) #This allows for text to not be text-wrapped.
ht_opt$message <- FALSE

server <- function(input, output, session, input_variable, generate_pseudo_sequence = FALSE, custom_schema = NULL, custom_signature = NULL, replace_schema = FALSE) {
#------------------State Check---------------------
  ## State container
  startup_done <- reactiveVal(FALSE)
  motif_plot_length <- 7:20
  netmhcpan_path <- "/mnt/c/Users/Yannic/netMHCpan-4.2/netMHCpan" # This is the absolute path in the WSL. Really want the system to read off .bashrc
  wsl_available <- reactiveVal(NULL)
  netmhcpan_available <- reactiveVal(NULL)
  
  observe({
    if (startup_done()) return() #Since there is no reactive dependency, this observe only runs once anyway. But just in case.
    #### --- WSL CHECK ---
    wsl_ok <- tryCatch({
      res <- system2("wsl", "--status", stdout = TRUE, stderr = TRUE)
      !is.null(res)
    }, error = function(e) FALSE)
    wsl_available(wsl_ok)
    message(paste0("WSL available: ", wsl_ok))
    
    #### --- netMHCpan CHECK (only if WSL exists) ---
    netmhcpan_ok <- FALSE
    if (wsl_ok) {
      netmhcpan_ok <- tryCatch({
        # Use -o (ignore output) and check exit status
        status <- system2(
          "wsl",
          c("test", "-x", shQuote(netmhcpan_path)),
          stdout = FALSE,
          stderr = FALSE
        )
        status == 0   # TRUE if executable exists
      }, error = function(e) FALSE)
    }
    
    netmhcpan_available(netmhcpan_ok)
    message(paste0("NetMHCpan available: ", netmhcpan_ok))
    
    startup_done(TRUE)
  })
  
  ##----Disable buttons------
  observe({
    req(!is.null(netmhcpan_available()))
    
    if (!netmhcpan_available()) {
      shinyjs::disable("run_netmhc")
    } else {
      shinyjs::enable("run_netmhc")
    }
  })
  
  output$netmhc_status <- renderText({
    # Check if WSL is available
    if (!wsl_available()) {
      return("Disabled: WSL is not available on this system.")
    }
    
    # Check if netMHCpan executable exists
    if (!netmhcpan_available()) {
      return(paste0("Disabled: netMHCpan not found at ", netmhcpan_path))
    }
    
    # Otherwise, no message
    ""
  })
  
#------------------Data management-----------------
  ## Reactive dataset container
  raw_list_r <- reactiveVal(NULL) #this is the list of raw data
  data_list_r <- reactiveVal(NULL) #this is the list of trimmed and filtered data
  data_info_r <- reactiveVal(NULL) #this is the info list to know which columns are used for what
  data_mod_map <- reactiveVal(NULL) #this is the conversion map when modified amino acids are given their own symbol. Only used in PTM analysis.
  prediction_cache <- reactiveVal(data.frame(Peptide = character())) #This is to save netMHCpan predictions
  binder_summary_all <- reactive(NULL)
  peptide_wide_unique <- reactive(NULL)
  annotation_provided <- reactiveVal(FALSE)
  annotation_df_r <- reactiveVal(NULL)
  condition_groups_r <- reactiveVal(list())  # list(name → list(samples, expr))
  active_expr_r      <- reactiveVal(list())  # expression being built
  editing_group_r    <- reactiveVal(NULL)    # name of group being edited, or NULL
  measurement_col_map_r <- reactiveVal(NULL)
  use_measurements <- reactiveVal(FALSE) # even when the measurement_col exist, sometimes we dont want to use the shortcut method without assigning groups by condition.
  
  #Preloaded Data
  observe({
    if (is.null(data_list_r())) {
      start <- Sys.time() #measure time
      
      register_custom_schema(custom_schema, custom_signature, replace_schema)
      
      if (!is.null(custom_signature)) {
        if (replace_schema) {
          signature <<- list()
          message(paste0("Replace default signature"))
        }
      }
      
      if (is.data.frame(input_variable)) {
        # Annotation table provided — load data from it
        loaded <- load_from_annotation(input_variable)
        annotation_provided(TRUE)
        raw_list_r(loaded)
        annotation_df_r(attr(loaded, "annotation"))
      } else if(is.list(input_variable)) {
        annotation_provided(FALSE)
        # must be named
        if (is.null(names(input_variable)) || any(names(input_variable) == "")) {
        }
        
        # all elements must be data.frames
        if (!all(vapply(input_variable, is.data.frame, logical(1)))) {
        }
        raw_list_r(input_variable)
      } else {
        stop("input must be either a data.frame or a named list of data.frames.")
      }
      
      processed <- lapply(raw_list_r(), normalize_df)
      print(Sys.time() - start)
      
      #Add the PTM_Pseudo sequence to each dataframe in the list. We do it outside of normalize_df() function because it unifies the mod_map across all dataframes in the list
      if (generate_pseudo_sequence) {
        start <- Sys.time()
        print("Generating PTM_Pseudo sequence")
        all_peptides <- unlist(lapply(processed, function(sample) {sample$df$PEPTIDE}), use.names = FALSE)
        all_tokens   <- extract_mod_tokens(all_peptides)
        global_mod_map <- build_mod_map(all_tokens)
        processed <- lapply(processed, function(sample) {
          df <- sample$df
          
          if (!is.null(df$PEPTIDE)) {
            df$PTM_Pseudo <- convert_with_mod_map(df$PEPTIDE, global_mod_map)
          } else {
            df$PTM_Pseudo <- NA_character_
          }
          
          sample$df <- df
          sample
        })
        print(Sys.time() - start)
      } else {
        print("SKIP: Generating PTM_Pseudo sequence")
        global_mod_map <- NA
      }

      start <- Sys.time()
    
      # Extract normalized dataframes
      dfs <- lapply(processed, function(x) x$df)
      
      # Add netMHCpan info to dataframes
      if (exists("netMHCpan")) {
        print("left_join netMHCpan pre-generated data.")
        
        # Collect all unique peptides across all dfs upfront
        all_peptides <- unique(unlist(lapply(dfs, `[[`, "STRIPPED")))
        
        # Filter netMHCpan once, before the loop
        netMHCpan_filtered <- netMHCpan[Peptide %in% all_peptides]
        
        dfs <- lapply(dfs, function(df) {
          df_dt <- as.data.table(df)
          result <- netMHCpan_filtered[df_dt, on = c(Peptide = "STRIPPED")]
          result <- as.data.frame(result)
          names(result)[names(result) == "Peptide"] <- "STRIPPED"
          result
        })
      } else {
        print("No netMHCpan pro-computed variable. SKIP binding prediction")
      }
      #The predicted_cache is initialized using all peptides in the input data. If the netMHCpan precomputed data has been left_joined, these will also be taken.
      prediction <- do.call(rbind, lapply(dfs, function(df) {
        hla_cols <- grep("^HLA", colnames(df), value = TRUE)
        
        if (length(hla_cols) == 0) {
          df_subset <- data.frame(Peptide = df$STRIPPED, stringsAsFactors = FALSE) #if only Peptide column exist, then R automatically formats to matrix. We need to enforce dataframe format.
        } else {
          df_subset <- df[, c("STRIPPED", hla_cols), drop = FALSE]
          colnames(df_subset)[1] <- "Peptide"
        }
        
        df_subset <- df_subset[!duplicated(df_subset$Peptide), , drop = FALSE]
        df_subset
      }))
      
      # Extract summary tables
      infos <- lapply(processed, function(x) x$table)
      merged_info <- infos %>%
        imap(~ .x %>% dplyr::rename(!!.y := coalesced_list)) %>%  # .y = sample name
        purrr::reduce(full_join, by = "final_name")

      print(Sys.time() - start)
      
      if(annotation_provided()) {
        quantity_cols <- get_quantity_cols(merged_info)
        
        dfs_subset <- lapply(dfs, function(df) {
          df[, intersect(colnames(df), quantity_cols), drop = FALSE]
        })
        
        measurement_col_map_r(
          build_measurement_col_map(dfs_subset, attr(loaded, "annotation"))
        )
      }
  
      prediction_cache(distinct(prediction))
      data_list_r(dfs)
      data_info_r(merged_info)
      data_mod_map(global_mod_map)
    }
  })
  
  ## Data upload from GUI
  observeEvent(input$show_transform, {
    req(input$files)
    
    df_list <- list()
    
    for (i in 1:nrow(input$files)) {
      file <- input$files[i, ]
      ext <- tools::file_ext(file$name)
      
      df <- switch(ext,
                   "csv" = read.csv(file$datapath),
                   "tsv" = read.delim(file$datapath),
                   read.delim(file$datapath))
      
      
      
      sample_name <- tools::file_path_sans_ext(file$name)
      df_list[[sample_name]] <- df
    }
    
    raw_list_r(df_list)
    
    processed <- lapply(df_list, normalize_df)
    dfs <- lapply(processed, function(x) x$df)
    infos <- lapply(processed, function(x) x$table)
    merged_info <- infos %>%
      imap(~ .x %>% dplyr::rename(!!.y := coalesced_list)) %>%  # rename value column to sample name
      purrr::reduce(full_join, by = "final_name")
    
    data_list_r(dfs)
    data_info_r(merged_info)
  })
  
  ## Allow selecting/deseleting data
  output$sample_selector <- renderUI({
    req(data_list_r())
    checkboxGroupInput(
      "selected_samples",
      tagList(icon("list-check"),"Select samples:"),
      choices = names(data_list_r()),
      selected = names(data_list_r())
    )
  })
  
  ## Active data (selected samples)
  active_data_list <- reactive({
    req(data_list_r(), input$selected_samples)
    data_list_r()[input$selected_samples]
  })
  
#-------------------Data Transformation tab-----------------------
  # Add delayed reaction
  filter_inputs <- reactive({
    list(
      samples               = input$selected_samples,
      length_range          = input$length_range,
      quantity_range        = input$quantity_range,
      score_range           = input$score_range,
      charge_range          = input$charge_range,
      mass_range            = input$mass_range,
      RT_range              = input$RT_range,
      upset_min_size    = input$upset_min_size,
      upset_min_degree  = input$upset_min_degree,
      upset_n_intersect = input$upset_n_intersect
    )
  })
  
  filter_inputs_d <- debounce(filter_inputs, millis = function() {
    isolate({
      delay_s <- input$debounce_delay_s
      if (is.null(delay_s) || is.na(delay_s) || delay_s < 0) return(2000)
      as.integer(delay_s * 1000)
    })
  })
  
  # Preprocessing data to align it with GUI
  processed_data_list <- reactive({
    
    lst <- active_data_list()  # Only selected samples
    req(lst)
    
    filters <- list(
      length_range   = filter_inputs_d()$length_range,
      quantity_range = filter_inputs_d()$quantity_range,
      score_range    = filter_inputs_d()$score_range,
      charge_range   = filter_inputs_d()$charge_range,
      mass_range     = filter_inputs_d()$mass_range,
      RT_range       = filter_inputs_d()$RT_range
    )
    apply_filters(lst, filters)
  })
  
  output$length_slider_ui   <- renderUI({ make_range_slider_ui(active_data_list(), "LENGTH",       "length_range",   tagList(icon("ruler-horizontal"), "Filter by Length:"),   step = 1) })
  output$quantity_slider_ui <- renderUI({ make_range_slider_ui(active_data_list(), "MAX_QUANTITY", "quantity_range", "Filter by Max Quantity:",                                step = 1) })
  output$score_slider_ui    <- renderUI({ make_range_slider_ui(active_data_list(), "SCORE",        "score_range",    "Filter by Score:") })
  output$charge_slider_ui   <- renderUI({ make_range_slider_ui(active_data_list(), "CHARGE",       "charge_range",   "Filter by Charge:",                                     step = 1) })
  output$mass_slider_ui     <- renderUI({ make_range_slider_ui(active_data_list(), "MASS",         "mass_range",     "Filter by Mass:",                                       digits = 2) })
  output$RT_slider_ui       <- renderUI({ make_range_slider_ui(active_data_list(), "RT",           "RT_range",       "Filter by RT:",                                         digits = 2) })
  
  
  observeEvent(input$generate_report, {
    
    tryCatch({
      # Choose an output file name
      out_file <- paste0("EpitoScope_Report_", Sys.Date(), ".html")
      
      # Use a temporary directory (required for Shiny Server / shinyapps.io)
      out_path <- file.path(tempdir(), out_file)
      
      # Show progress
      withProgress(message = "Generating report...", value = 0, {
        rmarkdown::render(
          input = "report.Rmd", #report.Rmd must be in the same folder. Otherwise change path here.
          output_file = out_path,
          params = list(
            data_info = data_info_r(),
            data_list = data_list_r(),
            mod_map_list = data_mod_map(),
            data_mod_map = data_mod_map(),
            processed_data_list = processed_data_list(),
            filters             = list(
              samples        = input$selected_samples,
              length_range   = input$length_range,
              quantity_range = input$quantity_range,
              score_range    = input$score_range,
              charge_range   = input$charge_range,
              mass_range     = input$mass_range,
              RT_range       = input$RT_range,
              min_presence_fraction = input$min_presence_fraction
            ),
            color_palette = input$color_palette,
            cluster_mode_ea     = input$cluster_mode_ea,
            upset_min_size    = input$upset_min_size,
            upset_min_degree  = input$upset_min_degree,
            upset_n_intersect = input$upset_n_intersect,
            binder_summary_all = safe_reactive(binder_summary_all),
            binder_unique = safe_reactive(peptide_wide_unique),
            group_list = safe_reactive(group_list),
            group_comp_data = safe_reactive(group_comp_data),
            col_map = safe_reactive(measurement_col_map_r),
            annotation_table = if (is.data.frame(input_variable)) input_variable else NULL,
            default_quantity_cols_r = safe_reactive(default_quantity_cols_r)
          ),
          envir = new.env(parent = globalenv())
        )
        
        incProgress(0.2, detail = "Finalizing report...")
        
        incProgress(0.0, detail = "Done!")
      })
      
      # Let the user download it
      showModal(modalDialog(
        title = "Report ready!",
        "Click below to download your HTML report.",
        downloadButton("download_report", "Download"),
        easyClose = TRUE
      ))
      
      # Store path so downloadHandler can access it
      output$download_report <- downloadHandler(
        filename = function() { out_file },
        content = function(file) {file.copy(out_path, file)}
      )
    }, error = function(e) {
      showNotification(paste("Report error:", conditionMessage(e)), type = "error", duration = 15)
    })
  })
  
  default_quantity_cols_r <- reactive({
    get_quantity_cols(data_info_r())
  })
  
#---------------------Summary Tab-------------------------
  ## ---------Column Mapping----------------
  output$summary_table <- DT::renderDT({
    req(data_info_r())
    DT::datatable(
      Colum_mapping_table(data_info_r()),
      options = list(pageLength = 10, autoWidth = TRUE),
      escape = FALSE,
      rownames = FALSE
    )
  })
  
  ## ---------Annotation Table----------------
  output$annotation_table <- DT::renderDT({
    if(!annotation_provided()){
      updatebs4Card(id = "annotation_card", session = session, action = "remove")
      shiny::validate(shiny::need(FALSE, "No annotation table provided."))
    }
    
    ann <- input_variable
    
    DT::datatable(
      ann,
      rownames = FALSE,
      options  = list(pageLength = 10, autoWidth = TRUE),
      escape = FALSE,
    )
  })
  
  ## ----Number of peptides and peptidoforms----
  output$summary_peptides_plot <- renderPlot({
    lst <- data_list_r()
    check_data_error(lst, na_policy = "ignore")
    plot_unique_counts(lst, column = "STRIPPED", y_label = "Number of unique peptides", color = input$color_palette)
  })
  
  output$summary_peptidoforms_plot <- renderPlot({
    lst <- data_list_r()
    check_data_error(lst, na_policy = "ignore")
    plot_unique_counts(lst, column = "PEPTIDE", y_label = "Number of unique peptidoforms", color = input$color_palette)
  })
  
  output$summary_proteins_plot <- renderPlot({
    lst <- data_list_r()
    check_data_error(lst, required_cols = "PROTEIN", na_policy = "ignore")
    # Use extract_protein_prefixes for proteins
    plot_unique_counts(lst, column = "PROTEIN", y_label = "Number of unique proteins",
                       transform_fn = extract_protein_prefixes, color = input$color_palette)
  })
  
  ## ----Charge / Mass / mz / RT / ppm----
  output$charge_plot <- renderPlot({
    lst <- data_list_r()
    check_data_error(lst, required_cols = "CHARGE", na_policy = "any") #Even if CHARGE is missing, preprocessing would add a charge of 1 to every row.
    plot_stacked_bar(lst, column = "CHARGE", fill_label = "Charge", percentage = FALSE, color = input$color_palette)
  })
  
  # Example usage for your Shiny outputs
  output$mass_plot <- plotly::renderPlotly({
    lst <- data_list_r()
    check_data_error(lst, required_cols = "MASS", na_policy = "any")
    plot_density(lst, column = "MASS", x_label = "Mass (Da)", color = input$color_palette)
  })
  
  output$mz_plot <- plotly::renderPlotly({
    lst <- data_list_r()
    check_data_error(lst, required_cols = "MZ", na_policy = "any")
    plot_density(lst, column = "MZ", x_label = "m/z", color = input$color_palette)
  })
  
  output$RT_plot <- renderUI({
    lst <- data_list_r()
    # Wrap plots in a grid (like motif plots)
    layout_column_wrap(
      width = "400px",  # each plot approx width
      !!!lapply(names(lst), function(sample_name) {
        plotOutput(paste0("RT_", sample_name), height = "300px")
      })
    )
  })
  
  observe({
    lst <- data_list_r()

    for (sample_name in names(lst)) {
      
      local({
        sample_local <- sample_name
        df_local <- lst[[sample_local]]
        
        output[[paste0("RT_", sample_local)]] <- renderPlot({
          check_data_error(df_local, required_cols = "RT", na_policy = "any")
          plot_histogram(df = df_local, column = "RT", x_label = "Retention Time", 
                         title_name = paste("RT Histogram –", sample_local), color = input$color_palette
          )
        })
      })
    }
  })
  
  output$ppm_plot <- plotly::renderPlotly({
    lst <- data_list_r()
    check_data_error(lst, required_cols = "PPM", na_policy = "any")
    plot_density(lst, column = "PPM", x_label = "ppm", color = input$color_palette)
  })
  
  ##----Score distribution----
  output$score_violin <- plotly::renderPlotly({
    lst <- data_list_r()
    check_data_error(lst, required_cols = "SCORE", na_policy = "any")
    plot_violin(lst, column = "SCORE", color = input$color_palette)
  })
  
  ##----Summary Table----
  output$RAW_summary_html <- renderUI({
    lst <- raw_list_r()
    check_data_error(lst, na_policy = "ignore")
    render_summary_pre(lst, font_size = "8px")
  })
#---------------------QC Tab-------------------------
  ## ---- length distribution ----
  output$length_plot <- renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "LENGTH", na_policy = "any")
    plot_length_distribution(lst, color = input$color_palette)
  })
  
  ## ---- length range percentage ----
  output$length_range_percentage <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "LENGTH", na_policy = "any")
    lst <- lapply(lst, function(df) {
      df$correct_range <- df$LENGTH >= 8 & df$LENGTH <= 13
      df
    })
    plot_stacked_bar(lst, column = "correct_range", fill_label = "8-13mer", percentage = TRUE, color = input$color_palette)
  })
  
  output$length_range_percentage_meas <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "LENGTH", na_policy = "any")
    shiny::validate(shiny::need(length(default_quantity_cols_r()) > 0, "No QUANTITY columns found."))
    plot_length_range_per_measurement(lst, default_quantity_cols_r(), color = input$color_palette)
  })
  
  ## ---- Other numeric columns on per measurement basis ----
  output$charge_plot_meas <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "CHARGE", na_policy = "all")
    plot_charge_per_measurement(lst, default_quantity_cols_r(), color = input$color_palette)
  })
  
  output$mass_plot_meas <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MASS", na_policy = "any")
    plot_density_envelope(lst, default_quantity_cols_r(), "MASS", "Mass (Da)", color = input$color_palette)
  })
  
  output$mz_plot_meas <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MZ", na_policy = "any")
    plot_density_envelope(lst, default_quantity_cols_r(), "MZ", "m/z", color = input$color_palette)
  })
  
  output$ppm_plot_meas <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PPM", na_policy = "any")
    plot_density_envelope(lst, default_quantity_cols_r(), column = "PPM", x_label = "ppm", color = input$color_palette)
  })
  
  output$score_violin_meas <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "SCORE", na_policy = "any")
    plot_violin_envelope(lst, default_quantity_cols_r(), column = "SCORE", x_label = "Score", color = input$color_palette)
  })
  
  output$RT_plot_meas <- renderUI({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "RT", na_policy = "all")
    qcols <- default_quantity_cols_r()
    layout_column_wrap(
      width = "400px",
      !!!lapply(names(lst), function(s) {
        pid <- paste0("rt_meas_", s)
        output[[pid]] <- renderPlot({
          plot_rt_histogram_range(lst[[s]], s, qcols)
        })
        plotOutput(pid, height = "300px")
      })
    )
  })
  
  ## ---- Motif Plot ----
  global_legend_plot <- suppressWarnings(ggseqlogo::ggseqlogo("ACDEFGHIKLMNPQRSTVWY") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Amino Acid Colors"))
  
  for (L in motif_plot_length) {
    local({
      L_local <- L
      output[[paste0("motif_legend_", L_local)]] <- renderPlot({
        global_legend_plot
      })
    })
  }
  
  ## Motif plots
  output$motif_tabs <- renderUI({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    
    # Build a list of tabPanels for each peptide length
    length_tabs <- lapply(motif_plot_length, function(L) {
      
      # For each length, build inner sample plots
      sample_plots <- lapply(names(lst), function(sample_name) {
        plotOutput(
          paste0("motif_", sample_name, "_", L),
          height = "180px"
        )
      })
      
      # Return the tabPanel for this length
      tabPanel(
        title = paste("Length", L),
        
        # Legend at the top
        tags$div(
          style = "text-align:center; margin-bottom:10px;",
          plotOutput(paste0("motif_legend_", L), height = "120px")
        ),
        
        # Layout for sample plots
        do.call(layout_column_wrap, c(width = "250px", sample_plots))
      )
    })
    
    # Pass the list of tabPanels to tabsetPanel
    do.call(tabsetPanel, c(id = "motif_length_tabs", length_tabs))
  })
  
  observe({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    
    for (L in motif_plot_length) {
      for (sample_name in names(lst)) {
        
        local({
          length_val <- L
          sample_val <- sample_name
          
          output[[paste0("motif_", sample_val, "_", length_val)]] <- renderPlot({
            
            df <- lst[[sample_val]]
            peptides <- unique(df$STRIPPED[df$LENGTH == length_val])
            
            shiny::validate(shiny::need(length(peptides) >= 5, "Not enough peptides. Need at least 5"))
            
            par(mar = c(1.5, 1.5, 2, 0.5))
            plot_seqlogo(peptides,title = paste(sample_val, "• Length", length_val))
          })
        })
      }
    }
  })
  
  ## ----Unique entries stats----
  output$summary_peptides_plot3 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    plot_unique_counts(lst, column = "STRIPPED", y_label = "Number of unique peptides",
                       color = input$color_palette, quantity_cols = default_quantity_cols_r())
  })
  
  output$summary_peptidoforms_plot3 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    plot_unique_counts(lst, column = "PEPTIDE", y_label = "Number of unique peptidoforms",
                       color = input$color_palette, quantity_cols = default_quantity_cols_r())
  })
  
  output$summary_proteins_plot3 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PROTEIN", na_policy = "all")
    plot_unique_counts(lst, column = "PROTEIN", y_label = "Number of unique proteins",
                       transform_fn = extract_protein_prefixes,
                       color = input$color_palette, quantity_cols = default_quantity_cols_r())
  })
  
  ## ----Dynamic Range plots (individual, one subplot figure)----
  output$dynrange_individual_ui <- renderUI({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
    n     <- sum(vapply(lst, function(df)
      nrow(df) >= 10 && "MAX_QUANTITY" %in% colnames(df), logical(1)))
    nrows <- max(1L, ceiling(n / 3))
    plotly::plotlyOutput("dynrange_grid", height = paste0(nrows * 320, "px"))
  })
  
  output$dynrange_grid <- renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
    dynamic_range_subplot(lst, data_col = "MAX_QUANTITY", ncol = 3)
  })
  
  observeEvent(list(input$dynrange_search, input$dynrange_pep_search), {
    lst        <- processed_data_list()
    prot_query <- trimws(if (is.null(input$dynrange_search))     "" else input$dynrange_search)
    pep_query  <- trimws(if (is.null(input$dynrange_pep_search)) "" else input$dynrange_pep_search)
    
    samples <- names(lst)[vapply(lst, function(df)
      nrow(df) >= 10 && "MAX_QUANTITY" %in% colnames(df), logical(1))]
    if (length(samples) == 0) return()
    
    match_idx_fn <- function(query, cols, n) {
      if (nchar(query) == 0) return(rep(FALSE, n))
      idx <- tryCatch(
        Reduce(`|`, lapply(cols, function(col)
          grepl(query, col, ignore.case = TRUE, perl = TRUE))),
        error = function(e) rep(FALSE, n))
      idx[is.na(idx)] <- FALSE
      idx
    }
    
    k <- length(samples)
    X <- vector("list", 3 * k); Y <- vector("list", 3 * k)
    T <- vector("list", 3 * k); H <- vector("list", 3 * k)
    
    for (i in seq_along(samples)) {
      df <- lst[[samples[i]]]
      df$MAX_QUANTITY[df$MAX_QUANTITY == 0] <- NA
      df <- df[!is.na(df$MAX_QUANTITY), , drop = FALSE]
      b <- (i - 1L) * 3L
      
      if (nrow(df) == 0) {
        X[[b+1]] <- list(); Y[[b+1]] <- list(); T[[b+1]] <- list(); H[[b+1]] <- "text"
        X[[b+2]] <- list(); Y[[b+2]] <- list(); T[[b+2]] <- list(); H[[b+2]] <- "text"
        X[[b+3]] <- list(); Y[[b+3]] <- list(); T[[b+3]] <- list(); H[[b+3]] <- "text"
        next
      }
      
      prot_vec <- if ("PROTEIN" %in% names(df)) df$PROTEIN else rep("", nrow(df))
      pep_vec  <- if ("PEPTIDE" %in% names(df)) df$PEPTIDE else rep("", nrow(df))
      df$Rank   <- rank(-df$MAX_QUANTITY, ties.method = "first")
      df$y_vals <- log2(df$MAX_QUANTITY)
      gene_nm   <- sub(".*?GN=([0-9A-Z/\\-]+).*", "\\1", prot_vec, perl = TRUE)
      disp_nm   <- ifelse(gene_nm == prot_vec, prot_vec, paste0(gene_nm, "<br>", prot_vec))
      df$hover_text <- paste0("<b>", disp_nm, "</b><br>", pep_vec, "<br>",
                              "log2(MAX_QUANTITY): ", round(df$y_vals, 2), "<br>Rank: ", df$Rank)
      
      n <- nrow(df)
      prot_idx <- match_idx_fn(prot_query, list(prot_vec), n)
      pep_idx  <- match_idx_fn(pep_query,  list(df$STRIPPED, pep_vec), n)
      pep_idx  <- pep_idx & !prot_idx
      base_idx <- !(prot_idx | pep_idx)
      
      base_hi <- if (sum(prot_idx) + sum(pep_idx) > 0) "skip" else "text"
      
      X[[b+1]] <- as.list(df$Rank[base_idx]); Y[[b+1]] <- as.list(df$y_vals[base_idx]); T[[b+1]] <- as.list(df$hover_text[base_idx]); H[[b+1]] <- base_hi
      X[[b+2]] <- as.list(df$Rank[prot_idx]); Y[[b+2]] <- as.list(df$y_vals[prot_idx]); T[[b+2]] <- as.list(df$hover_text[prot_idx]); H[[b+2]] <- "text"
      X[[b+3]] <- as.list(df$Rank[pep_idx]);  Y[[b+3]] <- as.list(df$y_vals[pep_idx]);  T[[b+3]] <- as.list(df$hover_text[pep_idx]);  H[[b+3]] <- "text"
    }
    
    idxs <- as.list(seq_len(3L * k) - 1L)
    plotly::plotlyProxy("dynrange_grid", session) %>%
      plotly::plotlyProxyInvoke("restyle",
                                list(x = X, y = Y, text = T, hoverinfo = H), idxs)
  })
  
  ## ----Dynamic Range plots combined----
  output$dynrange_combined <- renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
    dynamic_range_plot_combined(df_list = lst, data_col = "MAX_QUANTITY", color = input$color_palette)
  })
  
  observeEvent(list(input$dynrange_search, input$dynrange_pep_search), {
    lst        <- processed_data_list()
    prot_query <- trimws(if (is.null(input$dynrange_search))     "" else input$dynrange_search)
    pep_query  <- trimws(if (is.null(input$dynrange_pep_search)) "" else input$dynrange_pep_search)
    
    samples <- names(lst)[vapply(lst, function(df)
      nrow(df) >= 10 && "MAX_QUANTITY" %in% colnames(df), logical(1))]
    if (length(samples) == 0) return()
    
    match_idx_fn <- function(query, cols, n) {
      if (nchar(query) == 0) return(rep(FALSE, n))
      idx <- tryCatch(
        Reduce(`|`, lapply(cols, function(col)
          grepl(query, col, ignore.case = TRUE, perl = TRUE))),
        error = function(e) rep(FALSE, n))
      idx[is.na(idx)] <- FALSE
      idx
    }
    
    base_x <- base_y <- base_t <- vector("list", length(samples))
    prot_x <- numeric(0); prot_y <- numeric(0); prot_t <- character(0)
    pep_x  <- numeric(0); pep_y  <- numeric(0); pep_t  <- character(0)
    
    for (i in seq_along(samples)) {
      df <- lst[[samples[i]]]
      df$MAX_QUANTITY[df$MAX_QUANTITY == 0] <- NA
      df <- df[!is.na(df$MAX_QUANTITY), , drop = FALSE]
      if (nrow(df) == 0) { base_x[[i]] <- list(); base_y[[i]] <- list(); base_t[[i]] <- list(); next }
      
      prot_vec <- if ("PROTEIN" %in% names(df)) df$PROTEIN else rep("", nrow(df))
      pep_vec  <- if ("PEPTIDE" %in% names(df)) df$PEPTIDE else rep("", nrow(df))
      df$Rank   <- rank(-df$MAX_QUANTITY, ties.method = "first")
      df$y_vals <- log2(df$MAX_QUANTITY)
      gene_nm   <- sub(".*?GN=([0-9A-Z/\\-]+).*", "\\1", prot_vec, perl = TRUE)
      disp_nm   <- ifelse(gene_nm == prot_vec, prot_vec, paste0(gene_nm, "<br>", prot_vec))
      df$hover_text <- paste0("<b>", samples[i], "</b><br>", disp_nm, "<br>", pep_vec, "<br>",
                              "log2(MAX_QUANTITY): ", round(df$y_vals, 2), "<br>Rank: ", df$Rank)
      
      n <- nrow(df)
      prot_idx <- match_idx_fn(prot_query, list(prot_vec), n)
      pep_idx  <- match_idx_fn(pep_query,  list(df$STRIPPED, pep_vec), n)
      pep_idx  <- pep_idx & !prot_idx
      base_idx <- !(prot_idx | pep_idx)
      
      base_x[[i]] <- as.list(df$Rank[base_idx])
      base_y[[i]] <- as.list(df$y_vals[base_idx])
      base_t[[i]] <- as.list(df$hover_text[base_idx])
      
      prot_x <- c(prot_x, df$Rank[prot_idx]); prot_y <- c(prot_y, df$y_vals[prot_idx]); prot_t <- c(prot_t, df$hover_text[prot_idx])
      pep_x  <- c(pep_x,  df$Rank[pep_idx]);  pep_y  <- c(pep_y,  df$y_vals[pep_idx]);  pep_t  <- c(pep_t,  df$hover_text[pep_idx])
    }
    
    base_hoverinfo <- if (length(prot_x) + length(pep_x) > 0) "skip" else "text"
    
    x_vals  <- c(base_x, list(as.list(prot_x)), list(as.list(pep_x)))
    y_vals  <- c(base_y, list(as.list(prot_y)), list(as.list(pep_y)))
    t_vals  <- c(base_t, list(as.list(prot_t)), list(as.list(pep_t)))
    hi_vals <- c(rep(list(base_hoverinfo), length(samples)), list("text"), list("text"))
    idxs    <- as.list(seq_len(length(samples) + 2L) - 1L)
    
    plotly::plotlyProxy("dynrange_combined", session) %>%
      plotly::plotlyProxyInvoke("restyle",
                                list(x = x_vals, y = y_vals, text = t_vals, hoverinfo = hi_vals), idxs)
  })
  
  ##----1/k0 vs m/z----
  output$scatterplots_ui <- renderUI({
    lst <- processed_data_list()
    
    layout_column_wrap(
      width = "300px",
      !!!lapply(names(lst), function(sample_name) {
        plotOutput(paste0("scatter_", sample_name), height = "300px")
      })
    )
  })
  
  observe({
    lst <- processed_data_list()
    
    for (sample_name in names(lst)) {
      local({
        sample_local <- sample_name
        df_local     <- lst[[sample_local]]
        
        output[[paste0("scatter_", sample_local)]] <- renderPlot({
          check_data_error(df_local, required_cols = c("MZ", "K0"), na_policy = "all")
          generate_scatterplot(df_local, sample_local, color = input$color_palette)
        })
      })
    }
  })
  
  ## ---- measurement specific heatmap ----
  cor_mat_cache <- reactiveVal(matrix(NA, nrow = 5, ncol = 5)) #necessary as the space reserving rungs before the code.
  
  output$measurement_heatmap <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    
    shiny::validate(shiny::need(length(default_quantity_cols_r()) > 0, "No QUANTITY columns found."))
    
    pep_mat <- prepare_measurement_matrix(lst, default_quantity_cols_r())
    shiny::validate(shiny::need(nrow(pep_mat) > 0, "No peptides to plot."))
    
    # Drop columns with no data at all (would make cor() fail entirely)
    pep_mat <- pep_mat[, colSums(!is.na(pep_mat)) > 0, drop = FALSE]
    shiny::validate(shiny::need(ncol(pep_mat) >= 2, "Need at least 2 measurements with data for correlation."))
    
    cor_mat <- cor(pep_mat, method = "pearson", use = "pairwise.complete.obs")
    
    # Replace any remaining NAs (pairs with zero shared peptides) with 0
    cor_mat[is.na(cor_mat)] <- 0
    cor_mat_cache(cor_mat)
    
    group_names <- names(lst)
    get_group <- function(colname) {
      matched <- group_names[sapply(group_names, function(g) startsWith(colname, g))]
      if (length(matched) == 0) return(NA)
      matched[1]
    }
    groups <- sapply(colnames(cor_mat), get_group)
    
    #remove the sample prefix. However, if there will be duplicate names, keep the sample prefix
    colnames(cor_mat) <- sub("^.* \\| ", "", colnames(cor_mat))
    rownames(cor_mat) <- colnames(cor_mat)
    
    fs <- max(6, min(10, floor(800 / nrow(cor_mat))))
    
    stripped <- sub("^.* \\| ", "", colnames(cor_mat))
    if (anyDuplicated(stripped) == 0) {
      colnames(cor_mat) <- stripped
      rownames(cor_mat) <- stripped
    }
      
    if(input$cluster_mode_ea == "sample") {
      ht <- plot_heatmap(cor_mat,color = input$color_palette, cluster = "none", row_groups = groups, col_groups = groups, label = "Pearson", fontsize = fs)
    } else if (input$cluster_mode_ea == "mix") {
      ht <- plot_heatmap(cor_mat,color = input$color_palette, cluster = "rows",row_groups = NULL, col_groups = groups, label = "Pearson", fontsize = fs)
    } else {
      ht <- plot_heatmap(cor_mat,color = input$color_palette, cluster = input$cluster_mode_ea, label = "Pearson", fontsize = fs)
    }
    
    safe_draw(ht, heatmap_legend_side = "left"
    )
  }, width = "auto", 
  height = function() {
    n <- nrow(cor_mat_cache())
    min(500, n * 12 + 200)
  })
  
  ##----aa heatmap----
  output$aa_heatmap <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    
    shiny::validate(shiny::need(length(lst) >= 2, "need 2 or more samples to plot"))
    
    ht <- plot_aa_composition(
      lst,
      color = input$color_palette
    )
    
    safe_draw(ht)
    
  }, res = 144,
  height = function() {(length(processed_data_list()) * 50) + 500 })
#---------------------Results Tab-------------------------
  ##----Data Completeness----
  output$completeness_plot <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    spectra_cols <- data_info_r() %>%
      dplyr::filter(final_name == "SPECTRA") %>%          #Should be SPECTRA, but cannot currently do because 0 is missing in PEAKS, while 0 is existing in Fragpipe and NA is missing here.
      dplyr::select(-final_name) %>%   # all sample columns.
      unlist(recursive = TRUE, use.names = FALSE)
    
    plot_completeness(lst, spectra_cols, color = input$color_palette)
  })
  
  output$completeness_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    spectra_cols <- data_info_r() %>%
      dplyr::filter(final_name == "SPECTRA") %>%
      dplyr::select(-final_name) %>%   # all sample columns.
      unlist(recursive = TRUE, use.names = FALSE)
    
    plot_completeness(lst, spectra_cols, percent = TRUE, color = input$color_palette)
  })
  
  ##----Upset----
  output$upset_plot <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore") #By default plot_upset takes STRIPPED column. Need to adjust if we take PEPTIDE column instead.
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples to plot"))
    plot_upset(lst,
               min_size        = filter_inputs_d()$upset_min_size,
               min_degree      = filter_inputs_d()$upset_min_degree,
               n_intersections = filter_inputs_d()$upset_n_intersect)
  })
  
  ### ---- Intersection data ----
  upset_intersections <- reactive({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples"))
    compute_upset_intersections(
      lst,
      stripped        = TRUE,
      min_size        = input$upset_min_size,
      min_degree      = input$upset_min_degree,
      n_intersections = input$upset_n_intersect
    )
  })
  
  ### ---- Intersection table ----
  output$upset_intersection_table <- DT::renderDT({
    ints <- upset_intersections()
    df   <- data.frame(
      Intersection = names(ints),
      Count        = sapply(ints, length),
      stringsAsFactors = FALSE
    )
    DT::datatable(
      df,
      selection = "single",
      rownames  = FALSE,
      options   = list(pageLength = 10, dom = "frtip")
    )
  })
  
  ### ---- Modal on row click ----
  observeEvent(input$upset_intersection_table_rows_selected, {
    row  <- input$upset_intersection_table_rows_selected
    ints <- upset_intersections()
    nm   <- names(ints)[row]
    peps <- ints[[nm]]
    
    # Compute how many length panels we'll need for dynamic height
    lengths_present <- sort(unique(nchar(peps)))
    lengths_present <- lengths_present[lengths_present >= 7 & lengths_present <= 15]
    n_rows_motif    <- max(1, ceiling(length(lengths_present) / 3))
    plot_height     <- paste0(n_rows_motif * 220, "px")
    
    showModal(modalDialog(
      title      = paste0("Intersection: ", nm, "  (n = ", length(peps), ")"),
      size       = "xl",
      easyClose  = TRUE,
      
      h5("Sequence Motifs"),
      plotOutput("upset_modal_motif", height = plot_height),
      tags$hr(),
      downloadButton("upset_modal_download", "Download peptide list (.csv)"),
      footer = modalButton("Close")
    ))
    
    # Motif plot — one panel per length
    output$upset_modal_motif <- renderPlot({
      if (length(lengths_present) == 0) {
        return(
          ggplot() +
            annotate("text", x = .5, y = .5,
                     label = "No peptides with length 7–15 in this intersection.",
                     size = 5, color = "grey40") +
            theme_void()
        )
      }
      
      plots <- lapply(lengths_present, function(L) {
        p <- tryCatch(
          plot_seqlogo(peps[nchar(peps) == L], title = paste("Length", L)),
          error = function(e) NULL
        )
        if (is.null(p)) {
          ggplot() +
            annotate("text", x = .5, y = .5,
                     label = paste("Length", L, ": too few peptides (< 5)"),
                     size = 4, color = "grey50") +
            theme_void()
        } else {
          p
        }
      })
      
      patchwork::wrap_plots(plots, ncol = min(3, length(plots)))
    })
    
    # Download handler
    output$upset_modal_download <- downloadHandler(
      filename = function() {
        safe_nm <- gsub(" & ", "_", nm)
        safe_nm <- gsub("[^A-Za-z0-9_]", "", safe_nm)
        paste0("intersection_", safe_nm, ".csv")
      },
      content = function(file) {
        writeLines(peps, file)
      }
    )
  })
  
  ## ----Pairwise comparison of shared peptides----
  output$Pairwise_shared_peptide_matrix <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples to plot"))
    req(input$shared_mode)
    
    ht <- plot_shared_peptide(
      lst,
      mode  = input$shared_mode,
      color = input$color_palette,
      percent_type = "union"
    )
    
    safe_draw(ht)
  })
  
  ## ----Pairwise comparison of shared peptides quantity----
  output$pairwise_peptide_quant_correlation <- renderPlot({
    lst <- processed_data_list()
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples to plot"))
    check_data_error(lst, required_cols = "MAX_QUANTITY" , na_policy = "any")
    ht <- plot_pairwise_peptide_quant_correlation(lst, color = input$color_palette, cluster = input$cluster_mode)
    safe_draw(ht)
  })
  
  ## ----PCA plot----
  output$pca <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more sets to plot Upset"))
    plot_PCA(lst, color = input$color_palette)
  })
  
  ## ----Number of peptides and peptidoforms----
  output$summary_peptides_plot2 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    plot_unique_counts(lst, column = "STRIPPED", y_label = "Number of unique peptides", color = input$color_palette)
  })
  
  output$summary_peptidoforms_plot2 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    plot_unique_counts(lst, column = "PEPTIDE", y_label = "Number of unique peptidoforms", color = input$color_palette)
  })
  
  output$summary_proteins_plot2 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PROTEIN" , na_policy = "all")
    # Use extract_protein_prefixes for proteins
    plot_unique_counts(lst, column = "PROTEIN", y_label = "Number of unique proteins",
                       transform_fn = extract_protein_prefixes, color = input$color_palette)
  })
  
  ## ----Charge / Mass / mz / RT / ppm----
  output$charge_plot2 <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "CHARGE" , na_policy = "all")
    plot_stacked_bar(lst, column = "CHARGE", fill_label = "Charge", percentage = FALSE, color = input$color_palette)
  })
  
  # Example usage for your Shiny outputs
  output$mass_plot2 <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MASS" , na_policy = "all")
    plot_density(lst, column = "MASS", x_label = "Mass (Da)", color = input$color_palette)
  })
  
  output$mz_plot2 <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MZ" , na_policy = "all")
    plot_density(lst, column = "MZ", x_label = "m/z", color = input$color_palette)
  })
  
  output$RT_plot2 <- renderUI({
    lst <- processed_data_list()
    
    # Wrap plots in a grid (like motif plots)
    layout_column_wrap(
      width = "400px",  # each plot approx width
      !!!lapply(names(lst), function(sample_name) {
        plotOutput(paste0("RT2_", sample_name), height = "300px")
      })
    )
  })
  
  observe({
    lst <- processed_data_list()
    
    for (sample_name in names(lst)) {
      
      local({
        sample_local <- sample_name
        df_local <- lst[[sample_local]]
        
        output[[paste0("RT2_", sample_local)]] <- renderPlot({
          check_data_error(df_local, required_cols = "RT", na_policy = "any")
          plot_histogram(df = df_local, column = "RT", x_label = "Retention Time", 
                         title_name = paste("RT Histogram –", sample_local), color = input$color_palette
          )
        })
      })
    }
  })
  
  output$ppm_plot2 <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PPM" , na_policy = "all")
    plot_density(lst, column = "PPM", x_label = "ppm", color = input$color_palette)
  })
  
  ##----Score distribution----
  output$score_violin2 <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "SCORE" , na_policy = "all")
    plot_violin(lst, column = "SCORE", color = input$color_palette)
  })
  
#-------------------Group Comparison------------------------
  ## ----Render the UI for group assignment----
  output$group_assign_ui <- renderUI({
    req(names(active_data_list()))
    has_conditions <- annotation_provided() &&
      length(attr(annotation_df_r(), "condition_cols")) > 0
    
    tagList(
      tabsetPanel(
        id = "group_mode_tabs",
        
        tabPanel(
          title = "Manual (Samples)",
          value = "manual",
          tags$br(),
          uiOutput("manual_group_selector_ui")
        ),
        
        tabPanel(
          title = tagList(
            "Condition-based",
            if (!has_conditions)
              tags$small(" (no annotation/condition)", style = "color:#aaa; font-weight:normal;")
          ),
          value = "condition",
          tags$br(),
          if (!has_conditions) {
            tags$p(tags$em("Load an annotation table with condition columns to enable this feature."),
                   style = "color:grey;")
          } else {
            uiOutput("condition_group_builder_ui")
          }
        )
      ),
      
      # Grey out the condition tab when unavailable
      if (!has_conditions) {
        tags$script(HTML("
        setTimeout(function() {
          $('#group_mode_tabs a[data-value=\"condition\"]').css({
            'pointer-events': 'none',
            'color': '#aaa',
            'opacity': '0.55',
            'cursor': 'not-allowed'
          });
        }, 100);
      "))
      }
    )
  })
  
  #Manually select samples to group
  output$manual_group_selector_ui <- renderUI({
    req(names(active_data_list()))
    n <- input$n_groups
    tagList(lapply(seq_len(n), function(i) {
      tagList(
        textInput(paste0("group_name_", i), paste("Group", i, "name"), paste("Group", i)),
        selectInput(paste0("group_", i), paste("Select samples for Group", i),
                    choices = names(active_data_list()), multiple = TRUE),
        tags$hr()
      )
    }))
  })
  
  #build groups by condition
  output$condition_group_builder_ui <- renderUI({
    req(annotation_provided())
    cond_cols <- attr(annotation_df_r(), "condition_cols")
    req(length(cond_cols) > 0)
    
    tagList(
      # ── Expression builder ──────────────────────────────────────────
      tags$div(class = "well",
               tags$h6(tags$strong("Build expression")),
               fluidRow(
                 column(4, selectInput("cond_col", "Condition column:", choices = cond_cols)),
                 column(4, selectInput("cond_val", "Values:", choices = character(0), multiple = TRUE)),
                 column(4, selectInput("cond_op",  "Operator (2nd+ term):", choices = c("AND", "OR", "NOT")))
               ),
               actionButton("cond_add_term",   "Add term",  class = "btn-sm btn-primary mr-1"),
               actionButton("cond_undo_term",  "Undo last", class = "btn-sm btn-warning mr-1"),
               actionButton("cond_clear_expr", "Clear all", class = "btn-sm btn-danger"),
               tags$hr(),
               tags$h6(tags$strong("Current expression:")),
               verbatimTextOutput("cond_expr_text"),
               tags$h6(tags$strong("Sample mask preview:")),
               tableOutput("cond_mask_table")
      ),
      # ── Save / update group ─────────────────────────────────────────
      tags$div(class = "well",
               uiOutput("cond_edit_banner"),
               fluidRow(
                 column(8, textInput("cond_group_name", "Group name:", "")),
                 column(4, tags$br(), uiOutput("cond_save_btn_ui"))
               )
      ),
      # ── Saved groups ────────────────────────────────────────────────
      uiOutput("cond_saved_groups_ui"),
      tableOutput("cond_membership_table")
    )
  })
  
  observe({
    updateSelectizeInput(session, "HLA_alleles",
                         choices  = unname(unlist(hla_alleles)),
                         selected = "HLA-A02:01",
                         server   = TRUE)
  })
  
  # Populate value choices when column selection changes
  observeEvent(input$cond_col, {
    req(annotation_provided(), input$cond_col)
    vals <- sort(unique(as.character(annotation_df_r()[[input$cond_col]])))
    vals <- vals[!is.na(vals)]
    updateSelectInput(session, "cond_val", choices = vals, selected = character(0))
  })
  
  # Add a term to the expression
  observeEvent(input$cond_add_term, {
    req(input$cond_col, length(input$cond_val) > 0)
    expr <- active_expr_r()
    op   <- if (length(expr) == 0) NULL else input$cond_op
    expr[[length(expr) + 1]] <- list(col = input$cond_col, val = input$cond_val, op = op)
    active_expr_r(expr)
  })
  
  # Undo last term
  observeEvent(input$cond_undo_term, {
    expr <- active_expr_r()
    if (length(expr) > 0) active_expr_r(expr[-length(expr)])
  })
  
  # Clear expression
  observeEvent(input$cond_clear_expr, {
    active_expr_r(list())
    editing_group_r(NULL)
    updateTextInput(session, "cond_group_name", value = "")
  })
  
  # Cancel edit
  observeEvent(input$cond_cancel_edit, {
    active_expr_r(list())
    editing_group_r(NULL)
    updateTextInput(session, "cond_group_name", value = "")
  })
  
  # Expression text preview
  output$cond_expr_text <- renderText({
    expr <- active_expr_r()
    if (length(expr) == 0) return("(empty)")
    parts <- sapply(seq_along(expr), function(i) {
      t        <- expr[[i]]
      term_str <- paste0(t$col, " IN [", paste(t$val, collapse = ", "), "]")
      if (i == 1) term_str else paste(t$op, term_str)
    })
    paste(parts, collapse = " ")
  })
  
  # Sample mask preview
  output$cond_mask_table <- renderTable({
    expr <- active_expr_r()
    req(length(expr) > 0, annotation_provided())
    result_df        <- eval_condition_expr(expr, annotation_df_r())
    result_df$Match  <- ifelse(result_df$result, "Yes", "No")
    result_df$result <- NULL
    result_df
  }, striped = TRUE, bordered = TRUE, hover = TRUE)
  
  # Edit mode banner
  output$cond_edit_banner <- renderUI({
    nm <- editing_group_r()
    if (is.null(nm)) return(NULL)
    tags$div(class = "alert alert-info py-1 mb-2",
             tags$strong("Editing: "), nm,
             actionButton("cond_cancel_edit", "Cancel", class = "btn-sm btn-secondary ml-2")
    )
  })
  
  # Save/Update button label
  output$cond_save_btn_ui <- renderUI({
    if (is.null(editing_group_r())) {
      actionButton("cond_save_group", "Save as group", class = "btn-success btn-sm")
    } else {
      actionButton("cond_save_group", "Update group",  class = "btn-primary btn-sm")
    }
  })
  
  # Save or update group
  observeEvent(input$cond_save_group, {
    nm <- trimws(input$cond_group_name)
    req(nchar(nm) > 0)
    expr <- active_expr_r()
    req(length(expr) > 0, annotation_provided())
    
    groups <- condition_groups_r()
    
    if (is.null(groups[[nm]]) && length(groups) >= input$n_groups) {
      showNotification(paste0("Maximum of ", input$n_groups, " group(s) reached."), type = "warning")
      return()
    }
    
    result_df <- eval_condition_expr(expr, annotation_df_r())
    selected  <- if (isTRUE(attr(annotation_df_r(), "has_measurement"))) {
      list(name = result_df$name[result_df$result], measurement = result_df$measurement[result_df$result], expr = expr)
    } else {
      list(name = result_df$name[result_df$result], expr = expr)
    }
    
    if (length(selected$name) == 0) {
      showNotification("Expression matches no samples.", type = "warning")
      return()
    }
    
    groups[[nm]] <- selected
    condition_groups_r(groups)
    
    active_expr_r(list())
    editing_group_r(NULL)
    updateTextInput(session, "cond_group_name", value = "")
    showNotification(paste0("Group '", nm, "' saved (", length(selected$name), " samples)."), type = "message")
  })
  
  # Register per-group edit/delete observers dynamically
  observe({
    grps <- condition_groups_r()
    lapply(names(grps), function(nm) {
      local({
        local_nm <- nm
        observeEvent(input[[paste0("cond_delete_", local_nm)]], {
          g <- condition_groups_r()
          if (!is.null(g[[local_nm]])) {
            g[[local_nm]] <- NULL
            condition_groups_r(g)
            if (identical(editing_group_r(), local_nm)) {
              editing_group_r(NULL)
              active_expr_r(list())
              updateTextInput(session, "cond_group_name", value = "")
            }
          }
        }, ignoreInit = TRUE)
        
        observeEvent(input[[paste0("cond_edit_", local_nm)]], {
          g <- condition_groups_r()
          if (!is.null(g[[local_nm]])) {
            active_expr_r(g[[local_nm]]$expr)
            editing_group_r(local_nm)
            updateTextInput(session, "cond_group_name", value = local_nm)
          }
        }, ignoreInit = TRUE)
      })
    })
  })
  
  # Saved groups cards
  output$cond_saved_groups_ui <- renderUI({
    groups <- condition_groups_r()
    if (length(groups) == 0) return(tags$p(tags$em("No groups saved yet.")))
    
    cards <- lapply(names(groups), function(nm) {
      g          <- groups[[nm]]
      is_editing <- identical(editing_group_r(), nm)
      tags$div(
        class = paste("card mb-2", if (is_editing) "border-primary" else ""),
        tags$div(class = "card-body py-2",
                 tags$div(class = "d-flex justify-content-between align-items-center",
                          tags$strong(nm),
                          tags$span(class = "badge badge-secondary mr-auto ml-2",
                                    paste(length(g$names), "samples")),
                          tags$div(
                            actionButton(paste0("cond_edit_",   nm), icon("pencil-alt"),
                                         class = "btn-sm btn-outline-primary mr-1",   title = "Edit"),
                            actionButton(paste0("cond_delete_", nm), icon("trash"),
                                         class = "btn-sm btn-outline-danger",          title = "Delete")
                          )
                 )
        )
      )
    })
    
    tagList(tags$h6(tags$strong("Saved groups:")), tagList(cards))
  })
  
  # Sample membership summary table
  output$cond_membership_table <- renderTable({
    groups <- condition_groups_r()
    req(length(groups) > 0)
    
    ann <- annotation_df_r()
    colnames(ann) <- tolower(colnames(ann))
    has_meas <- isTRUE(attr(annotation_df_r(), "has_measurement"))
    
    if (has_meas) {
      mem <- ann[, c("name", "measurement")]
    } else {
      mem <- data.frame(name = names(active_data_list()), 
                        measurement = NA_character_, 
                        stringsAsFactors = FALSE)
    }
    
    for (nm in names(groups)) {
      grp <- groups[[nm]]
      if (has_meas) {
        mem[[nm]] <- mem$name %in% grp$name & mem$measurement %in% grp$measurement
      } else {
        mem[[nm]] <- mem$name %in% grp$name
      }
      mem[[nm]] <- ifelse(mem[[nm]], "X", "")
    }
    
    # Compute Status
    n_in <- rowSums(mem[, names(groups), drop = FALSE] == "X")
    mem$Status <- ifelse(n_in == 0, "excluded", ifelse(n_in > 1, "overlap", ""))
    
    mem
  }, striped = TRUE, bordered = TRUE, hover = TRUE)
  
  
  # Update group list
  group_list <- eventReactive(input$update_group_comp, {
    
    # Condition-based mode
    if (isTRUE(input$group_mode_tabs == "condition")) {
      groups_raw <- condition_groups_r()
      use_measurements(TRUE)
      if (length(groups_raw) < 2) {
        showNotification("Need at least 2 saved condition groups before updating.", type = "error")
        return(NULL)
      }
      showNotification("Groups updated successfully!", type = "message")
      has_meas <- isTRUE(attr(annotation_df_r(), "has_measurement"))
      if (has_meas) {
        return(lapply(groups_raw, function(x) {
          list(
            measurement = x$measurement,
            name = x$name
          )
        }))
      } else {
        use_measurements(FALSE)
        return(lapply(groups_raw, `[[`, "name"))
      }
    }
    
    # Manual mode (unchanged)
    n            <- input$n_groups
    groups       <- list()
    empty_groups <- c()
    
    use_measurements(FALSE)
    
    for (i in seq_len(n)) {
      custom_name <- input[[paste0("group_name_", i)]]
      samples     <- input[[paste0("group_", i)]]
      if (is.null(custom_name) || nchar(custom_name) == 0 ||
          is.null(samples)     || length(samples) == 0) {
        empty_groups <- c(empty_groups, paste("Group", i))
      } else {
        groups[[custom_name]] <- samples
      }
    }
    
    if (length(empty_groups) > 0) {
      showNotification(
        paste("The following group(s) are empty or have no name:",
              paste(empty_groups, collapse = ", ")),
        type = "error"
      )
      return(NULL)
    }
    
    showNotification("Groups updated successfully!", type = "message")
    groups
    
  }, ignoreNULL = TRUE)
  
  
  
  ## ----Group Comparison----
  group_peptide_sets <- reactive({
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    groups <- group_list()
    req(groups, length(groups) >= 2)
    col_map <- measurement_col_map_r()
    
    pep_col <- "STRIPPED"
    
    lapply(names(groups), function(g) {
      grp <- groups[[g]]
      
      if (!is.null(col_map) && is.list(grp) && !is.null(grp$measurement)) {
        # Measurement mode: grp = list(measurement = c(...), name = c(...))
        # A peptide counts for this measurement only if it has a non-NA quantity value
        measurements <- grp$measurement
        names_vec    <- grp$name
        
        peptide_counts <- table(unlist(lapply(seq_along(measurements), function(i) {
          entry <- col_map[[measurements[i]]]
          df <- lst[[names_vec[i]]]
          actual_col <- colnames(df)[tolower(colnames(df)) == tolower(entry$col)][1]
          df[[pep_col]][!is.na(df[[actual_col]])]
        })))
        
        n_items <- length(measurements)
      } else {
        # Sample mode: grp is a character vector of sample names
        peptide_counts <- table(unlist(lapply(grp, function(s) {
          df <- lst[[s]]
          if (is.null(df) || !pep_col %in% colnames(df)) return(character(0))
          unique(df[[pep_col]])
        })))
        
        n_items <- length(grp)
      }
      
      names(peptide_counts[
        peptide_counts / n_items >= input$min_presence_fraction
      ])
    }) |> setNames(names(groups))
  })
  
  output$group_venn_plot <- renderPlot({
    sets <- group_peptide_sets()
    shiny::validate(shiny::need(length(sets) >= 2, "Need 2 or more sets to compare"))
    
    group_venn_plotting(sets, color = input$color_palette)
  })
  
  heatmap_data_r  <- reactive({
    lst       <- processed_data_list()
    groups    <- group_list()
    data_info <- data_info_r()
    req(lst, groups, data_info, length(groups) > 0)
    quantity_cols <- data_info %>%
      dplyr::filter(final_name == "QUANTITY") %>%
      dplyr::select(-final_name) %>%
      unlist(recursive = TRUE, use.names = FALSE)
    shiny::validate(shiny::need(length(quantity_cols) > 0, "No QUANTITY columns found."))
    col_map_hm <- measurement_col_map_r()
    pep_mat <- prepare_peptide_matrix(lst, groups, quantity_cols, group_peptide_sets(),
                                      col_map = col_map_hm)
    shiny::validate(shiny::need(nrow(pep_mat) > 0, "No peptides to plot."))
    mat_full <- log10(pep_mat + 1)
    mat_full <- t(mat_full)
    clust_na0 <- function(x) { x2 <- x; x2[!is.finite(x2)] <- 0; dist(x2) }
    bin_map <- NULL; mat_display <- mat_full; col_ord <- seq_len(ncol(mat_full))
    if (ncol(mat_full) > 1000L) {
      mat_binned  <- bin_heatmap_columns(mat_full, max_cols = 1000L)
      bin_map     <- attr(mat_binned, "bin_map")
      mat_display <- mat_binned
      # bin_heatmap_columns already groups by presence pattern; no further column clustering
    } else {
      col_ord <- tryCatch(
        hclust(clust_na0(t(mat_full)), method = "complete")$order,
        error = function(e) seq_len(ncol(mat_full))
      )
      mat_display <- mat_full[, col_ord, drop = FALSE]
    }
    row_ord <- tryCatch(
      hclust(clust_na0(mat_display), method = "complete")$order,
      error = function(e) seq_len(nrow(mat_display))
    )
    list(mat_full = mat_full, mat_display = mat_display,
         bin_map = bin_map, row_ord = row_ord)
  })
  
  heatmap_mode_r  <- reactiveVal("binned")
  expanded_peps_r <- reactiveVal(NULL)
  zoom_info_r              <- reactiveVal(NULL)
  last_programmatic_t_r   <- reactiveVal(0)
  expand_counter_r <- reactiveVal(0L)
  
  observeEvent(plotly::event_data("plotly_relayout", source = "group_hm"), {
    ed <- plotly::event_data("plotly_relayout", source = "group_hm")
    if (!is.null(ed[["xaxis.range[0]"]])) {
      zoom_info_r(c(as.numeric(ed[["xaxis.range[0]"]]),
                    as.numeric(ed[["xaxis.range[1]"]])))
    } else if (isTRUE(ed[["xaxis.autorange"]])) {
      zoom_info_r(NULL)
      # Only treat as user double-click if it happened >1s after a programmatic render
      age <- as.numeric(Sys.time()) - isolate(last_programmatic_t_r())
      if (age > 1 && isolate(heatmap_mode_r()) == "expanded") {
        last_programmatic_t_r(as.numeric(Sys.time()))
        heatmap_mode_r("binned")
        expanded_peps_r(NULL)
      }
    }
  }, ignoreNULL = TRUE)
  
  zoom_debounced_r <- shiny::debounce(zoom_info_r, 500)
  
  # Tick-label proxy: only fires on zoom change, only in binned mode
  observeEvent(zoom_debounced_r(), {
    req(heatmap_mode_r() == "binned")
    zoom  <- zoom_debounced_r()
    hdata <- heatmap_data_r()
    req(hdata)
    
    n_visible <- if (!is.null(zoom)) {
      x0 <- max(1L, as.integer(round(zoom[1])) + 1L)
      x1 <- min(ncol(hdata$mat_display), as.integer(round(zoom[2])) + 1L)
      max(1L, x1 - x0 + 1L)
    } else {
      ncol(hdata$mat_display)
    }
    
    show <- n_visible <= 250
    plotly::plotlyProxy("group_peptide_heatmap_interactive", session) %>%
      plotly::plotlyProxyInvoke("relayout", list(
        "xaxis.showticklabels" = show,
        "margin.b"             = if (show) 120 else 30
      ))
  })
  
  # Dynamically show/hide x-axis tick labels based on visible column count
  # Uses plotlyProxy so no full re-render is triggered
  observe({
    zoom  <- zoom_debounced_r()
    hdata <- heatmap_data_r()
    mode  <- heatmap_mode_r()
    req(hdata)
    
    n_cols <- if (mode == "expanded" && !is.null(expanded_peps_r())) {
      length(expanded_peps_r())
    } else {
      ncol(hdata$mat_display)
    }
    
    n_visible <- if (!is.null(zoom)) {
      x0 <- max(1L, as.integer(round(zoom[1])) + 1L)
      x1 <- min(n_cols, as.integer(round(zoom[2])) + 1L)
      max(1L, x1 - x0 + 1L)
    } else {
      n_cols
    }
    
    show <- n_visible <= 250
    
    plotly::plotlyProxy("group_peptide_heatmap_interactive", session) %>%
      plotly::plotlyProxyInvoke("relayout", list(
        "xaxis.showticklabels" = show,
        "margin.b"             = if (show) 120 else 30
      ))
  })
  
  output$group_peptide_heatmap_interactive <- plotly::renderPlotly({
    hdata <- heatmap_data_r()
    color <- input$color_palette
    mode  <- heatmap_mode_r()
    
    if (mode == "expanded" && !is.null(expanded_peps_r())) {
      pep_names    <- expanded_peps_r()
      mat_to_show  <- hdata$mat_full[, pep_names, drop = FALSE]
      bin_map_show <- NULL
      ui_rev       <- paste0("exp_", expand_counter_r())
    } else {
      mat_to_show  <- hdata$mat_display
      bin_map_show <- hdata$bin_map
      ui_rev       <- "binned"
    }
    
    mat_to_show <- mat_to_show[hdata$row_ord, , drop = FALSE]
    
    hover_mat <- matrix("", nrow = nrow(mat_to_show), ncol = ncol(mat_to_show))
    for (j in seq_len(ncol(mat_to_show))) {
      cn <- colnames(mat_to_show)[j]
      if (!is.null(bin_map_show) && cn %in% names(bin_map_show)) {
        peps  <- bin_map_show[[cn]]
        n     <- length(peps)
        shown <- paste(head(peps, 30), collapse = "<br>")
        extra <- if (n > 30) paste0("<br><i>+", n - 30, " more</i>") else ""
        hover_mat[, j] <- paste0("<b>", n, " peptides in bin</b><br>", shown, extra)
      } else {
        hover_mat[, j] <- cn
      }
    }
    
    if (color == "default") {
      cscale <- list(list(0, "lightyellow"), list(1, "red"))
    } else {
      vcols  <- viridis::viridis(10, option = color)
      cscale <- lapply(seq_along(vcols) - 1,
                       function(i) list(i / (length(vcols) - 1), vcols[i + 1]))
    }
    
    show_ticks <- ncol(mat_to_show) <= 250
    
    plotly::plot_ly(
      source        = "group_hm",
      x             = colnames(mat_to_show),
      y             = rownames(mat_to_show),
      z             = mat_to_show,
      text          = hover_mat,
      type          = "heatmap",
      colorscale    = cscale,
      hovertemplate = "%{text}<extra></extra>",
      colorbar      = list(title = "log10(QUANTITY+1)")
    ) %>%
      plotly::layout(
        uirevision = ui_rev,
        xaxis = list(title = "", showticklabels = show_ticks,
                     tickfont = list(size = 8), tickangle = -45),
        yaxis = list(title = "", tickfont = list(size = 8), autorange = "reversed"),
        margin = list(l = 130, b = if (show_ticks) 120 else 30)
      ) %>%
      plotly::event_register("plotly_relayout")
  })
  
  # Manual expand: zoom in first, then click this button
  observeEvent(input$expand_heatmap_region, {
    zoom  <- zoom_info_r()
    hdata <- heatmap_data_r()
    req(hdata)
    
    if (is.null(hdata$bin_map)) {
      showNotification("Already showing individual peptides.", type = "message")
      return()
    }
    if (is.null(zoom)) {
      showNotification("Zoom into a region first, then click Expand.", type = "warning")
      return()
    }
    
    x0 <- max(1L, as.integer(round(zoom[1])) + 1L)
    x1 <- min(ncol(hdata$mat_display), as.integer(round(zoom[2])) + 1L)
    vis_bins  <- colnames(hdata$mat_display)[x0:x1]
    pep_names <- unlist(lapply(vis_bins, function(cn)
      if (cn %in% names(hdata$bin_map)) hdata$bin_map[[cn]] else cn))
    
    if (length(pep_names) > 2000) {
      showNotification(paste0(length(pep_names), " peptides in view — zoom in more."),
                       type = "warning", duration = 8)
      return()
    }
    
    expand_counter_r(isolate(expand_counter_r()) + 1L)   # always unique uirevision
    last_programmatic_t_r(as.numeric(Sys.time()))
    expanded_peps_r(pep_names)
    heatmap_mode_r("expanded")
    zoom_info_r(NULL)
  })
  
  observeEvent(input$reset_heatmap_view, {
    last_programmatic_t_r(as.numeric(Sys.time()))  # suppress re-render autorange
    heatmap_mode_r("binned")
    expanded_peps_r(NULL)
    zoom_info_r(NULL)
  })
  
  output$download_visible_peptides <- downloadHandler(
    filename = function() paste0("peptides_", Sys.Date(), ".txt"),
    content  = function(file) {
      hdata <- heatmap_data_r()
      zoom  <- zoom_info_r()
      mode  <- heatmap_mode_r()
      req(hdata)
      
      if (mode == "expanded" && !is.null(expanded_peps_r())) {
        peptides <- expanded_peps_r()
        if (!is.null(zoom)) {
          x0 <- max(1L, as.integer(round(zoom[1])) + 1L)
          x1 <- min(length(peptides), as.integer(round(zoom[2])) + 1L)
          if (x0 <= x1) peptides <- peptides[x0:x1]
        }
      } else {
        mat_cols <- colnames(hdata$mat_display)
        if (!is.null(zoom)) {
          x0 <- max(1L, as.integer(round(zoom[1])) + 1L)
          x1 <- min(length(mat_cols), as.integer(round(zoom[2])) + 1L)
          if (x0 <= x1) mat_cols <- mat_cols[x0:x1]
        }
        peptides <- unlist(lapply(mat_cols, function(cn)
          if (!is.null(hdata$bin_map) && cn %in% names(hdata$bin_map))
            hdata$bin_map[[cn]] else cn))
      }
      
      writeLines(unique(peptides), file)
    }
  )
  
  ## ----Group statistical analysis----
  group_comp_data <- eventReactive(input$update_group_comp, {

    lst <- processed_data_list()
    groups <- group_list()
    col_map <- measurement_col_map_r()
    
    groups <- groups[sapply(groups, function(g) {
      isTRUE(
        if (!is.null(col_map) && is.list(g) && !is.null(g$measurement)) {
          any(sapply(seq_along(g$measurement), function(i) {
            entry <- col_map[[g$measurement[i]]]
            !is.null(entry) && !is.null(lst[[g$name[i]]]) && nrow(lst[[g$name[i]]]) > 0
          }))
        } else {
          any(sapply(g, function(item) !is.null(lst[[item]]) && nrow(lst[[item]]) > 0))
        }
      )
    })]

    shiny::validate(shiny::need(length(groups) >= 2, "Need 2 or more sets to compare"))
    
    # Generate all unique pairwise combinations
    group_pairs <- combn(names(groups), 2, simplify = FALSE)
    pep_col <- "PEPTIDE" #WIP: need to think how to solve the PTM problem? Do I just add the same peptidoform together?
    # pep_col    <- if (!is.null(input$use_peptidoforms) && input$use_peptidoforms) "PEPTIDE" else "STRIPPED"
    
    # Compute volcano data for each pair
    group_comp_stats <- lapply(group_pairs, function(pair) {
      g1 <- pair[1]
      g2 <- pair[2]
      allowed_peptides_g1 <- group_peptide_sets()[[g1]]
      allowed_peptides_g2 <- group_peptide_sets()[[g2]]
      compute_group_comp_stats(lst, groups, allowed_peptides_g1, allowed_peptides_g2,
                               g1, g2, pep_col, default_quantity_cols_r(), col_map = col_map, use_measurements = use_measurements())
    })
    
    names(group_comp_stats) <- sapply(group_pairs, function(pair) paste(pair, collapse = "_vs_"))
    
    Filter(function(x) !is.null(x) && nrow(x) > 0, group_comp_stats)
  },ignoreNULL = TRUE)
  
  # Render the volcano tabset UI
  output$group_stats_tabs <- renderUI({
    volcano_list <- group_comp_data()
    req(volcano_list)
    shiny::validate(shiny::need(length(volcano_list) > 0, "Not enough data to compute any group comparisons."))
    
    output$volcano_comparison_tabs <- renderUI({
      
      # Build a list of tabPanels for each comparison
      comparison_tabs <- lapply(names(volcano_list), function(name) {
        tabPanel(
          title = name,
          fluidRow(
            column(6, plotlyOutput(paste0("volcano_", name))),
            column(6, plotlyOutput(paste0("ma_", name)))
          ),
          fluidRow(
            column(6, plotlyOutput(paste0("pval_hist_", name))),
            column(6, plotlyOutput(paste0("rank_fc_", name)))
          ),
          fluidRow(
            column(12, DT::DTOutput(paste0("peptide_table_", name)))
          ),
          fluidRow(
            column(6, h6("GO-term enrichment"), 
                   helpText("Only UniProt accession IDs (e.g. P04439) are supported. Other formats will be skipped."),
                   plotOutput(paste0("go_term_", name))),
            column(6, h6("STRING-DB network"), 
                   helpText("Only UniProt accession IDs (e.g. P04439) are supported. Other formats will be skipped."),
                   plotOutput(paste0("STRING_", name)))
          )
        )
      })
      
      # Generate the tabsetPanel from the list
      do.call(tabsetPanel, c(id = "volcano_comparison_tabs", comparison_tabs))
    })
    uiOutput("volcano_comparison_tabs")
  })
  
  # Render each volcano plot (optimized)
  observe({
    volcano_list <- group_comp_data()
    req(volcano_list)
    
    if (length(volcano_list) == 0) return()
    
    for (name in names(volcano_list)) {
      local({
        plot_name <- name
        df <- volcano_list[[plot_name]]
        
        #Significance calculation
        df$Significance <- ifelse(
          is.na(df$negLog10AdjP_BH) | is.na(df$log2FC), "Missing",
          ifelse(df$negLog10AdjP_BH > 1.3 & abs(df$log2FC) > 1,
                 "Significant", "Not significant")
        )
        ## ----Volcano plot----

        output[[paste0("volcano_", plot_name)]] <- renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          sel <- selected_peptide()
          group_volcano_plot(df,sel,plot_name)
        })
        
        ## ----MA plot----
        output[[paste0("ma_", plot_name)]] <- renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          sel <- selected_peptide()
          group_MA_plot(df, sel,plot_name)
        })
        
        ## ----P-value histogram----
        output[[paste0("pval_hist_", plot_name)]] <- renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          group_p_histogram(df,plot_name)
        })
        
        ## ----Ranked Fold Change----
        output[[paste0("rank_fc_", plot_name)]] <- renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          sel <- selected_peptide()
          group_rank_FC(df,plot_name,sel)
        })
        ## ----Peptide Fold Change table----
        output[[paste0("peptide_table_", plot_name)]] <- DT::renderDT({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          group_FC_table(df)
        })
        
        ## ----GO term----
        output[[paste0("go_term_", plot_name)]] <- renderPlot({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) != 0, "No data available."))
          run_go_enrichment(df)
        })
        
        ## ----STRING-DB----
        output[[paste0("STRING_", plot_name)]] <- renderPlot({
          shiny::validate(shiny::need(curl::has_internet(), "No Internet Connection."))
          shiny::validate(shiny::need(!is.null(df) && nrow(df) != 0, "No data available."))
          run_string(df)
        })
      })
    }
  })
  
  selected_peptide <- reactive({
    ed <- plotly::event_data("plotly_click", source = "group_diff")
    if (is.null(ed) || is.null(ed$key)) return(NULL)
    ed$key
  })
  
#-------------------PTM------------------------
  output$PTM_plot <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PTM" , na_policy = "ignore")
    
    lst <- lapply(lst, function(df) {
      
      df %>%
        dplyr::mutate(
          PTM = ifelse(is.na(PTM) | PTM == "", "Unmodified", PTM)
        ) %>%
        tidyr::separate_rows(PTM, sep = "\\s*[,;]\\s*") %>%
        dplyr::mutate(
          PTM = gsub("^\\d+", "", PTM) #This is specifically for fragpipe to remove the position information on PTMs
        )
    })
    
    lst <- lst[!vapply(lst, is.null, logical(1))]
    
    plot_stacked_bar(lst, column = "PTM", fill_label = "PTM", percentage = FALSE, rev_levels = FALSE)
  })
  
  ## Mapping of PTMs
  output$PTM_mod_map_table <- DT::renderDataTable({
    mod_map_list <- data_mod_map()   # reactive holding list of mod_maps per sample
    req(mod_map_list)
    
    if (is.null(mod_map_list) || length(mod_map_list) == 0 || all(is.na(mod_map_list))) {
      return(DT::datatable(data.frame(Message = "no PTM_pseudo sequence")))
    }
    
    # Flatten the list into a single data.frame
    mod_map_list <- purrr::imap_dfr(mod_map_list, ~ tibble(
      Sample = .y,
      Symbol = .x,
      Modification = names(.x)
    ))
    
    DT::datatable(mod_map_list, options = list(pageLength = 20))
  })
  
  ## PTM Psuedo Motif plots
  output$PTM_motif_tabs <- renderUI({
    lst <- processed_data_list()
    req(lst)
    
    mod_map_list <- data_mod_map()
    
    if (is.null(mod_map_list) || length(mod_map_list) == 0 || all(is.na(mod_map_list))) {
      return(tags$div("No PTM_pseudo sequence"))
    }
    
    length_tabs <- lapply(motif_plot_length, function(L) {
      
      sample_plots <- lapply(names(lst), function(sample_name) {
        plotOutput(
          paste0("PTM_motif_", sample_name, "_", L),
          height = "180px"
        )
      })
      
      tabPanel(
        paste("Length", L),
        do.call(layout_column_wrap, c(list(width = "250px"), sample_plots))
      )
    })
    
    do.call(tabsetPanel, c(list(id = "ptm_motif_length_tabs"), length_tabs))
  })
  
  observe({
    lst <- processed_data_list()
    req(lst)
    
    #Create custom namespace for ggseqplot
    AA_symbols <- LETTERS            # standard amino acids
    PTM_symbols <- c(as.character(1:9), letters)
    namespace <- c(AA_symbols, PTM_symbols)
    
    for (L in motif_plot_length) {
      for (sample_name in names(lst)) {
        
        local({
          length_val <- L
          sample_val <- sample_name
          
          output[[paste0("PTM_motif_", sample_val, "_", length_val)]] <- renderPlot({
            
            df <- lst[[sample_val]]
            req(df)
            
            # subset by length
            df_L <- df[df$LENGTH == length_val, ]
            df_L <- df_L[df_L$PTM != "", ]
            
            peptides <- unique(df_L[["PTM_Pseudo"]])
            peptides <- peptides[nchar(peptides) == length_val]
            
            if (length(peptides) < 5) {
              return(ggplot() + 
                       annotate("text", x = 0.5, y = 0.5,
                                label = "Not enough peptides",
                                size = 6, color = "red") +
                       theme_void())
            }
            
            par(mar = c(1.5, 1.5, 2, 0.5))
            
            plot_seqlogo(
              peptides,
              title = paste(sample_val, "• Length", length_val),
              namespace = namespace
            )
          })
        })
      }
    }
  })
  
#-------------------Binding prediction------------------------
  observeEvent(input$run_netmhc, {
    if (is.null(input$HLA_alleles) || length(input$HLA_alleles) == 0) {
      showNotification("Please select at least one HLA allele.", type = "error")
      return()
    }
    
    lst <- processed_data_list()
    cache <- prediction_cache()
    check_data_error(lst, required_cols = "LENGTH" ,na_policy = "any")
    req(cache)
    req(input$HLA_alleles)
    
    alleles_vec <- input$HLA_alleles
    
    # Extract peptides of correct length
    peptides <- unique(unlist(lapply(lst, `[[`, "STRIPPED"), use.names = FALSE))
    peptides <- peptides[nchar(peptides) %in% 8:11]
    
    shiny::validate(shiny::need(length(peptides) > 0, "No Peptides within the correct length."))
    
    # Path to netMHCpan
    withProgress(message = "Running netMHCpan predictions...", value = 0, {
      for (al in alleles_vec) {
        
        al_conversion <- sub("-", "\\.", al)
        al_conversion <- sub(":", "", al_conversion)
        
        # Determine which peptides need prediction
        if (!(al_conversion %in% colnames(cache)[-1])) {
          peptides_to_predict <- peptides       # new allele → predict all peptides
        } else {
          peptides_to_predict <- cache$Peptide[is.na(cache[[al_conversion]])]  # existing allele → only new peptides
          peptides_to_predict <- unique(peptides_to_predict[nchar(peptides_to_predict) %in% 8:11])
        }

        incProgress(1 / length(alleles_vec), detail = paste("Predicting for allele", al, " (# of peptides: ", length(peptides_to_predict), ")" ))
        
        if (length(peptides_to_predict) == 0) next
        
        res <- tryCatch({
          out <- run_netmhcpan(peptides_to_predict, al, netmhcpan_path)
          parse_netmhc_output(out)
        }, error = function(e) {
          showNotification(paste("netMHCpan failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        })
        if (is.null(res)) next
        
        if (!(al_conversion %in% colnames(cache)[-1])) {
          cache <- left_join(cache, res, by = "Peptide")
        } else {
          #avoid one-to-many.
          res   <- distinct(res,   Peptide, .keep_all = TRUE)
          cache <- distinct(cache, Peptide, .keep_all = TRUE)
          
          cache <- full_join(cache, res, by = "Peptide") %>%
            dplyr::mutate(
              !!al_conversion := coalesce(.data[[paste0(al_conversion, ".x")]],
                               .data[[paste0(al_conversion, ".y")]])
            ) %>%
            dplyr::select(-all_of(c(paste0(al_conversion, ".x"), paste0(al_conversion, ".y"))))
        }

      }
    })
    
    # Update reactive cache once at the end
    prediction_cache(cache)
  })
  
  output$allele_viz_selector_ui <- renderUI({
    df <- peptide_wide_unique()
    req(!is.null(df), ncol(df) > 1)
    allele_cols <- grep("^HLA", colnames(df), value = TRUE)
    req(length(allele_cols) > 0)
    selectInput(
      "allele_viz_select",
      tagList(icon("filter"), "Alleles to visualize:"),
      choices  = allele_cols,
      selected = allele_cols,   # all selected by default
      multiple = TRUE
    )
  })
  
  peptide_wide_all <- reactive({
      lst <- processed_data_list()
      cache <- prediction_cache()
      
      req(lst)
      
      v <- safe_validate(ncol(cache) > 1, "No prediction data")
      if (!is.null(v)) return(v)
      
      
      lst <- lapply(lst, function(df) {
        
        # Only keep necessary columns
        df <- distinct(df[, c("STRIPPED", "LENGTH")])
        
        # Left join with new cache
        df <- dplyr::left_join(df, cache, by = c("STRIPPED" = "Peptide"))

        df
      })
      
      out <- lapply(names(lst), function(nm) {
        df <- lst[[nm]]
        allele_cols <- grep("^HLA", colnames(df), value = TRUE)
        df <- df %>%
          dplyr::filter(LENGTH >= 8, LENGTH <= 11) %>%
          dplyr::select(STRIPPED, all_of(allele_cols)) %>%
          dplyr::mutate(Set = nm)
        
        df
      })
      
      dplyr::bind_rows(out)
  })
  
  peptide_wide_unique <- reactive({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    
    df <- peptide_wide_all()
    
    allele_cols <- grep("^HLA", colnames(df), value = TRUE)
    
    df %>%
      dplyr::group_by(Set, STRIPPED) %>%
      dplyr::summarise(dplyr::across(all_of(allele_cols), ~ if(all(is.na(.x))) {NA_real_} else min(.x, na.rm = TRUE)),
                       .groups = "drop")
  })
  
  binder_summary_all <- reactive({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    req(!is.null(peptide_wide_unique()))
    compute_binder_summary(peptide_wide_unique(), alleles = input$allele_viz_select)
  })
  
  output$binding_summary <- DT::renderDT({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    req(!is.null(peptide_wide_unique()))
    binder_summary_all()
  })
  
  output$binding_plot_percent <- renderPlot({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    req(!is.null(peptide_wide_unique()))
    plot_binders(peptide_wide_unique(), color = input$color_palette, percent = TRUE, alleles = input$allele_viz_select)
  })
  
  output$binding_plot_absolute <- renderPlot({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    req(!is.null(peptide_wide_unique()))
    plot_binders(peptide_wide_unique(), color = input$color_palette, percent = FALSE, alleles = input$allele_viz_select)
  })
  
  output$binding_table <- DT::renderDT({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    
    selected_alleles <- input$allele_viz_select
    
    df_collapsed <- peptide_wide_unique() %>%
      dplyr::select(-dplyr::any_of(
        setdiff(grep("^HLA", colnames(.), value = TRUE), selected_alleles)
      )) %>%
      dplyr::group_by(STRIPPED) %>%
      dplyr::summarise_all(~ {
        vals <- unique(.)
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) NA else paste(vals, collapse = ",")
      })
    
    DT::datatable(
      df_collapsed,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollX = TRUE
      )
    )
  })
  
  # All peptides
  output$dl_all <- downloadHandler(
    filename = function() { paste0("peptides_all_", Sys.Date(), ".csv") },
    content  = function(file) {
      df <- peptide_wide_unique()
      write.csv(df, file, row.names = FALSE, quote = FALSE)
    }
  )
  
  # Only binders (Strong or Weak)
  output$dl_binders <- downloadHandler(
    filename = function() { paste0("peptides_binders_", Sys.Date(), ".csv") },
    content = function(file) {
      df <- peptide_wide_unique()
      allele_cols <- grep("^HLA", colnames(df), value = TRUE)
      
      # Keep rows where at least one allele is <= 2 (Strong or Weak)
      binders <- df[rowSums(df[allele_cols] <= 2, na.rm = TRUE) > 0, ]
      write.csv(binders, file, row.names = FALSE, quote = FALSE)
    }
  )
  
  # Only non-binders (>2)
  output$dl_nonbinders <- downloadHandler(
    filename = function() { paste0("peptides_nonbinders_", Sys.Date(), ".csv") },
    content = function(file) {
      df <- peptide_wide_unique()
      allele_cols <- grep("^HLA", colnames(df), value = TRUE)
      
      # Keep rows where all non-NA alleles > 2
      nonbinders <- df[apply(df[allele_cols], 1, function(x) {
        non_na <- x[!is.na(x)]
        length(non_na) > 0 && all(non_na > 2)
      }), ]
      write.csv(nonbinders, file, row.names = FALSE, quote = FALSE)
    }
  )
  
  # Only missing (all NA)
  output$dl_missing <- downloadHandler(
    filename = function() { paste0("peptides_missing_", Sys.Date(), ".csv") },
    content = function(file) {
      df <- peptide_wide_unique()
      allele_cols <- grep("^HLA", colnames(df), value = TRUE)
      
      # Keep rows where all alleles are NA
      missing <- df[rowSums(!is.na(df[allele_cols])) == 0, ]
      write.csv(missing, file, row.names = FALSE, quote = FALSE)
    }
  )
  
#-------------------Peptide Search-------------
  output$peptide_table <- renderDT({
    lst <- data_list_r()
    
    keep_cols <- c("Sample", "PEPTIDE", "STRIPPED", "MAX_QUANTITY","LENGTH" ,"MASS" ,"CHARGE", "MZ" ,"K0", "RT", "PROTEIN")
    
    combined <- imap_dfr(lst, ~ dplyr::select(.x, dplyr::intersect(keep_cols, colnames(.x))) %>% dplyr::mutate(Sample = .y))
    
    peptide_counts <- combined %>% dplyr::distinct(STRIPPED, Sample) %>% dplyr::count(STRIPPED, name = "n_datasets")
    
    numeric_cols <- c("LENGTH", "MASS", "CHARGE", "MZ", "K0", "RT")
    
    combined <- combined %>% dplyr::inner_join(peptide_counts, by = "STRIPPED") %>%
      dplyr::arrange(Sample, STRIPPED) %>%
      dplyr::mutate(across(any_of(numeric_cols), as.numeric))
    
    combined <- combined[,c(1,2,3,4,5,6,7,8,9,11,12,10)] #rearramge column
    
    shiny::validate(shiny::need(nrow(combined) > 0, "No Peptide"))
    
    DT::datatable(
      combined,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollY = "60vh", #This should be good enough for the standard monitors. Otherwise need to make it responsive to browser window.
        scrollX = TRUE,
        dom = "Bfrtip"
      )
    )
  })
  
  #---------------------Dev Console-------------------------
  console_history <- reactiveVal("")
  
  observeEvent(input$console_run, {
    req(nchar(trimws(input$console_input)) > 0)
    
    result <- tryCatch(
      paste(capture.output(eval(parse(text = input$console_input),
                                envir = environment())),
            collapse = "\n"),
      error   = function(e) paste("Error:", conditionMessage(e)),
      warning = function(w) paste("Warning:", conditionMessage(w))
    )
    
    new_entry <- paste0(
      "> ", input$console_input, "\n",
      result, "\n",
      "---\n"
    )
    console_history(paste0(new_entry, console_history()))
  })
  
  observeEvent(input$console_clear, {
    console_history("")
  })
  
  output$console_output <- renderText({
    console_history()
  })
  
}
##---------------End------------

shinyApp(
  ui = ui,
  server = function(input, output, session) {
    server(input, output, session, 
           #input_variable = data_list,
           input_variable = preloaded_data[1:10,],
           generate_pseudo_sequence = FALSE, 
           custom_schema = NULL, 
           custom_signature = NULL, 
           replace_schema = FALSE)
  }
)