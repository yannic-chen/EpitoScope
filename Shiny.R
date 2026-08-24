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
  #motif_plot_length <- 7:20 # defined in global.R because ui.R uses it before server.
  netmhcpan_path <- "/mnt/c/Users/Yannic/netMHCpan-4.2/netMHCpan" # This is the absolute path in the WSL. Really want the system to read off .bashrc
  wsl_available <- reactiveVal(NULL)
  netmhcpan_available <- reactiveVal(NULL)
  software_r <- reactiveVal(NULL)
  
  observe({
    if (startup_done()) return() #Since there is no reactive dependency, this observe only runs once anyway. But just in case.
    
    os <- Sys.info()[["sysname"]] 
    
    if (os == "Windows") {
      #### --- WSL CHECK ---
      wsl_ok <- tryCatch(!is.null(system2("wsl", "--status", stdout = TRUE, stderr = TRUE)),
                         error = function(e) FALSE)
      wsl_available(wsl_ok)
      netmhcpan_use_wsl <<- TRUE
      message("WSL available: ", wsl_ok)
      
      #### --- netMHCpan CHECK (via WSL) ---
      netmhcpan_ok <- if (wsl_ok) tryCatch(
        system2("wsl", c("test", "-x", shQuote(netmhcpan_path)), stdout = FALSE, stderr = FALSE) == 0,
        error = function(e) FALSE) else FALSE
    } else {
      #### --- Linux / macOS: native, no WSL ---
      wsl_available(NA)
      netmhcpan_use_wsl <<- FALSE
      netmhcpan_ok <- nzchar(Sys.which(netmhcpan_path)) ||
        file.access(netmhcpan_path, mode = 1L) == 0
    }
    
    netmhcpan_available(netmhcpan_ok)
    message("Platform: ", os, " | netMHCpan available: ", netmhcpan_ok)
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
  
#------------------MHC setting---------------------  
  binder_thresholds <- reactive({
    s <- input$strong_cut; w <- input$weak_cut
    list(strong = if (is.null(s) || is.na(s)) 0.5 else s,
         weak   = if (is.null(w) || is.na(w)) 2   else w)
  })
  
  observeEvent(input$mhc_class, {
    req(input$mhc_class %in% c("I", "II"))
    thr <- if (input$mhc_class == "II") c(2, 10) else c(0.5, 2)
    updateNumericInput(session, "strong_cut", value = thr[1])
    updateNumericInput(session, "weak_cut",   value = thr[2])
    len <- if (input$mhc_class == "II") c(13, 25) else c(8, 11)
    updateSliderInput(session, "length_pct_range", value = len)
  }, ignoreInit = TRUE)
  
#---------------Visual setting-----------------
  plot_font_d <- debounce(reactive(if (is.null(input$plot_font)) 13 else input$plot_font), 400)
  
  plotly_plot_ids <- c(
    "Pairwise_shared_peptide_matrix","binding_plot","charge_plot_meas","dynrange_combined",
    "dynrange_grid","group_peptide_heatmap_interactive","length_plot","mass_plot","mass_plot2",
    "mass_plot_meas","measurement_heatmap","mz_plot","mz_plot2","mz_plot_meas",
    "pairwise_peptide_quant_correlation","pca","pca_variance","ppm_plot","ppm_plot2",
    "ppm_plot_meas","scatterplots","score_violin","score_violin2","score_violin_meas")
  
  # push font to every plotly WITHOUT re-rendering (client-side relayout, like a resize)
  observeEvent(plot_font_d(), {
    for (id in plotly_plot_ids)
      try(plotlyProxy(id, session) %>%
            plotlyProxyInvoke("relayout", list("font.size" = plot_font_d())), silent = TRUE)
  }, ignoreInit = TRUE)
  
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
      software_r(setNames(vapply(processed, function(x) x$software %||% NA_character_, character(1)),
                          names(processed)))
      
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
    req(data_list_r())
    sel <- intersect(input$selected_samples, names(data_list_r()))
    req(length(sel) > 0)                 # intersect guard against transient race for loading from SQL database.
    data_list_r()[sel]
  })
  
#-------------------Data Transformation tab-----------------------
  # Add delayed reaction
  filter_inputs <- reactive({
    list(
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
  
  output$length_slider_ui   <- renderUI({ make_range_slider_ui(data_list_r(), "LENGTH",       "length_range",   tagList(icon("ruler-horizontal"), "Filter by Length:"),   step = 1, current = isolate(input$length_range)) })
  output$quantity_slider_ui <- renderUI({ make_range_slider_ui(data_list_r(), "MAX_QUANTITY", "quantity_range", "Filter by Max Quantity:",                                step = 1, current = isolate(input$quantity_range)) })
  output$score_slider_ui    <- renderUI({ make_range_slider_ui(data_list_r(), "SCORE",        "score_range",    "Filter by Score:", current = isolate(input$score_range)) })
  output$charge_slider_ui   <- renderUI({ make_range_slider_ui(data_list_r(), "CHARGE",       "charge_range",   "Filter by Charge:",                                     step = 1, current = isolate(input$charge_range)) })
  output$mass_slider_ui     <- renderUI({ make_range_slider_ui(data_list_r(), "MASS",         "mass_range",     "Filter by Mass:",                                       digits = 2, current = isolate(input$mass_range)) })
  output$RT_slider_ui       <- renderUI({ make_range_slider_ui(data_list_r(), "RT",           "RT_range",       "Filter by RT:",                                         digits = 2, current = isolate(input$RT_range)) })
  
  
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
            mhc_length_range = input$mhc_length_range,
            binder_summary_all = safe_reactive(binder_summary_all),
            binder_unique = safe_reactive(peptide_wide_unique),
            binder_strong = binder_thresholds()$strong,
            binder_weak   = binder_thresholds()$weak,
            group_list = safe_reactive(group_list),
            group_comp_data = safe_reactive(group_comp_data),
            col_map = safe_reactive(measurement_col_map_r),
            annotation_table = if (is.data.frame(input_variable)) input_variable else NULL,
            default_quantity_cols_r = safe_reactive(default_quantity_cols_r),
            dynrange_prot_query = applied_search()$prot,
            dynrange_pep_query  = applied_search()$pep
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
    ann <- annotation_df_r()
    shiny::validate(shiny::need(is.data.frame(ann), "No annotation table provided."))
    
    DT::datatable(
      ann,
      rownames = FALSE,
      options  = list(pageLength = 10, autoWidth = TRUE),
      escape = FALSE,
    )
  })
  
  ## ----Number of peptides and peptidoforms----
  output$summary_peptides_plot <- renderPlot(
    unique_counts_plot(data_list_r(), "STRIPPED", "Number of unique peptides",plot_font_d(), input$color_palette))
  
  output$summary_peptidoforms_plot <- renderPlot(
    unique_counts_plot(data_list_r(), "PEPTIDE", "Number of unique peptidoforms",plot_font_d(), input$color_palette))
  
  output$summary_proteins_plot <- renderPlot(
    unique_counts_plot(data_list_r(), "PROTEIN", "Number of unique proteins",plot_font_d(), input$color_palette,
                       required_cols = "PROTEIN", transform_fn = extract_protein_prefixes))
  
  ## ----Charge / Mass / mz / RT / ppm----
  density_plot <- function(lst, column, x_label, font, palette, na_policy = "all"){
    check_data_error(lst, required_cols = column, na_policy = na_policy)
    scale_font(plot_density(lst, column = column, x_label = x_label, color = palette), font)
  }
  
  stacked_bar_plot <- function(lst, column, fill_label, font, palette,
                               percentage = FALSE, na_policy = "all"){
    check_data_error(lst, required_cols = column, na_policy = na_policy)
    scale_font(plot_stacked_bar(lst, column = column, fill_label = fill_label,
                                percentage = percentage, color = palette), font)
  }
  
  output$charge_plot <- renderPlot(stacked_bar_plot(data_list_r(), "CHARGE", "Charge", plot_font_d(), input$color_palette))
  output$mass_plot <- renderPlotly(density_plot(data_list_r(), "MASS", "Mass (Da)", isolate(plot_font_d()), input$color_palette))
  output$mz_plot   <- renderPlotly(density_plot(data_list_r(), "MZ", "m/z", isolate(plot_font_d()), input$color_palette))
  output$ppm_plot  <- renderPlotly(density_plot(data_list_r(), "PPM", "ppm", isolate(plot_font_d()), input$color_palette))
  
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
          scale_font(
            plot_histogram(df = df_local, column = "RT", x_label = "Retention Time", 
                         title_name = paste(sample_local), color = input$color_palette),
            plot_font_d())
        })
      })
    }
  })
  
  ##----Score distribution----
  build_score_violin <- function(font, palette) {
    lst <- data_list_r()
    check_data_error(lst, required_cols = "SCORE", na_policy = "any")
    scale_font(
      plot_violin(lst, column = "SCORE", color = palette),
      isolate(plot_font_d())) 
  }
  
  output$score_violin <- renderPlotly({build_score_violin(plot_font_d(), input$color_palette)})
  
  ##----Summary Table----
  output$RAW_summary_html <- renderUI({
    lst <- raw_list_r()
    check_data_error(lst, na_policy = "ignore")
    render_summary_pre(lst, font_size = "8px")
  })
#---------------------QC Tab-------------------------
  ## ---- length distribution ----
  build_length_plot <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "LENGTH", na_policy = "any")
    scale_font(plot_length_distribution(lst, color = palette), size = font)
  }
  output$length_plot <- renderPlotly({ build_length_plot(plot_font_d(),input$color_palette)})
  
  ## ---- length range percentage ----
  build_length_range_percentage <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "LENGTH", na_policy = "any")
    rng <- if (is.null(input$mhc_length_range)) c(8, 13) else input$mhc_length_range
    lst <- lapply(lst, function(df) {
      df$correct_range <- df$LENGTH >= rng[1] & df$LENGTH <= rng[2]
      df
    })
    scale_font(
      plot_stacked_bar(lst, column = "correct_range", fill_label = paste0(rng[1],"-",rng[2],"mer"), 
                       percentage = TRUE, color = palette, rev_levels = FALSE), font)
  }
  output$length_range_percentage <- renderPlot({build_length_range_percentage(plot_font_d(),input$color_palette)})
  
  build_length_range_percentage_meas <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "LENGTH", na_policy = "any")
    shiny::validate(shiny::need(length(default_quantity_cols_r()) > 0, "No QUANTITY columns found."))
    rng <- if (is.null(input$mhc_length_range)) c(8, 13) else input$mhc_length_range
    scale_font(
      plot_length_range_per_measurement(lst, default_quantity_cols_r(), color = palette,
                                        lo = rng[1], hi = rng[2]),font)
    }
  
  output$length_range_percentage_meas <- renderPlot({build_length_range_percentage_meas(plot_font_d(),input$color_palette)})
  
  ## ---- Other numeric columns on per measurement basis ----
  
  build_charge_plot_meas <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "CHARGE", na_policy = "all")
    scale_font(
      plot_charge_per_measurement(lst, default_quantity_cols_r(), color = palette, base_size = font),
      font)
    }
  output$charge_plot_meas <- plotly::renderPlotly({build_charge_plot_meas(plot_font_d(),input$color_palette)})
  
  build_mass_plot_meas <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MASS", na_policy = "any")
    scale_font(
      plot_density_envelope(lst, default_quantity_cols_r(), "MASS", "Mass (Da)", color = palette),
      font)
  }
  output$mass_plot_meas <- plotly::renderPlotly({build_mass_plot_meas(plot_font_d(),input$color_palette)})
  
  build_mz_plot_meas  <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MZ", na_policy = "any")
    scale_font(
      plot_density_envelope(lst, default_quantity_cols_r(), "MZ", "m/z", color = palette),
      font)
    }
  output$mz_plot_meas <- plotly::renderPlotly({build_mz_plot_meas(plot_font_d(),input$color_palette)})
  
  build_ppm_plot_meas <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PPM", na_policy = "any")
    scale_font(
      plot_density_envelope(lst, default_quantity_cols_r(), column = "PPM", x_label = "ppm", color = palette),
      font)
    }
  output$ppm_plot_meas <- plotly::renderPlotly({build_ppm_plot_meas(plot_font_d(),input$color_palette)})
  
  build_score_violin_meas <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "SCORE", na_policy = "any")
    scale_font(
      plot_violin_envelope(lst, default_quantity_cols_r(), column = "SCORE", x_label = "Score", color = palette),
      font)
    }
  output$score_violin_meas <- plotly::renderPlotly({build_score_violin_meas(plot_font_d(),input$color_palette)})
  
  build_RT_plot_meas <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "RT", na_policy = "all")
    qcols <- default_quantity_cols_r()
    layout_column_wrap(
      width = "400px",
      !!!lapply(names(lst), function(s) {
        pid <- paste0("rt_meas_", s)
        output[[pid]] <- renderPlot({
          scale_font(plot_rt_histogram_range(lst[[s]], s, qcols), font)
        })
        plotOutput(pid, height = "300px")
      })
    )
    }
  output$RT_plot_meas <- renderUI({build_RT_plot_meas(plot_font_d(),input$color_palette)})
  
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
    length_tabs <- lapply(motif_plot_length, function(L){
      sample_plots <- lapply(names(lst), function(sample_name)
        plotOutput(paste0("motif_", sample_name, "_", L), height = "180px"))
      tabPanel(
        title = paste("Length", L),
        div(`data-tab-key` = paste0("motif_len_", L),           # <- camera reads this
            tags$div(style = "text-align:center; margin-bottom:10px;",
                     plotOutput(paste0("motif_legend_", L), height = "120px")),
            do.call(layout_column_wrap, c(width = "250px", fixed_width = TRUE, sample_plots))
        )
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
  unique_counts_plot <- function(lst, column, y_label, font, palette,
                                 required_cols = NULL, na_policy = "ignore",
                                 transform_fn = identity, quantity_cols = NULL){
    check_data_error(lst, required_cols = required_cols, na_policy = na_policy)
    scale_font(
      plot_unique_counts(lst, column = column, y_label = y_label,
                         transform_fn = transform_fn,
                         color = palette, quantity_cols = quantity_cols),
      font)
  }
  
  output$summary_peptides_plot3 <- renderPlot(
    unique_counts_plot(processed_data_list(), "STRIPPED", "Number of unique peptides",plot_font_d(), input$color_palette,
                       quantity_cols = default_quantity_cols_r()))
  
  output$summary_peptidoforms_plot3 <- renderPlot(
    unique_counts_plot(processed_data_list(), "PEPTIDE", "Number of unique peptidoforms",plot_font_d(), input$color_palette,
                       quantity_cols = default_quantity_cols_r()))
  
  output$summary_proteins_plot3 <- renderPlot(
    unique_counts_plot(processed_data_list(), "PROTEIN", "Number of unique proteins",plot_font_d(), input$color_palette,
                       required_cols = "PROTEIN", na_policy = "all",transform_fn = extract_protein_prefixes,quantity_cols = default_quantity_cols_r()))
  
  ## ----Dynamic Range plots (individual, one subplot figure)----
  n_dynrange_panels <- reactive({
    lst <- processed_data_list()
    sum(vapply(lst, function(df) nrow(df) >= 10 && "MAX_QUANTITY" %in% colnames(df), logical(1)))
  })
  
  output$dynrange_individual_ui <- renderUI({
    check_data_error(processed_data_list(), required_cols = "MAX_QUANTITY", na_policy = "any")
    nrows <- max(1L, ceiling(n_dynrange_panels() / 3))
    plotly::plotlyOutput("dynrange_grid", height = paste0(nrows * 320, "px"))
  })
  
  output$dynrange_grid <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
    dynamic_range_subplot(lst, data_col = "MAX_QUANTITY", ncol = 3,
                          base_size = isolate(plot_font_d()))      # no scale_font wrapper
  })
  
  observeEvent(plot_font_d(), {
    s <- plot_font_d(); n <- n_dynrange_panels()
    ax <- unlist(lapply(seq_len(n), function(i) {
      sfx <- if (i == 1) "" else as.character(i)
      setNames(list(s, s), c(sprintf("xaxis%s.tickfont.size", sfx),
                             sprintf("yaxis%s.tickfont.size", sfx)))
    }), recursive = FALSE)
    ann <- if (n > 0) setNames(as.list(rep(s, n)), sprintf("annotations[%d].font.size", seq_len(n)-1)) else list()
    plotly::plotlyProxy("dynrange_grid", session) %>%
      plotly::plotlyProxyInvoke("relayout", c(list("font.size" = s), ax, ann))
  }, ignoreInit = TRUE)
  
  applied_search <- reactiveVal(list(prot = "", pep = ""))
  
  observeEvent(input$dynrange_go, {
    g <- function(x) if (is.null(x)) "" else trimws(x)
    applied_search(list(prot = g(input$dynrange_search),
                        pep  = g(input$dynrange_pep_search)))
  })
  
  observeEvent(list(applied_search()), {
    lst        <- processed_data_list()
    q          <- applied_search()
    prot_query <- q$prot
    pep_query  <- q$pep
    
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
  output$dynrange_combined <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
    rm <- isolate(if (is.null(input$dynrange_rank_mode)) "absolute" else input$dynrange_rank_mode)
    dynamic_range_plot_combined(df_list = lst, data_col = "MAX_QUANTITY",color = input$color_palette, rank_mode = rm)
  })
  
  observeEvent(list(applied_search(),input$dynrange_rank_mode), {
    lst        <- processed_data_list()
    q          <- applied_search()
    prot_query <- q$prot
    pep_query  <- q$pep
    rank_mode  <- if (is.null(input$dynrange_rank_mode)) "absolute" else input$dynrange_rank_mode
    
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
      n_pts     <- nrow(df)
      df$rank_x <- if (rank_mode == "relative") 100 * df$Rank / n_pts else df$Rank
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
      
      base_x[[i]] <- as.list(df$rank_x[base_idx]); base_y[[i]] <- as.list(df$y_vals[base_idx]); base_t[[i]] <- as.list(df$hover_text[base_idx])
      prot_x <- c(prot_x, df$rank_x[prot_idx]); prot_y <- c(prot_y, df$y_vals[prot_idx]); prot_t <- c(prot_t, df$hover_text[prot_idx])
      pep_x  <- c(pep_x,  df$rank_x[pep_idx]);  pep_y  <- c(pep_y,  df$y_vals[pep_idx]);  pep_t  <- c(pep_t,  df$hover_text[pep_idx])
    }
    
    base_hoverinfo <- if (length(prot_x) + length(pep_x) > 0) "skip" else "text"
    
    x_vals  <- c(base_x, list(as.list(prot_x)), list(as.list(pep_x)))
    y_vals  <- c(base_y, list(as.list(prot_y)), list(as.list(pep_y)))
    t_vals  <- c(base_t, list(as.list(prot_t)), list(as.list(pep_t)))
    hi_vals <- c(rep(list(base_hoverinfo), length(samples)), list("text"), list("text"))
    idxs    <- as.list(seq_len(length(samples) + 2L) - 1L)
    x_title <- if (rank_mode == "relative") "Rank (percentile %)" else "Rank"
    
    plotly::plotlyProxy("dynrange_combined", session) %>%
      plotly::plotlyProxyInvoke("restyle",
                                list(x = x_vals, y = y_vals, text = t_vals, hoverinfo = hi_vals), idxs) %>%
      plotly::plotlyProxyInvoke("relayout", list(xaxis = list(title = x_title)))
  })
  
  ##----1/k0 vs m/z----
  n_scatter_panels <- reactive({
    lst <- processed_data_list()
    sum(vapply(lst, function(df)
      !is.null(df) && nrow(df) > 0 && all(c("MZ","K0") %in% colnames(df)), logical(1)))
  })
  
  output$scatterplots_ui <- renderUI({
    nrows <- max(1L, ceiling(n_scatter_panels() / 3))
    plotly::plotlyOutput("scatterplots", height = paste0(nrows * 320, "px"))
  })
  
  output$scatterplots <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = c("MZ","K0"), na_policy = "all")
    plot_scatter_mz_k0_plotly(lst, color = input$color_palette, ncol = 3,
                              base_size = isolate(plot_font_d()))   # initial font, no font dependency
  })
  
  observeEvent(plot_font_d(), {
    s <- plot_font_d(); n <- n_scatter_panels()
    ax <- unlist(lapply(seq_len(n), function(i) {
      sfx <- if (i == 1) "" else as.character(i)
      setNames(list(s, s), c(sprintf("xaxis%s.tickfont.size", sfx),
                             sprintf("yaxis%s.tickfont.size", sfx)))
    }), recursive = FALSE)
    ann <- if (n > 0) setNames(as.list(rep(s, n)), sprintf("annotations[%d].font.size", seq_len(n)-1)) else list()
    plotly::plotlyProxy("scatterplots", session) %>%
      plotly::plotlyProxyInvoke("relayout", c(list("font.size" = s), ax, ann)) %>%
      plotly::plotlyProxyInvoke("restyle", list("textfont.size" = s))
  }, ignoreInit = TRUE)
  
  ## ---- measurement specific heatmap ----
  cor_mat_cache <- reactiveVal(matrix(NA, nrow = 5, ncol = 5)) #necessary as the space reserving rungs before the code.
  
  build_measurement_heatmap <- function(font, palette){
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    quantity_cols <- data_info_r() %>%
      dplyr::filter(final_name == "QUANTITY") %>% dplyr::select(-final_name) %>%
      unlist(recursive = TRUE, use.names = FALSE)
    shiny::validate(shiny::need(length(quantity_cols) > 0, "No QUANTITY columns found."))
    pep_mat <- prepare_measurement_matrix(lst, quantity_cols)
    shiny::validate(shiny::need(ncol(pep_mat) >= 2, "Need at least 2 groups for correlation."))
    shiny::validate(shiny::need(nrow(pep_mat) > 0, "No peptides to plot."))
    cor_mat <- cor(pep_mat, method = "pearson", use = "complete.obs")
    group_names <- names(lst)
    get_group <- function(colname){
      m <- group_names[sapply(group_names, function(g) startsWith(colname, g))]
      if (length(m) == 0) NA else m[1]
    }
    groups <- sapply(colnames(cor_mat), get_group)
    plot_correlation_heatmap_interactive(cor_mat, color = palette,
                                         cluster = input$cluster_mode_ea, groups = groups,
                                         label = "Pearson")
  }
  
  output$measurement_heatmap <- plotly::renderPlotly({build_measurement_heatmap(plot_font_d(), input$color_palette)})
  
#---------------------Results Tab-------------------------
  ##----Data Completeness----
  build_completeness_plot <- function(font, palette, percent = FALSE){
    lst <- processed_data_list()
    req(lst)
    spectra_cols <- data_info_r() %>%
      dplyr::filter(final_name == "SPECTRA") %>%
      dplyr::select(-final_name) %>%              # all sample columns
      unlist(recursive = TRUE, use.names = FALSE)
    scale_font(plot_completeness(lst, spectra_cols, percent = percent, color = palette), font)
  }
  
  output$completeness_plot  <- renderPlot(build_completeness_plot(plot_font_d(), input$color_palette))
  output$completeness_plot2 <- renderPlot(build_completeness_plot(plot_font_d(), input$color_palette, percent = TRUE))
  
  ##----Upset----
  build_upset_plot <- function(font){
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore") #By default plot_upset takes STRIPPED column. Need to adjust if we take PEPTIDE column instead.
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples to plot"))
    plot_upset(lst,
               min_size        = filter_inputs_d()$upset_min_size,
               min_degree      = filter_inputs_d()$upset_min_degree,
               n_intersections = filter_inputs_d()$upset_n_intersect,
               base_size = font)
  }
  
  output$upset_plot <- renderPlot({build_upset_plot(plot_font_d())})
  
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
  build_shared_peptide_matrix <- function(palette){
    lst <- processed_data_list()
    check_data_error(lst, na_policy = "ignore")
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples to plot"))
    req(input$shared_mode)
    plot_shared_peptide_interactive(lst, mode = input$shared_mode,
                                    color = palette, percent_type = "union")
  }
  
  build_pairwise_quant_correlation <- function(palette){
    lst <- processed_data_list()
    shiny::validate(shiny::need(length(lst) >= 2, "Need 2 or more samples to plot"))
    check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
    plot_pairwise_quant_correlation_interactive(lst, color = palette,
                                                cluster_mode = input$cluster_mode)
  }
  
  output$Pairwise_shared_peptide_matrix <- plotly::renderPlotly(build_shared_peptide_matrix(input$color_palette))
  
  output$pairwise_peptide_quant_correlation <- plotly::renderPlotly(build_pairwise_quant_correlation(input$color_palette))
  
  ## ----PCA plot----
  output$pca_group_ui <- renderUI({
    has_conditions <- annotation_provided() &&
      length(attr(annotation_df_r(), "condition_cols")) > 0
    cond_cols <- if (has_conditions) attr(annotation_df_r(), "condition_cols") else character(0)
    
    tagList(
      radioButtons("pca_group_mode", "Colour by",
                   choices = c("Sample" = "sample", "Manual groups" = "manual", "Condition" = "condition"),
                   selected = "sample", inline = TRUE),
      conditionalPanel("input.pca_group_mode == 'condition'",
                       selectInput("pca_condition_col", "Condition column", choices = cond_cols)),
      conditionalPanel("input.pca_group_mode == 'manual'",
                       helpText("Uses the Manual groups defined in the Group Comparison tab.")),
      if (!has_conditions) tags$script(HTML("
        setTimeout(function(){
          $('input[name=pca_group_mode][value=condition]').prop('disabled', true)
            .closest('label, .radio, .shiny-options-group > *')
            .css({'color':'#aaa','opacity':'0.55','cursor':'not-allowed'});
        }, 100);
      "))
    )
  })
  
  compute_sample_group <- function(mode, condition_col = NULL){
    samples <- names(processed_data_list())
    if (mode == "manual") {
      n <- if (is.null(input$n_groups)) 0 else input$n_groups
      m <- setNames(rep(NA_character_, length(samples)), samples)
      for (i in seq_len(n)) {
        gname <- input[[paste0("group_name_", i)]]
        samps <- input[[paste0("group_", i)]]
        if (!is.null(samps) && length(samps)) m[samps] <- gname
      }
      m
    } else if (mode == "condition" && annotation_provided()) {
      ann <- annotation_df_r()
      if (is.null(condition_col) || !condition_col %in% colnames(ann)) setNames(samples, samples)
      else setNames(as.character(ann[[condition_col]])[match(samples, ann$name)], samples)
    } else {
      setNames(samples, samples)
    }
  }
  
  pca_sample_group <- reactive({
    mode <- if (is.null(input$pca_group_mode)) "sample" else input$pca_group_mode
    compute_sample_group(mode, input$pca_condition_col)
  })
  
  compute_pca_fit <- function(lst, level, qc, group_names = names(lst)){
    shiny::validate(shiny::need(length(lst) >= 2, "Need \u22652 samples for PCA."))
    get_group <- function(cn){
      m <- group_names[sapply(group_names, function(g) startsWith(cn, g))]
      if (length(m) == 0) NA_character_ else m[1]
    }
    if (level == "measurement") {
      shiny::validate(shiny::need(length(qc) >= 3, "Need \u22653 measurements."))
      mat <- prepare_measurement_matrix(lst, qc)
      points <- colnames(mat); samp_of <- vapply(points, get_group, character(1))
    } else {
      peptides <- lapply(lst, function(df) df %>% dplyr::select(STRIPPED, MAX_QUANTITY) %>%
                           dplyr::group_by(STRIPPED) %>%
                           dplyr::summarise(MAX_QUANTITY = max(MAX_QUANTITY, na.rm = TRUE), .groups = "drop"))
      dfm <- Reduce(function(x, y) dplyr::full_join(x, y, by = "STRIPPED"), peptides)
      mat <- as.matrix(dfm[, -1]); colnames(mat) <- group_names; rownames(mat) <- dfm$STRIPPED
      points <- group_names; samp_of <- setNames(group_names, group_names)
    }
    mat[is.na(mat)] <- 0
    mat <- mat[apply(mat, 1, function(r) stats::sd(r) > 0), , drop = FALSE]
    shiny::validate(shiny::need(nrow(mat) >= 2 && ncol(mat) >= 3, "Not enough data/variance for PCA."))
    list(pca = prcomp(t(mat), scale. = TRUE), points = points, samp_of = samp_of)
  }
  
  pca_fit <- reactive({
    lst   <- processed_data_list()
    level <- if (is.null(input$pca_level)) "sample" else input$pca_level
    compute_pca_fit(lst, level, default_quantity_cols_r())
  })
  
  output$pca <- plotly::renderPlotly({
    fit <- pca_fit()
    s2g <- pca_sample_group()
    groups <- unname(s2g[fit$samp_of[fit$points]]); groups[is.na(groups)] <- "Ungrouped"
    dim <- if (is.null(input$pca_dim)) "2d" else input$pca_dim
    pca_scatter_plotly(fit$pca, labels = fit$points, groups = groups,
                       color = input$color_palette,
                       show_labels = length(fit$points) <= 30, dim = dim, base_size = isolate(plot_font_d()))
  })
  
  observeEvent(plot_font_d(), {
    s <- plot_font_d()
    plotly::plotlyProxy("pca", session) %>%
      plotly::plotlyProxyInvoke("relayout", list("font.size" = s)) %>%
      plotly::plotlyProxyInvoke("restyle", list("textfont.size" = s))
  }, ignoreInit = TRUE)
  
  output$pca_variance <- scale_font(plotly::renderPlotly({ plot_pca_variance(pca_fit()$pca) }),isolate(plot_font_d()))
  
  ## ----Number of peptides and peptidoforms----
  output$summary_peptides_plot2 <- renderPlot(
    unique_counts_plot(processed_data_list(), "STRIPPED", "Number of unique peptides",plot_font_d(), input$color_palette))
  
  output$summary_peptidoforms_plot2 <- renderPlot(
    unique_counts_plot(processed_data_list(), "PEPTIDE", "Number of unique peptidoforms",plot_font_d(), input$color_palette))
  
  output$summary_proteins_plot2 <- renderPlot(
    unique_counts_plot(processed_data_list(), "PROTEIN", "Number of unique proteins",plot_font_d(), input$color_palette,
                       required_cols = "PROTEIN", transform_fn = extract_protein_prefixes))
  
  ## ----Charge / Mass / mz / RT / ppm----
  output$charge_plot2 <- renderPlot(stacked_bar_plot(processed_data_list(), "CHARGE", "Charge", plot_font_d(), input$color_palette))
  output$mass_plot2 <- renderPlotly(density_plot(processed_data_list(), "MASS", "Mass (Da)", isolate(plot_font_d()), input$color_palette))
  output$mz_plot2   <- renderPlotly(density_plot(processed_data_list(), "MZ", "m/z", isolate(plot_font_d()), input$color_palette))
  output$ppm_plot2  <- renderPlotly(density_plot(processed_data_list(), "PPM", "ppm", isolate(plot_font_d()), input$color_palette))
  
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
          scale_font(plot_histogram(df = df_local, column = "RT", x_label = "Retention Time", 
                         title_name = paste(sample_local), color = input$color_palette),
                         plot_font_d())
        })
      })
    }
  })
  
  ##----Score distribution----
  output$score_violin2 <- plotly::renderPlotly({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "SCORE" , na_policy = "all")
    scale_font(
      plot_violin(lst, column = "SCORE", color = input$color_palette),
      isolate(plot_font_d()))
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
    n       <- input$n_groups
    samples <- names(active_data_list())
    
    tagList(lapply(seq_len(n), function(i) {
      prev_name <- isolate(input[[paste0("group_name_", i)]])
      prev_sel  <- isolate(input[[paste0("group_", i)]])
      
      tagList(
        textInput(paste0("group_name_", i), paste("Group", i, "name"),
                  value = if (is.null(prev_name) || !nzchar(prev_name)) paste("Group", i) else prev_name),
        selectInput(paste0("group_", i), paste("Select samples for Group", i),
                    choices  = samples,
                    selected = intersect(prev_sel, samples),
                    multiple = TRUE),
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
  
  ### ----Group Unique peptides----
  output$group_unique_bar <- renderPlot({
    sets <- group_peptide_sets()
    shiny::validate(shiny::need(length(sets) >= 1, "No groups defined"))
    lst <- lapply(sets, function(v) data.frame(STRIPPED = v, stringsAsFactors = FALSE))
    unique_counts_plot(lst, "STRIPPED", "Unique peptides", plot_font_d(), input$color_palette)
  })
  
  ### ----Group Euler----
  output$group_euler <- renderPlot({
    sets <- group_peptide_sets()
    shiny::validate(shiny::need(length(sets) >= 2, "Need 2 or more sets to compare"))
    
    group_euler_plot(sets, color = input$color_palette)
  })
  ### ----Group heatmap----
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
  
  
  observeEvent(
    plotly::event_data(
      "plotly_relayout",
      source = "group_hm",
      priority = "event"
    ),
    {
      ed <- plotly::event_data(
        "plotly_relayout",
        source = "group_hm",
        priority = "event"
      )
      
      req(ed)
      
      # ------------------------------------------------------------
      # 1. User zoom / pan
      # ------------------------------------------------------------
      if (!is.null(ed[["xaxis.range[0]"]]) &&
          !is.null(ed[["xaxis.range[1]"]])) {
        
        zoom_info_r(
          c(
            as.numeric(ed[["xaxis.range[0]"]]),
            as.numeric(ed[["xaxis.range[1]"]])
          )
        )
        
        return()
      }
      
      # ------------------------------------------------------------
      # 2. User double-click / explicit autorange
      # ------------------------------------------------------------
      if (isTRUE(ed[["xaxis.autorange"]])) {
        
        # Only reset if this wasn't caused by our own code
        age <- as.numeric(Sys.time()) - isolate(last_programmatic_t_r())
        
        if (age > 1) {
          zoom_info_r(NULL)
          
          if (isolate(heatmap_mode_r()) == "expanded") {
            expanded_peps_r(NULL)
            heatmap_mode_r("binned")
          }
        }
        
        return()
      }
    },
    ignoreInit = TRUE
  )

  zoom_debounced_r <- shiny::debounce(zoom_info_r, 500)

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
      
      x0 <- max(
        1L,
        as.integer(round(zoom[1])) + 1L
      )
      
      x1 <- min(
        n_cols,
        as.integer(round(zoom[2])) + 1L
      )
      
      max(1L, x1 - x0 + 1L)
      
    } else {
      
      n_cols
      
    }
    
    show <- n_visible <= 250
    
    plotly::plotlyProxy(
      "group_peptide_heatmap_interactive",
      session
    ) %>%
      plotly::plotlyProxyInvoke(
        "relayout",
        list(
          "xaxis.showticklabels" = show,
          "margin.b" = if (show) 120 else 30
        )
      )
  })
  
  heatmap_hover <- function(mat, bin_map){
    hv <- matrix("", nrow(mat), ncol(mat))
    for (j in seq_len(ncol(mat))) {
      cn <- colnames(mat)[j]
      if (!is.null(bin_map) && cn %in% names(bin_map)) {
        peps <- bin_map[[cn]]; n <- length(peps)
        shown <- paste(head(peps, 30), collapse = "<br>")
        extra <- if (n > 30) paste0("<br><i>+", n - 30, " more</i>") else ""
        hv[, j] <- paste0("<b>", n, " peptides in bin</b><br>", shown, extra)
      } else hv[, j] <- cn
    }
    hv
  }
  
  group_heatmap_figure <- function(mat, hover, palette, font = 8,
                                   source = NULL, ui_rev = NULL, register = FALSE){
    cscale <- if (palette == "default") list(list(0, "lightyellow"), list(1, "red"))
    else { vcols <- viridis::viridis(10, option = palette)
    lapply(seq_along(vcols) - 1, function(i) list(i/(length(vcols)-1), vcols[i+1])) }
    show_ticks <- ncol(mat) <= 250
    p <- plotly::plot_ly(
      source = source, x = colnames(mat), y = rownames(mat), z = mat, text = hover,
      type = "heatmap", colorscale = cscale, hovertemplate = "%{text}<extra></extra>",
      colorbar = list(title = "log10(QUANTITY+1)", tickfont = list(size = font))
    ) %>% plotly::layout(
      uirevision = ui_rev,
      xaxis  = list(title = "", showticklabels = show_ticks, tickfont = list(size = font), tickangle = -45),
      yaxis  = list(title = "", tickfont = list(size = font), autorange = "reversed"),
      margin = list(l = 130, b = if (show_ticks) 120 else 30),
      font   = list(size = font)
    )
    if (register) p <- plotly::event_register(p, "plotly_relayout")
    p
  }
  
  output$group_peptide_heatmap_interactive <- renderPlotly({
    hdata <- heatmap_data_r(); mode <- heatmap_mode_r()
    if (mode == "expanded" && !is.null(expanded_peps_r())) {
      mat <- hdata$mat_full[, expanded_peps_r(), drop = FALSE]; bin_map <- NULL
      ui_rev <- paste0("exp_", expand_counter_r())
    } else {
      mat <- hdata$mat_display; bin_map <- hdata$bin_map; ui_rev <- "binned"
    }
    mat   <- mat[hdata$row_ord, , drop = FALSE]
    hover <- heatmap_hover(mat, bin_map)
    group_heatmap_figure(mat, hover, input$color_palette,
                         font = isolate(plot_font_d()),
                         source = "group_hm", ui_rev = ui_rev, register = TRUE)
  })
  
  observeEvent(plot_font_d(), {
    plotly::plotlyProxy("group_peptide_heatmap_interactive", session) %>%
      plotly::plotlyProxyInvoke("relayout", list(
        "font.size" = plot_font_d(),
        "xaxis.tickfont.size" = plot_font_d(),
        "yaxis.tickfont.size" = plot_font_d()))
  }, ignoreInit = TRUE)
  
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
                   helpText("UniProt accession IDs (e.g. P04439) are first converted to genes. Remaining IDs are assumed to be protein IDs and converted via uniprot online request."),
                   plotlyOutput(paste0("go_term_", name))),
            column(6, h6("STRING-DB network"),
                   helpText("Protein IDs are directly send to STRING DB to let their side resolve the names."),
                   visNetwork::visNetworkOutput(paste0("STRING_", name), height = "500px"))
        ))
      })
      
      # Generate the tabsetPanel from the list
      do.call(tabsetPanel, c(id = "volcano_comparison_tabs", comparison_tabs))
    })
    uiOutput("volcano_comparison_tabs")
  })
  
  build_go_plot <- function(df, ont = NULL){
    bg <- if (is.null(input$go_background)) "genome" else input$go_background
    universe <- NULL; bg_note <- "whole genome"
    if (bg == "detected") { universe <- detected_universe(); bg_note <- "detected proteins" }
    else if (bg == "custom") {
      cust <- if (!is.null(input$go_custom_ids) && nzchar(input$go_custom_ids))
        unlist(strsplit(input$go_custom_ids, "[[:space:],;]+")) else NULL
      if (length(cust)) { universe <- .protein_to_entrez(cust); bg_note <- "custom list" }
    }
    ont_use <- if (!is.null(ont)) ont else if (is.null(input$go_ont)) "BP" else input$go_ont
    run_go_enrichment(df, universe = universe, bg_note = bg_note, ont = ont_use)
  }
  
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

        output[[paste0("volcano_", plot_name)]] <- plotly::renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          sel <- selected_peptide()
          scale_font(group_volcano_plot(df,sel,plot_name), isolate(plot_font_d()))
        })
        
        ## ----MA plot----
        output[[paste0("ma_", plot_name)]] <- plotly::renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          sel <- selected_peptide()
          scale_font(group_MA_plot(df, sel,plot_name), isolate(plot_font_d()))
        })
        
        ## ----P-value histogram----
        output[[paste0("pval_hist_", plot_name)]] <- plotly::renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          scale_font(group_p_histogram(df,plot_name), isolate(plot_font_d()))
        })
        
        ## ----Ranked Fold Change----
        output[[paste0("rank_fc_", plot_name)]] <- renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          sel <- selected_peptide()
          scale_font(group_rank_FC(df,plot_name,sel), isolate(plot_font_d()))
        })
        ## ----Peptide Fold Change table----
        output[[paste0("peptide_table_", plot_name)]] <- DT::renderDT({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
          group_FC_table(df)
        })
        
        ## ----GO term----
        output[[paste0("go_term_", plot_name)]] <- renderPlotly({
          shiny::validate(shiny::need(!is.null(df) && nrow(df) != 0, "No data available."))
          build_go_plot(df)
        })
        
        ## ----STRING-DB----
        output[[paste0("STRING_", plot_name)]] <- visNetwork::renderVisNetwork({
          shiny::validate(shiny::need(curl::has_internet(), "No Internet Connection."))
          shiny::validate(shiny::need(!is.null(df) && nrow(df) != 0, "No data available."))
          net <- run_string(df)
          shiny::validate(shiny::need(!is.null(net),
                                      "No STRING network (no significant hits, or none resolved by STRING)."))
          net
        })
      })
    }
  })
  
  observeEvent(plot_font_d(), {
    s  <- plot_font_d()
    vl <- isolate(group_comp_data())        # read names, but don't make it a trigger
    if (is.null(vl)) return()
    for (name in names(vl)) {
      for (id in paste0(c("volcano_", "ma_", "pval_hist_", "rank_fc_", "go_term_"), name)) {
        try({
          plotly::plotlyProxy(id, session) %>%
            plotly::plotlyProxyInvoke("relayout", list(
              "font.size"           = s,
              "xaxis.tickfont.size" = s,
              "yaxis.tickfont.size" = s)) %>%
            plotly::plotlyProxyInvoke("restyle", list("textfont.size" = s))
        }, silent = TRUE)
      }
    }
  }, ignoreInit = TRUE)
  
  detected_universe <- reactive({
    lst <- processed_data_list(); req(length(lst) > 0)
    all_prot <- unlist(lapply(lst, function(d) if ("PROTEIN" %in% names(d)) d$PROTEIN),
                       use.names = FALSE)
    .protein_to_entrez(unique(all_prot))
  })
  
  selected_peptide <- reactive({
    ed <- plotly::event_data("plotly_click", source = "group_diff")
    if (is.null(ed) || is.null(ed$key)) return(NULL)
    ed$key
  })
  
#-------------------PTM------------------------
  PTM_stacked_bar_plot <- function(lst, column, fill_label, font, palette,
                               percentage = FALSE, na_policy = "all", rev_levels = TRUE){
    check_data_error(lst, required_cols = column, na_policy = na_policy)
    scale_font(plot_stacked_bar(lst, column = column, fill_label = fill_label,
                                percentage = percentage, color = palette,
                                rev_levels = rev_levels), font)
  }
  
  prep_ptm_long <- function(lst){
    lapply(lst, function(df){
      df %>%
        dplyr::mutate(PTM = ifelse(is.na(PTM) | PTM == "", "Unmodified", PTM)) %>%
        tidyr::separate_rows(PTM, sep = "\\s*[,;]\\s*") %>%
        dplyr::mutate(PTM = gsub("^\\d+", "", PTM))   # strip FragPipe position prefix
    })
  }
  
  output$PTM_plot <- renderPlot({
    lst <- processed_data_list()
    check_data_error(lst, required_cols = "PTM", na_policy = "ignore")   # check before preprocessing
    PTM_stacked_bar_plot(prep_ptm_long(lst), "PTM", "PTM", plot_font_d(), input$color_palette,
                     na_policy = "ignore", rev_levels = FALSE)
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
        do.call(layout_column_wrap, c(list(width = "250px"), fixed_width = TRUE, sample_plots))
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
  
  observeEvent(list(input$strong_cut, input$weak_cut), {
    s <- input$strong_cut; w <- input$weak_cut
    sel <- if      (isTRUE(s == 0.5 && w == 2))  "I"
    else if (isTRUE(s == 2   && w == 10)) "II"
    else    character(0)                       # custom -> no radio ticked
    cur <- if (length(input$mhc_class)) input$mhc_class else character(0)
    if (!identical(cur, sel))
      updateRadioButtons(session, "mhc_class", selected = sel)
  }, ignoreInit = TRUE)
  
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
    set_order   <- unique(df$Set) 
    allele_cols <- grep("^HLA", colnames(df), value = TRUE)
    
    res <- df %>%
      dplyr::group_by(Set, STRIPPED) %>%
      dplyr::summarise(dplyr::across(all_of(allele_cols),
                                     ~ if (all(is.na(.x))) NA_real_ else min(.x, na.rm = TRUE)), .groups = "drop")
    res$Set <- factor(res$Set, levels = set_order)
    res
  })
  
  binder_summary_all <- reactive({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    req(!is.null(peptide_wide_unique()))
    compute_binder_summary(peptide_wide_unique(), alleles = input$allele_viz_select, 
                           strong = binder_thresholds()$strong, weak = binder_thresholds()$weak)
  })
  
  output$binding_summary <- DT::renderDT({
    cache <- prediction_cache()
    shiny::validate(shiny::need(ncol(cache) > 1, "No prediction data"))
    req(!is.null(peptide_wide_unique()))
    binder_summary_all()
  })
  
  output$binding_plot_ui <- renderUI({
    cache <- prediction_cache()
    if (ncol(cache) <= 1) return(helpText("No prediction data"))
    view <- if (is.null(input$binding_view)) "best" else input$binding_view
    
    if (view == "per_allele") {
      d      <- peptide_wide_unique()
      all_c  <- grep("^HLA", colnames(d), value = TRUE)
      sel    <- input$allele_viz_select
      n_all  <- if (is.null(sel) || !length(sel)) length(all_c) else length(intersect(all_c, sel))
      n_samp <- length(unique(d$Set))
      grp    <- if (is.null(input$binding_group)) "sample" else input$binding_group
      # panels = the facet dimension; rows within each panel = the other dimension
      if (grp == "sample") { n_panel <- n_samp; rows_each <- n_all  }
      else                 { n_panel <- n_all;  rows_each <- n_samp }
      per <- max(140, rows_each * 30 + 70)
      plotly::plotlyOutput("binding_plot", height = paste0(n_panel * per, "px"))
    } else {
      n_samp <- length(unique(peptide_wide_unique()$Set))
      plotly::plotlyOutput("binding_plot", height = paste0(max(300, n_samp * 55 + 150), "px"))
    }
  })
  
  build_binder_plot <- function(font, palette, view, percent, grp, alleles, orientation = "h"){
    shiny::validate(shiny::need(ncol(prediction_cache()) > 1, "No prediction data"))
    df <- peptide_wide_unique()
    shiny::validate(shiny::need(!is.null(df), "No prediction data"))
    if (view == "per_allele")
      plot_binders_per_allele_plotly(df, color = palette, percent = percent,
                                     alleles = alleles, facet_by = grp, orientation = orientation)
    else
      plot_binders_plotly(df, color = palette, percent = percent, alleles = alleles, orientation = orientation)
  }
  
  output$binding_plot <- renderPlotly({
    view    <- if (is.null(input$binding_view))  "best"     else input$binding_view
    percent <- (if (is.null(input$binding_scale)) "absolute" else input$binding_scale) == "percent"
    grp     <- if (is.null(input$binding_group)) "sample"   else input$binding_group
    build_binder_plot(plot_font_d(), input$color_palette, view, percent, grp, input$allele_viz_select)
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
  
  #---------------------SQL---------------------------------
  ## ----------------save SQL-------------------------------
  meta_table_rv <- reactiveVal(NULL)
  
  .meta_defaults <- function(annotation_df) {
    colget <- function(nm) {
      if (is.null(annotation_df)) return("")
      i <- which(tolower(names(annotation_df)) == nm)
      if (!length(i)) "" else paste(unique(annotation_df[[i[1]]]), collapse = "; ")
    }
    list(
      submitted_by = { op <- colget("operator"); if (nzchar(op)) op
      else tryCatch(unname(Sys.info()[["user"]]), error = function(e) "") },
      instrument   = colget("instrument"),
      conditions   = if (!is.null(annotation_df) && length(attr(annotation_df, "condition_cols")))
        paste(attr(annotation_df, "condition_cols"), collapse = ", ") else ""
    )
  }
  
  observeEvent(input$save_to_db, {
    lst <- data_list_r() 
    ann <- if (annotation_provided()) annotation_df_r() else NULL
    meta_table_rv(build_meta_table(lst, default_quantity_cols_r(),
                                   ann, measurement_col_map_r(), software_r()))
    showModal(modalDialog(
      title = "Save run to database", size = "xl",
      textInput("meta_user", "Submitted by:",
                value = tryCatch(unname(Sys.info()[["user"]]), error = function(e) "")),
      textAreaInput("meta_description", "Run description / notes:", rows = 2),
      tags$hr(),
      tags$b("Per-measurement metadata — double-click a cell to edit:"),
      rhandsontable::rHandsontableOutput("meta_edit_table"),
      footer = tagList(modalButton("Cancel"),
                       actionButton("confirm_save_db", "Save to database",
                                    class = "btn-primary", icon = icon("database"))),
      uiOutput("lowinfo_warn")
    ))
  })
  
  # px width for a column: wide enough for its longest single word (so wrapping
  # breaks cleanly on spaces) and its content, within sane bounds.
  .col_px <- function(nm, data = NULL, char_px = 8, pad = 26, min_px = 50, max_px = 240) {
    longest_word <- max(nchar(strsplit(nm, "\\s+")[[1]]), 0)          # header, per-word
    content_w    <- if (!is.null(data))
      suppressWarnings(max(nchar(as.character(data)), 0, na.rm = TRUE)) else 0
    w <- max(longest_word, min(content_w, 24)) * char_px + pad
    min(max(w, min_px), max_px)
  }
  
  output$meta_edit_table <- rhandsontable::renderRHandsontable({
    mt <- meta_table_rv(); req(mt)
    con <- .db_con(); on.exit(DBI::dbDisconnect(con))
    
    widths <- vapply(names(mt), function(cc) .col_px(cc, mt[[cc]]), numeric(1))
    
    ht <- rhandsontable::rhandsontable(mt, rowHeaders = NULL, height = 340,
                                       colWidths = unname(widths)) %>%
      rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
    
    for (cc in names(mt)) {
      if (is_vocab_field(cc)) {
        ht <- rhandsontable::hot_col(
          ht, col = cc, type = "dropdown",
          source = condition_term_source(cc, con),
          allowInvalid = TRUE, strict = FALSE)
      }
    }
    ht
  })
  
  db_refresh <- reactiveVal(0)
  
  # the modal's Save button does the actual write
  observeEvent(input$confirm_save_db, {
    edited <- if (is.null(input$meta_edit_table)) meta_table_rv()
    else rhandsontable::hot_to_r(input$meta_edit_table)
    spectra_cols <- data_info_r() %>%
      dplyr::filter(final_name == "SPECTRA") %>%
      dplyr::select(-final_name) %>%
      unlist(recursive = TRUE, use.names = FALSE)
    id <- tryCatch(
      save_analysis_to_db(lst = data_list_r(), meta_table = edited,
                          quantity_cols = default_quantity_cols_r(), col_map = measurement_col_map_r(),
                          spectra_cols = spectra_cols, data_info = data_info_r(), mod_map = data_mod_map(),
                          submitted_by = input$meta_user, description = input$meta_description),
      epito_duplicate = function(e) {
        showModal(modalDialog(
          title = "Possible duplicate",
          paste0("This analysis is ", conditionMessage(e),
                 ". Save it anyway as a new entry?"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("save_dup_anyway", "Save anyway", class = "btn-warning")
          )))
        NULL
      },
      error = function(e) { showNotification(paste("Save failed:", e$message), type = "error"); NULL })
    if (!is.null(id)) {
      removeModal()
      showNotification(paste("Saved as", id), type = "message")
      db_refresh(db_refresh() + 1)
    }
  })
  
  observeEvent(input$save_dup_anyway, {
    removeModal()
    edited <- rhandsontable::hot_to_r(input$meta_edit_table)
    spectra_cols <- data_info_r() %>% dplyr::filter(final_name == "SPECTRA") %>%
      dplyr::select(-final_name) %>% unlist(recursive = TRUE, use.names = FALSE)
    id <- save_analysis_to_db(lst = data_list_r(), meta_table = edited,
                              quantity_cols = default_quantity_cols_r(), col_map = measurement_col_map_r(),
                              spectra_cols = spectra_cols, data_info = data_info_r(), mod_map = data_mod_map(),
                              submitted_by = input$meta_user, description = input$meta_description,
                              allow_duplicate = TRUE)
    if (!is.null(id)) { showNotification(paste("Saved as", id), type = "message"); db_refresh(db_refresh() + 1) }
  })
  
  analyses_tbl <- reactive({ db_refresh(); list_analyses() })
  
  output$saved_runs_table <- DT::renderDT({
    DT::datatable(analyses_tbl(), rownames = FALSE, options = list(pageLength = 10), selection = "single")
  })
  
  output$meta_edit_table <- rhandsontable::renderRHandsontable({
    mt <- meta_table_rv(); req(mt)
    con <- .db_con(); on.exit(DBI::dbDisconnect(con))
    
    ht <- rhandsontable::rhandsontable(mt, rowHeaders = NULL, height = 340) %>%
      rhandsontable::hot_context_menu(allowRowEdit = FALSE, allowColEdit = FALSE)
    
    for (cc in names(mt)) {
      if (is_vocab_field(cc)) {
        ht <- rhandsontable::hot_col(
          ht, col = cc, type = "dropdown",
          source = condition_term_source(cc, con),
          allowInvalid = TRUE, strict = FALSE)   # permissive: new terms allowed & auto-register on save
      }
    }
    ht
  })
  
  output$lowinfo_warn <- renderUI({
    mt <- tryCatch(rhandsontable::hot_to_r(input$meta_edit_table),
                   error = function(e) meta_table_rv())
    req(mt)
    fl <- scan_lowinfo_conditions(mt)
    if (!nrow(fl)) return(NULL)
    items <- sprintf("<li><b>%s</b> &rarr; <code>%s = %s</code></li>",
                     fl$sample, fl$field, fl$value)
    HTML(paste0(
      "<div style='color:#8a5000;border:1px solid #e0b080;background:#fff8ee;",
      "padding:10px;border-radius:6px;margin-top:10px'>",
      "&#9888; These conditions look uninformative on their own &mdash; ",
      "will others know what it means?<ul style='margin:6px 0'>",
      paste(items, collapse = ""), "</ul>",
      "Prefer a descriptive state (e.g. <code>Tumor</code> / <code>Normal</code>) ",
      "over <code>yes</code> / <code>no</code>. You can still save as-is.",
      "</div>"))
  })
  
  ##-------------------load SQL------------------
  observeEvent(input$load_from_db, {
    sel <- input$saved_runs_table_rows_selected
    req(length(sel) == 1)
    aid <- analyses_tbl()$analysis_id[sel]
    
    loaded <- load_analysis_from_db(aid)
    dfs <- loaded$data
    di <- loaded$data_info
    dmm <- loaded$data_mod_map
    
    # core data reactiveVals
    raw_list_r(dfs)                              # note: normalized data, not original raw format
    data_list_r(dfs)
    data_info_r(di)
    data_mod_map(dmm)
    software_r(loaded$software)
    
    # annotation / grouping — route through the SAME validator as an upload
    ann <- if (!is.null(loaded$annotation))
      tryCatch(check_annotation_table(loaded$annotation), error = function(e) {
        showNotification(paste("Annotation rebuild skipped:", conditionMessage(e)),
                         type = "warning"); NULL
      }) else NULL
    
    if (!is.null(ann)) {
      annotation_df_r(ann)
      annotation_provided(TRUE)
      measurement_col_map_r(build_measurement_col_map(dfs, ann))
    } else {
      annotation_provided(FALSE)
      measurement_col_map_r(NULL)
    }
    
    condition_groups_r(list())                   # clear any groups from the prior session
    showNotification(sprintf("Loaded %s (%d samples)", aid, length(dfs)), type = "message")
  })
  
  ##--------------------query SQL-------------------
  # picking a table pre-fills a SELECT
  observeEvent(input$sql_table_pick, {
    req(input$sql_table_pick)
    updateTextAreaInput(session, "sql_query",
                        value = paste0("SELECT * FROM ", input$sql_table_pick, " LIMIT 50;"))
  })
  
  sql_result_rv <- reactiveVal(NULL)
  observeEvent(input$sql_run, {
    q <- trimws(input$sql_query)
    # read-only guard: only allow SELECT / WITH / PRAGMA
    if (!grepl("^(select|with|pragma)\\b", tolower(q))) {
      showNotification("Only SELECT / WITH / PRAGMA queries are allowed here.", type = "error"); return()
    }
    res <- tryCatch({
      con <- .db_con(); on.exit(DBI::dbDisconnect(con))
      DBI::dbGetQuery(con, q)
    }, error = function(e) { showNotification(paste("Query error:", conditionMessage(e)), type = "error"); NULL })
    sql_result_rv(res)
  })
  
  output$sql_result <- DT::renderDT({
    req(sql_result_rv())
    DT::datatable(sql_result_rv(), rownames = FALSE, filter = "top",
                  options = list(pageLength = 25, scrollX = TRUE))
  })
  
  ##--------------------browse SQL-------------------
  observe({ db_refresh(); updateSelectInput(session, "browse_table", choices = db_tables()) })
  
  output$browse_result <- DT::renderDT({
    req(input$browse_table)
    con <- .db_con(); on.exit(DBI::dbDisconnect(con))
    tbl <- DBI::dbQuoteIdentifier(con, input$browse_table)
    n   <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) n FROM %s", tbl))$n
    df  <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s LIMIT 5000", tbl))
    DT::datatable(df, filter = "top", rownames = FALSE,
                  caption = if (n > 5000) sprintf("Showing first 5,000 of %s rows.", n) else NULL,
                  options = list(pageLength = 25, scrollX = TRUE))
  })
  
  #---------------------Plot publication modal------------
  export_catalog <- list(
    peptides_unique = list(
      title = "Unique peptides", engine = "ggplot",
      build = function(font, palette)
        unique_counts_plot(processed_data_list(),"STRIPPED", "Number of unique peptides", font, palette)
    ),
    peptides_unique_meas = list(
      title = "Unique peptides distribution", engine = "ggplot",
      build = function(font, palette)
        unique_counts_plot(processed_data_list(),"STRIPPED", "Number of unique peptides", font, palette, quantity_cols = default_quantity_cols_r())
    ),
    peptidoforms_unique = list(
      title = "Unique peptidoforms", engine = "ggplot",
      build = function(font, palette)
        unique_counts_plot(processed_data_list(),"PEPTIDE", "Number of unique peptidoforms", font, palette)
    ),
    peptidoforms_unique_meas = list(
      title = "Unique peptidoforms distribution", engine = "ggplot",
      build = function(font, palette)
        unique_counts_plot(processed_data_list(),"PEPTIDE", "Number of unique peptidoforms", font, palette, quantity_cols = default_quantity_cols_r())
    ),
    proteins_unique = list(
      title = "Unique proteins", engine = "ggplot",
      build = function(font, palette)
        unique_counts_plot(processed_data_list(),"PROTEIN", "Number of unique proteins", font, palette,required_cols = "PROTEIN", na_policy = "all", 
                           transform_fn = extract_protein_prefixes)
    ),
    proteins_unique_meas = list(
      title = "Unique proteins distribution", engine = "ggplot",
      build = function(font, palette)
        unique_counts_plot(processed_data_list(),"PROTEIN", "Number of unique proteins", font, palette,required_cols = "PROTEIN", na_policy = "all", 
                           transform_fn = extract_protein_prefixes, quantity_cols = default_quantity_cols_r())
    ),
    length_distribution = list(
      title = "Peptide Length Distribution", engine = "plotly",
      build = function(font, palette) build_length_plot(font, palette)
    ),
    length_range_distribution = list(
      title = "Peptide Length Range %", engine = "ggplot",
      build = function(font, palette) build_length_range_percentage(font, palette)
    ),
    length_range_distribution_meas = list(
      title = "Peptide Length Range % Distribution", engine = "ggplot",
      build = function(font, palette) build_length_range_percentage_meas(font, palette)
    ),
    charge_plot = list(
      title = "Charge", engine = "ggplot",
      build = function(font, palette) stacked_bar_plot(processed_data_list(), "CHARGE", "Charge", font, palette)
    ),
    charge_plot_meas = list(
      title = "Charge (measurement)", engine = "plotly",
      build = function(font, palette) build_charge_plot_meas(font, palette)
    ),
    mass_plot = list(
      title = "Mass", engine = "plotly",
      build = function(font, palette) density_plot(processed_data_list(), "MASS", "Mass (Da)", font, palette)
      ),
    mass_plot_meas = list(
      title = "Mass (measurement)", engine = "plotly",
      build = function(font, palette) build_mass_plot_meas(font, palette)
    ),
    mz_plot = list(
      title = "M/Z", engine = "plotly",
      build = function(font, palette) density_plot(processed_data_list(), "MZ", "m/z", font, palette)
      ),
    mz_plot_meas = list(
      title = "M/Z (measurement)", engine = "plotly",
      build = function(font, palette) build_mz_plot_meas(font, palette)
    ),
    ppm_plot  = list(
      title = "Mass Uncertainty", engine = "plotly",
      build = function(font, palette) density_plot(processed_data_list(), "PPM", "ppm", font, palette)
      ),
    ppm_plot_meas = list(
      title = "Mass Uncertainty (measurement)", engine = "plotly",
      build = function(font, palette) build_ppm_plot_meas(font, palette)
    ),
    score_plot = list(
      title = "Score", engine = "plotly",
      build = function(font, palette) build_score_violin(font, palette)
    ),
    score_plot_meas = list(
      title = "Score (measurement)", engine = "plotly",
      build = function(font, palette) build_score_violin_meas(font, palette)
    ),
    rt_histograms = list(
      title = "Retention Time", engine = "ggplot",
      build = function(font, palette)
        panel_grid(processed_data_list(), function(df, s){
          check_data_error(df, required_cols = "RT", na_policy = "any")
          scale_font(plot_histogram(df, column = "RT", x_label = "Retention Time",
                                    title_name = s, color = palette), font)
        })
    ),
    rt_per_measurement = list(
      title = "Retention Time (measurement)", engine = "ggplot",
      build = function(font, palette){
        lst   <- processed_data_list()
        check_data_error(lst, required_cols = "RT", na_policy = "all")
        qcols <- default_quantity_cols_r()
        panel_grid(lst, function(df, s) scale_font(plot_rt_histogram_range(df, s, qcols), font))
      }
    ),
    mz_k0 = list(
      title = "1/k0 vs m/z", engine = "plotly",
      build = function(font, palette){
        lst <- processed_data_list()
        check_data_error(lst, required_cols = c("MZ","K0"), na_policy = "all")
        plot_scatter_mz_k0_plotly(lst, color = palette, ncol = 3, base_size = font)
      }
    ),
    dynrange_individual = list(
      title = "Dynamic range (per sample)", engine = "plotly",
      build = function(font, palette){
        lst <- processed_data_list()
        check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
        dynamic_range_subplot(lst, data_col = "MAX_QUANTITY", ncol = 3, base_size = font)
      }
    ),
    dynrange_combined = list(
      title = "Dynamic range (combined)", engine = "plotly",
      controls = radioButtons("exp_dynrange_rank_mode", "Rank mode",
                              c("Absolute" = "absolute", "Relative" = "relative"),
                              selected = "absolute", inline = TRUE),
      build = function(font, palette){
        lst <- processed_data_list()
        check_data_error(lst, required_cols = "MAX_QUANTITY", na_policy = "any")
        rm  <- if (is.null(input$exp_dynrange_rank_mode)) "absolute" else input$exp_dynrange_rank_mode
        dynamic_range_plot_combined(lst, "MAX_QUANTITY", color = palette, rank_mode = rm)
      }
    ),
    measurement_heatmap = list(
      title = "Measurement correlation", engine = "plotly",
      build = function(font, palette) build_measurement_heatmap(font, palette)
      ),
    completeness_count = list(
      title = "Data completeness (count)", engine = "ggplot",
      build = function(font, palette) build_completeness_plot(font, palette, percent = FALSE)
    ),
    completeness_percent = list(
      title = "Data completeness (%)", engine = "ggplot",
      build = function(font, palette) build_completeness_plot(font, palette, percent = TRUE)
    ),
    upset_plot = list(
      title = "Upset plot", engine = "ggplot",
      build = function(font, palette) build_upset_plot(font)
    ),
    shared_peptide_matrix = list(
      title = "Pairwise shared peptides", engine = "plotly",
      build = function(font, palette) build_shared_peptide_matrix(palette)
    ),
    pairwise_quant_correlation = list(
      title = "Pairwise quant correlation", engine = "plotly",
      build = function(font, palette) build_pairwise_quant_correlation(palette)
    ),
    pca_scatter = list(
      title = "PCA", engine = "plotly",
      controls = function(){
        cond_cols <- if (annotation_provided()) attr(annotation_df_r(), "condition_cols") else character(0)
        tagList(
          radioButtons("exp_pca_level",  "Level", c("Sample"="sample","Measurement"="measurement"),
                       selected = "sample", inline = TRUE),
          radioButtons("exp_pca_dim",    "View",  c("2D"="2d","3D"="3d"), selected = "2d", inline = TRUE),
          radioButtons("exp_pca_colour", "Colour by",
                       c("Sample"="sample","Manual groups"="manual","Condition"="condition"),
                       selected = "sample", inline = TRUE),
          conditionalPanel("input.exp_pca_colour == 'condition'",
                           selectInput("exp_pca_condition_col", "Condition column", choices = cond_cols))
        )
      },
      build = function(font, palette){
        lst   <- processed_data_list()
        level <- if (is.null(input$exp_pca_level))  "sample" else input$exp_pca_level
        mode  <- if (is.null(input$exp_pca_colour)) "sample" else input$exp_pca_colour
        dim   <- if (is.null(input$exp_pca_dim))    "2d"     else input$exp_pca_dim
        
        fit    <- compute_pca_fit(lst, level, default_quantity_cols_r())
        s2g    <- compute_sample_group(mode, input$exp_pca_condition_col)
        groups <- unname(s2g[fit$samp_of[fit$points]]); groups[is.na(groups)] <- "Ungrouped"
        
        pca_scatter_plotly(fit$pca, labels = fit$points, groups = groups, color = palette,
                           show_labels = length(fit$points) <= 30, dim = dim, base_size = font)
      }
    ),
    group_unique_bar = list(
      title = "Unique peptides per group", engine = "ggplot",
      build = function(font, palette){
        sets <- group_peptide_sets()
        shiny::validate(shiny::need(length(sets) >= 1, "No groups defined"))
        lst  <- lapply(sets, function(v) data.frame(STRIPPED = v, stringsAsFactors = FALSE))
        unique_counts_plot(lst, "STRIPPED", "Unique peptides", font, palette)
      }
    ),
    group_euler = list(
      title = "Group Euler", engine = "ggplot",
      build = function(font, palette){
        sets <- group_peptide_sets()
        shiny::validate(shiny::need(length(sets) >= 2, "Need 2 or more sets to compare"))
        scale_font(group_euler_plot(sets, color = palette), font)
      }
    ),
    group_heatmap = list(
      title = "Group peptide heatmap", engine = "plotly",
      build = function(font, palette){
        hdata <- heatmap_data_r()
        mat   <- hdata$mat_display[hdata$row_ord, , drop = FALSE]
        hover <- heatmap_hover(mat, hdata$bin_map)
        group_heatmap_figure(mat, hover, palette, font = font)     # static: source/ui_rev/register omitted
      }
    ),
    group_stats = list(
      title = "Group comparison stats", engine = "plotly",
      controls = function(){
        comps <- tryCatch(names(group_comp_data()), error = function(e) character(0))
        tagList(
          selectInput("exp_stat_comp", "Comparison", choices = comps),
          radioButtons("exp_stat_type", "Plot",
                       c("Volcano"="volcano", "MA"="ma", "P-value hist"="pval",
                         "Rank-FC"="rankfc", "GO-term"="go"),          # <- added
                       selected = "volcano", inline = TRUE),
          conditionalPanel("input.exp_stat_type == 'go'",              # <- GO-only control
                           radioButtons("exp_go_ont", "GO ontology",
                                        c("BP"="BP","CC"="CC","MF"="MF"),
                                        selected = isolate(if (is.null(input$go_ont)) "BP" else input$go_ont),
                                        inline = TRUE))
        )
      },
      build = function(font, palette){
        vl   <- group_comp_data()
        comp <- input$exp_stat_comp
        shiny::validate(shiny::need(!is.null(comp) && comp %in% names(vl), "Select a comparison"))
        df   <- vl[[comp]]
        shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No data to plot"))
        df$Significance <- ifelse(
          is.na(df$negLog10AdjP_BH) | is.na(df$log2FC), "Missing",
          ifelse(df$negLog10AdjP_BH > 1.3 & abs(df$log2FC) > 1, "Significant", "Not significant"))
        
        sel <- NULL
        switch(if (is.null(input$exp_stat_type)) "volcano" else input$exp_stat_type,
               volcano = group_volcano_plot(df, sel, comp),
               ma      = group_MA_plot(df, sel, comp),
               pval    = group_p_histogram(df, comp),
               rankfc  = group_rank_FC(df, comp, sel),
               go      = { g <- build_go_plot(df, input$exp_go_ont)
               shiny::validate(shiny::need(!is.null(g), "No enriched GO terms for this comparison"))
               g }
               )
      }
    ),
    ptm_distribution = list(
      title = "PTM distribution", engine = "ggplot",
      build = function(font, palette){
        lst <- processed_data_list()
        check_data_error(lst, required_cols = "PTM", na_policy = "ignore")
        PTM_stacked_bar_plot(prep_ptm_long(lst), "PTM", "PTM", font, palette,
                         na_policy = "ignore", rev_levels = FALSE)
      }
    ),
    predicted_binders = list(
      title = "Predicted binders", engine = "plotly",
      controls = function(){
        d     <- tryCatch(peptide_wide_unique(), error = function(e) NULL)
        all_c <- if (!is.null(d)) grep("^HLA", colnames(d), value = TRUE) else character(0)
        tagList(
          radioButtons("exp_binding_orient", "Bars",
                       c("Horizontal" = "h", "Vertical" = "v"),
                       selected = "h", inline = TRUE),
          radioButtons("exp_binding_view",  "View",
                       c("Best per peptide" = "best", "Per allele" = "per_allele"),
                       selected = "best", inline = TRUE),
          radioButtons("exp_binding_scale", "Scale",
                       c("Absolute" = "absolute", "Percent" = "percent"),
                       selected = "absolute", inline = TRUE),
          conditionalPanel("input.exp_binding_view == 'per_allele'",
                           radioButtons("exp_binding_group", "Facet by",
                                        c("Sample" = "sample", "Allele" = "allele"),
                                        selected = "sample", inline = TRUE)),
          selectInput("exp_allele_viz_select", "Alleles",
                      choices = all_c, selected = all_c, multiple = TRUE)
        )
      },
      build = function(font, palette){
        view    <- if (is.null(input$exp_binding_view))  "best"     else input$exp_binding_view
        percent <- (if (is.null(input$exp_binding_scale)) "absolute" else input$exp_binding_scale) == "percent"
        grp     <- if (is.null(input$exp_binding_group)) "sample"   else input$exp_binding_group
        orient <- if (is.null(input$exp_binding_orient)) "h" else input$exp_binding_orient
        build_binder_plot(font, palette, view, percent, grp, input$exp_allele_viz_select, orientation = orient)
      }
    )
  )
  
  motif_entries <- setNames(
    lapply(motif_plot_length, function(L){
      force(L)                                   # capture L per closure
      list(
        title  = paste("Sequence motif — length", L), engine = "ggplot",
        build  = function(font, palette){
          lst <- processed_data_list()
          panel_grid(lst, function(df, s){
            peps <- unique(df$STRIPPED[df$LENGTH == L])
            if (length(peps) < 5) return(NULL)   # panel_grid drops NULLs
            scale_font(plot_seqlogo(peps, title = paste(s, "• Length", L)), font)
          })
        }
      )
    }),
    paste0("motif_len_", motif_plot_length)      # keys: motif_len_7 … motif_len_20
  )
  export_catalog <- c(export_catalog, motif_entries)
  
  open_export_modal <- function(initial){
    showModal(modalDialog(
      title = "Publication export", size = "l", easyClose = TRUE, footer = modalButton("Close"),
      
      # (1) near-fullscreen sizing + preview area/paper styling (auto-removed on close)
      tags$style(HTML("
  .modal-dialog { max-width: 94vw !important; width: 94vw; margin: 3vh auto; }
  .modal-content { height: 94vh; background: #fff; }
  .modal-body   { height: calc(94vh - 58px); overflow: hidden; display: flex; flex-direction: column; }
  .exp-content-row { flex: 1 1 auto; min-height: 0; }
  .exp-controls { height: 100%; overflow-y: auto; min-height: 0; }

  /* flex-centered paper — matches the JS (scale about center, no translate) */
  .exp-preview-area {
    height: 100%; position: relative; overflow: hidden;
    display: flex; align-items: flex-start; justify-content: flex-start;
  }
  .exp-paper {
    flex: 0 0 auto;
    background: #fff; box-shadow: 0 0 0 1px rgba(0,0,0,.18);
  }

  /* height chain so the plot fills the paper */
  #exp_preview_ui { height: 100%; width: 100%; display: block; }
  #exp_preview_ui > div, #exp_preview_ui .plotly,
  #exp_preview_ui .shiny-plot-output, #exp_preview_ui img {
    height: 100% !important; width: 100% !important;
  }

  .exp-controls { height: 100%; overflow-y: auto; }

  .modal.fade .modal-dialog { transition: none !important; transform: none !important; }
  .modal.fade            { transition: none !important; }
  .modal-backdrop.fade   { transition: none !important; }
  .modal-backdrop.show   { opacity: .5; }
")),
      
      fluidRow(
        column(12,
               selectInput("exp_which", "Figure",
                           choices  = setNames(names(export_catalog),
                                               vapply(export_catalog, `[[`, "", "title")),
                           selected = initial, width = "100%")
        )
      ),
      fluidRow(class = "exp-content-row",
        column(2, class = "exp-controls",
               uiOutput("exp_extra_controls"),
               sliderInput ("exp_font", "Font size", 8, 28, 13),
               selectInput ("exp_palette", "Palette", c("default", "viridis", "magma", "inferno", "plasma", "cividis", "mako", "rocket", "turbo"),
                            selected = isolate(input$color_palette)),
               textInput("exp_title", "Title",        placeholder = "(keep default)"),
               textInput("exp_xlab",  "X-axis label", placeholder = "(keep default)"),
               textInput("exp_ylab",  "Y-axis label", placeholder = "(keep default)"),
               fluidRow(column(6, numericInput("exp_w","Width (in)", 8, 1, 40, 0.5)),
                        column(6, numericInput("exp_h","Height (in)",6, 1, 40, 0.5))),
               fluidRow(column(6, numericInput("exp_dpi","DPI", 300, 72, 600, 1)),
                        column(6, selectInput ("exp_fmt","Format", c("png","svg","jpeg","pdf")))),
               uiOutput("exp_dl_ui")
        ),
        # (2) wrap the preview in area + paper so fitExportPaper has targets
        column(10,
               div(class = "exp-preview-area",
                   div(id = "exp_paper", class = "exp-paper",
                       uiOutput("exp_preview_ui"))))
      )
    ))
    shinyjs::runjs("setTimeout(window.fitExportPaper, 120);")
  }
  
  observeEvent(input$open_export, open_export_modal(names(export_catalog)[1]))
  
  exp_entry <- reactive({ req(input$exp_which); export_catalog[[input$exp_which]] })
  
  exp_plot <- reactive({
    e <- exp_entry(); req(e)
    p <- e$build(input$exp_font, input$exp_palette)
    p <- apply_export_labels(p, e$engine, input$exp_title, input$exp_xlab, input$exp_ylab)
    if (identical(e$engine, "plotly")) {
      if (is.null(p$x$layout)) p$x$layout <- list()
      p$x$layout$width  <- input$exp_w * 96
      p$x$layout$height <- input$exp_h * 96
      p <- plotly::config(p, responsive = FALSE)
    }
    p
  }) %>% debounce(300)
  
  output$exp_preview_ui <- renderUI({
    e <- exp_entry(); req(e)
    if (e$engine == "plotly") plotlyOutput("exp_preview_plotly", height = "100%", width = "100%")
    else                      plotOutput ("exp_preview_ggplot",  height = "100%", width = "100%")
  })
  
  output$exp_extra_controls <- renderUI({
    e <- exp_entry(); req(e)
    if (is.null(e$controls)) return(NULL)
    if (is.function(e$controls)) e$controls() else e$controls
  })
  
  output$exp_preview_plotly <- renderPlotly({ req(exp_entry()$engine == "plotly"); exp_plot() })
  output$exp_preview_ggplot <- renderPlot(
    { req(exp_entry()$engine == "ggplot"); exp_plot() },
    width  = function() input$exp_w * 96,
    height = function() input$exp_h * 96,
    res = 96
  )
  
  output$exp_dl_ui <- renderUI({
    e <- exp_entry(); req(e)
    if (e$engine == "plotly")
      actionButton("exp_dl_plotly", "Download", icon = icon("download"), class = "btn-primary")
    else
      downloadButton("exp_dl_gg", "Download", class = "btn-primary")
  })
  
  # ggplot -> ggsave at exact inches/DPI
  output$exp_dl_gg <- downloadHandler(
    filename = function() paste0(input$exp_which, ".", input$exp_fmt),
    content  = function(file)
      ggsave(file, exp_plot(), width = input$exp_w, height = input$exp_h,
             dpi = input$exp_dpi, units = "in", device = input$exp_fmt, limitsize = FALSE)
  )
  
  # plotly -> client-side capture (reuses your existing downloadPlotly handler)
  observeEvent(input$exp_dl_plotly, {
    session$sendCustomMessage("downloadPlotly", list(
      id = "exp_preview_plotly", format = input$exp_fmt,
      width = round(input$exp_w * 96), height = round(input$exp_h * 96),
      scale = input$exp_dpi / 96, filename = input$exp_which))
  })
  
  # Adjust modal window on window size change
  observe({
    req(input$exp_w, input$exp_h, input$exp_which)
    shinyjs::runjs(sprintf(
      "var p=document.getElementById('exp_paper');
     if(p){ p.dataset.pxw='%f'; p.dataset.pxh='%f'; setTimeout(window.fitExportPaper, 30); }",
      input$exp_w * 96, input$exp_h * 96))
  })
  
  # per-card camera buttons  -> open on the clicked plot
  observeEvent(input$export_open, open_export_modal(input$export_open))
  
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
           input_variable = preloaded_data,
           generate_pseudo_sequence = FALSE, 
           custom_schema = NULL, 
           custom_signature = NULL, 
           replace_schema = FALSE)
  }
)
