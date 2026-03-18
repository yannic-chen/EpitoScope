# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Script Name:        Shiny.R
# Purpose:            This script creates a shiny app for easy immunopeptidomics analysis
# Author:             Yannic Chen
# Date Created:       2025-11-19
# Last Modified:      2026-02-02
# Version:            0.2
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
#
# Version:
# --- v0.2 ---
# switched from base shiny to bs4Dash UI
# now also works with PeaksXPro peptide.csv. -> missing charge column now temporarily assigns charge 0 to everything.
# changed RT density plot to histogram. Density plot still exists.
#
# --- v0.1 ---
# Initial version

source("global.R") #global.R must be in the same folder. Otherwise change this path.
source("ui.R") #ui.R must be in the same folder. Otherwise change this path.

options(shiny.maxRequestSize = 5*1024^3) #Increase upload limit (in bytes) if needed. 1024^3 = 1 GB
options(width=10000) #This allows for text to not be text-wrapped.

server <- function(input, output, session, preloaded_data = NULL, generate_pseudo_sequence = FALSE) {
  
#------------------Data management-----------------  
  ## Reactive dataset container
  raw_list_r <- reactiveVal(NULL) #this is the list of raw data
  data_list_r <- reactiveVal(NULL) #this is the list of trimmed and filtered data
  data_info_r <- reactiveVal(NULL) #this is the info list to know which columns are used for what
  data_mod_map <- reactiveVal(NULL) #this is the conversion map when modified amino acids are given their own symbol. Only used in PTM analysis.
  prediction_cache <- reactiveVal(data.frame(Peptide = character())) #This is to save netMHCpan predictions

  #Check if we preload_data
  observe({
    if (is.null(data_list_r())) {
      start <- Sys.time() #measure time
      
      raw_list_r(preloaded_data)
      
      # Process but split df and table
      processed <- lapply(preloaded_data, normalize_df)
      
      print(Sys.time() - start)
      
      
      #Add the PTM_Pseudo sequence to each dataframe in the list. We do it outside of normalize_df() function because it unifies the mod_map across all dataframes in the list
      if (generate_pseudo_sequence) {
        start <- Sys.time()
        print("Generating PTM_Pseudo sequence")
        all_peptides <- unlist(lapply(processed, function(sample) {
          sample$df$PEPTIDE
        }), use.names = FALSE)
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
        dfs <- lapply(dfs, function(df) {
          df_dt <- as.data.table(df)   # temporary conversion
          result <- netMHCpan[df_dt, on = c(Peptide = "STRIPPED")]  # fast left join
          
          result <- as.data.frame(result)        # convert back to data.frame
          names(result)[names(result) == "Peptide"] <- "STRIPPED" #Rename back
          result
        })
      } else {
        print("No netMHCpan pro-computed variable. SKIP binding prediction")
      }
      #The predicted_cache is initialized using all peptides in the input data. If the netMHCpan precomputed data has been left_joined, these will also be taken.
      prediction <- do.call(rbind, lapply(dfs, function(df) {
          # pick STRIPPED + all columns starting with HLA
          hla_cols <- grep("^HLA", colnames(df), value = TRUE)
          df_subset <- df[, c("STRIPPED", hla_cols), drop = FALSE]
          colnames(df_subset)[1] <- "Peptide"
          df_subset <- df_subset[!duplicated(df_subset$Peptide), ]
          df_subset
        }))
      
      # Extract summary tables
      infos <- lapply(processed, function(x) x$table)
      merged_info <- infos %>%
        imap(~ .x %>% dplyr::rename(!!.y := coalesced_list)) %>%  # .y = sample name
        purrr::reduce(full_join, by = "final_name")

      print(Sys.time() - start)
  
      prediction_cache(prediction)
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
    req(input$selected_samples)
    data_list_r()[input$selected_samples]
  })
  
  PTM_Pseudo_sequence <- reactive({
    
  })
  
#-------------------Data Transformation tab-----------------------
  
  # Preprocessing data to align it with GUI
  processed_data_list <- reactive({
    lst <- active_data_list()  # Only selected samples
    req(lst)
    
    # Apply filters
    lst_transformed <- lapply(lst, function(df) {
      # Length filter
      if(!is.null(input$length_range)) {
        df <- df[df$LENGTH >= input$length_range[1] & df$LENGTH <= input$length_range[2], ]
      }
      
      # Quantity filter
      if ("MAX_QUANTITY" %in% names(df)) {
        
        req(input$quantity_range)
        req(length(input$quantity_range) == 2)
        
        df <- df %>%
          dplyr::filter(
            MAX_QUANTITY >= input$quantity_range[1],
            MAX_QUANTITY <= input$quantity_range[2]
          )
      }
      
      # Score filter
      if("SCORE" %in% colnames(df) && !is.null(input$score_range)) {
        df <- df[df[["SCORE"]] >= input$score_range[1] &
                   df[["SCORE"]] <= input$score_range[2], ]
      }
      
      # charge filter
      if(!is.null(input$charge_range)) {
        df <- df[df$CHARGE >= input$charge_range[1] & df$CHARGE <= input$charge_range[2], ]
      }
      
      # mass filter
      if(!is.null(input$mass_range)) {
        df <- df[df$MASS >= input$mass_range[1] & df$MASS <= input$mass_range[2], ]
      }
      
      # RT filter
      if(!is.null(input$RT_range)) {
        df <- df[df$RT >= input$RT_range[1] & df$RT <= input$RT_range[2], ]
      }
      
      df <- df[rowSums(!is.na(df)) > 0, , drop = FALSE] #Remove rows that are all NA
      
      df
    })
    
    lst_transformed
  })
  
  output$length_slider_ui <- renderUI({
    lst <- active_data_list()
    req(lst)
    
    if (!any(sapply(lst, function(df) "LENGTH" %in% colnames(df)))) {
      return(tags$div(style = "color: #b30000; font-style: italic;","No LENGTH column in data."))
    }
    
    # Compute min/max across all selected samples
    min_val <- min(sapply(lst, function(df) min(df[["LENGTH"]], na.rm = TRUE)), na.rm = TRUE)
    max_val <- max(sapply(lst, function(df) max(df[["LENGTH"]], na.rm = TRUE)), na.rm = TRUE)
    
    sliderInput(
      "length_range",
      tagList(icon("ruler-horizontal"),"Filter by Length:"),
      min = min_val,
      max = max_val,
      value = c(min_val, max_val),
      step = 1
    )
  })
  
  output$quantity_slider_ui <- renderUI({
    lst <- active_data_list()
    req(lst)
    
    if (!any(sapply(lst, function(df) "MAX_QUANTITY" %in% colnames(df)))) {
      return(tags$div(style = "color: #b30000; font-style: italic;","No MAX_QUANTITY column in data."))
    }
    
    # Compute min/max across all selected samples
    min_val <- min(sapply(lst, function(df) min(df[["MAX_QUANTITY"]], na.rm = TRUE)), na.rm = TRUE)
    max_val <- max(sapply(lst, function(df) max(df[["MAX_QUANTITY"]], na.rm = TRUE)), na.rm = TRUE)
    
    # Ensure valid slider range
    if (min_val == max_val) {
      max_val <- min_val + 1
    }
    
    sliderInput(
      "quantity_range",
      "Filter by Max Quantity:",
      min = min_val,
      max = max_val,
      value = c(min_val, max_val),
      step = 1
    )
  })

  output$score_slider_ui <- renderUI({
    lst <- active_data_list()  # already standardized/filtered datasets
    req(lst)
    
    if (!any(sapply(lst, function(df) "SCORE" %in% colnames(df)))) {
      return(tags$div(style = "color: #b30000; font-style: italic;","No SCORE column in data."))
    }
    
    # Compute min/max across all selected samples
    min_val <- min(sapply(lst, function(df) min(df[["SCORE"]], na.rm = TRUE)), na.rm = TRUE)
    max_val <- max(sapply(lst, function(df) max(df[["SCORE"]], na.rm = TRUE)), na.rm = TRUE)
    
    sliderInput("score_range",
                "Filter by Score:",
                min = min_val,
                max = max_val,
                value = c(min_val, max_val)
    )
  })
  
  output$charge_slider_ui <- renderUI({
    lst <- active_data_list()
    req(lst)
    
    if (!any(sapply(lst, function(df) "CHARGE" %in% colnames(df)))) {
      return(tags$div(style = "color: #b30000; font-style: italic;","No CHARGE column in data."))
    }

    # Compute min/max across all selected samples
    min_val <- min(sapply(lst, function(df) min(df[["CHARGE"]], na.rm = TRUE)), na.rm = TRUE)
    max_val <- max(sapply(lst, function(df) max(df[["CHARGE"]], na.rm = TRUE)), na.rm = TRUE)
    
    sliderInput(
      "charge_range",
      "Filter by Charge:",
      min = min_val,
      max = max_val,
      value = c(min_val, max_val),
      step = 1
    )
  })
  
  output$mass_slider_ui <- renderUI({
    lst <- active_data_list()
    req(lst)
    
    if (!any(sapply(lst, function(df) "MASS" %in% colnames(df)))) {
      return(tags$div(style = "color: #b30000; font-style: italic;","No MASS column in data."))
    }
    
    # Compute min/max across all selected samples
    min_val <- min(sapply(lst, function(df) min(df[["MASS"]], na.rm = TRUE)), na.rm = TRUE)
    max_val <- max(sapply(lst, function(df) max(df[["MASS"]], na.rm = TRUE)), na.rm = TRUE)
    
    sliderInput(
      "mass_range",
      "Filter by Mass:",
      min = min_val,
      max = max_val,
      value = c(min_val, max_val),
      step = 1
    )
  })
  
  output$RT_slider_ui <- renderUI({
    lst <- active_data_list()
    req(lst)
    
    if (!any(sapply(lst, function(df) "RT" %in% colnames(df)))) {
      return(tags$div(style = "color: #b30000; font-style: italic;","No RT column in data."))
    }
    
    # Compute min/max across all selected samples
    min_val <- min(sapply(lst, function(df) min(df[["RT"]], na.rm = TRUE)), na.rm = TRUE)
    max_val <- max(sapply(lst, function(df) max(df[["RT"]], na.rm = TRUE)), na.rm = TRUE)
    
    sliderInput(
      "RT_range",
      "Filter by RT:",
      min = min_val,
      max = max_val,
      value = c(min_val, max_val),
      step = 1
    )
  })
  
  observeEvent(input$generate_report, {
    
    # Choose an output file name
    out_file <- paste0("EpitoScope_Report_", Sys.Date(), ".html")
    
    # Use a temporary directory (required for Shiny Server / shinyapps.io)
    out_path <- file.path(tempdir(), out_file)
    
    # Show progress
    withProgress(message = "Generating report...", value = 0, {
      
      incProgress(0.1, detail = "Preparing data...")
      Sys.sleep(0.1)  # optional, simulate preprocessing
      
      incProgress(0.2, detail = "Processing tables and plots...")
      Sys.sleep(0.1)  # optional, simulate heavy processing
      
      # Render the R Markdown report
      incProgress(0.5, detail = "Rendering R Markdown...")
      rmarkdown::render(
        input = "report.Rmd", #report.Rmd must be in the same folder. Otherwise change path here.
        output_file = out_path,
        params = list(
          data_info = data_info_r(),
          data_list = data_list_r(),
          processed_data_list = processed_data_list()
        ),
        envir = new.env(parent = globalenv())
      )
      
      incProgress(0.2, detail = "Finalizing report...")
      Sys.sleep(0.1)
      
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
      content = function(file) {
        file.copy(out_path, file)
      }
    )
  })
  
#---------------------Summary Tab-------------------------
  ## ---------Column Mapping----------------
  output$summary_table <- DT::renderDT({
    req(data_info_r())
    DT::datatable(
      Colum_mapping_table(data_info_r()),
      options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
      escape = FALSE,
      rownames = FALSE
    )
  })
  
  ## ----Number of peptides and peptidoforms----
  output$summary_peptides_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "STRIPPED" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No STRIPPED column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_unique_counts(lst, column = "STRIPPED", y_label = "Number of unique peptides", color = input$color_palette)
  })
  
  output$summary_peptidoforms_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "PEPTIDE" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No PEPTIDE column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_unique_counts(lst, column = "PEPTIDE", y_label = "Number of unique peptidoforms", color = input$color_palette)
  })
  
  output$summary_proteins_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "PROTEIN" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No PROTEIN column in data.", cex = 1.2)
      return(invisible())
    }
    
    # Use extract_protein_prefixes for proteins
    plot_unique_counts(lst, column = "PROTEIN", y_label = "Number of unique proteins",
                       transform_fn = extract_protein_prefixes, color = input$color_palette)
  })
  
  ## ----Charge / Mass / mz / RT / ppm----
  output$charge_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "CHARGE" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No CHARGE column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_stacked_bar(lst, column = "CHARGE", fill_label = "Charge", percentage = FALSE, color = input$color_palette)
  })
  
  # Example usage for your Shiny outputs
  output$mass_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MASS" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MASS column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_density(lst, column = "MASS", x_label = "Mass (Da)", color = input$color_palette)
  })
  
  output$mz_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MZ" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MZ column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_density(lst, column = "MZ", x_label = "m/z", color = input$color_palette)
  })
  
  output$RT_plot <- renderUI({
    lst <- data_list_r()
    req(lst)
    
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
    req(lst)
    
    for (sample_name in names(lst)) {
      
      local({
        sample_local <- sample_name
        df_local <- lst[[sample_local]]
        
        output[[paste0("RT_", sample_local)]] <- renderPlot({
          plot_histogram(df = df_local, column = "RT", x_label = "Retention Time (min)", 
                         title_name = paste("RT Histogram –", sample_local), color = input$color_palette
          )
        })
      })
    }
  })
  
  output$ppm_plot <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "PPM" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No PPM column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_density(lst, column = "PPM", x_label = "ppm", color = input$color_palette)
  })
  
  ##----Score distribution----
  output$score_violin <- renderPlot({
    lst <- data_list_r()
    req(lst)
    
    if (!any(vapply(lst, function(df) "SCORE" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No SCORE column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_violin(lst, column = "SCORE", color = input$color_palette)
  })
  
  ##----Summary Table----
  output$RAW_summary_html <- renderUI({
    lst <- raw_list_r()
    req(lst)
    
    render_summary_pre(lst, font_size = "8px")
  })
#---------------------QC Tab-------------------------
  ## ---- length distribution ----
  output$length_plot <- renderPlotly({
    lst <- processed_data_list()
    req(lst)
    
    plot_length_distribution(lst, color = input$color_palette)
  })
  
  ## ---- length range percentage ----
  output$length_range_percentage <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    lst <- lapply(lst, function(df) {
      df$correct_range <- df$LENGTH >= 8 & df$LENGTH <= 13
      df
    })
    
    plot_stacked_bar(lst, column = "correct_range", fill_label = "8-13mer", percentage = TRUE, color = input$color_palette)
  })
  
  ## ---- Motif Plot ----
  global_legend_plot <- ggseqlogo::ggseqlogo("ACDEFGHIKLMNPQRSTVWY") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Amino Acid Colors")
  
  lengths <- 7:20
  for (L in lengths) {
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
    req(lst)
    
    lengths <- 7:20
    
    # Build a list of tabPanels for each peptide length
    length_tabs <- lapply(lengths, function(L) {
      
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
    req(lst)
    
    lengths <- 7:20
    
    for (L in lengths) {
      for (sample_name in names(lst)) {
        
        local({
          length_val <- L
          sample_val <- sample_name
          
          output[[paste0("motif_", sample_val, "_", length_val)]] <- renderPlot({
            
            df <- lst[[sample_val]]
            req(df)
            
            # subset by length
            df_L <- df[df$LENGTH == length_val, ]
            
            peptides <- unique(df_L[["STRIPPED"]])
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
              title = paste(sample_val, "• Length", length_val)
            )
          })
        })
      }
    }
  })
  
  ## ----Dynamic Range plots----
  output$dynrange_individual_ui <- renderUI({
    lst <- processed_data_list()
    req(lst)
    
    # Wrap plots in a grid (like motif plots)
    layout_column_wrap(
      width = "400px",  # each plot approx width
      !!!lapply(names(lst), function(sample_name) {
        plotOutput(paste0("dynrange_", sample_name), height = "300px")
      })
    )
  })
  
  output$dynrange_combined <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MAX_QUANTITY" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MAX_QUANTITY column in data.", cex = 1.2)
      return(invisible())
    }
    
    dynamic_range_plot_combined(df_list = lst, data_col = "MAX_QUANTITY", color = input$color_palette)
  })
  
  observe({
    lst <- processed_data_list()
    req(lst)
    
    for (sample_name in names(lst)) {
      
      local({
        sample_local <- sample_name
        df_local <- lst[[sample_local]]
        
        output[[paste0("dynrange_", sample_local)]] <- renderPlot({
          req(df_local)
          
          if (nrow(df_local) < 10) {
            plot.new()
            text(0.5, 0.5, "Not enough peptides", cex = 1.4)
            return()
          }
          
          dynamic_range_plot(
            df = df_local,
            data_col = "MAX_QUANTITY",
            title_name = paste("Dynamic Range –", sample_local),
            name_col = "PROTEIN",
            gene = "HLA",               # Default highlight can be HLA
            gene_regex = "HLA[A-C]+" #Only HLA-A,B and C
          )
        })
      })
    }
  })
  
  ##----1/k0 vs m/z----
  output$scatterplots_ui <- renderUI({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MZ" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MZ column in data.", cex = 1.2)
      return(invisible())
    }
    
    if (!any(vapply(lst, function(df) "K0" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No K0 column in data.", cex = 1.2)
      return(invisible())
    }
    
    layout_column_wrap(
      width = "300px",
      !!!lapply(names(lst), function(sample_name) {
        plotOutput(paste0("scatter_", sample_name), height = "300px")
      })
    )
  })
  
  observe({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MZ" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MZ column in data.", cex = 1.2)
      return(invisible())
    }
    
    if (!any(vapply(lst, function(df) "K0" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No K0 column in data.", cex = 1.2)
      return(invisible())
    }
    
    plots <- generate_scatterplots(lst, color = input$color_palette)
    
    for (i in seq_along(names(lst))) {
      local({
        sample_name_local <- names(lst)[i]
        plot_local <- plots[[i]]
        
        output[[paste0("scatter_", sample_name_local)]] <- renderPlot({
          plot_local
        })
      })
    }
  })
  
  ##----aa heatmap----
  output$aa_heatmap <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if(length(lst) < 2) {
      plot.new()
      text(0.5, 0.5,"At least 2 samples required.", cex = 1.2)
      return(invisible())
    }
    
    ht <- plot_aa_composition(
      lst,
      color = input$color_palette
    )
    
    draw(ht)
    
  }, res = 144,
  height = function() {
    lst <- processed_data_list()
    n_rows <- length(lst)
    row_height_px <- 50
    (n_rows * row_height_px) + 500
  })
#---------------------Results Tab-------------------------
  ##----Data Completeness----
  output$completeness_plot <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    spectra_cols <- data_info_r() %>%
      filter(final_name == "SPECTRA") %>%          #Should be SPECTRA, but cannot currently do because 0 is missing in PEAKS, while 0 is existing in Fragpipe and NA is missing here.
      dplyr::select(-final_name) %>%   # all sample columns.
      unlist(recursive = TRUE, use.names = FALSE)
    
    plot_completeness(lst, spectra_cols, color = input$color_palette)
  })
  
  output$completeness_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    spectra_cols <- data_info_r() %>%
      filter(final_name == "SPECTRA") %>%
      dplyr::select(-final_name) %>%   # all sample columns.
      unlist(recursive = TRUE, use.names = FALSE)
    
    plot_completeness(lst, spectra_cols, percent = TRUE, color = input$color_palette)
  })
  
  ##----Upset----
  output$upset_plot <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    plot_upset(lst)
  })
  
  ## ----Pairwise comparison of shared peptides----
  output$Pairwise_shared_peptide_matrix <- renderPlot({
    lst <- processed_data_list()
    req(lst, input$shared_mode)
    
    ht <- plot_shared_peptide(
      lst,
      mode  = input$shared_mode,
      color = input$color_palette,
      percent_type = "union"
    )
    
    draw(ht)
  })
  
  ## ----Pairwise comparison of shared peptides quantity----
  output$pairwise_peptide_quant_correlation <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    ht <- plot_pairwise_peptide_quant_correlation(lst, color = input$color_palette, cluster = input$cluster_mode)
    ComplexHeatmap::draw(ht)
  })
  
  ## ----PCA plot----
  output$pca <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    plot_PCA(lst, color = input$color_palette)
  })
  
  ## ----Number of peptides and peptidoforms----
  output$summary_peptides_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    plot_unique_counts(lst, column = "STRIPPED", y_label = "Number of unique peptides", color = input$color_palette)
  })
  
  output$summary_peptidoforms_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    plot_unique_counts(lst, column = "PEPTIDE", y_label = "Number of unique peptidoforms", color = input$color_palette)
  })
  
  output$summary_proteins_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    # Use extract_protein_prefixes for proteins
    plot_unique_counts(lst, column = "PROTEIN", y_label = "Number of unique proteins",
                       transform_fn = extract_protein_prefixes, color = input$color_palette)
  })
  
  ## ----Charge / Mass / mz / RT / ppm----
  output$charge_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "CHARGE" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No CHARGE column in data.", cex = 1.2)
      return(invisible())
    }
    
    
    plot_stacked_bar(lst, column = "CHARGE", fill_label = "Charge", percentage = FALSE, color = input$color_palette)
  })
  
  # Example usage for your Shiny outputs
  output$mass_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MASS" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MASS column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_density(lst, column = "MASS", x_label = "Mass (Da)", color = input$color_palette)
  })
  
  output$mz_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "MZ" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No MZ column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_density(lst, column = "MZ", x_label = "m/z", color = input$color_palette)
  })
  
  output$RT_plot2 <- renderUI({
    lst <- processed_data_list()
    req(lst)
    
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
    req(lst)
    
    for (sample_name in names(lst)) {
      
      local({
        sample_local <- sample_name
        df_local <- lst[[sample_local]]
        
        output[[paste0("RT2_", sample_local)]] <- renderPlot({
          plot_histogram(df = df_local, column = "RT", x_label = "Retention Time (min)", 
                         title_name = paste("RT Histogram –", sample_local), color = input$color_palette
          )
        })
      })
    }
  })
  
  output$ppm_plot2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "PPM" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No PPM column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_density(lst, column = "PPM", x_label = "ppm", color = input$color_palette)
  })
  
  ##----Score distribution----
  output$score_violin2 <- renderPlot({
    lst <- processed_data_list()
    req(lst)
    
    if (!any(vapply(lst, function(df) "SCORE" %in% colnames(df), logical(1)))) {
      plot.new()
      text(0.5, 0.5,"No SCORE column in data.", cex = 1.2)
      return(invisible())
    }
    
    plot_violin(lst, column = "SCORE", color = input$color_palette)
  })
  
#-------------------Group Comparison------------------------
  ## ----Render the UI for group assignment----
  output$group_assign_ui <- renderUI({
    req(names(active_data_list()))
    n <- input$n_groups
    
    # Wrap everything in a tagList so Shiny renders the list properly
    tagList(
      lapply(seq_len(n), function(i) {
        tagList(
          textInput(
            inputId = paste0("group_name_", i),
            label   = paste("Group", i, "name"),
            value   = paste("Group", i)  # default name
          ),
          selectInput(
            inputId = paste0("group_", i),
            label   = paste("Select samples for", paste0("Group ", i)),
            choices = names(active_data_list()),
            multiple = TRUE
          ),
          tags$hr()
        )
      })
    )
  })
  
  # Update group list
  group_list <- eventReactive(input$update_group_comp, {
    
    n <- input$n_groups
    groups <- list()
    empty_groups <- c()
    
    for (i in seq_len(n)) {
      custom_name <- input[[paste0("group_name_", i)]]
      samples     <- input[[paste0("group_", i)]]
      
      # Check if the group is empty or has no name
      if (is.null(custom_name) || nchar(custom_name) == 0 || is.null(samples) || length(samples) == 0) {
        empty_groups <- c(empty_groups, paste("Group", i))
      } else {
        groups[[custom_name]] <- samples
      }
    }
    
    # If any group is empty → show error and return NULL
    if (length(empty_groups) > 0) {
      showNotification(
        paste(
          "The following group(s) are empty or have no name and must be filled before continuing:", 
          paste(empty_groups, collapse = ", ")
        ),
        type = "error",
        duration = NULL
      )
      return(NULL)
    }
    
    # If all groups are valid → success message
    showNotification("Groups updated successfully!", type = "message")
    
    groups
  })
  
  ## ----Group Comparison----
  group_peptide_sets <- reactive({
    lst <- processed_data_list()
    groups <- group_list()
    pep_col <- "PEPTIDE"
    
    lapply(names(groups), function(g) {
      
      sample_names <- groups[[g]]
      
      peptide_counts <- table(unlist(lapply(sample_names, function(s) {
        df <- lst[[s]]
        if (!pep_col %in% colnames(df)) return(character(0))
        unique(df[[pep_col]])
      })))
      
      n_samples <- length(sample_names)
      
      names(peptide_counts[
        peptide_counts / n_samples >= input$min_presence_fraction
      ])
    }) |> setNames(names(groups))
  })
  
  output$group_venn_plot <- renderPlot({
    
    sets <- group_peptide_sets()
    req(length(sets) >= 2)
    
    n_groups <- length(sets)
    
    p <- ggVennDiagram::ggVennDiagram(
      sets,
      label_alpha = 0
    ) +
      ggplot2::theme_void()  +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = margin(10,10,10,10)
      )
    
    if (input$color_palette != "default") {
      p <- p + ggplot2::scale_fill_viridis_c(
        option = input$color_palette
      )  
    } else {
      p <- p + scale_fill_distiller(palette = "RdBu")
    }
    
    p
  })
  
  output$group_peptide_heatmap <- renderPlot({
    lst <- processed_data_list()
    groups <- group_list()
    
    quantity_cols <- data_info_r() %>%
      filter(final_name == "QUANTITY") %>%
      dplyr::select(-final_name) %>% 
      unlist(recursive = TRUE, use.names = FALSE)
    
    pep_mat <- prepare_peptide_matrix(lst, groups, quantity_cols, group_peptide_sets())
    
    req(nrow(pep_mat) > 0)
    
    plot_heatmap(pep_mat, color = input$color_palette, transpose = TRUE, log_transform = TRUE)
  })
  
  ## ----Group statistical analysis----
  group_comp_data <- eventReactive(input$update_group_comp, {
    lst <- processed_data_list()
    groups <- group_list()
    #Get all columns containing the quantity info.
    quantity_cols <- data_info_r() %>%
      filter(final_name == "QUANTITY") %>%
      dplyr::select(-final_name) %>%   # all sample columns.
      unlist(recursive = TRUE, use.names = FALSE)

    # Keep only non-empty groups
    groups <- groups[sapply(groups, function(g) {
      any(sapply(g, function(s) nrow(lst[[s]]) > 0))
    })]
    req(length(groups) > 1) # at least 2 groups needed
    
    # Generate all unique pairwise combinations
    group_pairs <- combn(names(groups), 2, simplify = FALSE)
    
    # Compute volcano data for each pair
    pairwise_volcano <- lapply(group_pairs, function(pair) {
      g1 <- pair[1]
      g2 <- pair[2]
      
      # Peptide column
      pep_col <- "PEPTIDE" #WIP: need to think how to solve the PTM problem? Do I just add the same peptidoform together?
      keep_cols <- c(pep_col, quantity_cols, "PROTEIN")
      
      #Here we filter based on union-intersect criteria
      allowed_peptides_g1 <- group_peptide_sets()[[g1]]
      allowed_peptides_g2 <- group_peptide_sets()[[g2]]
      
      # Combine data for each group
      df_g1 <- dplyr::bind_rows(lapply(groups[[g1]], function(s) {
        df <- lst[[s]]
        cols <- intersect(keep_cols, colnames(df))
        if (length(cols) < 2) return(NULL)  # need peptide column + ≥1 quantity column
        df[, cols, drop = FALSE]
        df[df[["PEPTIDE"]] %in% allowed_peptides_g1, , drop = FALSE]
      }))
      
      df_g2 <- dplyr::bind_rows(lapply(groups[[g2]], function(s) {
        df <- lst[[s]]
        cols <- intersect(keep_cols, colnames(df))
        if (length(cols) < 2) return(NULL)  # need peptide column + ≥1 quantity column
        df[, cols, drop = FALSE]
        df[df[["PEPTIDE"]] %in% allowed_peptides_g2, , drop = FALSE]
      }))
      
      # Skip pair if either group is empty
      if (is.null(df_g1) || is.null(df_g2) ||
          nrow(df_g1) == 0 || nrow(df_g2) == 0) {
        return(NULL)
      }
      
      df_long <- bind_rows(
        df_g1 %>%
          pivot_longer(
            cols = any_of(quantity_cols),
            names_to = "Sample",
            values_to = "Quantity"
          ) %>%
          mutate(Group = g1),
        
        df_g2 %>%
          pivot_longer(
            cols = any_of(quantity_cols),
            names_to = "Sample",
            values_to = "Quantity"
          ) %>%
          mutate(Group = g2)
      )
      
      #safe_mean <- function(x) if(length(x) > 0) mean(x, na.rm = TRUE) else NA_real_
      safe_ttest <- function(x, y) {
        x <- x[!is.na(x)]
        y <- y[!is.na(y)]
        
        # Not enough data
        if (length(x) < 2 || length(y) < 2) {
          return(NA_real_)
        }
        
        sx <- sd(x)
        sy <- sd(y)
        
        # Zero variance
        if (is.na(sx) || is.na(sy) || (sx == 0 && sy == 0)) {
          return(NA_real_)
        }
        
        tryCatch(
          t.test(x, y)$p.value,
          error = function(e) NA_real_
        )
      }

      volcano_df <- df_long %>%
        group_by(.data[[pep_col]]) %>%
        summarise(
          PROTEIN = dplyr::first(PROTEIN),
          Mean_G1 = mean(Quantity[Group == g1], na.rm = TRUE),
          Mean_G2 = mean(Quantity[Group == g2], na.rm = TRUE),
          log2FC = log2(Mean_G2 + 1) - log2(Mean_G1 + 1),
          pval = safe_ttest(
            Quantity[Group == g2],
            Quantity[Group == g1]
          ),
          .groups = "drop"
        ) %>%
        ungroup() %>%  # important before applying p.adjust
        mutate(
          adj_pval_BH = p.adjust(pval, method = "BH"),          # Benjamini-Hochberg FDR
          adj_pval_Bonf = p.adjust(pval, method = "bonferroni"),# Bonferroni
          negLog10P = -log10(pval),
          negLog10AdjP_BH = -log10(adj_pval_BH),
          negLog10AdjP_Bonf = -log10(adj_pval_Bonf)
        ) %>% #this is for MAplot
        mutate(
          A = 0.5 * (log2(Mean_G1 + 1) + log2(Mean_G2 + 1)),
        )
      volcano_df
    })
    
    names(pairwise_volcano) <- sapply(group_pairs, function(pair) paste(pair, collapse = "_vs_"))
    
    #for developmental purpose
    temp <<- pairwise_volcano[sapply(pairwise_volcano, nrow) > 0]
    
    pairwise_volcano[sapply(pairwise_volcano, nrow) > 0]
  })
  
  # Render the volcano tabset UI
  output$volcano_tabs <- renderUI({
    volcano_list <- group_comp_data()
    req(volcano_list)
    
    if (length(volcano_list) == 0) {
      return(tags$div(
        style = "color:red; font-weight:bold; padding:20px;",
        "Not enough data to compute any group comparisons."
      ))
    }
    
    output$volcano_comparison_tabs <- renderUI({
      req(volcano_list)  # make sure the list exists
      
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
            column(6, h6("GO-term enrichment"), plotOutput(paste0("go_term_", name))),
            column(6, h6("STRING-DB network"), plotOutput(paste0("STRING_", name)))
          )
        )
      })
      
      # Generate the tabsetPanel from the list
      do.call(tabsetPanel, c(id = "volcano_comparison_tabs", comparison_tabs))
    })
    
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
          
          # Check if groups have data
          if (is.null(df) || nrow(df) == 0 ||
              all(is.na(df$log2FC)) || all(is.na(df$negLog10AdjP_BH))) {
            
            return(
              plot_ly() %>%
                layout(
                  annotations = list(
                    x = 0.5, y = 0.5, xref = "paper", yref = "paper",
                    showarrow = FALSE,
                    text = paste("Not enough data for comparison:", plot_name),
                    font = list(size = 16, color = "red")
                  ),
                  xaxis = list(visible = FALSE),
                  yaxis = list(visible = FALSE)
                )
            )
          }
          
          sel <- selected_peptide()
          if (is.null(sel)) sel <- NA_character_  # <- avoids length 0
          
          # Add 'selected' column
          volc_df <- df %>%
            mutate(selected = !is.na(sel) & PEPTIDE == sel)
          
          # Split by significance for plotting
          ns_df        <- volc_df %>% filter(Significance == "Not significant")
          sig_df       <- volc_df %>% filter(Significance == "Significant", !selected)
          selected_df  <- volc_df %>% filter(selected)
          
          #Axis setting
          safe_max_neglogp <- suppressWarnings(max(df$negLog10AdjP_BH, na.rm = TRUE))
          safe_min_fc      <- suppressWarnings(min(df$log2FC,   na.rm = TRUE))
          safe_max_fc      <- suppressWarnings(max(df$log2FC,   na.rm = TRUE))
          
          if (!is.finite(safe_max_neglogp)) safe_max_neglogp <- 1
          if (!is.finite(safe_min_fc))      safe_min_fc <- -1
          if (!is.finite(safe_max_fc))      safe_max_fc <- 1
          
          #Plotting
          plot_ly(source = "group_diff") %>%
            add_trace(
              data = ns_df, x = ~log2FC, y = ~negLog10AdjP_BH,
              type = "scatter", mode = "markers",
              name = "Not significant",
              marker = list(color = "grey40", size = 6),
              hoverinfo = "none",
              inherit = FALSE
            ) %>%
            add_trace(
              data = sig_df, x = ~log2FC, y = ~negLog10AdjP_BH,
              type = "scatter", mode = "markers",
              name = "Significant",
              marker = list(color = "red", size = 8),
              text = ~PEPTIDE,
              key = ~PEPTIDE,
              hoverinfo = "text+x+y"
            ) %>%
            add_trace(
              data = selected_df,
              x = ~log2FC, y = ~negLog10AdjP_BH,
              type = "scatter", mode = "markers",
              name = "Selected peptide",
              marker = list(color = "lightgreen", size = 12),
              text = ~PEPTIDE,
              key = ~PEPTIDE,
              hoverinfo = "text+x+y"
            ) %>%
            layout(
              title = paste0("volcano plot: ",plot_name),
              xaxis = list(title = "log2 Fold Change"),
              yaxis = list(title = "-log10 p-value"),
              shapes = list(
                # fold change vertical lines
                list(type = "line", x0 = -1, x1 = -1, y0 = 0, y1 = safe_max_neglogp,
                     line = list(dash = "dash", color = "grey")),
                list(type = "line", x0 =  1, x1 =  1, y0 = 0, y1 = safe_max_neglogp,
                     line = list(dash = "dash", color = "grey")),
                # p-value cutoff
                list(type = "line", x0 = safe_min_fc, x1 = safe_max_fc,
                     y0 = 1.3, y1 = 1.3,
                     line = list(dash = "dash", color = "grey"))
              )
            )
        })
        
        ## ----MA plot----
        output[[paste0("ma_", plot_name)]] <- renderPlotly({
          
          if (is.null(df) || nrow(df) == 0 || all(is.na(df$A))) {
            return(plot_ly())
          }
          
          sel <- selected_peptide()
          if (is.null(sel)) sel <- NA_character_  # <- avoids length 0
          
          # Add 'selected' column
          ma_df <- df %>%
            mutate(selected = !is.na(sel) & PEPTIDE == sel)
          
          # Split by significance for plotting
          ns_df        <- ma_df %>% filter(Significance == "Not significant")
          sig_df       <- ma_df %>% filter(Significance == "Significant", !selected)
          selected_df  <- ma_df %>% filter(selected)
          
          plot_ly(source = "group_diff") %>%
            add_trace(
              data = ns_df,
              x = ~A, y = ~log2FC,
              type = "scatter", mode = "markers",
              name = "Not significant",
              marker = list(color = "grey40", size = 6),
              hoverinfo = "none",
              inherit = FALSE
            ) %>%
            add_trace(
              data = sig_df,
              x = ~A, y = ~log2FC,
              type = "scatter", mode = "markers",
              name = "Significant",
              marker = list(color = "red", size = 8),
              text = ~PEPTIDE,
              key = ~PEPTIDE,
              hoverinfo = "text+x+y"
            ) %>%
            add_trace(
              data = selected_df,
              x = ~A, y = ~log2FC,
              type = "scatter", mode = "markers",
              name = "Selected peptide",
              marker = list(color = "lightgreen", size = 12),
              text = ~PEPTIDE,
              key = ~PEPTIDE,
              hoverinfo = "text+x+y"
            ) %>%
            layout(
              title = paste("MA plot:", plot_name),
              xaxis = list(title = "Mean abundance (A)"),
              yaxis = list(title = "log2 Fold Change (M)"),
              shapes = list(
                list(type = "line", x0 = min(df$A, na.rm = TRUE),
                     x1 = max(df$A, na.rm = TRUE),
                     y0 = 0, y1 = 0,
                     line = list(dash = "dash", color = "grey"))
              )
            )
        })
        
        
        ## ----P-value histogram----
        output[[paste0("pval_hist_", plot_name)]] <- renderPlotly({
          
          if (is.null(df) || nrow(df) == 0 || all(is.na(df$negLog10AdjP_BH))) {
            return(plot_ly())
          }
          
          pvals <- 10^(-df$negLog10AdjP_BH)
          pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
          
          if (length(pvals) == 0) {
            return(plot_ly())
          }
          
          plot_ly(
            x = pvals,
            type = "histogram",
            nbinsx = 50,
            marker = list(color = "grey40"),
            hoverinfo = "x+y"
          ) %>%
            layout(
              title = paste("P-value distribution:", plot_name),
              xaxis = list(title = "Adjusted p-value", range = c(0, 1)),
              yaxis = list(title = "Count", type = "log"),
              shapes = list(
                list(
                  type = "line",
                  x0 = 0.05, x1 = 0.05,
                  y0 = 0, y1 = 1,
                  xref = "x",
                  yref = "paper",
                  line = list(color = "red", dash = "dash")
                )
              ),
              annotations = list(
                list(
                  x = 0.05, y = 1,
                  xref = "x", yref = "paper",
                  text = "0.05",
                  showarrow = FALSE,
                  xanchor = "left",
                  font = list(color = "red")
                )
              )
            )
          
        })
        
        ## ----Ranked Fold Change----
        output[[paste0("rank_fc_", plot_name)]] <- renderPlotly({
          
          if (is.null(df) || nrow(df) == 0 || all(is.na(df$log2FC))) {
            return(plot_ly())
          }
          
          sel <- selected_peptide()
          if (is.null(sel)) sel <- NA_character_
          
          rank_df <- df %>%
            filter(!is.na(log2FC)) %>%
            arrange(log2FC) %>%          # ascending; use desc(log2FC) if you prefer
            mutate(
              rank = row_number(),
              selected = !is.na(sel) & PEPTIDE == sel
            )
          
          selected_df <- rank_df %>% filter(selected)
          rest_df     <- rank_df %>% filter(!selected)
          
          plot_ly(source = "group_diff") %>%
            
            # All peptides
            add_trace(
              data = rest_df,
              x = ~rank,
              y = ~log2FC,
              type = "scatter",
              mode = "markers",
              name = "Peptides",
              marker = list(color = "grey40", size = 6),
              hoverinfo = "none"
            ) %>%
            
            # Selected peptide
            add_trace(
              data = selected_df,
              x = ~rank,
              y = ~log2FC,
              type = "scatter",
              mode = "markers",
              name = "Selected peptide",
              marker = list(color = "lightgreen", size = 12),
              text = ~PEPTIDE,
              key = ~PEPTIDE,
              hoverinfo = "text+x+y"
            ) %>%
            
            layout(
              title = paste("Ranked fold change:", plot_name),
              xaxis = list(title = "Rank (based on Fold Change)"),
              yaxis = list(title = "log2 Fold Change"),
              shapes = list(
                list(
                  type = "line",
                  x0 = 0,
                  x1 = max(rank_df$rank),
                  y0 = 0,
                  y1 = 0,
                  line = list(dash = "dash", color = "grey")
                )
              )
            )
        })
        ## ----Peptide Fold Change table----
        output[[paste0("peptide_table_", plot_name)]] <- DT::renderDT({
          
          if (is.null(df) || nrow(df) == 0) {
            return(
              DT::datatable(
                data.frame(Message = "No data available"),
                options = list(dom = "t")
              )
            )
          }
          
          table_df <- df %>%
            filter(!is.na(log2FC)) %>%
            dplyr::select(
              PEPTIDE,
              log2FC,
              adj_pval_BH,
              adj_pval_Bonf,
              Significance
            ) %>%
            arrange(adj_pval_BH)
          
          DT::datatable(
            table_df,
            rownames = FALSE,
            selection = "none",
            options = list(
              pageLength = 10,
              scrollX = TRUE,
              order = list(list(2, "asc"))
            )
          ) %>%
            DT::formatRound(c("log2FC", "adj_pval_BH","adj_pval_Bonf"), digits = 3)
        })
        
        ## ----GO term----
        output[[paste0("go_term_", plot_name)]] <- renderPlot({
          if (is.null(df) || nrow(df) == 0) {
            plot.new()
            text(0.5, 0.5, "No data available")
            return()
          }
          
          df <- df %>%
            filter(!is.na(log2FC)) %>%
            filter(Significance == "Significant") %>%
            dplyr::select(PEPTIDE, Significance, PROTEIN)
          
          uni_ids <- df$PROTEIN %>%
            strsplit(";") %>%                # split multiple proteins
            lapply(function(x) sapply(strsplit(x, "\\|"), `[`, 1)) %>%  # take first part of each
            unlist() %>%
            unique()
          
          if (is.null(uni_ids) || nrow(as.data.frame(uni_ids)) == 0) {
            plot.new()
            text(0.5, 0.5, "No significant IDs")
            return()
          }
          
          gene_map <- bitr(uni_ids, fromType="UNIPROT", toType="ENTREZID", OrgDb=org.Hs.eg.db)
          if (nrow(gene_map) == 0) {
            plot.new()
            text(0.5, 0.5, "No valid UniProt->Entrez mapping")
            return()
          }
          
          ego <- enrichGO(
            gene = gene_map$ENTREZID,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            universe = NULL,
            readable = TRUE
          )
          
          if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
            plot.new()
            text(0.5, 0.5, "No significant GO terms")
            return()
          }
          
          barplot(ego, showCategory = 10)
        })
        
        ## ----STRING-DB----
        output[[paste0("STRING_", plot_name)]] <- renderPlot({
          if (is.null(df) || nrow(df) == 0) {
            plot.new()
            text(0.5, 0.5, "No data available")
            return()
          }
          
          df <- df %>%
            filter(!is.na(log2FC)) %>%
            filter(Significance == "Significant") %>%
            dplyr::select(PEPTIDE, Significance, PROTEIN)
          
          uni_ids <- df$PROTEIN %>%
            strsplit(";") %>%                # split multiple proteins
            lapply(function(x) sapply(strsplit(x, "\\|"), `[`, 1)) %>%  # take first part of each
            unlist() %>%
            unique()
          
          
          if (is.null(uni_ids) || nrow(as.data.frame(uni_ids)) == 0) {
            plot.new()
            text(0.5, 0.5, "No significant IDs")
            return()
          }
          
          url <- paste0(
            "https://string-db.org/api/json/network?",
            "identifiers=", paste(uni_ids, collapse = "%0d"),
            "&species=", 9606 # human
          )
          
          res <- GET(url)
          
          if (http_status(res)$category != "Success") {
            plot.new()
            text(0.5, 0.5,
                 paste0("STRING request failed (HTTP ", res$status_code, ")"),
                 cex = 1.2)
            return()
          }
          
          # 2) Try to parse JSON safely
          data <- tryCatch(
            {
              fromJSON(content(res, "text", encoding = "UTF-8"))
            },
            error = function(e) {
              plot.new()
              text(0.5, 0.5,
                   paste0("JSON parse error:\n", e$message),
                   cex = 0.9)
              return(NULL)
            }
          )
          
          # 3) If parsing failed, stop here
          if (is.null(data)) {
            plot.new()
            text(0.5, 0.5,
                 paste0("no STRING-DB result", res$status_code, ")"),
                 cex = 1.2)
            return()
          }
          
          g <- graph_from_data_frame(
            data[, c("preferredName_A", "preferredName_B", "score")],
            directed = FALSE
          )
          
          ggraph(g, layout = "fr") +
            geom_edge_link(aes(width = score), alpha = 0.8) +
            geom_node_point(size = 5, color = "steelblue") +
            geom_node_text(aes(label = name), repel = TRUE) +
            theme_void()
          
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
    req(lst)
    
    lst <- lapply(lst, function(df) {
      
      if (!"PTM" %in% colnames(df)) {
        return(NULL)
      }
      
      df %>%
        mutate(
          PTM = ifelse(is.na(PTM) | PTM == "", "Unmodified", PTM)
        ) %>%
        tidyr::separate_rows(PTM, sep = "\\s*[,;]\\s*") %>%
        mutate(
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
    
    # Check if mod_map is NULL or empty
    if (is.null(mod_map_list) || length(mod_map_list) == 0 || all(is.na(mod_map_list))) {
      # Return a simple message in the UI
      return(tags$div("No PTM_pseudo sequence"))
    }
    
    lengths <- 7:11
    
    tabsetPanel(
      id = "motif_length_tabs",
      !!!lapply(lengths, function(L) {
        
        tabPanel(
          paste("Length", L),
          
          layout_column_wrap(
            width = "250px",  # each plot gets ~250px width
            !!!lapply(names(lst), function(sample_name) {
              plotOutput(
                paste0("PTM_motif_", sample_name, "_", L),
                height = "180px"
              )
            })
          )
        )
        
      })
    )
  })
  
  observe({
    lst <- processed_data_list()
    req(lst)
    
    #Create custom namespace for ggseqplot
    AA_symbols <- LETTERS            # standard amino acids
    PTM_symbols <- c(as.character(1:9), letters)
    namespace <- c(AA_symbols, PTM_symbols)
    
    lengths <- 7:11
    
    for (L in lengths) {
      for (sample_name in names(lst)) {
        
        local({
          length_val <- L
          sample_val <- sample_name
          
          output[[paste0("PTM_motif_", sample_val, "_", length_val)]] <- renderPlot({
            
            df <- lst[[sample_val]]
            req(df)
            
            # subset by length
            df_L <- df[df$LENGTH == length_val, ]
            df_L <- df[df$PTM != "", ]
            
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
    req(lst)
    req(cache)
    req(input$HLA_alleles)
    
    alleles_vec <- input$HLA_alleles
    
    # Extract peptides of correct length
    peptides <- unique(unlist(lapply(lst, `[[`, "STRIPPED"), use.names = FALSE))
    peptides <- peptides[nchar(peptides) %in% 8:11]
    req(length(peptides) > 0)
    
    # Path to netMHCpan
    netmhcpan_path <- "/mnt/c/Users/Yannic/netMHCpan-4.2/netMHCpan" # user-defined
    
    withProgress(message = "Running netMHCpan predictions...", value = 0, {
      for (al in alleles_vec) {
        
        al_conversion <- sub("-", "\\.", al)
        al_conversion <- sub(":", "", al_conversion)
        print(al_conversion)
        
        # Determine which peptides need prediction
        if (!(al_conversion %in% colnames(cache)[-1])) {
          peptides_to_predict <- peptides       # new allele → predict all peptides
        } else {
          peptides_to_predict <- cache$Peptide[is.na(cache[[al_conversion]])]  # existing allele → only new peptides
          peptides_to_predict <- peptides_to_predict[nchar(peptides_to_predict) %in% 8:11]
        }
        incProgress(1 / length(alleles_vec), detail = paste("Predicting for allele", al, " (# of peptides: ", length(peptides_to_predict), ")" ))
        
        if (length(peptides_to_predict) == 0) next
        
        # Temp files
        peptide_file <- tempfile(fileext = ".txt")
        output_file  <- tempfile(fileext = ".txt")
        writeLines(peptides_to_predict, peptide_file)
        
        # Convert Windows paths to WSL paths
        peptide_wsl <- trimws(system2("wsl", c("wslpath", "-a", shQuote(peptide_file)), stdout = TRUE))
        out_wsl     <- trimws(system2("wsl", c("wslpath", "-a", shQuote(output_file)), stdout = TRUE))
        
        # Build and run netMHCpan command
        cmd <- paste(
          shQuote(netmhcpan_path),
          "-p", shQuote(peptide_wsl),
          "-a", shQuote(al),
          "-l 8,9,10,11",
          "-xls",
          "-xlsfile", shQuote(out_wsl)
        )
        system2("wsl", c("bash", "--login", "-c", shQuote(cmd)),stdout = NULL)
        
        # Read and clean output
        res <- read.table(output_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        
        # --- Rename columns ---
        header1 <- colnames(res)
        header2 <- as.character(unlist(res[1, ]))
        colnames_new <- header1
        colnames_new[2] <- "Peptide"
        
        current_hla <- NULL
        for (i in seq_along(colnames_new)) {
          if (grepl("^HLA", header1[i])) current_hla <- header1[i]
          if (!is.null(current_hla) && header2[i] != "") colnames_new[i] <- paste0(current_hla, "_", header2[i])
        }
        colnames(res) <- colnames_new
        res <- res[-1, ]  # remove header row
        
        # Keep only Peptide + Rank columns
        keep_cols <- c("Peptide", grep("_Rank$", colnames(res), value = TRUE))
        res <- res[, keep_cols, drop = FALSE]
        
        # Clean column names
        colnames(res) <- gsub("_Rank$", "", colnames(res))
        colnames(res)[-1] <- sub("^([^.]+\\.[^.]+)\\.", "\\1", colnames(res)[-1])
        
        # Convert numeric columns
        res[-1] <- lapply(res[-1], as.numeric)
        res <<- res
        
        
        if (!(al_conversion %in% colnames(cache)[-1])) {
          cache <- left_join(cache, res, by = "Peptide")
        } else {
          cache <- full_join(cache, res, by = "Peptide") %>%
            mutate(
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
  
  
  peptide_wide_all <- reactive({
      lst <- processed_data_list()
      cache <- prediction_cache()
      req(lst)
      
      # Check if "netMHCpan" exists in the list
      if (ncol(cache) <= 1) {
        data.frame(Message = "No binding predictions available. No netMHCpan precomputed data available and netMHCpan has not ran yet")
        req(FALSE)  # Stops this reactive, downstream reactives won't run
      }
      
      lst <- lapply(lst, function(df) {
        
        # Remove existing HLA columns
        allele_cols <- grep("^HLA", colnames(df), value = TRUE)
        cols_to_keep <- setdiff(colnames(df), allele_cols)
        df <- df[, cols_to_keep]
        
        # Left join with new predictions
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
      
      temp <<- out
      
      dplyr::bind_rows(out)
  })
  
  peptide_wide_unique <- reactive({
    df <- peptide_wide_all()
    req(df)
    
    allele_cols <- grep("^HLA", colnames(df), value = TRUE)
    
    df %>%
      dplyr::group_by(Set, STRIPPED) %>%
      dplyr::summarise(dplyr::across(all_of(allele_cols), ~ if(all(is.na(.x))) {NA_real_} else min(.x, na.rm = TRUE)),
                       .groups = "drop")
  })
  
  binder_summary_all <- reactive({
    df <- peptide_wide_unique()
    req(df)
    
    allele_cols <- grep("^HLA", colnames(df), value = TRUE)
    
    df %>%
      tidyr::pivot_longer(
        cols = all_of(allele_cols),
        names_to = "Allele",
        values_to = "Rank"
      ) %>%
      dplyr::mutate(
        Class = dplyr::case_when(
          is.na(Rank) ~ "Missing",
          Rank <= 0.5 ~ "Strong",
          Rank <= 2   ~ "Weak",
          TRUE        ~ "Non"
        ),
        Class = factor(Class, levels = c("Strong", "Weak", "Non", "Missing"))
      ) %>%
      dplyr::count(Set, Allele, Class) %>%
      tidyr::complete(
        Set,
        Allele,
        Class,
        fill = list(n = 0)
      ) %>%
      tidyr::pivot_wider(
        names_from = Class,
        values_from = n
      ) %>%
      dplyr::mutate(
        Total = Strong + Weak + Non + Missing,
        Strong_pct = round(100 * Strong / Total, 2),
        Weak_pct   = round(100 * Weak   / Total, 2),
        Non_pct    = round(100 * Non    / Total, 2),
        NA_pct     = round(100 * Missing / Total, 2)
      )
  })
  
  output$binding_summary <- DT::renderDT({
    if (!exists("netMHCpan")) {
      # Return a small placeholder table with a message
      data.frame(Message = "netMHCpan pre-generated data missing; analysis skipped.")
    } else {
      binder_summary_all()
    }
  })
  
  output$binding_plot_percent <- renderPlot({
    df <- peptide_wide_unique()
    req(df)
    
    if (!exists("netMHCpan")) {
      # Return a small placeholder table with a message
      return(data.frame(Message = "netMHCpan pre-generated data missing; analysis skipped."))
    } 
    
    plot_binders(df, color = input$color_palette, percent = TRUE)
  })
  
  output$binding_plot_absolute <- renderPlot({
    df <- peptide_wide_unique()
    req(df)
    
    if (!exists("netMHCpan")) {
      # Return a small placeholder table with a message
      return(data.frame(Message = "netMHCpan pre-generated data missing; analysis skipped."))
    } 
    
    plot_binders(df, color = input$color_palette, percent = FALSE)
  })
  
  output$binding_table <- DT::renderDT({
    df <- peptide_wide_unique()
    req(df)
    
    if (!exists("netMHCpan")) {
      # Return a small placeholder table with a message
      return(data.frame(Message = "netMHCpan pre-generated data missing; analysis skipped."))
    } 
    
    # Collapse rows by peptide
    df_collapsed <- df %>%
      group_by(STRIPPED) %>%
      summarise_all(~ {
        vals <- unique(.)
        vals <- vals[!is.na(vals)]
        if(length(vals) == 0) NA else paste(vals, collapse = ",")
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
    content = function(file) {
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
  search_results <- eventReactive(input$search_peptide, {
    req(input$peptide_query)
    query <- trimws(input$peptide_query)
    req(nchar(query) > 0)
    
    lst <- data_list_r()
    req(lst)
    
    keep_cols <- c("Sample", "PEPTIDE", "STRIPPED", "LENGTH" ,"MASS" ,"CHARGE", "MZ" ,"K0", "RT", "PROTEIN")
    
    filtered_list <- lapply(lst, function(df) {
      df_sub <- df[, intersect(keep_cols, colnames(df)), drop = FALSE]
      if (input$exact) df_sub[df_sub$STRIPPED == query, , drop = FALSE] else df_sub[grepl(query, df_sub$STRIPPED, ignore.case = TRUE), , drop = FALSE]
    })
    
    # combine all samples (optional – remove bind_rows if per-sample)
    df_combined <- dplyr::bind_rows(filtered_list, .id = "Sample")
    head(df_combined, 100)
  })
  
  output$peptide_table <- renderDT({
    df <- search_results()
    shiny::validate(shiny::need(nrow(df) > 0, "No matching peptides found"))
    
    DT::datatable(
      df,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollY = "60vh", #This should be good enough for the standard monitors. Otherwise need to make it responsive to browser window.
        scrollX = TRUE,
        dom = "Bfrtip"
      )
    )
  })
  
  unique_peptides_table <- reactive({
    lst <- data_list_r()
    req(lst)
    
    # Handle case when there is only one or no dataset
    if (length(lst) < 2) {
      return(
        tibble(
          Message = paste0(
            "Analysis is redundant, since less than 2 datasets are given."
          )
        )
      )
    }
    
    keep_cols <- c("Sample", "PEPTIDE", "STRIPPED", "LENGTH" ,"MASS" ,"CHARGE", "MZ" ,"K0", "RT", "PROTEIN")
    
    combined <- imap_dfr(lst, ~ dplyr::select(.x, intersect(keep_cols, colnames(.x))) %>%
                           mutate(Sample = .y))
    
    peptide_counts <- combined %>%
      count(STRIPPED, name = "n_datasets")
    
    combined %>%
      inner_join(
        peptide_counts %>% filter(n_datasets == 1),
        by = "STRIPPED"
      ) %>%
      dplyr::select(-n_datasets) %>%
      arrange(Sample, STRIPPED)
  })
  
  output$unique_peptide_table <- renderDT({
    df <- unique_peptides_table()
    shiny::validate(shiny::need(nrow(df) > 0, "No matching peptides found"))
    
    DT::datatable(
      df,
      filter = "top",
      options = list(
        pageLength = 25,
        scrollY = "60vh", #This should be good enough for the standard monitors. Otherwise need to make it responsive to browser window.
        scrollX = TRUE,
        dom = "Bfrtip"
      )
    )
  })
  
}
##---------------End------------

shinyApp(
  ui = ui,
  server = function(input, output, session) {
    server(input, output, session, preloaded_data = preloaded_data, generate_pseudo_sequence = FALSE)
  }
)