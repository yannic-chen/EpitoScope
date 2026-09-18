app_ui <- function() {bs4DashPage(
  title = "EpitoScope",
  help = TRUE, # automatically enable/disable all bs4Dash::tooltips and popover that are present in the shiny app
  fullscreen = FALSE, # an icon is displayed in the navbar to switch to full screen mode
  scrollToTop = FALSE, # allows to toggle the scroll to top button shown in the bottom right corner

  ## ---- Sidebar ----
  sidebar = bs4DashSidebar(
    title = "Filters",
    collapsed = F,
    minified = F, # Minified means a little part of the sidebar is still visible.
    expandOnHover = T, # when minified is TRUE, if this property is TRUE, the sidebar opens when hovering but re-collapses as soon as the focus is lost.
    elevation = 3,
    uiOutput("sample_selector"),
    selectInput("na_policy", "Missing values (NA) in filters:",
                choices = c("Include" = "include",
                            "Exclude" = "exclude",
                            "Treat as 0" = "zero"),
                selected = "exclude"),
    uiOutput("length_slider_ui"),
    uiOutput("quantity_slider_ui"),
    uiOutput("score_slider_ui"),
    uiOutput("charge_slider_ui"),
    uiOutput("mass_slider_ui"),
    uiOutput("RT_slider_ui"),
    actionButton("generate_report", tagList(icon("file-arrow-down"),"Generate HTML Report")),
    actionButton("open_export", "Export figure", icon = icon("camera"))
  ),

  ## ---- Header ----
  header = bs4DashNavbar(
    fixed = TRUE,
    title = dashboardBrand(
      title = "EpitoScope",
      color = "primary",
      href = "https://github.com/yannic-chen/EpitoScope",
      image = "www/Epitoscope.png"
    ),
    rightUi = tags$li(
      class = "nav-item dropdown d-flex align-items-center px-3",
      tags$div(
        style = "display: flex; align-items: center; gap: 6px;",
        tags$span(
          "",
          bs4Dash::tooltip(
            icon("info-circle"),
            title     = "How long (seconds) to wait after the last filter change before recomputing plots. Higher = less interruptions.",
            placement = "bottom"
          ),
          style = "color: rgba(0,0,0,1); font-size: 11px; white-space: nowrap;"
        ),
        numericInput(
          inputId = "debounce_delay_s",
          label   = NULL,
          value   = 2,
          min     = 0,
          max     = 30,
          step    = 0.5,
          width   = "70px"
        )
      )
    ),
    navbarMenu(
      id = "navmenu",
      navbarTab(tabName = "raw_summary", text = "RAW"),
      navbarTab(tabName = "qc", text = "QC"),
      navbarTab(tabName = "results", text = "Results"),
      navbarTab(tabName = "group_comp", text = "Groups"),
      navbarTab(tabName = "ptm", text = "PTM"),
      navbarTab(tabName = "binding_pred", text = "Binding"),
      navbarTab(tabName = "peptide_lookup", text = "Lookup"),
      navbarTab(
        text = "Consoles",
        dropdownHeader("Advanced functions"),
        navbarTab(tabName = "advanced_settings", text = "Settings"),
        navbarTab(tabName = "dev_console", text = "Console"),
        dropdownDivider(),
        navbarTab(
          text = "SQL",
          navbarTab(tabName = "SQL_console", text = "Save/Load"),
          navbarTab(tabName = "SQL_query",   text = "Query SQL"),
          navbarTab(tabName = "SQL_browser", text = "Browse SQL")
          )
        )
      )
  ),

  ## ---- Control Bar ----
  controlbar = bs4DashControlbar( #this is just an extra sidebar on the right
    #skinSelector(),
    #tags$hr(),
    selectInput(
      inputId = "color_palette",
      label = tagList(icon("palette"), "Color palette"),
      choices = c("default", "viridis", "magma", "inferno", "plasma", "cividis", "mako", "rocket", "turbo"),
      selected = "default"
    ),
    sliderInput("plot_font", "Plot font size", min = 8, max = 28, value = 13, step = 1),
    pinned = FALSE, value = ""
    ),

  ## ---- Main ----
  body = bs4DashBody(
    shinyjs::useShinyjs(), #needed to make button grey out
    tags$head(tags$style(HTML("
  html { scrollbar-gutter: stable; }
  .main-header.split-mode { background: #6f42c1 !important; }
"))),
    tags$script(HTML("
      document.addEventListener('keydown', function(e){
        if (e.key === 'Enter' &&
            (e.target.id === 'dynrange_search' || e.target.id === 'dynrange_pep_search')) {
          e.preventDefault();
          document.getElementById('dynrange_go').click();
        }
      });
    ")),
    plotExportJS <- tags$head(tags$script(HTML("
  // ---- existing downloadPlotly handler stays here ----
  Shiny.addCustomMessageHandler('downloadPlotly', function(msg){
    var gd = document.getElementById(msg.id);
    if(!gd) return;
    Plotly.downloadImage(gd, {
      format: msg.format, width: msg.width, height: msg.height,
      scale: msg.scale, filename: msg.filename
    });
  });

  // ---- new: resize plotly when its card maximizes/restores ----
  function resizeCardPlots(card){
    card.querySelectorAll('.plotly.html-widget').forEach(function(gd){
      // defer one frame so the card has its final size first
      window.requestAnimationFrame(function(){ Plotly.Plots.resize(gd); });
    });
  }
  // bs4Dash toggles the .maximized-card / .card-maximized class on the card element.
  var mo = new MutationObserver(function(muts){
    muts.forEach(function(m){
      if(m.attributeName === 'class'){
        resizeCardPlots(m.target);
      }
    });
  });
  document.addEventListener('DOMContentLoaded', function(){
    document.querySelectorAll('.card').forEach(function(c){
      mo.observe(c, { attributes: true });
    });
    // cards rendered later (renderUI) — observe on the fly
    var bodyMo = new MutationObserver(function(muts){
      muts.forEach(function(m){
        m.addedNodes.forEach(function(n){
          if(n.nodeType===1){
            if(n.classList && n.classList.contains('card')) mo.observe(n,{attributes:true});
            n.querySelectorAll && n.querySelectorAll('.card').forEach(function(c){ mo.observe(c,{attributes:true}); });
          }
        });
      });
    });
    bodyMo.observe(document.body, { childList: true, subtree: true });
  });
window.fitExportPaper = function(){
  var area  = document.querySelector('.exp-preview-area');
  var paper = document.getElementById('exp_paper');
  if(!area || !paper) return;
  var pxW = parseFloat(paper.dataset.pxw || '768');
  var pxH = parseFloat(paper.dataset.pxh || '576');
  paper.style.width  = pxW + 'px';
  paper.style.height = pxH + 'px';
  var s = Math.min(area.clientWidth / pxW, area.clientHeight / pxH);
  paper.style.transformOrigin = 'top left';
  paper.style.transform = 'scale(' + s + ')';
};
"))),
    bs4TabItems(
      ### ---- RAW summary ----
      bs4TabItem(
        tabName = "raw_summary",  # must match menuItem
        fluidRow(
          ### ---- Column Map ----
          bs4Card(title = tagList("Column Map",
                                  span(bs4Dash::tooltip(icon("info-circle"), title = "This table maps the columns of your dataset to the expected schema.
                                                    Note: m/z values are taken from the report and not calculated from the mass and charge column. These two values do differ.", placement = "right")
                                  )), width = 12, maximizable = TRUE,
                  #DT::DTOutput("summary_table")
                  div(style = 'overflow-x: auto;', DT::DTOutput("summary_table"))
                  ),

          ### ---- Annotation Table (if given) ----
          bs4Card(inputId = "annotation_card", title = tagList("Annotation Table", bs4Dash::tooltip(icon("info-circle"),"This is the original input annotation table.", placement = "right")
          ), width = 12, maximizable = TRUE,
          #DT::DTOutput("annotation_table")
          div(style = 'overflow-x: auto;', DT::DTOutput("annotation_table"))
          ),


          ### ---- Unique Entries ----
          bs4Card(title = "Unique entries", width = 12, maximizable = TRUE,
                  tabsetPanel(
                    tabPanel("Peptides",
                             h6("Unique Peptides (no PTMs)"),
                             plotOutput("summary_peptides_plot")
                             ),
                    tabPanel("Peptidoforms",
                             h6("Peptidoforms (including PTMs)"),
                             plotOutput("summary_peptidoforms_plot")
                             ),
                    tabPanel("Proteins",
                             h6("Proteins", bs4Dash::tooltip(icon("info-circle"),"Protein names is obtained from the Accession column and is truncated to the first space.", placement = "right")),
                             plotOutput("summary_proteins_plot")
                             )
                    )
                  ),

          ### ---- Other Numeric Columns ----
          bs4Card(title = "Other Numeric Columns", width = 12, maximizable = TRUE,
                  tabsetPanel(
                    tabPanel("Charge",
                             bs4Dash::tooltip(icon("info-circle"),"In case multiple charges are given, the minimum is taken.", placement = "right"),
                             plotOutput("charge_plot")
                             ),
                    tabPanel("Mass",
                             plotlyOutput("mass_plot")
                             ),
                    tabPanel("m/z",
                             plotlyOutput("mz_plot")
                             ),
                    tabPanel("RT",
                             uiOutput("RT_plot")
                             ),
                    tabPanel("Mass Error",
                             bs4Dash::tooltip(icon("info-circle"),"Either ppm (PEAKS) or delta Mass (Fragpipe). Ignore the x-axis label.", placement = "right"),
                             plotlyOutput("ppm_plot")
                             ),
                    tabPanel("score",
                             plotlyOutput("score_violin")
                             )
                    )
                  ),

          ### ---- Summary Table ----
          bs4Card(title = tagList("Summary Table", bs4Dash::tooltip(icon("info-circle"),"This is simply the summary() output. Each sample occupies one row.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  uiOutput("RAW_summary_html")
                  )
          )
        ),

      ## ---- QC ----
      bs4TabItem(
        tabName = "qc",  # must match menuItem
        fluidRow(

          ### ---- Peptide Length Distribution ----
          with_export(
            bs4Card(title = "Peptide Length Distribution", width = 12, maximizable = TRUE,
                    plotlyOutput("length_plot")),
            "length_distribution"
          ),
          ### ---- Length Range Percentage ----
          with_export(
            bs4Card(title = "Length Range Percentage", width = 12, maximizable = TRUE,
                    sliderInput("mhc_length_range", "length window:", min = 5, max = 30,
                                value = c(8, 13), step = 1),
                    tabsetPanel(
                      tabPanel("Per Sample",
                               div(`data-tab-key` = "length_range_distribution",
                                   plotOutput("length_range_percentage"))),
                    tabPanel("Per Measurement",
                             div(`data-tab-key` = "length_range_distribution_meas",
                                 plotOutput("length_range_percentage_meas")))
                    )),
            "length_range_distribution"
            ),

          ### ---- Other Numeric Columns ----
          with_export(
            bs4Card(title = "Other Numeric Columns", width = 12, maximizable = TRUE,
                    tabsetPanel(
                      tabPanel("Charge",
                               tabsetPanel(
                                 tabPanel("Per Sample",
                                          div(`data-tab-key` = "charge_plot",
                                              plotOutput("charge_plot2"))),
                                 tabPanel("Across Measurement",
                                          div(`data-tab-key` = "charge_plot_meas",
                                              plotlyOutput("charge_plot_meas"))
                                 )
                               )
                      ),

                      tabPanel("Mass",
                               tabsetPanel(
                                 tabPanel("Per Sample",
                                          div(`data-tab-key` = "mass_plot",
                                              plotlyOutput("mass_plot2"))),
                                 tabPanel("Across Measurement",
                                          h6("Density Curve with the line being the mean and the shaded band representing the range between measurements."),
                                          div(`data-tab-key` = "mass_plot_meas",
                                              plotlyOutput("mass_plot_meas"))
                                 )
                               )
                      ),
                      tabPanel("m/z",
                               tabsetPanel(
                                 tabPanel("Per Sample",
                                          div(`data-tab-key` = "mz_plot",
                                              plotlyOutput("mz_plot2"))),
                                 tabPanel("Across Measurement",
                                          h6("Density Curve with the line being the mean and the shaded band representing the range between measurements."),
                                          div(`data-tab-key` = "mz_plot_meas",
                                              plotlyOutput("mz_plot_meas"))
                                 )
                               )
                      ),
                      tabPanel("RT",
                               tabsetPanel(
                                 tabPanel("Per Sample",
                                          div(`data-tab-key` = "rt_histograms",
                                              uiOutput("RT_plot2"))),
                                 tabPanel("Across Measurement",
                                          h6("The light shaded band across the histogram represents the range between measurements."),
                                          div(`data-tab-key` = "rt_per_measurement",
                                              uiOutput("RT_plot_meas"))
                                 )
                               )
                      ),
                      tabPanel("Mass Error",
                               bs4Dash::tooltip(icon("info-circle"), "Either ppm (PEAKS) or delta Mass (Fragpipe).", placement = "right"),
                               tabsetPanel(
                                 tabPanel("Per Sample",
                                          div(`data-tab-key` = "ppm_plot",
                                              plotlyOutput("ppm_plot2"))),
                                 tabPanel("Across Measurement",
                                          h6("Density Curve with the line being the mean and the shaded band representing the range between measurements."),
                                          div(`data-tab-key` = "ppm_plot_meas",
                                              plotlyOutput("ppm_plot_meas"))
                                 )
                               )
                      ),
                      tabPanel("score",
                               tabsetPanel(
                                 tabPanel("Per Sample",
                                          div(`data-tab-key` = "score_plot",
                                              plotlyOutput("score_violin2"))),
                                 tabPanel("Across Measurement",
                                          h6("Line is mean distribution. THe dashed line is min and the shaded area is max."),
                                          div(`data-tab-key` = "score_plot_meas",
                                              plotlyOutput("score_violin_meas"))
                                 )
                               )
                      )
                    )),
            "charge_plot"
          ),
          ### ---- Motif Plot ----
          with_export(
            bs4Card(title = tagList("Motif Plot", bs4Dash::tooltip(icon("info-circle"),"Minimum of 5 sequences are required for Motif generation.", placement = "right")
                                    ), width = 12, maximizable = TRUE,
                    uiOutput("motif_tabs")),
            paste0("motif_len_", motif_plot_length[1])
            ),
          ### ---- Unique Entries ----
          with_export(
            bs4Card(title = "Identification distribution", width = 12, maximizable = TRUE,
                    tabsetPanel(
                      tabPanel("Peptides",
                               h6("Unique Peptides (no PTMs)", bs4Dash::tooltip(icon("info-circle"),"0s are considered identified but not quantified if NA also exist. Otherwise 0 is considered not identified.", placement = "right")),
                               div(`data-tab-key` = "peptides_unique_meas",
                                   plotOutput("summary_peptides_plot3"))
                      ),
                      tabPanel("Peptidoforms",
                               h6("Peptidoforms (including PTMs)", bs4Dash::tooltip(icon("info-circle"),"0s are considered identified but not quantified if NA also exist. Otherwise 0 is considered not identified.", placement = "right")),
                               div(`data-tab-key` = "peptidoforms_unique_meas",
                                   plotOutput("summary_peptidoforms_plot3"))
                      ),
                      tabPanel("Proteins",
                               h6("Proteins", bs4Dash::tooltip(icon("info-circle"),"Protein names is obtained from the Accession column and is truncated to the first space.", placement = "right")),
                               div(`data-tab-key` = "proteins_unique_meas",
                                   plotOutput("summary_proteins_plot3"))
                      )
                    )),
            "peptides_unique_meas"
          ),

          ### ---- Dynamic Range plot ----
          with_export(
          bs4Card(title = tagList("Dynamic Rang", bs4Dash::tooltip(icon("info-circle"),"Plot generation and highlighting can take long time with large number of samples", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  textInput("dynrange_search", "Highlight protein (regex supported):", placeholder = "HLA[ABC]"),
                  textInput("dynrange_pep_search","Highlight peptide (regex, stripped or peptidoform):",placeholder = "e.g. SLLQHLIGL|SINFKL"),
                  actionButton("dynrange_go", "Highlight", icon = icon("magnifying-glass")),
                  h6("Plots with Protein matches have hoverinfo for non-matches deactivated to allow better hovering over matches."),
                  tabsetPanel(
                    tabPanel("Individual",
                             div(`data-tab-key` = "dynrange_individual",
                             uiOutput("dynrange_individual_ui"))),
                    tabPanel("Combined",
                             radioButtons("dynrange_rank_mode", "Rank scale:",
                                          choices  = c("Absolute" = "absolute", "Relative (%)" = "relative"),
                                          selected = "absolute", inline = TRUE),
                             div(`data-tab-key` = "dynrange_combined",
                             plotlyOutput("dynrange_combined")))
                    )
                  ),
          "dynrange_individual"
          ),
          ### ---- 1/k0 vs mz ----
          with_export(
            bs4Card(title = "k0 vs mz", width = 12, maximizable = TRUE,
                    uiOutput("scatterplots_ui")),
            "mz_k0"
          ),
          ### ---- measurement specific heatmap ----
          with_export(
          bs4Card(title = "Measurement Specific Heatmap", width = 12, maximizable = TRUE,
                  selectInput("cluster_mode_ea", "Clustering:",
                              choices = c(
                                "Sample" = "sample",
                                "Rows only" = "rows",
                                "Columns only" = "columns",
                                "Rows + columns" = "both"
                              ), selected = "sample"

                  ),
                  plotlyOutput("measurement_heatmap", width = "auto", height = "650px")
                  ),
          "measurement_heatmap"
          )
          )
        ),

      ## ---- Results ----
      bs4TabItem(
        tabName = "results",  # must match menuItem
        fluidRow(
          ### ---- Unique Entries ----
          with_export(
            bs4Card(title = "Unique entries", width = 12, maximizable = TRUE,
                    tabsetPanel(
                      tabPanel("Peptides",
                               h6("Unique Peptides (no PTMs)"),
                               div(`data-tab-key` = "peptides_unique",
                                   plotOutput("summary_peptides_plot2"))
                      ),
                      tabPanel("Peptidoforms",
                               h6("Peptidoforms (including PTMs)"),
                               div(`data-tab-key` = "peptidoforms_unique",
                                   plotOutput("summary_peptidoforms_plot2"))
                      ),
                      tabPanel("Proteins",
                               h6("Proteins", bs4Dash::tooltip(icon("info-circle"),"Protein names is obtained from the Accession column and is truncated to the first space.", placement = "right")),
                               div(`data-tab-key` = "proteins_unique",
                                   plotOutput("summary_proteins_plot2"))
                      )
                    )),
            "peptides_unique"
          ),
          ### ---- Data Completeness ----
          with_export(
            bs4Card(title = tagList("Data Completeness", bs4Dash::tooltip(icon("info-circle"),"NA is used for missing/not identified. If no NA exist, then 0 will be used for missing/not identified", placement = "right")
            ), width = 12, maximizable = TRUE,
            h6("WARNING: For PEAKS 12 Studio, the column X.Spec has been replaced with X.Feature. X.Feature returns 0 even when a peptide has been identified but could not be quantified. X.Spec on the other hand only returns 0 if it is not identified at all."),
            tabsetPanel(
              # --- absolute ---
              tabPanel(
                "Absolute",
                div(`data-tab-key` = "completeness_count",
                    plotOutput("completeness_plot"))
              ),

              # --- percentage ---
              tabPanel(
                "Percentage",
                div(`data-tab-key` = "completeness_percent",
                    plotOutput("completeness_plot2"))
              )
            )),
            "completeness_count"
          ),
          ### ---- Upset Plot of Peptides ----
          with_export(
            bs4Card(title = tagList("Upset Plot of Peptides", bs4Dash::tooltip(icon("info-circle"),"Currently intersection min_size is 2% of combined number of unique peptides.", placement = "right")
            ), width = 12, maximizable = TRUE,
            fluidRow(
              column(4, numericInput("upset_min_size",    "Min. intersection size (percent):", value = 2,  min = 0, step = 0.1)),
              column(4, numericInput("upset_min_degree",  "Min. degree:",            value = 1,  min = 1, step = 1)),
              column(4, numericInput("upset_n_intersect", "Max. intersections:",     value = 40, min = 3, step = 1))
            ),
            plotOutput("upset_plot"),
            tags$hr(),
            helpText("Click a row below to explore the peptides in that intersection."),
            DT::DTOutput("upset_intersection_table")
            ),
            "upset_plot"
          ),
          ### ---- Pairwise shared peptide matrix ----
          with_export(
            bs4Card(title = tagList("Pairwise shared peptide matrix", bs4Dash::tooltip(icon("info-circle"),"Pairwise comparison of number of shared peptides. For Percent visualization union is used (i.e. jaccard style).", placement = "right")
            ), width = 12, maximizable = TRUE,
            selectInput("shared_mode", "Visualization:", choices = c("Count" = "count", "Percent" = "percent")),
            plotlyOutput("Pairwise_shared_peptide_matrix")
            ),
            "shared_peptide_matrix"
          ),
          ### ---- Pairwise shared peptide quantity comparison matrix ----
          with_export(
            bs4Card(title = tagList("Pairwise shared peptide quantity comparison matrix", bs4Dash::tooltip(icon("info-circle"),"Pairwise comparison of peptide max quantity of shared peptides, using pearson correlation. Clustering distance is 'euclidean' and method is 'complete'.
                                                                                                         Peptide label is removed if # > 200. NA is coloured Gray. Rastering is used for large data.", placement = "right")
            ), width = 12, maximizable = TRUE,
            selectInput("cluster_mode", "Clustering:",
                        choices = c(
                          "None" = "none",
                          "Rows only" = "rows",
                          "Columns only" = "columns",
                          "Rows + columns" = "both"
                        ), selected = "none"

            ),
            plotlyOutput("pairwise_peptide_quant_correlation")
            ),
            "pairwise_quant_correlation"
          ),
          ### ---- PCA plot ----
          with_export(
            bs4Card(title = "PCA plot", width = 12, maximizable = TRUE,
                    radioButtons("pca_level", "PCA level",
                                 c("Sample-level" = "sample", "Measurement-level" = "measurement"),
                                 selected = "sample", inline = TRUE),
                    radioButtons("pca_dim", "View", c("2D" = "2d", "3D" = "3d"),
                                 selected = "2d", inline = TRUE),
                    uiOutput("pca_group_ui"),
                    plotlyOutput("pca", height = "600px"),
                    tags$hr(),
                    h5("Explained variance"),
                    h6("This plot shows how many principle components (PC) are required to explain how much of the variance
                     and how much each PC controbutes. The more PC are required to explain a given variance, the more complex the data is to group."),
                    plotlyOutput("pca_variance", height = "300px")
            ),
            "pca_scatter"
          )
        )
      ),

      ## ---- Group Comparison ----
      bs4TabItem(
        tabName = "group_comp",  # must match menuItem
        fluidRow(
          h6("Note: No normalization is being done. Max Quantity values are used as is to calculate the group mean which is then compared."),
          ### ---- Create Groups ----
          bs4Card(title = tagList("Create Groups", bs4Dash::tooltip(icon("info-circle"),"Here we can create groups out of one or more samples to do group based comparison.
                                                       The groups are only valid for this page. Pairwise comparisons will be done for all possible group pairings.
                                                       This can be time consuming, as such the maximum number of groups is limited to 5.
                                                       Group data are combined via unionized, instead of intersected.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  sliderInput(
                    "min_presence_fraction",
                    tagList("Minimum fraction the peptide has to be present per group", bs4Dash::tooltip(icon("info-circle"),"0 = union, 1 = intersect", placement = "right")
                    ),
                    min = 0,
                    max = 1,
                    value = 0.7,
                    step = 0.05
                  ),
                  numericInput("n_groups", "Number of groups (2-5):", 2, min = 2, max = 5),
                  uiOutput("group_assign_ui"),
                  actionButton("update_group_comp", "Update Groups"),
                  radioButtons("go_background", "Enrichment background:",
                               c("Whole genome" = "genome", "Detected proteins" = "detected", "Custom list" = "custom"),
                               selected = "genome", inline = TRUE),
                  conditionalPanel("input.go_background == 'custom'",
                                   textAreaInput("go_custom_ids", "Custom background — UniProt IDs (space/comma/newline separated):",
                                                 rows = 4, placeholder = "P04439\nP01889\n...")),

                  # --- warning box: only visible for detected/custom ---
                  conditionalPanel(
                    condition = "input.go_background == 'detected' || input.go_background == 'custom'",
                    div(class = "alert alert-warning", style = "margin-top:8px; padding:8px 12px;",
                        icon("exclamation-triangle"),
                        tags$b(" Online lookup required. "),
                        "Protein symbols through UniProt web requests. The first run on a dataset can be slow \u2014. Results are cached afterward, so later ",
                        "changes are fast.")
                  ),

                  selectInput("go_ont", "GO ontology:",
                              c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"),
                              selected = "BP")
          ),

          ### ---- Group Stuff ----
          with_export(
            bs4Card(title = tagList("Group based Unique Peptides", bs4Dash::tooltip(icon("info-circle"),"Violin-plot showing distribution of peptides quantities within groups. All measurements from all samples are used.", placement = "right")
            ), width = 6, maximizable = TRUE,
            plotOutput("group_unique_bar")
            ),
            "group_unique_bar"
          ),
          with_export(
            bs4Card(title = tagList("Group based Euler Diagram", bs4Dash::tooltip(icon("info-circle"),"Venn Diagram with circle size representing group size", placement = "right")
            ), width = 6, maximizable = TRUE,
            plotOutput("group_euler")
            ),
            "group_euler"
          ),
          with_export(
            bs4Card(title = tagList("Peptide based heatmap", bs4Dash::tooltip(icon("info-circle"),"Cluster peptides based on differential pattern across groups. Clustering distance is 'euclidean' and method is 'complete'.
                                                                            The MAX quantity for the peptides is used. Peptides are binned for over 1000+ peptides. Peptide Llabels are remoed when 250+ peptides are visible. Zoom to adjust. NA is coloured Gray.", placement = "right")
            ), width = 12, maximizable = TRUE,
            plotlyOutput("group_peptide_heatmap_interactive", height = "600px"),
            fluidRow(
              column(4,tagList(actionButton("expand_heatmap_region", "Expand Visible Bins",icon = icon("search-plus"),class = "btn-sm btn-outline-info mt-2"),
                               bs4Dash::tooltip(
                                 icon("info-circle"),
                                 title = "Replaces the heatmap with individual peptides from the currently zoomed region. Peptides outside this regions are removed and needs to be reset via Full View button.",
                                 placement = "top")
                               )
                     ),
              column(4,tagList(actionButton("reset_heatmap_view", "Full View",icon = icon("compress"),class = "btn-sm btn-outline-secondary mt-2"),
                               bs4Dash::tooltip(
                                 icon("info-circle"),
                                 title = "Returns to the full binned heatmap with all peptides. Use this to restore the full data that was cut by the Expand Visible Bins.",
                                 placement = "top")
                               )
                     ),
              column(4,downloadButton("download_visible_peptides", "Download Visible",class = "btn-sm btn-outline-secondary mt-2")
                     )
              )
            ),
            "group_heatmap"
          ),
          ### ---- Statistical Plots ----
          with_export(
          bs4Card(title = tagList("Statistical Plots", bs4Dash::tooltip(icon("info-circle"),"Only significant datapoints are interactable (to improve speed).
                                                       Both Bonferroni and Benjamin-hochberg adjusted p-value are calculated. Only Benjamin-hochberg used for now.
                                                      p-value can only be calculated when more than 2 datapoints/measurements per peptide in each group exist. Otherwise, no p-value is calculated, which means no visualization.
                                                      Significant here means a corrected p-value of less than 0.05 and a log2FC larger than 1.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  uiOutput("group_stats_tabs")
                  ),
          "group_stats"
          )
        )
      ),

      ## ---- PTMs ----
      bs4TabItem(
        tabName = "ptm",  # must match menuItem
        fluidRow(
          ### ---- PTM distribution ----
          with_export(
            bs4Card(title = "PTM distribution", width = 12, maximizable = TRUE,
                    plotOutput("PTM_plot")
            ),
            "ptm_distribution"
          ),

          ### ---- PTM sequence motif ----
          bs4Card(title = tagList("PTM sequence motif", bs4Dash::tooltip(icon("info-circle"),"The below table and motif plots are experimental analysis to identify the significance of PTMs on the motif.
                                                       To do this, modified amino acids are treated as their own unique amino acid with a given symbol as mapped on the table.
                                                       Note: Since binding predictions do not account for PTMs, it is not possible to determine whether these peptidoforms are predicted binders or only their native form.
                                                       Only modified peptides are used to plot the motifs. Otherwise signals will be buried under the sheer quantity of non-modified peptides", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  DT::dataTableOutput("PTM_mod_map_table"),
                  uiOutput("PTM_motif_tabs")
          )
        )
      ),

      ## ---- Binding Predictions ----
      bs4TabItem(
        tabName = "binding_pred",  # must match menuItem
        fluidRow(
          ### ---- Call netMHCpan from windows subsystem for linux ----
          bs4Card(title = "run netMHCpan", width = 12, maximizable = TRUE,
                  selectizeInput(
                    "HLA_alleles",
                    "Select HLA allele(s):",
                    choices  = NULL,
                    selected = NULL,
                    multiple = TRUE,
                    options  = list(placeholder = "Type to search alleles...")
                  ),
                  actionButton("run_netmhc", "Run netMHCpan"),
                  selectizeInput(
                    "HLA_alleles_II",
                    "Select MHC-II allele(s):",
                    choices = NULL,
                    multiple = TRUE,
                    options = list(placeholder = "Type to search DR / DQ / DP alleles...")
                    ),
                  actionButton("run_netmhcII", "Run netMHCIIpan"),
                  textOutput("netmhc_status"),
                  radioButtons("mhc_class", "Default thresholds:", c("MHC-I" = "I", "MHC-II" = "II"), selected = "I", inline = TRUE),
                  h6("For custom cutoffs, go to Console \u2192 Settings.")
          ),

          ### ---- Allele selection ----
          uiOutput("allele_viz_selector_ui"),

          ### ---- Summary Table ----
          bs4Card(title = tagList("Summary Table", bs4Dash::tooltip(icon("info-circle"),"On default the threshold for weak binder is 2.0 and for strong binder is 0.5 Rank_EL. (Hard coded).", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  DT::DTOutput("binding_summary")
          ),

          ### ---- Summary barchart ----
          with_export(
            bs4Card(title = tagList("Binding summary", bs4Dash::tooltip(icon("info-circle"),
                                                                        "Best = best class across selected alleles (min rank). Per allele = each allele's own count breakdown, one panel per sample.",
                                                                        placement = "right")),
                    width = 12, maximizable = TRUE,
                    radioButtons("binding_view", "View:",
                                 c("Best (any allele)" = "best", "Per allele" = "per_allele"),
                                 selected = "best", inline = TRUE),
                    radioButtons("binding_scale", "Scale:",
                                 c("Absolute" = "absolute", "Percentage" = "percent"),
                                 selected = "absolute", inline = TRUE),
                    conditionalPanel("input.binding_view == 'per_allele'",
                                     radioButtons("binding_group", "Group by:",
                                                  c("Sample" = "sample", "Allele" = "allele"),
                                                  selected = "sample", inline = TRUE)),
                    uiOutput("binding_plot_ui")
            ),
            "predicted_binders"
          ),
          ### ---- Peptide Table ----
          bs4Card(title = "Peptide Table", width = 12, maximizable = TRUE,
                  DT::DTOutput("binding_table")
          ),
          fluidRow(
            column(3, downloadButton("dl_all", "Download All Peptides")),
            column(3, downloadButton("dl_binders", "Download Binders")),
            column(3, downloadButton("dl_nonbinders", "Download Non-Binders")),
            column(3, downloadButton("dl_missing", "Download Missing"))
          )
        )
      ),

      ## ---- Peptide Lookup ----
      bs4TabItem(
        tabName = "peptide_lookup",  # must match menuItem
        fluidRow(
          ### ---- Peptide Lookup ----
          h6("If you want to find peptides only found in one sample, set n_datasets = 1"),
          bs4Card(title = tagList("Peptide Lookup", bs4Dash::tooltip(icon("info-circle"),"The stripped peptides rather than peptidoforms are used as identifier.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  DT::DTOutput("peptide_table")
                  )
          )
        ),
      ## ---- SQL ----
      bs4TabItem(
        tabName = "SQL_console",
        bs4Card(title = "Save to database", width = 12,
                h6("The normalized raw data is saved, without any filtering."),
                actionButton("save_to_db", "Save results to database", icon = icon("database")),
                tags$hr(),
                tags$style(HTML("
                      .handsontable thead th .colHeader {
                        white-space: normal;
                        word-break: break-word;
                        line-height: 1.15;
                        display: inline-block;
                        width: 100%;
                      }
                      .handsontable thead th {
                        white-space: normal;
                        vertical-align: middle;
                        height: auto;
                      }
                    ")),
                h6("Saved runs"),
                DT::DTOutput("saved_runs_table"),
                actionButton("load_from_db", "Load selected analysis", icon = icon("upload"))
        )
      ),

      bs4TabItem(tabName = "SQL_query",
                 bs4Card(title = "SQL help — click to expand", width = 12, collapsible = TRUE, collapsed = TRUE,
                         tags$p("A query reads: ", tags$code("SELECT columns FROM table WHERE conditions"), "."),
                         tags$b("Tables:"),
                         tags$ul(
                           tags$li(tags$code("analyses"), " — one row per saved run."),
                           tags$li(tags$code("sample_metadata"), " — per-measurement info, long form (field_name / field_value)."),
                           tags$li(tags$code("peptides"), " — the data: STRIPPED, PROTEIN, quantity, SCORE, RT, MZ, LENGTH…")),
                         tags$b("Try these:"),
                         tags$pre("SELECT * FROM analyses LIMIT 20;"),
                         tags$pre("SELECT DISTINCT field_name FROM sample_metadata;")
                 ),
                 bs4Card(title = "Query", width = 12,
                         textAreaInput("sql_query", NULL, rows = 4, width = "100%",
                                       value = "SELECT * FROM analyses ORDER BY timestamp DESC LIMIT 50;"),
                         actionButton("sql_run", "Run", icon = icon("play"), class = "btn-primary"),
                         tags$hr(),
                         DT::DTOutput("sql_result")
                 )
      ),
      bs4TabItem(tabName = "SQL_browser",
                 bs4Card(title = "Browse tables", width = 12,
                         helpText("Pick a table to view its rows — filter and search the columns directly. ",
                                  "For questions that combine tables, use ", tags$b("query SQL"), "."),
                         selectInput("browse_table", "Table:", choices = NULL),
                         DT::DTOutput("browse_result")
                 )
      ),

      ## ---- Dev Console ----
      bs4TabItem(
        tabName = "dev_console",
        fluidRow(
          bs4Card(title = "R Console", width = 12,
                  textAreaInput("console_input", NULL,
                                value    = "",
                                rows     = 6,
                                width    = "100%",
                                placeholder = "Type R code here..."),
                  actionButton("console_run", "Run", class = "btn-primary"),
                  actionButton("console_clear", "Clear", class = "btn-secondary ml-2"),
                  tags$hr(),
                  verbatimTextOutput("console_output")
          )
        )
      ),

      ## ---- Advanced settings ----
      bs4TabItem(
        tabName = "advanced_settings",
        checkboxInput("split_measurements",
                      "Treat each measurement as a separate sample",
                      value = FALSE),
        helpText("For datasets where one file contains multiple conditions as measurement columns. ",
                 "Splits each measurement into its own sample. Original data is preserved for saving. ",
                 "Changing this resets sample selection and manual groups."),
        fluidRow(
          bs4Card(title = "Binding thresholds", width = 6,
                  numericInput("strong_cut", "Strong binder  %Rank \u2264", value = 0.5, min = 0, step = 0.1),
                  numericInput("weak_cut",   "Weak binder  %Rank \u2264",   value = 2.0, min = 0, step = 0.1)
          ),
          bs4Card(title = "netMHCpan raw results", width = 6,
                  p("Raw netMHCpan outputs are written to R's session temp directory:"),
                  verbatimTextOutput("netmhc_tempdir"),
                  actionButton("open_netmhc_folder", "Open folder", icon = icon("folder-open"), class = "btn-sm")
          )
        )
      )
      )
    )
  )}
