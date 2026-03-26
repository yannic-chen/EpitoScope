ui <- bs4DashPage(
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
    selectInput(
      inputId = "color_palette",
      label = tagList(icon("palette"), "Color palette"),
      choices = c("default", "viridis", "magma", "inferno", "plasma", "cividis", "mako", "rocket", "turbo"),
      selected = "default"
    ),
    # Your other UI elements here
    uiOutput("sample_selector"),
    uiOutput("length_slider_ui"),
    uiOutput("quantity_slider_ui"),
    uiOutput("score_slider_ui"),
    uiOutput("charge_slider_ui"),
    uiOutput("mass_slider_ui"),
    uiOutput("RT_slider_ui"),
    actionButton("generate_report", tagList(icon("file-arrow-down"),"Generate HTML Report"))
  ),
  
  ## ---- Header ----
  header = bs4DashNavbar(
    fixed = TRUE,
    title = dashboardBrand(
      title = "EpitoScope",
      color = "primary",
      #href = "https://adminlte.io/themes/v3", #host link
      #image = "https://adminlte.io/themes/v3/dist/img/AdminLTELogo.png"#add image
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
        text = "ExtraMenu",
        dropdownHeader("Dropdown header"),
        navbarTab(tabName = "Tab3", text = "Tab 3"),
        dropdownDivider(),
        navbarTab(
          text = "Sub menu",
          dropdownHeader("Another header"),
          navbarTab(tabName = "Tab4", text = "Tab 4"),
          dropdownHeader("Yet another header"),
          navbarTab(tabName = "Tab5", text = "Tab 5"),
          navbarTab(
            text = "Sub sub menu",
            navbarTab(tabName = "Tab6", text = "Tab 6"),
            navbarTab(tabName = "Tab7", text = "Tab 7")
            )
          )
        )
      )
  ),
  
  ## ---- Control Bar ----
  controlbar = bs4DashControlbar( #this is just an extra sidebar on the right
    skinSelector(), pinned = FALSE, value = ""
    ),
  
  ## ---- Main ----
  body = bs4DashBody(
    shinyjs::useShinyjs(), #needed to make button grey out
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
                  p(DT::DTOutput("summary_table"))
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
                             h6("Proteins", bs4Dash::tooltip(icon("info-circle"),"This depends heavily on the style of the header and how the software identifies it.", placement = "right")),
                             plotOutput("summary_proteins_plot")
                             )
                    )
                  ),
          
          ### ---- Unique Entries ----
          bs4Card(title = "Other Numeric Columns", width = 12, maximizable = TRUE, 
                  tabsetPanel(
                    tabPanel("Charge",
                             plotOutput("charge_plot")                               
                             ),
                    tabPanel("Mass",
                             plotOutput("mass_plot")
                             ),
                    tabPanel("m/z",
                             plotOutput("mz_plot")
                             ),
                    tabPanel("RT",
                             uiOutput("RT_plot")
                             ),
                    tabPanel("Mass Error",
                             bs4Dash::tooltip(icon("info-circle"),"Either ppm (PEAKS) or delta Mass (Fragpipe). Ignore the x-axis label.", placement = "right"),
                             plotOutput("ppm_plot")
                             ),
                    tabPanel("score",
                             plotOutput("score_violin")
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
          bs4Card(title = "Peptide Length Distribution", width = 12, maximizable = TRUE, 
                  plotlyOutput("length_plot")
                  ),
          
          ### ---- Length Range Percentage ----
          bs4Card(title = "Length Range Percentage", width = 12, maximizable = TRUE,
                  plotOutput("length_range_percentage")
                  ),
          
          ### ---- Motif Plot ---
          bs4Card(title = tagList("Motif Plot", bs4Dash::tooltip(icon("info-circle"),"Minimum of 5 sequences are required for Motif generation.", placement = "right")
                                  ), width = 12, maximizable = TRUE, 
                  uiOutput("motif_tabs")
                  ),
          
          ### ---- Dynamic Range plot ----
          bs4Card(title = tagList("Dynamic Rang", bs4Dash::tooltip(icon("info-circle"),"Minimum of 10 sequences are required for Motif generation.", placement = "right")
                                  ), width = 12, maximizable = TRUE, 
                  tabsetPanel(
                    tabPanel("Individual", h6("HLA[A-C] are highlighted"),uiOutput("dynrange_individual_ui")),
                    tabPanel("Combined", plotOutput("dynrange_combined"))
                    )
                  ),
          
          ### ---- 1/k0 vs mz ----
          bs4Card(title = "k0 vs mz", width = 12, maximizable = TRUE,
                  uiOutput("scatterplots_ui")
                  ),
          
          ### ---- measurement specific heatmap ----
          bs4Card(title = "Measurement Specific Heatmap", width = 12, maximizable = TRUE,
                  selectInput("cluster_mode_ea", "Clustering:",
                              choices = c(
                                "Sample" = "sample",
                                "Rows only" = "rows",
                                "Columns only" = "columns",
                                "Rows + columns" = "both",
                                "Sample + rows" = "mix"
                              ), selected = "sample"
                              
                  ),
                  plotOutput("measurement_heatmap", width= "auto", height = "auto")
                  ),
          
          ### ---- Composition profiling ----
          bs4Card(title = tagList("Composition profiling", bs4Dash::tooltip(icon("info-circle"),"Ideally I want to integrate the whole C.profiler from Vacic et al. 2007, but that is only written in python.", placement = "right")
                                  ), width = 12, maximizable = TRUE, solidHeader = TRUE, status = "warning", collapsed = TRUE,
                  h6("WIP: If we want to compute enrichment, then we need to do it compared to a background. Either the reference proteome or the sum of peptides of all samples can be used. For reference proteome, probably just hardcode the info.", style = "color: red;"),
                  h6("WIP: Colour or number in the heatmap should represent difference to background. Is it possible to get numbers on the bargraph?", style = "color: red;"),
                  plotOutput("aa_heatmap", width = "1200", height = "auto")
                  )
          )
        ),
      
      ## ---- Results ----
      bs4TabItem(
        tabName = "results",  # must match menuItem
        fluidRow(
          
          ### ---- Data Completeness ----
          bs4Card(title = "Data Completeness", width = 12, maximizable = TRUE,
                  h6("WARNING: For PEAKS 12 Studio, the column X.Spec has been replaced with X.Feature. X.Feature returns 0 even when a peptide has been identified but could not be quantified. X.Spec on the other hand only returns 0 if it is not identified at all."),
                  tabsetPanel(
                    # --- absolute ---
                    tabPanel(
                      "Absolute",
                      plotOutput("completeness_plot")
                    ),
                    
                    # --- percentage ---
                    tabPanel(
                      "Percentage",
                      plotOutput("completeness_plot2")
                    )
                  )
          ),
          
          ### ---- Upset Plot of Peptides ----
          bs4Card(title = tagList("Upset Plot of Peptides", bs4Dash::tooltip(icon("info-circle"),"Currently intersection min_size is 2% of combined number of unique peptides.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  plotOutput("upset_plot")
                  ),
          
          ### ---- Pairwise shared peptide matrix ----
          bs4Card(title = tagList("Pairwise shared peptide matrix", bs4Dash::tooltip(icon("info-circle"),"Pairwise comparison of number of shared peptides. For Percent visualization union is used (i.e. jaccard style).", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  selectInput("shared_mode", "Visualization:", choices = c("Count" = "count", "Percent" = "percent")),
                  plotOutput("Pairwise_shared_peptide_matrix")
                  ),
          
          ### ---- Pairwise shared peptide quantity comparison matrix ----
          bs4Card(title = tagList("Pairwise shared peptide quantity comparison matrix", bs4Dash::tooltip(icon("info-circle"),"Pairwise comparison of peptide max quantity of shared peptides, using pearson correlation. Clustering distance is 'euclidean' and method is 'complete'.", placement = "right")
          ), width = 12, maximizable = TRUE,
          selectInput("cluster_mode", "Clustering:",
                      choices = c(
                        "None" = "none",
                        "Rows only" = "rows",
                        "Columns only" = "columns",
                        "Rows + columns" = "both"
                      ), selected = "none"
                      
          ),
          plotOutput("pairwise_peptide_quant_correlation")
          ),
          
          ### ---- PCA plot ----
          bs4Card(title = "PCA plot", width = 12, maximizable = TRUE,
                  plotOutput("pca")
          ),
          
          ### ---- Unique Entries ----
          bs4Card(title = "Unique entries", width = 12, maximizable = TRUE, 
                  tabsetPanel(
                    tabPanel("Peptides",
                             h6("Unique Peptides (no PTMs)"),
                             plotOutput("summary_peptides_plot2")
                    ),
                    tabPanel("Peptidoforms", 
                             h6("Peptidoforms (including PTMs)"),
                             plotOutput("summary_peptidoforms_plot2")
                    ),
                    tabPanel("Proteins",
                             h6("Proteins", bs4Dash::tooltip(icon("info-circle"),"This depends heavily on the style of the header and how the software identifies it.", placement = "right")),
                             plotOutput("summary_proteins_plot2")
                    )
                  )
          ),
          
          ### ---- Unique Entries ----
          bs4Card(title = "Other Numeric Columns", width = 12, maximizable = TRUE, 
                  tabsetPanel(
                    tabPanel("Charge",
                             plotOutput("charge_plot2")                               
                    ),
                    tabPanel("Mass",
                             plotOutput("mass_plot2")
                    ),
                    tabPanel("m/z",
                             plotOutput("mz_plot2")
                    ),
                    tabPanel("RT",
                             uiOutput("RT_plot2")
                    ),
                    tabPanel("Mass Error",
                             bs4Dash::tooltip(icon("info-circle"),"Either ppm (PEAKS) or delta Mass (Fragpipe). Ignore the x-axis label.", placement = "right"),
                             plotOutput("ppm_plot2")
                    ),
                    tabPanel("score",
                             plotOutput("score_violin2")
                    )
                  )
          )
        )
      ),
      
      ## ---- Group Comparison ----
      bs4TabItem(
        tabName = "group_comp",  # must match menuItem
        fluidRow(
          
          ### ---- Create Groups ----
          bs4Card(title = tagList("Create Groups", bs4Dash::tooltip(icon("info-circle"),"Here we can create groups out of one or more samples to do group based comparison. 
                                                       The groups are only valid for this page. Pairwise comparisons will be done for all possible group pairings.
                                                       This can be time consuming, as such the maximum number of groups is limited to 5.
                                                       Group data are combined via unionized, instead of intersected.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  #numericInput(
                  #  inputId = "max_missing_per_group",
                  #  label   = "Max missing samples per peptide (per group)",
                  #  value   = 0,
                  #  min     = 0,
                  #  step    = 1
                  #),
                  sliderInput(
                    "min_presence_fraction",
                    tagList("Minimum fraction of samples the peptide has to be present per group", bs4Dash::tooltip(icon("info-circle"),"0 = union, 1 = intersect", placement = "right")
                    ),
                    min = 0,
                    max = 1,
                    value = 0.7,
                    step = 0.05
                  ),
                  numericInput("n_groups", "Number of groups (2-5):", 2, min = 2, max = 5),
                  uiOutput("group_assign_ui"),
                  actionButton("update_group_comp", "Update Groups"),
                  h6("WIP: currently this checkbox doesnt do anything. Need to solve how to combine info same peptide with different peptidoforms"),
                  checkboxInput(
                    inputId = "use_peptidoforms",
                    label   = tagList(
                      "Compare at peptidoform level ",
                      tags$small(
                        style = "color: grey; font-weight: normal;",
                        "(includes PTMs)"
                      )
                    ),
                    value = FALSE  # default: compare stripped sequences
                  )
          ),
          
          ### ---- Group Stuff ----
          bs4Card(title = tagList("Quantity variance", bs4Dash::tooltip(icon("info-circle"),"Violin-plot showing distribution of peptides quantities within groups. All measurements from all samples are used.", placement = "right")
                                  ), width = 6, maximizable = TRUE
          ),
          bs4Card(title = tagList("Group based Venn Diagram", bs4Dash::tooltip(icon("info-circle"),"This is nicer to identify biological differences given multiple biological samples.", placement = "right")
                                  ), width = 6, maximizable = TRUE,
                  plotOutput("group_venn_plot")
          ),
          bs4Card(title = tagList("Peptide based heatmap", bs4Dash::tooltip(icon("info-circle"),"Cluster peptides based on differential pattern across groups. The MAX quantity for the peptides is used.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  plotOutput("group_peptide_heatmap")
          ),
          
          ### ---- Statistical Plots ----
          bs4Card(title = tagList("Statistical Plots", bs4Dash::tooltip(icon("info-circle"),"Only significant datapoints are interactable (to improve speed). 
                                                       Both Bonferroni and Benjamin-hochberg adjusted p-value are calculated. Only Benjamin-hochberg used for now. Could include a switch.
                                                      Values for p-value and logFC calculations come from all columns designated as QUANTITY information from the 'Column Map' info.
                                                      MA plot is helpful to determine whether the interesting hits are low or high expressed peptides 
                                                      (while volcano plot gives you highly differential expressed peptides that are also significant).
                                                      Significant here means a corrected p-value of less than 0.05 and a log2FC larger than 1.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  uiOutput("group_stats_tabs")
          ),
          
          ### ---- GO-terms ----
          bs4Card(title = "GO-term and STRING-DB analysis. Proof of work", width = 6, maximizable = TRUE,
                  h6("WIP: This is wrong, since we would need to translate peptide difference to gene/protein level difference. And accessions need to be correct. 
                     Currently we only pick the significant peptides and only get their first accession. These accessions are then used for GO-term. Also we do not differentiate between negative and positive change.
                     We use bitr() to convert from uniprot to entrezID, so if the protein is not uniprot ID, we will likely get an error.
                     we also only keep the first mapping per UniProt.", style = "color: red;")
          )
        )
      ),
      
      ## ---- PTMs ----
      bs4TabItem(
        tabName = "ptm",  # must match menuItem
        fluidRow(
          ### ---- PTM distribution ----
          bs4Card(title = "PTM distribution", width = 12, maximizable = TRUE,
                  plotOutput("PTM_plot")
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
          h6("Ensure correct format of the netMHCpan data. THis is now important with the implementation of netMHCpan call."),
          h6("For the sake of filesize, Rank_EL has been limited to 2 decimal places and number have a ceiling of 9.99."),
          h6("When 9.99 -> '', one allele adds around 48MB data (12MB per length)"),
          
          ### ---- Call netMHCpan from windows subsystem for linux ----
          bs4Card(title = "run netMHCpan", width = 12, maximizable = TRUE,
                  selectInput(
                    "HLA_alleles",
                    "Select HLA allele(s):",
                    choices = c(
                      "HLA-A02:01",
                      "HLA-A01:01",
                      "HLA-A03:01",
                      "HLA-A24:02",
                      "HLA-B07:02",
                      "HLA-B08:01",
                      "HLA-B15:01",
                      "HLA-C07:01",
                      "HLA-C07:02"
                    ),
                    selected = "HLA-A02:01",
                    multiple = TRUE
                  ),
                  actionButton("run_netmhc", "Run netMHCpan"),
                  textOutput("netmhc_status")
          ),
          ### ---- Summary Table ----
          bs4Card(title = tagList("Summary Table", bs4Dash::tooltip(icon("info-circle"),"On default the threshold for weak binder is 2.0 and for strong binder is 0.5 Rank_EL. (Hard coded).", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  DT::DTOutput("binding_summary")
          ),
          
          ### ---- Summary barchart ----
          bs4Card(title = tagList("Summary barchart", bs4Dash::tooltip(icon("info-circle"),"The binders is the best category for all alleles.", placement = "right")
                                  ), width = 12, maximizable = TRUE,
                  tabsetPanel(
                    # --- absolute ---
                    tabPanel(
                      "Absolute",
                      plotOutput("binding_plot_absolute")
                    ),
                    
                    # --- percentage ---
                    tabPanel(
                      "Percentage",
                      plotOutput("binding_plot_percent")
                    )
                  )
          ),
          
          ### ---- Peptide Table ----
          bs4Card(title = "Peptide Table", width = 12, maximizable = TRUE,
                  DT::DTOutput("binding_table")
          ),
          h6("WIP: Maybe remove everything (sample and EL_rank) and only keep the peptide sequence? Makes copy-pasting easier.", style = "color: red;"),
          h6("WIP: With this we would be able to unique by peptide to remove duplicates. Since comma separated list (Set) destroys the csv style file format.", style = "color: red;"),
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
        )
      )
    )
  )
