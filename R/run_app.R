#' Launch the EpitoScope Shiny application
#'
#' @param input_variable Optional pre-loaded annotation table. Default `NULL`.
#' @param generate_pseudo_sequence Logical, passed to the server. Default `FALSE`.
#' @param custom_schema,custom_signature Optional named lists to extend/replace the schema.
#' @param replace_schema Logical; replace vs. extend the built-in schema.
#' @param db_path SQLite results DB path. Defaults under [tools::R_user_dir()].
#' @param ... Passed to [shiny::shinyApp()].
#' @return A Shiny app object
#' @export
run_app <- function(input_variable = NULL,
                    generate_pseudo_sequence = FALSE,
                    custom_schema = NULL,
                    custom_signature = NULL,
                    replace_schema = FALSE,
                    db_path = NULL,
                    ...) {
  options(shiny.maxRequestSize = 5 * 1024^3, width = 10000)
  ComplexHeatmap::ht_opt(message = FALSE)

  if (is.null(db_path)) {
    dir <- tools::R_user_dir("EpitoScope", which = "data")
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    db_path <- file.path(dir, "epitoscope_results.sqlite")
  }
  options(epitoscope.db = db_path)
  seed_db_once()

  shiny::shinyApp(
    ui = app_ui(),
    server = function(input, output, session) {
      app_server(input, output, session,
                 input_variable           = input_variable,
                 generate_pseudo_sequence = generate_pseudo_sequence,
                 custom_schema            = custom_schema,
                 custom_signature         = custom_signature,
                 replace_schema           = replace_schema)
    },
    ...
  )
}
