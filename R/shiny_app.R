#' Launch TuNeS Interactive Polygon Selector
#'
#' Opens a Shiny app for interactively drawing polygons on spatial data
#'
#' @return Launches Shiny application
#' @export
#' @examples
#' \dontrun{
#' launch_tunes()
#' }
launch_tunes <- function() {
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("TuNeS: Tumor Nest Selector"),
    
    shiny::fluidRow(
      shiny::column(12, 
                    shiny::textInput("seurat_name", 
                                     "Seurat object name:", 
                                     value = "seurat_obj"),
                    shiny::textInput("fov_name", 
                                     "FOV name:", 
                                     value = "fov"),
                    shiny::actionButton("load_data", 
                                        "Load Data", 
                                        style = "margin-bottom: 20px;")
      )
    ),
    
    plotly::plotlyOutput("plot", height = "600px"),
    
    shiny::fluidRow(
      shiny::column(6, 
                    shiny::textInput("polygon_name", 
                                     "Polygon object name:", 
                                     value = "drawn_polygons")),
      shiny::column(6, 
                    shiny::actionButton("extract", 
                                        "Extract Polygons", 
                                        style = "margin-top: 25px;"))
    ),
    
    shiny::verbatimTextOutput("shapes")
  )
  
  server <- function(input, output, session) {
    
    rv <- shiny::reactiveValues(
      shapes = NULL,
      coords = NULL,
      data_loaded = FALSE
    )
    
    shiny::observeEvent(input$load_data, {
      tryCatch({
        seurat_obj <- get(input$seurat_name, envir = .GlobalEnv)
        coords <- Seurat::GetTissueCoordinates(seurat_obj, image = input$fov_name)
        coords$is_cancer <- seurat_obj$is_cancer
        
        rv$coords <- coords
        rv$data_loaded <- TRUE
        
        output$shapes <- shiny::renderPrint({
          cat(paste0("✓ Loaded data from '", input$seurat_name, 
                     "' with FOV '", input$fov_name, "'\n"))
          cat(paste0("  ", nrow(coords), " cells loaded\n"))
        })
        
      }, error = function(e) {
        output$shapes <- shiny::renderPrint({
          cat("ERROR loading data:\n")
          cat(paste0("  ", e$message, "\n"))
          cat("\nMake sure the Seurat object name and FOV name are correct.\n")
        })
        rv$data_loaded <- FALSE
      })
    })
    
    output$plot <- plotly::renderPlotly({
      shiny::req(rv$data_loaded)
      
      coords <- rv$coords
      
      plotly::plot_ly(data = coords, 
                      x = ~x, 
                      y = ~y,
                      type = 'scatter', 
                      mode = 'markers',
                      color = ~is_cancer,
                      colors = c("TRUE" = "red", "FALSE" = "lightgray"),
                      marker = list(size = 3, opacity = 0.5),
                      source = "plot") %>%
        plotly::layout(
          dragmode = 'drawclosedpath',
          newshape = list(
            line = list(color = 'blue', width = 2),
            fillcolor = 'rgba(0,0,255,0.2)'
          ),
          shapes = list(),
          title = "Draw tumor regions"
        ) %>%
        plotly::config(
          modeBarButtonsToAdd = c('drawclosedpath', 'drawrect', 'drawcircle', 'eraseshape'),
          displayModeBar = TRUE,
          editable = TRUE
        )
    })
    
    shiny::observe({
      event_data <- plotly::event_data("plotly_relayout", source = "plot")
      if (!is.null(event_data)) {
        rv$shapes <- event_data
      }
    })
    
    shiny::observeEvent(input$extract, {
      
      if (!rv$data_loaded) {
        output$shapes <- shiny::renderPrint({
          cat("ERROR: Please load data first by clicking 'Load Data'.\n")
        })
        return()
      }
      
      shapes_data <- rv$shapes
      obj_name <- input$polygon_name
      
      if (obj_name == "" || !grepl("^[a-zA-Z][a-zA-Z0-9._]*$", obj_name)) {
        output$shapes <- shiny::renderPrint({
          cat("ERROR: Invalid object name.\n")
          cat("Object name must start with a letter and contain only letters, numbers, dots, or underscores.\n")
        })
        return()
      }
      
      output$shapes <- shiny::renderPrint({
        cat("Attempting to extract shapes...\n\n")
        
        if (is.null(shapes_data)) {
          cat("No shapes data captured yet.\n")
          return()
        }
        
        if ("shapes" %in% names(shapes_data) && is.data.frame(shapes_data$shapes)) {
          shapes_df <- shapes_data$shapes
          
          cat(paste("Found", nrow(shapes_df), "shape(s)\n\n"))
          
          polygons_list <- list()
          
          for (i in 1:nrow(shapes_df)) {
            path <- shapes_df$path[i]
            
            cat(paste("Processing shape", i, "\n"))
            
            coords_parsed <- parse_svg_path(path)
            
            if (!is.null(coords_parsed)) {
              polygon_sf <- sf::st_polygon(list(as.matrix(coords_parsed))) %>%
                sf::st_sfc() %>%
                sf::st_sf(region_id = i)
              
              polygons_list[[i]] <- polygon_sf
              cat(paste("  ✓ Shape", i, "parsed successfully\n"))
            } else {
              cat(paste("  ✗ Failed to parse shape", i, "\n"))
            }
          }
          
          if (length(polygons_list) > 0) {
            all_polygons <- do.call(rbind, polygons_list)
            
            assign(obj_name, all_polygons, envir = .GlobalEnv)
            
            cat(paste0("\n✓ SUCCESS! Polygons created as '", obj_name, 
                       "' in your environment\n\n"))
            print(all_polygons)
          } else {
            cat("\nFailed to parse any polygons from the shape data.\n")
          }
        } else {
          cat("No shapes found in the expected format.\n")
        }
      })
    })
  }
  
  shiny::shinyApp(ui, server)
}
