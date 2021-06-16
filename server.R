server <- function(input, output) {
  #Initiate variables
  typefile<-reactiveValues(df=NULL)
  data<-reactiveValues(df=NULL)
  shinyDirChoose(input, 'dataset', filetypes=c("csv", "tsv"), roots=c(wd = normalizePath("~")))
  data_clean<-reactiveValues(df=NULL)
  abundance<-reactiveValues(df=NULL)
  printer<-reactiveValues(df=NULL)
  #Initiate placeholders
  output$imported_table <- DT::renderDataTable(
    {
      showData <- tibble(FileNames = vector("character"))
      isolate(showData)
    },
    rownames = FALSE,
    options = (list(scrollX = TRUE))
  )
  output$dir <- renderPrint({
    cat(c("No", "directory", "selected"), sep = " ")
  })
  
  #Initiate observeEvents
  observeEvent(input$dataset_import, { #button to upload data
    req(input$dataset)
    
    if (input$dataset_type=="vConTACT2 (csv)"){
      typefile$df<-"csv"
    } else { #DRAM-v filetype
      typefile$df<-"tsv"
      
    }
    x<<-input$dataset
    wd<<-file.path(normalizePath("~"), paste(unlist(x$path[-1]), collapse = .Platform$file.sep))
    output$dir<-renderText({wd})
    
    data$df<-MetaViral::import_files(path=wd, filetype = typefile$df)
    output$imported_table<-DT::renderDataTable({
      data$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    shinyjs::show("downloadRaw")
    shinyjs::show("downTypeRaw")
    printer$df<-data$df
    
  })#end observeEvente import
  
  observeEvent(input$dataset_clean, { #button to clean data
    if ( typefile$df == "csv"){
      output_from<-"vcontact2"
    }else{
      output_from<-"DRAMv"
    }
    data_clean$df<- MetaViral::cleaning(viral_data = data$df, output_from = output_from)
    output$imported_table<-DT::renderDataTable({
      data_clean$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    if (output_from == "vcontact2"){
      shinyjs::show("dataset_matrix")
      shinyjs::show("taxa_selection")
      shinyjs::show("dataset_abundance")
    }else{
      shinyjs::show("database")
      shinyjs::show("database_explore")
    }
    printer$df<-data_clean$df
  })#end observeEvent clean
  observeEvent(input$dataset_abundance, {
    showNotification("Wait patiently...")
    abundance$df<-MetaViral::abundance_vc(viral_vc = data_clean$df,
                                          taxa = input$taxa_selection,
                                          abuntype = input$dataset_matrix,
                                          output_type = "matrix")
    
    output$imported_table<-DT::renderDataTable({
      abundance$df
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
    
    if (input$dataset_matrix=="absolute"){
      shinyjs::show("bio_stats")
    }
    printer$df<- abundance$df
    
  })#end of observeEvent abundance
  observeEvent(input$bio_stats, {
    biostats<-MetaViral::bio_statistics(abundance$df)
    
    output$imported_table<-DT::renderDataTable({
      biostats
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
  })#end of observeEvent bio_stats
  
  observeEvent(input$database_explore, {
    data_base<-MetaViral::db_exploring(viral_annotations = data_clean$df, database = input$database )
    
    output$imported_table<-DT::renderDataTable({
      data_base
    },
    options=(list(pageLength=10,scrollX=TRUE))
    )
    showNotification("Done")
  })
  
  
  
  # Download Raw Data
  output$downloadRaw <- downloadHandler(
    filename = function() {
      today <- toupper(format(Sys.Date(), format = "%d%b%y"))
      paste("ViralData_", today, ".", input$downTypeRaw, sep = "")
    },
    content = function(file) {
      tmp <- printer$df
      
      if (input$downTypeRaw == "rds") {
        saveRDS(object = tmp, file = file)
      } else {
        
        write.csv(x = tmp, file = file, quote = F, row.names = F, col.names = T)
      }
    }
  )
  
}
shinyApp(ui = ui, server = server)
