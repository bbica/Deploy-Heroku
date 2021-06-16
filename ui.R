library(shiny)
library(DT)
library(shinyjs)
library(tibble)
library(dplyr)
library(KEGGREST)
library(readr)
library(purrr)
library(tidyr)
library(gtools)
library(shinyFiles)

list_csv <- function(flnm) {
  '%>%' <- tidyr::`%>%`
  return_namefile <-
    readr::read_csv(flnm) %>%
    dplyr::mutate(filename = flnm)
  return(return_namefile)
}

list_tsv <- function(flnm) {
  '%>%' <- tidyr::`%>%`
  return_namefile <-
    readr::read_tsv(flnm) %>%
    dplyr::mutate(filename = flnm)
  return(return_namefile)
}

import_files <- function(path, filetype="tsv"){
  '%>%' <- tidyr::`%>%`
  wd_dram<-path
  setwd(wd_dram)
  if (filetype=="tsv"){
    export_files <-
      list.files(pattern = paste0("*.", filetype),
                 full.names = T) %>%
      purrr::map_df(~list_tsv(.))
  } else if (filetype=="csv"){
    export_files <-
      list.files(pattern = paste0("*.", filetype),
                 full.names = T) %>%
      purrr::map_df(~list_csv(.))
  } else {
    warning("Filetype must be 'csv' or 'tsv'")
  }

  return(export_files)
}

cleaning<-function(viral_data, output_from="vcontact2", ...){
  '%>%' <- tidyr::`%>%`
  if (output_from=="vcontact2"){
    viral_genome <- tidyr::separate(viral_data, col=filename, into=c("Biome", "Category"), sep="-(?=[^-])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Biome, into=c("junk", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("junk", "Category"), sep="(?=[0-9])", remove = FALSE)
    viral_genome <- tidyr::separate(viral_genome, col=Category, into=c("Category", "junk"), sep="[-]", extra = "merge")
    viral_genome$junk <- NULL
    viral_genome<-viral_genome[,-1]

    viral_vc<-viral_genome %>%
      dplyr::group_by(`VC Subcluster`) %>%
      dplyr::arrange(Genome, Family, Order, Genus, Biome)

    viral_vc<-viral_vc[!is.na(viral_vc$VC), ]

    viral_vc <- tidyr::separate(viral_vc, col=VC, into=c("VC", "VC2"), sep="\\D", extra = "merge")

    viral_vc<-transform(viral_vc, VC = as.numeric(VC))
    viral_vc<-transform(viral_vc, VC2 = as.numeric(VC2))
    viral_vc <- viral_vc[order(viral_vc$VC, viral_vc$VC2, decreasing = FALSE), ]

    return(viral_vc)



  }else if (output_from=="DRAMv") {
    #extract only the name of the biome from the filename column (FAZER pipes!!!)
    DRAMv_with_sources <- tidyr::separate(viral_data, col=filename, into=c("Biome", "junk"), sep="_(?=[^_]+$)", extra = "merge")
    DRAMv_with_sources <- tidyr::separate(DRAMv_with_sources, col=Biome, into=c("junk2", "Biome"), sep="(?=[A-z0-9])", extra = "merge")
    DRAMv_with_sources$junk <- NULL
    DRAMv_with_sources$junk2 <- NULL

    DRAMv_with_sources2 <- DRAMv_with_sources%>%
      tidyr::separate(col=viral_hit, into=c("viral_function", "protein_and_vname"), sep="\\s", extra="merge")%>%
      tidyr::separate(col=protein_and_vname, into=c("viral_function", "viral_tax"), sep="(?=\\[)")%>% #separates the column by the fist  "["
      #kegg
      tidyr::separate(col=kegg_hit, into=c("kegg_hit", "kegg_hit2"), sep=";")%>%
      tidyr::separate(col=kegg_hit, into=c("kegg_hit", "junk"), sep="(?=\\[)")%>%
      #vogdb
      tidyr::separate(col=vogdb_description, into=c("vogdb_description", "junk"), sep=";")%>%
      tidyr::separate(col=vogdb_description, into=c("junk2", "vogdb_description"), sep="(?=\\s)", extra = "merge")%>%
      #pfam
      tidyr::separate(col=pfam_hits, into=c("pfam_hits", "pfam_id"), sep="(?=\\[)")%>%
      #peptidase
      tidyr::separate(col=peptidase_hit, into=c("junk", "peptidase_hit"), sep="-", extra = "merge")%>%
      tidyr::separate(col=peptidase_hit, into=c("peptidase_hit", "junk2"), sep="(?≥\\))", extra = "merge")%>%
      tidyr::separate(col=peptidase_hit, into=c("peptidase_hit", "peptidase_tax"), sep="(?=\\()", extra = "drop")%>%
      tidyr::separate(col=peptidase_tax, into=c("junk", "peptidase_tax"), sep="(?=\\w)", extra = "merge")

    DRAMv_with_sources2$junk <- NULL
    DRAMv_with_sources2$junk2 <- NULL

    #Extracts the taxonomic predition from the square brackets
    DRAMv_with_sources2$viral_tax = sub(".*\\[([^][]+)].*", "\\1", DRAMv_with_sources2$viral_tax)

    #removes rows without any predictions
    DRAMv_with_sources3 <- DRAMv_with_sources2 %>%
      dplyr::filter_at( dplyr::vars(kegg_hit,viral_function,pfam_hits,cazy_hits,vogdb_description,peptidase_hit),dplyr::any_vars(!is.na(.)))

    #removes rows with NA in the amg_flags
    DRAMv_with_sources4 <- DRAMv_with_sources3[!is.na(DRAMv_with_sources3$amg_flags),]
    #filter for only auxiliary scores inferior to 4
    DRAMv_with_sources4<-dplyr::filter(DRAMv_with_sources4, auxiliary_score<4)
    #filter out flags related with viral function
    unwanted_flags <- c("F", "T", "P", "A")

    for (i in 1:max(nchar(DRAMv_with_sources4$amg_flags))) {
      prm<-gtools::permutations(n=4,r=i,v=unwanted_flags, repeats.allowed = TRUE)
      prm<-apply(prm, 1, function(x) paste0(x, collapse=''))
      DRAMv_with_sources4<-DRAMv_with_sources4[ ! DRAMv_with_sources4$amg_flags %in% prm, ]
    }

    viral_annotations<-DRAMv_with_sources4
    return(viral_annotations)

  }else {
    warning("output_from must be vcontact2 or DRAMv")
  }
}#end of function

abundance_vc<-function(viral_vc, taxa="Family",output_type="matrix", abuntype="relative"){
  '%>%' <- tidyr::`%>%`
  #Majority rules
  viral_test1<-viral_vc %>%
    # add a column n with count by categories
    dplyr::add_count(VC.Subcluster, Genus, Biome, Category) %>%
    # select max or first occurrence
    dplyr::group_by(VC.Subcluster) %>%
    # keep only first TRUE
    dplyr::mutate(VC_Genus = Genus[n == max(n)][1]) %>%
    # do not keep temp var
    dplyr::select(-n)

  viral_test1<-viral_test1 %>%
    dplyr::add_count(VC.Subcluster, Family, Biome, Category) %>%
    dplyr::group_by(VC.Subcluster) %>%
    dplyr::mutate(VC_Family = Family[n == max(n)][1]) %>%
    dplyr::select(-n)

  viral_test1<-viral_test1 %>%
    dplyr::add_count(VC.Subcluster, Genome, Biome, Category) %>%
    dplyr::group_by(VC.Subcluster) %>%
    dplyr::mutate(VC_Genome = Genome[n == max(n)][1]) %>%
    dplyr::select(-n)

  viral_test1<-viral_test1 %>%
    dplyr::add_count(VC.Subcluster, Order, Biome, Category) %>%
    dplyr::group_by(VC.Subcluster) %>%
    dplyr::mutate(VC_Order = Order[n == max(n)][1]) %>%
    dplyr::select(-n)


  biomas<-unique(viral_vc$Biome)
  categories<-unique(viral_vc$Category)
  viralclusters<-unique(viral_vc$VC.Subcluster)
  vec<-NULL
  for (i in 1:length(biomas)){
    for (j in 1:length(categories)){
      for (k in 1:length(viralclusters)){
        rowselected<-which(viral_test1$Biome==biomas[i]&viral_test1$Category==categories[j]&viral_test1$VC.Subcluster==viralclusters[k])
        vector_id<-paste0(biomas[i], "-", categories[j], "-", viralclusters[k])
        occurencies<-length(rowselected)
        if (taxa=="Family"){
          Majority_taxa<-"VC_Family"
        }else if (taxa=="Order"){
          Majority_taxa<-"VC_Order"
        }else if (taxa=="Genus"){
          Majority_taxa<-"VC_Genus"
        }else if (taxa=="Genome"){
          Majority_taxa<-"VC_Genome"
        }else {
          Majority_taxa<-"VC_Family"
          print("Family selected as the taxa")
        }

        taxa_column<-which(colnames(viral_test1)==Majority_taxa)
        correspondent_taxa<-viral_test1[rowselected[1], taxa_column]
        vec<-rbind(vec, cbind(vector_id, correspondent_taxa, occurencies))
      }#loop k
    }#loop j
    print(biomas[i])
  }#loop i
  vec2<-vec[!is.na(vec[,2]), ]
  vec_biomas<-c()
  for (i in 1:nrow(vec2)){
    vec_biomas<-c(vec_biomas,strsplit(vec2[i,1], split = "-")[[1]][1])
  }
  vec_biomas_mx<-c()
  for (i in 1:length(biomas)){
    assign(paste0(biomas[i],"_data"),vec2[which(vec_biomas==biomas[i]),])
    vec_biomas_mx<-c(vec_biomas_mx, paste0(biomas[i],"_data"))
  }
  vec4<-list()
  for (i in 1:length(vec_biomas_mx)){
    bioma_mx<-get(vec_biomas_mx[i])
    unique_taxa<-unique(bioma_mx[,2])
    vec3<-NULL
    for (j in 1:length(unique_taxa)){
      if (abuntype=="absolute"){
        sum_abs<-sum(bioma_mx[which(bioma_mx[,2]==unique_taxa[j]), 3])
        vec3<-rbind(vec3,cbind(unique_taxa[j], sum_abs))

      }else{
        sum_abs<-sum(bioma_mx[which(bioma_mx[,2]==unique_taxa[j]), 3])
        sum_rel<-sum_abs/sum(bioma_mx[,3])
        vec3<-rbind(vec3,cbind(unique_taxa[j], sum_rel))
      }
      #dividor<-nrow(vec3)/length(vec_biomas_mx)
      #vec4<-vec3[1:dividor,]
    }
    vec4[[i]]<-vec3
  }
  names(vec4)<-biomas
  if (output_type=="list"){
    return(vec4)

  }else{
    vecall<-c()

    for (i in 1:length(vec4)){
      vectmp<-vec4[[i]][,1]
      vecall<-c(vecall,vectmp)
    }
    vecuniquetaxa<-unique(vecall)
    matrix_final<-matrix(data=NA, nrow=length(biomas), ncol=length(vecuniquetaxa))
    rownames(matrix_final)<-biomas
    colnames(matrix_final)<-vecuniquetaxa
    for (i in 1:length(biomas)){
      for (j in 1:length(vecuniquetaxa)){
        corresp<-which(vec4[[i]][,1]==vecuniquetaxa[j])
        insert_in_matrix<-c()
        if (rlang::is_empty(corresp)==TRUE){
          insert_in_matrix<-0
        }else{
          insert_in_matrix<-as.numeric(vec4[[i]][corresp,2])
        }
        matrix_final[i,j]<-insert_in_matrix
      }
    }
    return(matrix_final)
  }#end else
}#end of function

bio_statistics<-function(abundance_mx){
  Biomes<-rownames(abundance_mx)
  df_stats <- as.data.frame(matrix(0, ncol = 1, nrow = length(Biomes)))
  df_stats[,1]<-Biomes
  colnames(df_stats)<-"Biomes"

  df_stats$Abundance<-apply(abundance_mx,1,sum)#Total number of individuals(Abundance)
  df_stats$Richness<-apply(abundance_mx>0,1,sum)#Taxa Richness
  df_stats$Menhinick<-df_stats$Richness/sqrt(df_stats$Abundance)#Menhinick's index
  df_stats$Margalef<-(df_stats$Richness-1)/log(df_stats$Abundance)#Margalef's index
  df_stats$Shannon<-vegan::diversity(abundance_mx, index="shannon")#Shannon-Wiener Index
  df_stats$Simpson<-vegan::diversity(abundance_mx, index="simpson")#Simpson's Index

  df_stats$Pilou_evenness<-df_stats$Simpson/log(df_stats$Richness)#Pilou evenness
  df_stats$Simpson_evenness<-exp(df_stats$Simpson)/df_stats$Richness#Hill's evenness ratio for Simpson's Index
  df_stats$Shannon_evenness<-exp(df_stats$Shannon)/df_stats$Richness#Hill's evenness ratio for Shannon-Wiener Index

  df_stats$True_Shannon<-exp(df_stats$Shannon)#True diversity for Shannon-Wiener Index
  df_stats$True_Simpson<-1/(df_stats$Simpson)#True diversity for Simpson's Index

  return(df_stats)
}#end of function

db_exploring<-function(viral_annotations, database="all"){
  '%>%' <- tidyr::`%>%`
  if (database=="peptidase"){
    peptidase<-dplyr::select(viral_annotations, peptidase_hit, peptidase_id, peptidase_family, peptidase_RBH, peptidase_eVal, Biome, fasta)
    peptidase<-na.omit(peptidase)
    peptidase<-peptidase %>%
      dplyr::add_count(peptidase_hit,name="peptidase_counts", Biome)
    return(peptidase)

  }else if (database=="vog"){
    vog<-dplyr::select(viral_annotations, vogdb_description, Biome, fasta)
    vog<-na.omit(vog)
    vog<-vog %>%
      dplyr::add_count(vogdb_description,name="vog_counts", Biome)
    return(vog)

  }else if (database=="kegg"){
    kegg<-dplyr::select(viral_annotations, Biome, kegg_hit, kegg_id, fasta)
    kegg<-na.omit(kegg)
    kegg<-kegg %>%
      dplyr::add_count(kegg_hit,name="kegg_counts", Biome)

    kegg$Pathway<-NA
    kegg$Pathway2<-NA
    kegg$Pathway3<-NA

    for (i in 1:nrow(kegg)) {
      protein_id<-c(kegg$kegg_id[i])
      query<-tryCatch(KEGGREST::keggGet(protein_id), error=function(e) NULL)

      if (length(query[[1]]$PATHWAY)==1) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY

      } else if (length(query[[1]]$PATHWAY)==2) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]

      } else if (length(query[[1]]$PATHWAY)==3) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]
        kegg$Pathway3[i]<-query[[1]]$PATHWAY[3]
      }
    }#end of for loop

    kegg<-kegg %>%
      dplyr::add_count(Pathway,name="pathway_counts", Biome)

    return(kegg)

  }else if (database=="pfam"){
    pfam<-dplyr::select(viral_annotations, Biome, pfam_hits, pfam_id, fasta)
    pfam<-na.omit(pfam)
    pfam<-pfam %>%
      dplyr::add_count(pfam_hits,name="pfam_counts", Biome)
    return(pfam)

  }else if (database=="viral"){
    viraldb<-dplyr::select(viral_annotations, Biome, viral_function, viral_id, viral_tax, viral_identity,
                           viral_bitScore, viral_eVal, fasta)
    viraldb<-na.omit(viraldb)
    viraldb<-viraldb %>%
      dplyr::add_count(viral_function,name="viral_counts", Biome)
    return(viraldb)

  }else if (database=="all"){
    peptidase<-dplyr::select(viral_annotations, peptidase_hit, peptidase_id, peptidase_family, peptidase_RBH, peptidase_eVal, Biome)
    peptidase<-na.omit(peptidase)
    peptidase<-peptidase %>%
      dplyr::add_count(peptidase_hit,name="peptidase_counts", Biome)

    vog<-dplyr::select(viral_annotations, vogdb_description, Biome)
    vog<-na.omit(vog)
    vog<-vog %>%
      dplyr::add_count(vogdb_description,name="vog_counts", Biome)

    kegg<-dplyr::select(viral_annotations, Biome, kegg_hit, kegg_id)
    kegg<-na.omit(kegg)
    kegg<-kegg %>%
      dplyr::add_count(kegg_hit,name="kegg_counts", Biome)

    kegg$Pathway<-NA
    kegg$Pathway2<-NA
    kegg$Pathway3<-NA

    for (i in 1:nrow(kegg)) {
      protein_id<-c(kegg$kegg_id[i])
      query<-tryCatch(KEGGREST::keggGet(protein_id), error=function(e) NULL)

      if (length(query[[1]]$PATHWAY)==1) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY

      } else if (length(query[[1]]$PATHWAY)==2) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]

      } else if (length(query[[1]]$PATHWAY)==3) {
        kegg$Pathway[i]<-query[[1]]$PATHWAY[1]
        kegg$Pathway2[i]<-query[[1]]$PATHWAY[2]
        kegg$Pathway3[i]<-query[[1]]$PATHWAY[3]
      }
    }#end of for loop

    kegg<-kegg %>%
      dplyr::add_count(Pathway,name="pathway_counts", Biome)

    pfam<-dplyr::select(viral_annotations, Biome, pfam_hits, pfam_id)
    pfam<-na.omit(pfam)
    pfam<-pfam %>%
      dplyr::add_count(pfam_hits,name="pfam_counts", Biome)

    viraldb<-dplyr::select(viral_annotations, Biome, viral_function, viral_id, viral_tax, viral_identity,
                           viral_bitScore, viral_eVal)
    viraldb<-na.omit(viraldb)
    viraldb<-viraldb %>%
      dplyr::add_count(viral_function,name="viral_counts", Biome)
    #merge all db
    database_all<-Reduce(function(x,y) merge(x,y,by="Biome",all=TRUE) ,list(kegg, pfam, vog, viraldb, peptidase))
    return(database_all)

  }else {
    warning("choose one of the following: kegg, pfam, vog, viraldb, peptidase or all")
  }
}#end of function

######
# UI #
######
ui <- shinyUI(
  fluidPage(
    tags$head(


      # Note the wrapping of the string in HTML()
      tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Rajdhani:wght@600&display=swap');
      body {
        background-color: lightgrey;
        color: black;
      }
      h2 {
        font-family: 'Rajdhani', sans-serif;
      }
      .shiny-input-container {
        color: #474747;
      }"))
    ),
    shinyjs::useShinyjs(),
    div(
      id="myapp", #Start of header
      shiny::titlePanel(title = div(
        shiny::splitLayout(
          cellWidths = NULL,
          h2("MetaViral app", align = "left"),

          div(
            style = "position:absolute;top:1em; right:1em;",
            a(img(height = 50, width = 50, src="logo.png"), href="https://github.com/bbica/MetaViral") ,
          ))
      ))

    ),#end of div

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        h4("Choose the directory with your data"),

        shinyFiles::shinyDirButton(id = "dataset", label = "Select path", title="Select path"),
        #shiny::selectInput(inputId = "dataset_type", label = "File type", choices = c("DRAM-v", "vConTACT2")),
        br(),
        br(),
        shiny::radioButtons(inputId = "dataset_type", label = "File type", choices = c("DRAM-v (tsv)", "vConTACT2 (csv)")),
        shiny::actionButton(inputId = "dataset_import", label = "Upload", icon = icon("log-out", lib = "glyphicon")),
        br(),
        br(),
        shiny::actionButton(inputId = "dataset_clean", label = "Clean", icon = icon("clean", lib = "glyphicon")),
        br(),
        br(),
        hr(style = "border-color: black"),
        br(),
        shinyjs::hidden(shiny::radioButtons(inputId = "dataset_matrix", label = "Abundance type to generate", choices = c("absolute", "relative"))),
        shinyjs::hidden(shiny::selectInput(inputId = "taxa_selection", label = "Taxa to analyse", choices = c("Order", "Family", "Genus", "Genome"))),
        shinyjs::hidden(shiny::actionButton(inputId = "dataset_abundance", label = "Calculate", icon = icon("thumbnails-small", lib = "glyphicon"))),
        br(),
        br(),
        shinyjs::hidden(shiny::actionButton(inputId = "bio_stats", label = "Determine Biodiversity Indices", icon = icon("stats-circle", lib="glyphicon"))),
        #shiny::selectizeInput(inputId = "database", label = "Database exploration",
        #choices=c("kegg", "pfam", "vog", "viraldb", "peptidase", "all"), multiple=T),
        br(),
        shinyjs::hidden(shiny::selectInput(inputId = "database", label = "Database exploration",
                                           choices=c("kegg", "pfam", "vog", "viraldb", "peptidase", "all"))),
        shinyjs::hidden(shiny::actionButton(inputId = "database_explore", label = "Explore"))



      ), #end of sidebarPanel
      shiny::mainPanel(fluidRow(
        h2("Dataset"),

        div(
          style = "position:absolute;top:1em; right:1em;",
          shinyjs::hidden(selectInput(inputId ="downTypeRaw", label = "Save as:", choices=c("rds", "csv"))),
          shinyjs::hidden(downloadButton(outputId = "downloadRaw", label = "Download"))
        ),
        br(),
        br(),
        h4("Filepath"),
        verbatimTextOutput("dir", placeholder = TRUE),
        uiOutput(outputId = "output_dataset"),
        h4("Data"),
        DT::dataTableOutput(outputId = "imported_table")
      ))
    )#end of sidebarLayout
  )#end of fluidPage
)#end of shintUI
