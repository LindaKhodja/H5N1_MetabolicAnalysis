# Charger les packages nécessaires
library(shiny)
library(shinydashboard)
library(KEGGREST)
library(visNetwork)

# Interface utilisateur de l'application
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "AnaMetabo™"),
  
  # Barre latérale pour la navigation
  dashboardSidebar(
    sidebarMenu(
      menuItem("Accueil", tabName = "home", icon = icon("home")),
      menuItem("Importation des Pathways", tabName = "import_pathway", icon = icon("database")),
      menuItem("Visualisation de Réseaux", tabName = "visualisation", icon = icon("project-diagram")),
      menuItem("Analyse Qualitative", tabName = "qualitative_analysis", icon = icon("edit"))
    )
  ),
  
  # Contenu principal
  dashboardBody(
    tabItems(
      # Page d'accueil
      tabItem(tabName = "home",
              h2("Bienvenue dans l'application AnaMetabo"),
              p("Un outil pour visualiser et effectuer des analyses qualitatives des pathways du virus H5N1."),
              # Contenu principal de la page d'accueil
              tags$div(style = "flex: 1"),
              tags$footer(
                p("Créé par Céline Hosteins, Linda Khodja, Franck Sanchez, Maroa Alani, Ilham Bedraoui"),
                style = "position: fixed; bottom: 0; width: 100%; text-align: center; background-color: #f8f9fa; padding: 10px; font-size: 12px; color: #555;"
              )
      ),
      
      # Onglet Importation des Pathways
      tabItem(tabName = "import_pathway",
              fluidRow(
                box(width = 4,
                    title = "Importer un Pathway", status = "primary",
                    selectInput("pathway_id", "Choisir un Pathway KEGG:", choices = c("hsa05164"), selected = "hsa05164"),
                    actionButton("import_kegg", "Importer le Pathway depuis KEGG")
                ),
                box(width = 8,
                    title = "Aperçu du Pathway", status = "info",
                    verbatimTextOutput("xml_preview")
                )
              )
      ),
      
      # Onglet Visualisation de Réseaux
      tabItem(tabName = "visualisation",
              fluidRow(
                box(width = 12,
                    title = "Visualisation du Pathway", status = "info",
                    visNetworkOutput("network", height = "600px")
                )
              )
      ),
      
      # Onglet Analyse Qualitative
      tabItem(tabName = "qualitative_analysis",
              h2("Analyse Qualitative du Pathway"),
              fluidRow(
                box(width = 4, 
                    title = "Options d'Annotation", status = "primary",
                    textInput("node_annotation", "Annoter un nœud (ID):", ""),
                    textAreaInput("annotation_text", "Texte de l'annotation:", "", rows = 3),
                    actionButton("add_annotation", "Ajouter l'annotation")
                ),
                box(width = 8,
                    title = "Annotations", status = "info",
                    tableOutput("annotation_table")
                )
              )
      )
    )
  )
)

# Serveur de l'application
server <- function(input, output, session) {
  # Stocker les annotations
  annotations <- reactiveVal(data.frame(Node = character(), Annotation = character(), stringsAsFactors = FALSE))
  
  # Importer et afficher un aperçu du fichier XML depuis KEGG
  observeEvent(input$import_kegg, {
    req(input$pathway_id)
    pathway_id <- input$pathway_id
    kgml_url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
    kgml_file <- tempfile(fileext = ".xml")
    download.file(kgml_url, destfile = kgml_file, mode = "wb")
    xml_content <- readLines(kgml_file, n = 10)  # Lire les 10 premières lignes pour un aperçu
    output$xml_preview <- renderPrint({
      xml_content
    })
  })
  
  # Visualisation du pathway avec visNetwork
  observeEvent(input$import_kegg, {
    req(input$pathway_id)
    pathway_id <- input$pathway_id
    kgml_url <- paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
    kgml_file <- tempfile(fileext = ".xml")
    download.file(kgml_url, destfile = kgml_file, mode = "wb")
    
    # Parser le fichier XML pour extraire les nœuds et les arêtes (simplification)
    nodes <- data.frame(id = 1:5, label = paste("Nœud", 1:5))  # Exemple de nœuds
    edges <- data.frame(from = c(1, 2, 3), to = c(2, 3, 4))  # Exemple d'arêtes
    
    output$network <- renderVisNetwork({
      visNetwork(nodes, edges) %>%
        visNodes(shape = "dot", size = 15) %>%
        visEdges(arrows = "to") %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
    })
  })
  
  # Ajouter des annotations
  observeEvent(input$add_annotation, {
    req(input$node_annotation, input$annotation_text)
    new_annotation <- data.frame(Node = input$node_annotation, Annotation = input$annotation_text, stringsAsFactors = FALSE)
    annotations(rbind(annotations(), new_annotation))
  })
  
  # Afficher les annotations
  output$annotation_table <- renderTable({
    annotations()
  })
}

# Créer l'application Shiny
shinyApp(ui, server)
