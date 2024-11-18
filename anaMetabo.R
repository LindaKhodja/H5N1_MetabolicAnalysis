# Charger les packages nécessaires
library(shiny)
library(shinydashboard)
library(KEGGREST)
library(visNetwork)
library(xml2)
library(dplyr)
library(DT)

# Fonction pour charger et traiter le fichier SBML
load_sbml_data <- function(sbml_file) {
  doc <- read_xml(sbml_file)
  ns <- xml_ns_rename(xml_ns(doc), d1 = "sbml")
  
  # Extraire les métabolites (nœuds)
  species_nodes <- xml_find_all(doc, ".//sbml:listOfSpecies/sbml:species", ns)
  nodes <- data.frame(
    id = xml_attr(species_nodes, "id"),
    label = xml_attr(species_nodes, "name"),
    stringsAsFactors = FALSE
  )
  
  # Extraire les réactions (arêtes)
  reaction_nodes <- xml_find_all(doc, ".//sbml:listOfReactions/sbml:reaction", ns)
  edges <- data.frame()
  
  for (reaction in reaction_nodes) {
    reaction_id <- xml_attr(reaction, "id")
    reaction_name <- xml_attr(reaction, "name")
    
    # Réactifs et produits
    reactants <- xml_find_all(reaction, ".//sbml:listOfReactants/sbml:speciesReference", ns)
    products <- xml_find_all(reaction, ".//sbml:listOfProducts/sbml:speciesReference", ns)
    
    reactant_ids <- xml_attr(reactants, "species")
    product_ids <- xml_attr(products, "species")
    
    # Créer des liens pour chaque combinaison réactif-produit
    for (reactant in reactant_ids) {
      for (product in product_ids) {
        edges <- rbind(edges, data.frame(from = reactant, to = product, label = reaction_name, stringsAsFactors = FALSE))
      }
    }
  }
  
  list(nodes = nodes, edges = edges)
}

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
              tags$footer(
                p("Créé par Céline Hosteins, Linda Khodja, Franck Sanchez, Maroa Alani, Ilham Bedraoui"),
                style = "position: fixed; bottom: 0; width: 100%; text-align: center; background-color: #f8f9fa; padding: 10px; font-size: 12px; color: #555;"
              )
      ),
      
      # Onglet Visualisation de Réseaux
      tabItem(tabName = "visualisation",
              fluidRow(
                box(width = 4,
                    title = "Charger le graphe",
                    status = "primary",
                    h4("Depuis un fichier SBML"),
                    fileInput("sbmlFile", "Sélectionner un fichier SBML", accept = c(".sbml", ".xml")),
                    actionButton("loadGraphFromSBML", "Charger depuis SBML"),
                    hr(),
                    h4("Depuis des fichiers CSV"),
                    fileInput("nodesCSV", "Fichier des nœuds (nodes.csv)", accept = ".csv"),
                    fileInput("edgesCSV", "Fichier des arêtes (edges.csv)", accept = ".csv"),
                    actionButton("loadGraphFromCSV", "Charger depuis CSV"),
                    hr(),
                    actionButton("saveGraph", "Sauvegarder le graphe")
                ),
                box(width = 8,
                    title = "Visualisation du Réseau",
                    status = "info",
                    visNetworkOutput("network", height = "600px"),
                    textOutput("graph_stats"),
                    hr(),
                    actionButton("addNode", "Ajouter un nœud"),
                    actionButton("deleteNode", "Supprimer un nœud sélectionné"),
                    actionButton("addEdge", "Ajouter une arête"),
                    actionButton("deleteEdge", "Supprimer une arête sélectionnée")
                )
              ),
              fluidRow(
                box(width = 6, title = "Table des Nœuds", status = "info", DTOutput("nodes_table")),
                box(width = 6, title = "Table des Arêtes", status = "info", DTOutput("edges_table"))
              )
      )
    )
  )
)

# Serveur de l'application
server <- function(input, output, session) {
  
  # Réactif pour stocker les données du graphe
  graph_data <- reactiveVal(list(nodes = data.frame(id = character(), label = character(), stringsAsFactors = FALSE),
                                 edges = data.frame(from = character(), to = character(), label = character(), stringsAsFactors = FALSE)))
  
  # Fonction pour mettre à jour les statistiques du graphe
  update_graph_stats <- function(nodes, edges) {
    paste("Nombre de nœuds :", nrow(nodes), "| Nombre d'arêtes :", nrow(edges))
  }
  
  # Charger un fichier SBML
  observeEvent(input$loadGraphFromSBML, {
    req(input$sbmlFile)
    sbml_data <- load_sbml_data(input$sbmlFile$datapath)
    graph_data(sbml_data)
    output$graph_stats <- renderText({
      update_graph_stats(sbml_data$nodes, sbml_data$edges)
    })
    output$network <- renderVisNetwork({
      visNetwork(sbml_data$nodes, sbml_data$edges) %>%
        visEdges(arrows = "to") %>%
        visNodes(color = list(border = "black", background = "#97C2FC")) %>%
        visLayout(randomSeed = 42) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
    })
  })
  
  # Charger depuis des fichiers CSV
  observeEvent(input$loadGraphFromCSV, {
    req(input$nodesCSV, input$edgesCSV)
    nodes_data <- read.csv(input$nodesCSV$datapath, stringsAsFactors = FALSE)
    edges_data <- read.csv(input$edgesCSV$datapath, stringsAsFactors = FALSE)
    graph_data(list(nodes = nodes_data, edges = edges_data))
    output$graph_stats <- renderText({
      update_graph_stats(nodes_data, edges_data)
    })
    output$network <- renderVisNetwork({
      visNetwork(nodes_data, edges_data) %>%
        visEdges(arrows = "to") %>%
        visNodes(color = list(border = "black", background = "#97C2FC")) %>%
        visLayout(randomSeed = 42) %>%
        visInteraction(navigationButtons = TRUE, nodesIdSelection = TRUE)
    })
  })
  
  # Affichage des tables
  output$nodes_table <- renderDT({
    datatable(graph_data()$nodes, options = list(pageLength = 10), selection = "single")
  })
  output$edges_table <- renderDT({
    datatable(graph_data()$edges, options = list(pageLength = 10), selection = "single")
  })
}

# Lancer l'application
shinyApp(ui, server)