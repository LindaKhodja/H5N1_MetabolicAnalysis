# Charger les packages nécessaires
library(shiny)
library(shinydashboard)
library(KEGGREST)
library(visNetwork)
library(xml2)
library(dplyr)
library(DT)

# Définir des couleurs de bordures spécifiques aux compartiments
compartment_colors <- c(
  "cytosol" = "blue",
  "endoplasmic reticulum membrane" = "green",
  "golgi membrane" = "orange",
  "endocytic vesicle membrane" = "purple",
  "nucleoplasm" = "red",
  "endosome lumen" = "pink",
  "plasma membrane" = "yellow",
  "unknown" = "black"
)

# Fonction pour charger et traiter le fichier SBML
load_sbml_data <- function(sbml_file) {
  doc <- read_xml(sbml_file)
  ns <- xml_ns_rename(xml_ns(doc), d1 = "sbml")
  
  # Extraire les métabolites (nœuds)
  species_nodes <- xml_find_all(doc, ".//sbml:listOfSpecies/sbml:species", ns)
  nodes <- data.frame(
    id = xml_attr(species_nodes, "id"),
    label = xml_attr(species_nodes, "name"),
    shape = "dot",  # Mettre les espèces sous forme de cercle
    color.background = "lightgrey", # Couleur par défaut des espèces
    size = 25,
    font.size = 20,
    stringsAsFactors = FALSE
  )
  
  # Extraire le compartiment à partir du label du nœud
  nodes$compartment <- ifelse(
    grepl("\\[.*\\]$", nodes$label),
    gsub(".*\\[(.*)\\]$", "\\1", nodes$label),
    "unknown"
  )
  
  nodes$color.border <- ifelse(
    nodes$compartment %in% names(compartment_colors),
    compartment_colors[nodes$compartment],
    "black"
  )
  
  # Extraire les réactions (liens)
  reaction_nodes <- xml_find_all(doc, ".//sbml:listOfReactions/sbml:reaction", ns)
  
  # Ajouter les nœuds de réaction (carrés rouges)
  for (reaction in reaction_nodes) {
    reaction_id <- xml_attr(reaction, "id")
    reaction_name <- xml_attr(reaction, "name")
    
    new_node <- data.frame(
      id = reaction_id,
      label = ifelse(!is.na(reaction_name), reaction_name, reaction_id),
      shape = "square",  # Forme carrée pour les réactions
      color.background = "red",     # Couleur rouge pour les réactions
      size = 10,  # Taille réduite
      font.size = 12,  # Taille de la police du label
      color.border = "black",
      compartment = "reaction",
      stringsAsFactors = FALSE
    )
    nodes <- rbind(nodes, new_node)
  }
  
  # Extraire les liens entre espèces et réactions (arêtes)
  edges <- data.frame()
  
  for (reaction in reaction_nodes) {
    reaction_id <- xml_attr(reaction, "id")
    
    # Réactifs
    reactants <- xml_find_all(reaction, ".//sbml:listOfReactants/sbml:speciesReference", ns)
    reactant_ids <- xml_attr(reactants, "species")
    for (reactant in reactant_ids) {
      edges <- rbind(edges, data.frame(
        from = reactant,
        to = reaction_id,
        arrows = "to",
        dashes = FALSE,
        color = "black",
        stringsAsFactors = FALSE
      ))
    }
    
    # Produits
    products <- xml_find_all(reaction, ".//sbml:listOfProducts/sbml:speciesReference", ns)
    product_ids <- xml_attr(products, "species")
    for (product in product_ids) {
      edges <- rbind(edges, data.frame(
        from = reaction_id,
        to = product,
        arrows = "to",
        dashes = FALSE,
        color = "black",
        stringsAsFactors = FALSE
      ))
    }
    # Régulateurs (modulateurs)
    modifiers <- xml_find_all(reaction, ".//sbml:listOfModifiers/sbml:modifierSpeciesReference", ns)
    modifier_ids <- xml_attr(modifiers, "species")
    for (modifier in modifiers) {
      sbo_term <- xml_attr(modifier, "sboTerm")  # Extraire le terme SBO pour déterminer le type de régulation
      
      # Déterminer la couleur et le type de ligne pour les régulations
      if (!is.na(sbo_term)) {
        if (sbo_term == "SBO:0000013") {  # Activation (positive)
          color <- "green"
          dashes <- TRUE
        } else if (sbo_term == "SBO:0000020") {  # Inhibition (négative)
          color <- "red"
          dashes <- TRUE
        } else {
          color <- "black"
          dashes <- FALSE
        }
      } else {
        color <- "black"
        dashes <- FALSE
      }
      
      edges <- rbind(edges, data.frame(
        from = modifier_ids,
        to = reaction_id,
        arrows = "to",
        dashes = dashes,
        color = color,
        stringsAsFactors = FALSE
      ))
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
      menuItem("Table des Nœuds et Arêtes", tabName = "table_nodes_edges", icon = icon("table"))
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
                    title = "Charger et Modifier le Graphe",
                    status = "primary",
                    fileInput("sbml_file", "Sélectionner un fichier SBML", accept = c(".sbml", ".xml")),
                    actionButton("generate_graph", "Charger le Graphe"),
                    hr(),
                    selectInput("layout_choice", "Choisir un layout:", choices = c(
                      "Force Atlas 2 Based" = "forceAtlas2Based",
                      "Barnes Hut" = "barnesHut",
                      "Hierarchical" = "hierarchical",
                      "Circular" = "circular"
                    )),
                    hr(),
                    h4("Modifier le Graphe"),
                    textInput("add_node_label", "Label du Nouveau Nœud"),
                    actionButton("add_node", "Ajouter un Nœud"),
                    actionButton("delete_node", "Supprimer le Nœud Sélectionné"),
                    hr(),
                    h4("Ajouter une Arête"),
                    selectInput("edge_from", "Depuis le Nœud (Label)", choices = NULL),
                    selectInput("edge_to", "Vers le Nœud (Label)", choices = NULL),
                    textInput("edge_label", "Label de l'Arête"),
                    actionButton("add_edge", "Ajouter une Arête"),
                    hr(),
                    h4("Supprimer une Arête"),
                    selectInput("delete_edge_from", "Depuis le Nœud (Label)", choices = NULL),
                    selectInput("delete_edge_to", "Vers le Nœud (Label)", choices = NULL),
                    actionButton("delete_edge", "Supprimer l'Arête"),
                    hr(),
                    actionButton("saveGraph", "Sauvegarder le Graphe")
                ),
                box(width = 8,
                    title = "Visualisation du Réseau",
                    status = "info",
                    visNetworkOutput("network", height = "600px"),
                    textOutput("graph_stats")
                )
              )
      ),
      
      # Onglet Table des Nœuds et Arêtes
      tabItem(tabName = "table_nodes_edges",
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
  graph_data <- reactiveValues(
    nodes = data.frame(id = character(), label = character(), color = character(), compartment = character(), stringsAsFactors = FALSE),
    edges = data.frame(from = character(), to = character(), label = character(), stringsAsFactors = FALSE)
  )
  
  # Charger un fichier SBML
  observeEvent(input$generate_graph, {
    req(input$sbml_file)
    sbml_data <- load_sbml_data(input$sbml_file$datapath)
    graph_data$nodes <- sbml_data$nodes
    graph_data$edges <- sbml_data$edges
    output$graph_stats <- renderText({
      paste("Nombre de nœuds :", nrow(sbml_data$nodes), "| Nombre d'arêtes :", nrow(sbml_data$edges))
    })
    
    # Mettre à jour les choix pour les sélecteurs
    updateSelectInput(session, "edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "delete_edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "delete_edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
  })
  
  # Ajouter un nœud
  observeEvent(input$add_node, {
    req(input$add_node_label)
    new_id <- paste0("node_", nrow(graph_data$nodes) + 1)
    new_node <- data.frame(
      id = new_id,
      label = input$add_node_label,
      stringsAsFactors = FALSE
    )
    # Extraire le compartiment à partir du label du nœud
    new_node$compartment <- ifelse(
      grepl("\\[.*\\]$", new_node$label),
      gsub(".*\\[(.*)\\]$", "\\1", new_node$label),
      "unknown"
    )
    new_node$color.border <- ifelse(
      new_node$compartment %in% names(compartment_colors),
      compartment_colors[new_node$compartment],
      "black"
    )
    # Harmonisation des colonnes
    missing_cols <- setdiff(names(graph_data$nodes), names(new_node))
    new_node[missing_cols] <- NA  # Ajouter les colonnes manquantes
    new_node <- new_node[, names(graph_data$nodes)]  # Réordonner les colonnes
    
    # Vérification et ajout
    if (!new_id %in% graph_data$nodes$id) {
      graph_data$nodes <- rbind(graph_data$nodes, new_node)
      # Mettre à jour les choix pour les sélecteurs
      updateSelectInput(session, "edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "delete_edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "delete_edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      showNotification("Nœud ajouté avec succès!", type = "message")
    } else {
      showNotification("L'ID du nœud existe déjà!", type = "error")
    }
  })
  
  # Ajouter une arête
  observeEvent(input$add_edge, {
    req(input$edge_from, input$edge_to)
    new_edge <- data.frame(
      from = input$edge_from,
      to = input$edge_to,
      label = input$edge_label,
      stringsAsFactors = FALSE
    )
    # Vérification et ajout
    if (!any(graph_data$edges$from == input$edge_from & graph_data$edges$to == input$edge_to)) {
      graph_data$edges <- rbind(graph_data$edges, new_edge)
      showNotification("Arête ajoutée avec succès!", type = "message")
    } else {
      showNotification("L'arête existe déjà!", type = "error")
    }
  })
  
  # Supprimer une arête
  observeEvent(input$delete_edge, {
    req(input$delete_edge_from, input$delete_edge_to)
    edge_exists <- graph_data$edges$from == input$delete_edge_from & graph_data$edges$to == input$delete_edge_to
    if (any(edge_exists)) {
      graph_data$edges <- graph_data$edges[!edge_exists, ]
      showNotification("Arête supprimée avec succès!", type = "message")
    } else {
      showNotification("L'arête n'existe pas!", type = "error")
    }
  })
  
  # Graphe interactif avec choix de layout
  output$network <- renderVisNetwork({
    req(graph_data$nodes, graph_data$edges)
    
    # Sélection du layout en fonction de l'entrée utilisateur
    layout_choice <- input$layout_choice
    
    vis_net <- visNetwork(graph_data$nodes, graph_data$edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visEdges(arrows = "to") %>%
      visInteraction(zoomView = TRUE, dragView = TRUE, multiselect = TRUE)
    
    if (layout_choice == "forceAtlas2Based") {
      vis_net <- vis_net %>%
        visPhysics(
          solver = "forceAtlas2Based",
          stabilization = TRUE,
          forceAtlas2Based = list(
            gravitationalConstant = -50,
            centralGravity = 0.01,
            springLength = 100,
            springConstant = 0.08,
            damping = 0.4
          )
        )
    } else if (layout_choice == "barnesHut") {
      vis_net <- vis_net %>%
        visPhysics(
          solver = "barnesHut",
          stabilization = TRUE
        )
    } else if (layout_choice == "hierarchical") {
      vis_net <- vis_net %>%
        visHierarchicalLayout()
    } else if (layout_choice == "circular") {
      vis_net <- vis_net %>%
        visIgraphLayout(layout = "layout_in_circle")
    }
    
    vis_net
  })
  
  # Affichage des tables dans l'onglet Table des Nœuds et Arêtes
  output$nodes_table <- renderDT({
    datatable(graph_data$nodes, options = list(pageLength = 10))
  })
  output$edges_table <- renderDT({
    datatable(graph_data$edges, options = list(pageLength = 10))
  })
}

# Lancer l'application
shinyApp(ui, server)