library(shiny)
library(visNetwork)
library(xml2)
library(RColorBrewer)
library(DT)

# Fonction pour générer des couleurs dynamiques
generate_colors <- function(n) {
  max_colors <- 12
  if (n <= max_colors) {
    brewer.pal(n, "Set3")
  } else {
    colorRampPalette(brewer.pal(max_colors, "Set3"))(n)
  }
}

# Fonction pour créer un fichier SBML à partir des données du graphe
save_sbml <- function(nodes, edges, file_path) {
  sbml_doc <- xml_new_root("sbml", version = "1.2", xmlns = "http://www.sbml.org/sbml/level2/version4")
  
  # Ajouter la liste des espèces (nœuds)
  list_of_species <- xml_add_child(sbml_doc, "listOfSpecies")
  for (i in 1:nrow(nodes)) {
    species <- xml_add_child(list_of_species, "species",
                             id = nodes$id[i],
                             name = nodes$label[i])
  }
  
  # Ajouter la liste des réactions (arêtes)
  list_of_reactions <- xml_add_child(sbml_doc, "listOfReactions")
  for (i in 1:nrow(edges)) {
    reaction <- xml_add_child(list_of_reactions, "reaction",
                              id = paste0("reaction_", i),
                              name = ifelse(!is.na(edges$label[i]), edges$label[i], paste0("reaction_", i)))
    list_of_reactants <- xml_add_child(reaction, "listOfReactants")
    xml_add_child(list_of_reactants, "speciesReference", species = edges$from[i])
    list_of_products <- xml_add_child(reaction, "listOfProducts")
    xml_add_child(list_of_products, "speciesReference", species = edges$to[i])
  }
  
  # Sauvegarder le fichier
  write_xml(sbml_doc, file_path)
}

# Fonction pour charger et traiter le fichier SBML
load_sbml_data <- function(sbml_file) {
  doc <- read_xml(sbml_file)
  ns <- xml_ns_rename(xml_ns(doc), d1 = "sbml")
  
  # Extraire les métabolites (nœuds) sans les compartiments
  species_nodes <- xml_find_all(doc, ".//sbml:listOfSpecies/sbml:species", ns)
  nodes <- data.frame(
    id = xml_attr(species_nodes, "id"),
    label = xml_attr(species_nodes, "name"),
    shape = "dot",  # Mettre les espèces sous forme de cercle
    color = "blue", # Couleur par défaut des espèces
    stringsAsFactors = FALSE
  )
  
  # Extraire le compartiment à partir de l'ID
  nodes$group <- ifelse(
    grepl("\\[.*\\]$", nodes$label),
    gsub(".*\\[(.*)\\]$", "\\1", nodes$label),
    "unknown"
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
      color = "red",     # Couleur rouge pour les réactions
      group = "reaction",
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
        label = "reactant",
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
        label = "product",
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
        if (sbo_term == "SBO:0000170") {  # Activation (positive)
          color <- "blue"
          dashes <- TRUE
        } else if (sbo_term == "SBO:0000169") {  # Inhibition (négative)
          color <- "orange"
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
        label = "modulator",
        arrows = "to",
        dashes = dashes,
        color = color,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  list(nodes = nodes, edges = edges)
}

# Interface utilisateur (UI)
ui <- fluidPage(
  titlePanel("Reactome SBML Graph Visualization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("sbml_file", "Upload an SBML File (.xml)", accept = c(".xml", ".sbml")),
      actionButton("generate_graph", "Generate Graph"),
      tags$hr(),
      h4("Modify Graph"),
      textInput("add_node_label", "New Node Label"),
      textInput("add_node_group", "New Node Group"),
      actionButton("add_node", "Add Node"),
      actionButton("delete_node", "Delete Selected Node"),
      tags$hr(),
      h4("Add Edge"),
      selectInput("edge_from", "From Node (Label)", choices = NULL),
      selectInput("edge_to", "To Node (Label)", choices = NULL),
      textInput("edge_label", "Edge Label"),
      actionButton("add_edge", "Add Edge"),
      tags$hr(),
      h4("Delete Edge"),
      selectInput("delete_edge_from", "From Node (Label)", choices = NULL),
      selectInput("delete_edge_to", "To Node (Label)", choices = NULL),
      actionButton("delete_edge", "Delete Edge"),
      tags$hr(),
      h4("Save Graph"),
      textInput("save_filename", "File Name", value = "modified_graph.xml"),
      downloadButton("download_sbml", "Download SBML File"),
      tags$hr(),
      h4("Graph Summary"),
      textOutput("graph_summary"),
      tags$hr(),
      h4("Cellular Compartments"),
      uiOutput("legend")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Graph Visualization",
                 visNetworkOutput("network", height = "700px")),
        tabPanel("Nodes and Edges",
                 h4("Nodes Table"),
                 DTOutput("nodes_table"),
                 h4("Edges Table"),
                 DTOutput("edges_table"))
      )
    )
  )
)

# Serveur (backend principal)
server <- function(input, output, session) {
  # Réactif pour stocker les données du graphe
  graph_data <- reactiveValues(
    nodes = data.frame(id = character(), label = character(), group = character(), stringsAsFactors = FALSE),
    edges = data.frame(from = character(), to = character(), label = character(), stringsAsFactors = FALSE)
  )
  
  # Charger les données initiales du SBML
  observeEvent(input$generate_graph, {
    req(input$sbml_file)
    file_path <- input$sbml_file$datapath
    data <- load_sbml_data(file_path)
    graph_data$nodes <- data$nodes
    graph_data$edges <- data$edges
    
    # Mettre à jour les choix pour les sélecteurs en utilisant les labels
    updateSelectInput(session, "edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "delete_edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "delete_edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
  })
  
  # Ajouter un nœud
  observeEvent(input$add_node, {
    req(input$add_node_label, input$add_node_group)
    new_id <- paste0("node_", nrow(graph_data$nodes) + 1)
    new_node <- data.frame(
      id = new_id,
      label = input$add_node_label,
      group = input$add_node_group,
      stringsAsFactors = FALSE
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
      showNotification("Node successfully added!", type = "message")
    } else {
      showNotification("Node ID already exists!", type = "error")
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
      showNotification("Edge successfully added!", type = "message")
    } else {
      showNotification("Edge already exists!", type = "error")
    }
  })
  
  # Supprimer une arête
  observeEvent(input$delete_edge, {
    req(input$delete_edge_from, input$delete_edge_to)
    edge_exists <- graph_data$edges$from == input$delete_edge_from & graph_data$edges$to == input$delete_edge_to
    if (any(edge_exists)) {
      graph_data$edges <- graph_data$edges[!edge_exists, ]
      showNotification("Edge successfully deleted!", type = "message")
    } else {
      showNotification("Edge does not exist!", type = "error")
    }
  })
  
  # Supprimer un nœud sélectionné
  observeEvent(input$delete_node, {
    selected_node <- input$network_selectedNodes
    if (!is.null(selected_node)) {
      graph_data$nodes <- graph_data$nodes[graph_data$nodes$id != selected_node, ]
      graph_data$edges <- graph_data$edges[graph_data$edges$from != selected_node & graph_data$edges$to != selected_node, ]
      # Mettre à jour les choix pour les sélecteurs
      updateSelectInput(session, "edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "delete_edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "delete_edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      showNotification("Node successfully deleted!", type = "message")
    } else {
      showNotification("No node selected!", type = "error")
    }
  })
  
  # Télécharger le fichier SBML
  output$download_sbml <- downloadHandler(
    filename = function() {
      input$save_filename
    },
    content = function(file) {
      save_sbml(graph_data$nodes, graph_data$edges, file)
    }
  )
  
  # Résumé du graphe
  output$graph_summary <- renderText({
    num_nodes <- nrow(graph_data$nodes)
    num_edges <- nrow(graph_data$edges)
    paste("Nodes:", num_nodes, "| Edges:", num_edges)
  })
  
  # Légende
  output$legend <- renderUI({
    unique_groups <- unique(graph_data$nodes$group)
    colors <- generate_colors(length(unique_groups))
    prefix_colors <- setNames(colors, unique_groups)
    
    legend_items <- lapply(unique_groups, function(group) {
      div(
        style = sprintf("margin-bottom: 5px; padding: 5px; background-color: %s; color: black; border-radius: 5px;", prefix_colors[group]),
        group
      )
    })
    
    do.call(tagList, legend_items)
  })
  
  # Tableaux interactifs
  output$nodes_table <- renderDT({
    datatable(graph_data$nodes, options = list(pageLength = 10))
  })
  
  output$edges_table <- renderDT({
    datatable(graph_data$edges, options = list(pageLength = 10))
  })
  
  # Graphe interactif
  output$network <- renderVisNetwork({
    req(graph_data$nodes, graph_data$edges)
    unique_groups <- unique(graph_data$nodes$group)
    colors <- generate_colors(length(unique_groups))
    prefix_colors <- setNames(colors, unique_groups)
    graph_data$nodes$color <- prefix_colors[graph_data$nodes$group]
    
    visNetwork(graph_data$nodes, graph_data$edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visEdges(arrows = "to") %>%
      visInteraction(zoomView = TRUE, dragView = TRUE, multiselect = TRUE) %>%
      visEvents(select = "function(properties) {
                   Shiny.setInputValue('network_selectedNodes', properties.nodes[0], {priority: 'event'});
                 }")
  })
}

# Exécuter l'application Shiny
shinyApp(ui = ui, server = server)
