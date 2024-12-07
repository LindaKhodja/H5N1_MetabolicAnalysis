#server.R

# Serveur (backend principal)
server <- function(input, output, session) {
  # Réactif pour stocker les données du graphe
  graph_data <- reactiveValues(
    nodes = data.frame(id = character(), label = character(), color = character(), compartment = character(), stringsAsFactors = FALSE),
    edges = data.frame(from = character(), to = character(), label = character(), stringsAsFactors = FALSE)
  )
  
  # Réactif pour suivre si le graphe est généré
  graph_generated <- reactiveVal(FALSE)
  
  # Gestion du clic sur le bouton "Add Node"
  observeEvent(input$add_node, {
    showModal(modalDialog(
      title = "Ajouter un nœud",
      textInput("node_label", "Nom du nœud :", value = ""),
      selectInput("node_color", "Couleur du nœud :", 
                  choices = c("Rouge", "Vert", "Bleu"), selected = "Rouge"),
      numericInput("node_size", "Taille du nœud :", value = 10, min = 1, max = 100),
      footer = tagList(
        modalButton("Annuler"),
        actionButton("confirm_add_node", "Valider")
      )
    ))
  })
  
  # Gestion du clic sur le bouton "Valider" dans la modale
  observeEvent(input$confirm_add_node, {
    req(input$node_label, input$node_color, input$node_size)
    
    new_node <- data.frame(
      id = if (nrow(graph_data$nodes) == 0) 1 else max(as.numeric(graph_data$nodes$id)) + 1,
      label = input$node_label,
      color = input$node_color,
      size = input$node_size,
      compartments = input$compartment,
      stringsAsFactors = FALSE
    )
    
    if (!new_node$id %in% graph_data$nodes$id) {
      graph_data$nodes <- rbind(graph_data$nodes, new_node)
      
      updateSelectInput(session, "edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "delete_edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      updateSelectInput(session, "delete_edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
      
      showNotification("Node successfully added!", type = "message")
      removeModal()
    } else {
      showNotification("Node ID already exists!", type = "error")
    }
  })
  
  # Gestion du clic sur le bouton "Add Edge"
  observeEvent(input$add_edge, {
    showModal(modalDialog(
      title = "Ajouter une arête",
      selectInput("edge_from_modal", "Depuis le nœud (Label) :", choices = setNames(graph_data$nodes$id, graph_data$nodes$label)),
      selectInput("edge_to_modal", "Vers le nœud (Label) :", choices = setNames(graph_data$nodes$id, graph_data$nodes$label)),
      textInput("edge_label_modal", "Label de l'arête :", value = ""),
      footer = tagList(
        modalButton("Annuler"),
        actionButton("confirm_add_edge", "Valider")
      )
    ))
  })
  
  observeEvent(input$confirm_add_edge, {
    req(input$edge_from_modal, input$edge_to_modal)
    
    new_edge <- data.frame(
      from = input$edge_from_modal,
      to = input$edge_to_modal,
      label = input$edge_label_modal,
      stringsAsFactors = FALSE
    )
    
    if (!any(graph_data$edges$from == input$edge_from_modal & graph_data$edges$to == input$edge_to_modal)) {
      graph_data$edges <- rbind(graph_data$edges, new_edge)
      showNotification("Edge successfully added!", type = "message")
      removeModal()
    } else {
      showNotification("Edge already exists!", type = "error")
    }
  })
  
  # Gestion du clic sur le bouton "Delete Edge"
  observeEvent(input$delete_edge, {
    showModal(modalDialog(
      title = "Supprimer une arête",
      selectInput("delete_edge_from_modal", "Depuis le nœud (Label) :", 
                  choices = setNames(graph_data$nodes$id, graph_data$nodes$label)),
      selectInput("delete_edge_to_modal", "Vers le nœud (Label) :", 
                  choices = setNames(graph_data$nodes$id, graph_data$nodes$label)),
      footer = tagList(
        modalButton("Annuler"),
        actionButton("confirm_delete_edge", "Valider")
      )
    ))
  })
  
  observeEvent(input$confirm_delete_edge, {
    req(input$delete_edge_from_modal, input$delete_edge_to_modal)
    
    edge_exists <- graph_data$edges$from == input$delete_edge_from_modal & 
      graph_data$edges$to == input$delete_edge_to_modal
    
    if (any(edge_exists)) {
      graph_data$edges <- graph_data$edges[!edge_exists, ]
      showNotification("Edge successfully deleted!", type = "message")
      removeModal()
    } else {
      showNotification("Edge does not exist!", type = "error")
    }
  })
  
  # Charger un fichier SBML
  observeEvent(input$generate_graph, {
    req(input$sbml_file)
    sbml_data <- load_sbml_data(input$sbml_file$datapath)
    graph_data$nodes <- sbml_data$nodes
    graph_data$edges <- sbml_data$edges
    graph_generated(TRUE)
    output$graph_stats <- renderText({
      paste("Nombre de nœuds :", nrow(sbml_data$nodes), "| Nombre d'arêtes :", nrow(sbml_data$edges))
    })
    
    # Mettre à jour les choix pour les sélecteurs
    updateSelectInput(session, "edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "delete_edge_from", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
    updateSelectInput(session, "delete_edge_to", choices = setNames(graph_data$nodes$id, graph_data$nodes$label))
  })
  
  
  # Légende dynamique dans la sidebar
  output$sidebar_legend <- renderUI({
    if (graph_generated()) {
      # Couleurs uniques des compartiments
      compartment_legend <- unique(graph_data$nodes[, c("compartment", "color.border"), drop = FALSE])
      
      div(
        h4("Legend"),
        h5("Compartments:"),
        tags$ul(
          lapply(1:nrow(compartment_legend), function(i) {
            tags$li(
              style = "list-style-type: none; margin-bottom: 5px; display: flex; align-items: center;",
              div(
                style = paste0(
                  "width: 15px; height: 15px; border: 1px solid ", compartment_legend$color.border[i], 
                  "; border-radius: 50%; margin-right: 10px; background-color: transparent;"
                )
              ),
              span(
                style = "color: black; font-weight: bold; font-size: 14px;",
                paste0(compartment_legend$compartment[i])
              )
            )
          })
        )
      )
    } else {
      NULL  # Rien n'est affiché si le graphe n'est pas généré
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
  
  # Tableaux interactifs
  output$nodes_table <- renderDT({
    datatable(graph_data$nodes, options = list(pageLength = 10))
  })
  
  output$edges_table <- renderDT({
    datatable(graph_data$edges, options = list(pageLength = 10))
  })
  
  output$network <- renderVisNetwork({
    req(graph_data$nodes, graph_data$edges)
    
    # Sélection du layout en fonction de l'entrée utilisateur
    layout_choice <- input$layout_choice
    
    vis_net <- visNetwork(graph_data$nodes, graph_data$edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visEdges(arrows = "to") %>%
      visInteraction(zoomView = TRUE, dragView = TRUE, multiselect = TRUE)
    
    # Ajouter le layout choisi
    if (layout_choice == "forceAtlas2Based") {
      vis_net <- vis_net %>%
        visPhysics(
          solver = "forceAtlas2Based",
          stabilization = TRUE,
          forceAtlas2Based = list(
            gravitationalConstant = -100,
            centralGravity = 0.01,
            springLength = 150,
            springConstant = 0.08,
            avoidOverlap = 0.2,
            damping = 0.4
          )
        )
    } else if (layout_choice == "barnesHut") {
      vis_net <- vis_net %>%
        visPhysics(
          solver = "barnesHut",
          stabilization = TRUE,
          barnesHut = list(
            gravitationalConstant = -20000,
            centralGravity = 0.3,
            springLength = 200,
            springConstant = 0.05,
            damping = 0.09
          )
        )
    } else if (layout_choice == "hierarchical") {
      vis_net <- vis_net %>%
        visHierarchicalLayout(
          direction = "UD",  # "UD" signifie de haut en bas (Up-Down)
          levelSeparation = 150,  # Augmenter la séparation entre les niveaux
          nodeSpacing = 150,  # Augmenter l'espacement entre les nœuds
          treeSpacing = 200,  # Espacement entre les sous-arbres
          blockShifting = TRUE,  # Ajuster les blocs pour éviter les chevauchements
          edgeMinimization = TRUE,
          parentCentralization = TRUE
        )
    } else if (layout_choice == "circular") {
      vis_net <- vis_net %>%
        visIgraphLayout(layout = "layout_in_circle")
    }
    
    # Ajouter l'option pour ajuster la position du texte des labels
    vis_net <- vis_net %>%
      visNodes(font = list(vadjust = -20))  # Ajuster la position du texte des labels des nœuds
    
    vis_net
  })
}