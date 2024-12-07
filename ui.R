#ui.R

# Interface Utilisateur (UI)
ui <- dashboardPage(
  skin = "purple",
  
  dashboardHeader(title = "AnaMetabo™"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Accueil", tabName = "home", icon = icon("home")),
      menuItem("Visualisation de Réseaux", tabName = "network_vis", icon = icon("project-diagram"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Onglet Accueil
      tabItem(tabName = "home",
              h2("Bienvenue dans l'application AnaMetabo"),
              p("Un outil pour visualiser et effectuer des analyses qualitatives des pathways du virus H5N1."),
              tags$footer(
                p("Créé par Céline Hosteins, Linda Khodja, Franck Sanchez, Maroa Alani, Ilham Bedraoui"),
                style = "position: fixed; bottom: 0; width: 100%; text-align: center; background-color: #f8f9fa; padding: 10px; font-size: 12px; color: #555;"
              )
      ),
      
      # Onglet Visualisation de Réseaux
      tabItem(tabName = "network_vis",
              fluidPage(
                titlePanel("Visualisation de Réseaux SBML"),
                
                # Corps principal
                sidebarLayout(
                  sidebarPanel(
                    fileInput("sbml_file", "Importer un fichier SBML (.xml)", accept = c(".xml", ".sbml")),
                    actionButton("generate_graph", "Générer le Graphe"),
                    tags$hr(),
                    
                    selectInput("layout_choice", "Choisir un Layout :", 
                                choices = c(
                                  "Force Atlas 2 Based" = "forceAtlas2Based", 
                                  "Barnes Hut" = "barnesHut", 
                                  "Hiérarchique" = "hierarchical", 
                                  "Circulaire" = "circular"),
                                selected = "barnesHut"  # Définir Barnes Hut comme valeur par défaut
                    ),
                    tags$hr(),
                    
                    # Légende conditionnelle
                    uiOutput("sidebar_legend"),  # Conteneur pour la légende
                    width = 3  # Ajuste la largeur de la sidebar
                  ),
                  
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Graph Visualization",
                               div(style = "position: relative;",
                                   downloadButton("download_sbml", label = NULL, 
                                                  icon = icon("download"), 
                                                  style = "position: absolute; top: 10px; right: 10px; z-index: 1000;"),
                                   
                                   div(style = "text-align: center; font-size: 16px; margin-top: 10px;", 
                                       textOutput("graph_summary")),
                                   
                                   visNetworkOutput("network", height = "580px")
                               )
                      ),
                      tabPanel("Nodes and Edges",
                               h4("Nodes Table"),
                               DTOutput("nodes_table"),
                               h4("Edges Table"),
                               DTOutput("edges_table"))
                    ),
                    width = 9
                  )
                ),
                
                # Footer
                tags$footer(
                  style = "background-color: #f8f9fa; padding: 15px; border-top: 1px solid #dee2e6; text-align: center;",
                  h4("Graph Modification"),
                  div(
                    actionButton("add_node", "Add Node", style = "margin-right: 10px;"),
                    actionButton("delete_node", "Delete Selected Node", style = "margin-right: 10px;"),
                    actionButton("add_edge", "Add Edge", style = "margin-right: 10px;"),
                    actionButton("delete_edge", "Delete Edge"),
                    style = "margin-top: 10px;"
                  )
                )
              )
      )
    )
  )
)