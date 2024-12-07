#utils.R

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

# Définir une palette Spectral avec des couleurs natives de R
compartment_colors <- c(
  "cytosol" = "darkred", # Rouge vif
  "extracellular region" = "purple",
  "endoplasmic reticulum membrane" = "darkorange", # Orange foncé
  "endoplasmic reticulum lumen" = "magenta",
  "Golgi membrane" = "gold", # Jaune doré
  "endocytic vesicle membrane" = "yellowgreen", # Vert jaunâtre
  "nucleoplasm" = "pink", # Rose vif
  "endosome lumen" = "cyan", # Cyan
  "plasma membrane" = "skyblue", # Bleu ciel
  "unknown" = "blue", # Bleu foncé
  "nuclear envelope" = "darkgreen", # Violet
  "mitochondrial inner membrane" = "darkslategray", # Gris-vert sombre
  "mitochondrial intermembrane space" = "lightcoral" # Rouge clair
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
    color.background = "#F0F0F0", # Couleur par défaut des espèces
    size = 50,
    font.size = 20,
    stringsAsFactors = FALSE
  )
  
  # Vérifier le sboTerm pour chaque nœud et modifier la forme en fonction
  sbo_terms <- xml_attr(species_nodes, "sboTerm")
  
  # Modifier les nœuds ayant le terme SBO:0000278 en "diamond"
  nodes$shape <- ifelse(sbo_terms == "SBO:0000278", "diamond", nodes$shape) #rna
  
  nodes$shape <- ifelse(sbo_terms == "SBO:0000252", "dot", nodes$shape) #gene
  
  nodes$shape <- ifelse(sbo_terms == "SBO:000013", "triangleDown", nodes$shape) 
  
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
  
  # Modifier le label pour supprimer la partie entre crochets (ex. "name [compartment]" devient "name")
  nodes$label <- gsub("\\[.*\\]$", "", nodes$label)
  nodes$label <- trimws(nodes$label)  # Supprimer les espaces blancs au début ou à la fin si nécessaire
  
  # Extraire les réactions (liens)
  reaction_nodes <- xml_find_all(doc, ".//sbml:listOfReactions/sbml:reaction", ns)
  
  # Ajouter les nœuds de réaction (carrés rouges)
  for (reaction in reaction_nodes) {
    reaction_id <- xml_attr(reaction, "id")
    reaction_name <- xml_attr(reaction, "name")
    
    #Noeud des réactions
    new_node <- data.frame(
      id = reaction_id,
      label = ifelse(!is.na(reaction_name), reaction_name, reaction_id),
      shape = "square",  # Forme carrée pour les réactions
      color.background = "red",     # Couleur rouge pour les réactions
      size = 10,  # Taille réduite
      font.size = 20,  # Taille de la police du label
      color.border = "black",
      compartment = "reaction", #par défaut nouveau noeud dans cytosol
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
        if (sbo_term == "SBO:0000459") {  # Activation (positive)
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

