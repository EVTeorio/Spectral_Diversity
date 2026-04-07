
# ============================================================
# Interactive Taxonomic Tree: Order → Family → Genus → Species
# ============================================================

# 1. Install & load packages

install.packages("igraph")

library(tidyverse)
library(dplyr)
library(collapsibleTree)
library(ggraph)
library(igraph)
library(data.tree)
library(visNetwork)

# 1. Load CSV with UTF-8 encoding
heirachy <- read.csv("SOFOR/big_data/heirchy table prfdp.csv")
heirachy <- heirachy[,-1]
# unique(heirachy$Class)
#heirachy <- heirachy[!is.na(heirachy$Species) & heirachy$Species != "", ]


# 2. Remove/replace non-ASCII characters
heirachy[] <- lapply(heirachy, function(x) iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT"))

#collapsible hierarchical tree
collapsibleTree(
  heirachy,
  hierarchy = c("Class", "Order", "Family", "Genus", "Species"),
  root = "All Tree Species",
  width = 900,
  height = 700,
  zoomable = TRUE,
  collapsed = TRUE,      
  fill = "lightgreen",
)


collapsibleTree(
  heirachy,
  hierarchy = c("Class", "Order", "Family", "Genus", "sp_code"),
  root = "All Tree Species",
  width = 500,
  height = 1000,
  zoomable = TRUE,
  collapsed = FALSE,
  fill = "lightgreen",
  fontSize = 12,
  direction = "v"   # 'v' for vertical (top-down), 'h' for horizontal
)

# ggraph ########################
# -----------------------------
# 2. Create root node
# -----------------------------
tax_tree <- Node$new("All Tree Species")

# -----------------------------
# 3. Function to get existing child or create new one
# -----------------------------
get_or_add_child <- function(parent, name) {
  existing <- parent$children[names(parent$children) == name]
  if (length(existing) > 0) {
    return(existing[[1]])
  } else {
    return(parent$AddChild(name))
  }
}

# -----------------------------
# 4. Build tree 
# -----------------------------
for (i in 1:nrow(heirachy)) {
  class_node  <- get_or_add_child(tax_tree, heirachy$Class[i])
  order_node  <- get_or_add_child(class_node, heirachy$Order[i])
  family_node <- get_or_add_child(order_node, heirachy$Family[i])
  genus_node  <- get_or_add_child(family_node, heirachy$Genus[i])
  species_node<- get_or_add_child(genus_node, heirachy$Species[i])
  #code_node<- get_or_add_child(genus_node, heirachy$sp_code[i])
}

# -----------------------------
# 5. Convert tree to network for visNetwork
# -----------------------------
network <- ToDataFrameNetwork(tax_tree, nameCol = "name")

nodes <- data.frame(
  id = unique(c(network$from, network$to)),
  label = unique(c(network$from, network$to)),
  title = unique(c(network$from, network$to))
)

nodes$label <- sub(".*/", "", nodes$label)


edges <- network %>% rename(from = from, to = to) %>% distinct()

# -----------------------------
# 6. Draw interactive top-down tree
# -----------------------------
visNetwork(nodes, edges, width = "100%", height = "800px") %>%
  visHierarchicalLayout(direction = "UD", sortMethod = "directed") %>%
  visNodes(shape = "box", color = list(background = "lightgreen")) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)








