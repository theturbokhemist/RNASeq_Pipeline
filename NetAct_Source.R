plot_network_v = function(tf_links = tf_links){
  require(visNetwork)
  topology=data.frame(as.matrix(tf_links), stringsAsFactors = F)
  
  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list, label = node_list, font.size =30, value=c(rep(1,length(node_list))))
  
  
  #nodes <- data.frame(id = node_list, label = node_list, font.size =30,shape='circle',value=c(rep(1,length(node_list))))
  edge_col <- data.frame(c(1,2),c("blue","darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1,2),c("arrow","circle"))
  colnames(arrow_type) <- c("type", "color")
  edges <- data.frame(from =c(topology[,1]), to = c(topology[,2])
                      , arrows.to.type	=arrow_type$color[c(as.numeric(topology[,3]))]
                      , width = 3
                      , color = edge_col$color[c(as.numeric(topology[,3]))]
  )
  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  #  file  <- paste("network_",file,".html",sep="")
  #  visSave(network, file = file, selfcontained = F)
}

####EDIT TO PLOT_NETWORK_V function####
plot_network_v_html = function(tf_links = tf_links){
  require(visNetwork)
  topology=data.frame(as.matrix(tf_links), stringsAsFactors = F)
  
  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list, label = node_list, font.size =30, value=c(rep(1,length(node_list))))
  
  
  #nodes <- data.frame(id = node_list, label = node_list, font.size =30,shape='circle',value=c(rep(1,length(node_list))))
  edge_col <- data.frame(c(1,2),c("blue","darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1,2),c("arrow","circle"))
  colnames(arrow_type) <- c("type", "color")
  edges <- data.frame(from =c(topology[,1]), to = c(topology[,2])
                      , arrows.to.type	=arrow_type$color[c(as.numeric(topology[,3]))]
                      , width = 3
                      , color = edge_col$color[c(as.numeric(topology[,3]))]
  )
  network <- visNetwork(nodes, edges, height = "1000px", width = "100%") %>% 
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  network
  #  file  <- paste("network_",file,".html",sep="")
  file <- paste(Directory_Plots_NetAct_RACIPE, "/TFNetwork_Lucas_html", "_lfc", gsub(".", "", as.character(lfc), fixed = TRUE),
                "_qval", gsub(".", "", as.character(qval), fixed = TRUE), "_miTh", gsub(".", "", as.character(miTh), fixed = TRUE),filter_type, gsub(".*filter", "", Directory_CSVs_NetAct_RACIPE),
                if (outliers == TRUE) {print("_outliers-removed")}, ".html", sep = "")
  visSave(network, file = file, selfcontained = F)
}

