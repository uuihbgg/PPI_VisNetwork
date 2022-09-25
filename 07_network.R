rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options(BioC_mirror="https://anaconda.org/bioconda/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(download.file.method = 'libcurl')
# options(url.method='libcurl')

getwd()
setwd('C:/Users/LWH/Downloads/bioinfo/IPF&CHD&SSc')

source("C:/Users/LWH/Downloads/bioinfo/AD/AD_functions.R")
source("C:/Users/LWH/Downloads/bioinfo/AD/violin.R")
library(dplyr)
library(plyr)
library(tidyr)
library(showtext)
library(visNetwork)
library(ggplot2)
library(STRINGdb)
library(igraph)
library(ggraph)
library(graphlayouts)
library(networkD3)
library(clustAnalytics)
# devtools::install_github("martirm/clustAnalytics")
# load(file = "genes_mapped2.Rdata")
# font_all <- font_files()
# windowsFonts()
font_add("Calibri", "calibri.ttf",
         bold = NULL, italic = NULL,
         bolditalic = NULL, symbol = NULL)
font_add("constan", regular = "constan.ttf",
         italic = "constani.ttf")
showtext_auto()
vars1 = c("input_CHP_down","input_CHP_up" ,"input_CTD_down","input_CTD_up" , 
          "input_IPF_down","input_IPF_up")
vars2 = paste0(vars1, "_node")

filenames2 <- list.files("04_PPI/", pattern = "*node.csv")
for (i in 1:length(filenames2)) {
  node_table <- read.csv(paste0("04_PPI/",filenames2[i]),header = TRUE) %>% 
    .[.$selected == "true",]
  assign(vars2[i], node_table)
}

filenames1 <- list.files("04_PPI/", pattern = "^input*") %>% 
  .[-grep("node",.)]
for (i in 1:length(filenames1)) {
  edge_table <- read.csv(paste0("04_PPI/",filenames1[i]),header = TRUE) %>% 
    .[.$source %in% 
        (vars2[i] %>% get() %>% .$shared.name)|.$target %in% 
        (vars2[i] %>% get() %>% .$shared.name),]
  assign(vars1[i], edge_table)
}
# (vars2[1] %>% get() %>% .$shared.name)

vars3 <- vars1 %>% strsplit("_") %>% lapply(function(x) paste0(x[2],"_",x[3])) %>%
  unlist() %>% paste("edges",.,sep = "_")
for (i in 1:length(vars3)) {
  assign(vars3[i],vars1[i] %>% get() %>% .[,c(1,2,3)])
  # colnames(vars1[i] %>% get()) <- c("from", "to", "weight")
}

get_network = function(edges) {
  tmp <- graph.data.frame(edges, directed = T)
  tmp1 <- graph.data.frame(edges, directed = F)
  igraph::V(tmp)$deg <- igraph::degree(tmp) # 每个节点连接的节点数
  igraph::V(tmp)$size <- igraph::degree(tmp)/5 #
  igraph::E(tmp)$width <- igraph::E(tmp)$combined_score
  igraph::V(tmp)$grp <- membership(cluster_leiden(tmp1, 
                                                  resolution_parameter=0.01,
                                                  weights=E(tmp)$combined_score ))
  return(tmp)
}
# V(vars4[1] %>% get())$grp
vars4 <- paste0("net","_",vars1 %>% 
                  lapply(function(x) strsplit(x,"_") %>% unlist() %>% 
                           .[2:3] %>% paste0(collapse = "_")) %>% unlist())
for (i in 1:length(vars4)) {
  assign(vars4[i], get_network(vars3[i] %>% get()))
}

get_ggraph <- function(net,UP){
  deg_cut <- V(net)$deg %>% .[order(.,decreasing = T)] %>% .[20]
  if(UP == TRUE){
    plot_ggraph = ggraph(net,layout = 'centrality',cent=deg)+
      geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
      geom_node_point(aes(size=size,color=deg), alpha=0.6)+
      geom_node_text(aes(filter=deg>deg_cut, size = deg/2,label=name),  
                     repel = T,family="serif")+
      scale_edge_width(range = c(0.2,1))+
      scale_size_continuous(range = c(1,10) )+
      guides(size="none")+ 
      scale_color_continuous(name="Degree",low="red",high ="blue",
                             guide = guide_legend(title.theme = element_text(
                               size = 15,face = "bold",colour = "red"
                             ), label.theme = element_text(size = 13,face = "bold"),
                             direction = "horizontal",title.position = "top",
                             title.vjust = 0.5))+ 
      theme_graph()+
      theme(legend.position=c(0.35,0))}
  else{
    plot_ggraph = ggraph(net,layout = 'centrality',cent=deg)+
      geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
      geom_node_point(aes(size=size,color=deg), alpha=0.6)+
      geom_node_text(aes(filter=deg>deg_cut, size = deg/2,label=name),  
                     repel = T,family="serif") + 
      scale_edge_width(range = c(0.2,1))+
      scale_size_continuous(range = c(1,10) )+
      guides(size="none")+ 
      scale_color_continuous(name="Degree",low="lightgreen",high ="darkblue",
                             guide = guide_legend(title.theme = element_text(
                               size = 15,face = "bold",colour = "lightgreen"
                             ), label.theme = element_text(size = 13,face = "bold"),
                             direction = "horizontal",title.position = "top",
                             title.vjust = 0.5))+ 
      theme_graph()+
      theme(legend.position=c(0.35,0))}
  return(plot_ggraph)
}
vars5 <- paste0("graph","_",vars1 %>% 
                  lapply(function(x) strsplit(x,"_") %>% unlist() %>% 
                           .[2:3] %>% paste0(collapse = "_")) %>% unlist())
list_UP <- rep(c(F,T),3)
for (i in 1:length(vars5)) {
  assign(vars5[i],get_ggraph(vars4[i] %>% get(), list_UP[i]))
}
library(patchwork)
load(file = "vars_go.Rdata")
# (vars5[2] %>% 
#     get())/(vars5[1] %>% get())|(goplot_up_CHP/goplot_down_CHP)+
#   plot_layout(widths = c(5,2))
design <- "
  13
  ##
  24
"
vars5[2] %>% get() + vars5[1] %>% get() + 
  goplot_up_CHP + goplot_down_CHP +
  plot_layout(design=design,widths=c(3,2),heights=c(6,-1.1,6))

(vars5[4] %>% get())/(vars5[3] %>% get()) + 
  goplot_up_CTD + goplot_down_CTD +
  plot_layout(design=design,widths=c(3,2),heights=c(6,-1.1,6))

(vars5[6] %>% get())/(vars5[5] %>% get()) + 
  goplot_up_IPF + goplot_down_IPF +
  plot_layout(design=design,widths=c(3,2),heights=c(6,-1.1,6))

# colnames(edges_IPF_down)[3] = "weight"
# tmp <- graph.data.frame(edges_IPF_up, directed = F)
# membership(cluster_leiden(tmp, resolution_parameter=0.01,
#                           weights=E(tmp)$combined_score ))
# plot(tmp, edge.width=E(tmp)$combined_score, 
#      edge.arrow.size=.2,edge.curved=0,
#      vertex.color="orange", vertex.frame.color="#555555",
#      vertex.label.color="black",
#      vertex.label.cex=.7, layout = layout_in_circle)
# degree(tmp,normalized = T)
# vcount(tmp)/ecount(tmp) #返回图g的定点数、边数
# is.connected(tmp) #图g是否连通
# clusters(tmp) #图g有多少分支
# edge_density(tmp)
# transitivity(tmp, type="global")
# transitivity(tmp, type="local")
# cluster_fast_greedy(tmp, merges = TRUE, modularity = TRUE,
#                     membership = TRUE, weights = E(tmp)$combined_score)


# plot(tmp, vertex.shape="none",edge.arrow.size=.1,edge.width=E(tmp)$combined_score,
#      vertex.label.font=2, vertex.label.color="gray40",
#      vertex.label.cex=.7, edge.color="gray85")
pattern <- "name|IsSingleNode|selected|SelfLoops|PartnerOfMultiEdgedNodePairs|NumberOfDirectedEdges"
cytohub <- vars2[2] %>% get() %>% 
  .[,-grep(pattern, colnames(.))]
cytohub = cytohub[,-13] %>%  #MRPL24 RPS15 RPL3 RPL12 RPS3 RPL36 CASP3 GNL2 PSMB6 FBL RPL18
  mutate(., score = apply(., 1, function(x) prod(x+1))) %>% 
  .[order(.$score, decreasing = T),] %>% rownames() %>% .[1:20]
p = lapply(colnames(cytohub)[-c(1,13)], 
           function(x) rownames(cytohub[order(cytohub[,x], decreasing = T),])[1:10])

############################################################################

rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

setwd('C:/Users/LWH/Downloads/bioinfo/IPF&CHD&SSc')
source("C:/Users/LWH/Downloads/bioinfo/AD/AD_functions.R")
load(file = "Drug.Rdata")
library(tidyverse)
library(dplyr)
library(visNetwork)
library(ggplot2)
library(igraph)

rownames(proteins) <- proteins$parent_key
drug$Symbol <-  proteins[drug$id,]$Symbol
edges_drug <- drug[,c(6,7)] %>% as.data.frame() %>% mutate_all(unlist)
# edges_drug$name = unlist(lapply(c(1:316), 
#                             function(i) ifelse(grepl("[",edges_drug$name[i],fixed = T), 
#                                                edges_drug$name[i] %>% 
#                                                  str_extract("\\[.+\\]") %>% 
#                                                  str_sub(2,-2), 
#                                                edges_drug$name[i])))
tmp <- graph.data.frame(edges_drug, directed = T)
# E(tmp)$weight = count.multiple(tmp)
# tmp<-igraph::simplify(tmp,remove.multiple=TRUE,remove.loops=TRUE,
#               edge.attr.comb="mean")
# igraph::V(tmp)$size <- igraph::degree(tmp)/5
# plot(tmp,edge.width=E(tmp)$weight %>% sqrt(),main="有向图 g.undir",
#      edge.label=E(tmp)$weight,edge.arrow.size=.05)
links <- get.edgelist(tmp) %>% as.data.frame() %>% .[-c(312,315),]
links[,3] <- count.multiple(tmp) %>% .[-c(312,315)]
colnames(links) = c("from","to","weight")
GENES <- c("SLC6A4","CHRM1","CP","SERPIND1","MMP13","MMP10",
           "MMP11","MMP7","MMP1") #"COL1A1","COL3A1",
DRUGS <- drug1 %>% .[.$Drug_score>22|.$Gene_number>1,] %>% rownames()
links_2 <- links %>% mutate(from_c = count(., "from")[match(.$from, count(., "from")$from),2]) %>%
  mutate(to_c = count(., "to")[match(.$to, count(., "to")$to),2]) %>%
  filter(!(from_c == 1 & to_c == 1)) %>%
  dplyr::select(1,2,3) %>% .[.$to %in% GENES & .$from %in% DRUGS,]
tmp_id <- links_2 %>% { data.frame(id = c(.$from, .$to)) } %>% distinct() %>% 
  .$id
links_2$from <- drug1[links_2$from,]$Name
nodes_2 <- links_2 %>% { data.frame(id = c(.$from, .$to)) } %>% distinct()
nodes_2$group <- lapply(nodes_2$id,function(x) ifelse(isTRUE(x %in% DF[DF$logFC>0,]$SYMBOL),
                                                      "UP—Regulated",
                                                      ifelse(isTRUE(x %in% DF[DF$logFC<0,]$SYMBOL),
                                                             "DOWN—Regulated",
                                                             "Drug"))) %>% 
  unlist()
nodes_2$value <- igraph::degree(tmp) %>% as.data.frame() %>% .[tmp_id,] %>%
  .^(1/2)
simple_id1 <- nodes_2[nodes_2$id %>% grep("-",.),]$id %>% lapply(function(x) strsplit(x,"-") %>% 
                                                     unlist() %>% tail(1)) %>% 
  unlist() %>% str_to_title()
nodes_2[nodes_2$id %>% grep("-",.),]$id <- simple_id1
simple_id2 <- links_2[links_2$from %>% grep("-",.),]$from %>% 
  lapply(function(x) strsplit(x,"-") %>% unlist() %>% tail(1)) %>% 
  unlist() %>% str_to_title()
links_2[links_2$from %>% grep("-",.),]$from <- simple_id2

# nodes_2$shape <- lapply(nodes_2$id,function(x) ifelse(grepl("^DB",x),"dot","circle"))
library(qgraph) 
g<-graph_from_data_frame(links_2)
coords <- qgraph.layout.fruchtermanreingold(
  as_edgelist(g, names = F), 
  weights=E(g)$weight, 
  vcount=vcount(g), 
  area=vcount(g)^2/2, 
  repulse.rad=vcount(g)^3/250
) 
addNodes <- data.frame(label = c("Drug","UP","DOWN"), 
                       shape = "icon",
                       icon.code = c("f111", "f139","f13a"), 
                       icon.size = c(10,20,20), 
                       icon.color = c("skyblue","#B84D64","#B6CCD7")
                       # font = list("face"= "Times",size = 10)
                       )
visNetwork(nodes_2, links_2,height = "2000px") %>% visInteraction(hover = T) %>% 
  visOptions(highlightNearest = list(enabled = T, degree = 2, hover = T), 
             nodesIdSelection = TRUE)%>% 
  visEdges(shadow = TRUE,
           arrows =list(to = list(enabled = TRUE, scaleFactor = 0.6)),
           color = list(color = "lightblue", highlight = "red")) %>%
  visNodes(font = list("face"= "Times",size = 40),physics=T) %>%
  visGroups(groupname = "Drug",shape = "icon", 
            icon = list(code = "f111", color = "skyblue",size = 25),
            font = list(size = 30,face="bold"))%>%
  visGroups(groupname = "UP—Regulated", shape = "icon", 
            icon = list(code = "f139", color = "#B84D64",size = 75)) %>% 
  visGroups(groupname = "DOWN—Regulated", shape = "icon", 
            icon = list(code = "f13a", color = "#B6CCD7",size = 75)) %>%
  addFontAwesome(name="font-awesome")%>%
  visIgraphLayout("layout.norm", layoutMatrix = coords)  %>%
  visPhysics(enabled = FALSE) %>% 
  visLegend(position = "left",useGroups = FALSE,addNodes = addNodes,
            stepY=50)
                                            
  # visPhysics(forceAtlas2Based = list(gravitationalConstant = -50, 
  #                                    springConstant = 0.2,
  #                                    "centralGravity" = 1),
  #            repulsion=list(nodeDistance=200,springConstant = 0.2))


  
#   nodes <- data.frame(id = 1:3, 
#                       shape = "icon", 
#                       icon.face = "FontAwesome",
#                       icon.color = c("#800000", "#0000ff", "#ffa500"), # doesn't have any effect on icon color
#                       icon.code = c("f102", "f0c3", "f0f0"))
#   edges <- data.frame(from = c(1,2), to = c(2,3))
# #f139|f102:up;f13a|f103:down
#   visNetwork(nodes, edges) %>%
#     addFontAwesome()
