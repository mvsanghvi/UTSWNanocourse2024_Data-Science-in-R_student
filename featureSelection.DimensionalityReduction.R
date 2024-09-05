install.packages("factoextra")
install.packages("umap")
install.packages("Rtsne")
install.packages("fpc")
install.packages("ggbiplot")
library(ggplot2)
library(cowplot)
library(factoextra)
library(corrplot)
library(Rtsne)
library(umap)
library(cluster)
library(patchwork)
library(fpc)
library(ggbiplot)
iris_filepath <- "C:/Users/S438187/Desktop/Wu Lab/Other Projects/Comp Bio/UTSWNanocourse2024_Data Science in R_student/data/iris.csv"
iris1 <- read_csv(iris_filepath)
unique.iris1 <- iris1[!duplicated(iris1),]
features <- subset(unique.iris1, select = -c(Species))
labels <- unique.iris1$Species
head(as.data.frame(features))
scaled.features <- scale(features)
corrplot(cor(scaled.features), method = 'number')
data.pca <- princomp(scaled.features)
summary(data.pca)
#Run PCA
fviz_pca_var(data.pca, col.var = "black", repel = T)
fviz_cos2(data.pca, choice = "var", axes = 1:2)
ggbiplot(data.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = TRUE, ellipse.prob = 0.68)
#Run t-SNE
perplexity.value <- c(2,5,30,45)
plot.list <- list()
tsne.list <- list()
set.seed(0) #for reproducible results
for (i in 1:length(perplexity.value)){
  tsne.list[[i]] <- Rtsne(scaled.features, dims = 2, perplexity = perplexity.value[i])
  tmp.df <- data.frame(tsne_1 = tsne.list[[i]]$Y[,1], tsne_2 = tsne.list[[i]]$Y[,2])
  plot.list[[i]] <- ggplot(tmp.df, aes(tsne_1,tsne_2)) + geom_point() + 
    ggtitle(paste0("Perplexity: ", perplexity.value[i]))
}
plot_grid(plotlist = plot.list)
#Make the UMAP
min.dist <- c(0.01, 0.05, 0.1, 0.5)
n.neighbors <- c(5, 15, 30, 50)
plot.list <- list()
umap.list <- list()
set.seed(0) #for reproducible results
ct <- 1
for (i in 1:length(min.dist)){
  for (j in 1:length(n.neighbors)){
    umap.list[[ct]] <- umap(features, min_dist = min.dist[i], n_neighbors = n.neighbors[j])
    tmp.df <- data.frame(umap_1 = umap.list[[ct]]$layout[,1], umap_2 = umap.list[[ct]]$layout[,2])
    plot.list[[ct]] <- ggplot(tmp.df, aes(umap_1,umap_2)) + geom_point() + 
      ggtitle(paste0("Min_dist: ", min.dist[i], ", n_neighbors: ", n.neighbors[j]))
    ct <- ct + 1
  }
}
plot_grid(plotlist = plot.list, nrow = length(min.dist))

#for BCA
ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = unique.iris1$Species, ellipse = TRUE)

#for tsne and umap
tmp.df <- data.frame(umap_1 = umap.list[[ct]]$layout[,1], umap_2=umap.list[[ct]]$layout[,2],
                     Species= unique.iris1$Species)
plot.list[[ct]] <- ggplot(tmp.df)
#k-means clustering
scaled.features <- scale(features)
k <- seq(1,10) #list of values from 1 to 10
km.list <- list()
sil.list <- list() 
wcss <- c()
avg.score <- c() #average silhouette score 
set.seed(0)
for (i in 1:length(k)){
  print(i)
  km.list[[i]] <- kmeans(scaled.features, k[i], iter.max = 10, nstart = 1)
  wcss <- c(wcss, km.list[[i]]$tot.withinss)
  if (i > 1){
    sil.list[[i]] <- silhouette(km.list[[i]]$cluster, dist(scaled.features))
    avg.score <- c(avg.score, mean(sil.list[[i]][,3]))
  }
}
plot(k, wcss, type = 'b', main = paste('The Elbow Method'), xlab = 'Number of clusters', ylab = 'WCSS')
plot(seq(2,10), avg.score, type = 'b', main = paste('Silhouette scores'),xlab = 'Number of clusters',
     ylab = 'Average silhouette scores')
#Comparing Elbow and silhouette
p1 <- fviz_nbclust(scaled.features, kmeans, method = "wss") + geom_vline(xintercept = 7, linetype = 2)
p2 <- fviz_nbclust(scaled.features, kmeans, method = "silhouette")
p1 + p2
#Hierarchical clustering
scaled.features <- scale(features)
set.seed(0)
dist_mat <- dist(scaled.features, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

cut_avg <- cutree(hclust_avg, k = 2) # k – specifies the desired number of clusters
plot(hclust_avg)
rect.hclust(hclust_avg , k = 2, border = 2:6)


fviz_nbclust(scaled.features, hcut, method = "silhouette")
## testing with other linkages
set.seed(0)
hclus_single <- hclust(dist_mat, method = 'single')
hclus_complete <- hclust(dist_mat, method = 'complete')
hclus_centroid <- hclust(dist_mat, method = 'centroid') 
hclus_ward <- hclust(dist_mat, method = 'ward.D')
gen_silhouette <- function(h.clus, dist_mat){
  cut_avg.list <- list()
  sil.list <- list()
  avg.score <- c() #average silhouette score for each k
  for (i in 2:10){
    cut_avg.list[[i]] <- cutree(h.clus, k = i) # k - specifies the desired number of clusters
    sil.list[[i]] <- silhouette(cut_avg.list[[i]], dist_mat)
    avg.score <- c(avg.score, mean(sil.list[[i]][,3]))
  }
  return (avg.score)
}
sil_single <- gen_silhouette(hclus_single, dist_mat)
sil_complete <- gen_silhouette(hclus_complete, dist_mat)
sil_centroid <- gen_silhouette(hclus_centroid, dist_mat)
sil_ward <- gen_silhouette(hclus_ward, dist_mat)
sil_avg <- gen_silhouette(hclust_avg, dist_mat)

plot.data <- data.frame(k=seq(2,10), single=sil_single, complete=sil_complete, centroid=sil_centroid, ward=sil_ward, average=sil_avg)
data_long <- reshape::melt(plot.data, id = "k")
ggplot(data_long, aes(x = k, y = value, color = variable)) + geom_line() + xlab('Number of clusters') + 
  ylab('Average silhouette scores') + labs(color='Linkage') 
cc <- data.frame(labels = labels, k.means=km2.list[[7]]$cluster)
res <- cc %>% group_by(across(everything())) %>%
  mutate(Count = n()) %>%
  ungroup() %>%
  distinct()
#DBSCAN
data("multishapes")
df <- multishapes[, 1:2]
ggplot(df, aes(x,y)) + geom_point()
set.seed(123)
km.res <- kmeans(df, 5, nstart = 25)
fviz_cluster(km.res, df, frame = FALSE, geom = "point")


