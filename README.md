# ClusteringVariableSelectionUsingR

These are R codes for implementing variable (feature) selection in model-based clustering 

Function "gmmVarSel" for selecting clustering relevant variables. Function "gmmEM" for identifying optimal clusters using these clustering relevant variables

#examples
## Variable selection
s <- gmmVarSel(iris[,1:4])

## Identify optimal clusters using these relevant variables (s)
m <- gmmEM(iris[,s], c=3, initialize = "kmeans")

summary(m)

table(m$class, iris[,5])


