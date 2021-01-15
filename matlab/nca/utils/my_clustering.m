function [order,h] = my_clustering(X,distance_measure)

% [order,h] = my_clustering(X,distance_measure)
%
%'euclidean' 

if ~exist('distance_measure','var'), distance_measure = 'euclidean'; end

% see pdist.m

% cluster the rows of a matrix 
l=linkage(pdist(X,distance_measure),'ward');
[h,t,order]=dendrogram(l,0,'orientation','left');
