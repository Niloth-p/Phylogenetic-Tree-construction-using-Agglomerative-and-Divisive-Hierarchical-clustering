N IMAYAVALLI
2015A8PS0444H
Assignment 2
Phylogenetic Tree using Hierarchical clustering

Description:
The prog is written to construct the phylogenetic tree (dendrogram) based on DNA/Protein sequences of species 
given in a dataset using Agglomerative and Divisive Hierarchical Clustering and to compare Agglomerative and Divisive
methods. I have used min linkage to calculate the proximities.

Dataset used:
Human Gene DNA Sequences
http://genome.crg.es/datasets/ggalhsapgenes2005/hg16.311.putative.cds.fa

How to run:
python Clusterer.py

For using some other dataset, change the global variable 'filename'
with the name of your dataset.
And change the global variable 'first_time' to create and store the pickle file.

PseudoCode/Explanation of the algo:
Read data, preprocessing (also store the data in an ordered dictionary,
    so that a list of exclusive keys can be obtained - to get
    all pairs of pts for calculating proximity matrix)
Needleman Wunsch algo - for distance calculation
Calc Z - proximity matrix, dump pickle

Clustering algo - Agglomerative
    Set numOfClusters = N
    Set clusters = list of keys
    Loop function call to clustering method until numOfClusters = 1
        Choose min value in Z -> get i and j indices
        Append to matrix Z
            #1 - Cluster A
            #2 - CLuster B
            #3 - Proximity value
            #4 - numOfpts in new cluster formed i.e Total pts in clusters A & B
            Append 1 to 4 to a list, and append that list to Z
        Update numOfClusters -> --1
        Update Clusters -> remove A or B
        Update Proximity Matrix by Ward's method using Lance Williams Formula
Draw dendrogram

Divisive clustering algo
    Set numOfClusters = 1
    Set clusters = all points as 1 cluster
    Loop function call to clustering method until numOfClusters = N
        Choose the cluster to split
            Get the max dissimilarity for each cluster i.e max distance between 2 points in cluster
            Find the cluster with max dissimilarity - this is the cluster to split first
        Find ele to form the new cluster i.e find the elemenet in the cluster
            whose sum of distances to all other points in the cluster is maximum
        Remove the ele to new cluster A from old cluster B
        Split eles in old original cluster B between the 2 clusters A and B 
            Find ele with max(avdisttoB - avdisttoA)
                max - the max value 
                avdist - average distance (sum of distances / num of points)
                If avdisttoB < avdisttoA, skip to next step
            Add the ele to cluster A and remove from cluster B
            Keep finding eles 1 by 1 and adding them to cluster A, and removing from cluster B
                Do it 1 by 1, bcoz the value of avdist will vary with addition/removal of each point
        Update clusters by updating the old cluster B, and adding new cluster A
        Update numOfClusters += 1
        Add [clustA, clustB] to a temp list which will be needed for computing linkage matrix
    Get the linkage matrix X such that
        #1 - num
        #2 - num + 1
        #3 - Proximity value
        #4 - [len(clustA), len(clustB)]
        Get #4 from the temp list after reversing its order
    Draw dendrogram

