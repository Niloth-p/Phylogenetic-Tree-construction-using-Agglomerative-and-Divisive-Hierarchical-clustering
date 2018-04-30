'''
The prog is written to construct the phylogenetic tree (dendrogram)
based on DNA/Protein sequences of species given in a dataset using
Agglomerative and Divisive Hierarchical Clustering and to compare the 2 methods
'''
# pylint: disable=invalid-name
import itertools
import collections
from math import inf
import pickle
import numpy as np
#import plotly.plotly as py
#import plotly.figure_factory as ff
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt
#import fastcluster
#import pandas


filename = "human_gene.data"
first_time = False


def ReadData():
    '''
    Reads the gene sequences from the file
    and stores as a dictionary
    '''
    f = open(filename, "r") 
    data = {}
    label = ''
    genes = ''
    count = 0
    if f.mode == 'r':
        fl = f.readlines()
        for x in fl:
            #print(x)
            if x[0] == '>':
                count = count + 1
                if label != '':
                    data[label] = genes
                    genes = ''
                label = x[4:].rstrip()
            else:
                genes = genes + x.rstrip()
        data[label] = genes
    print("Read the data")
    return data


def NeedlemanWunsch(seqA, seqB):
    '''
    Generating the initial proximity matrix
    '''
    #Scoring Scheme
    MatchScore = 1
    MismatchScore = -1
    GapPenalty = -1

    sizeA = len(seqA)
    sizeB = len(seqB)
    a = np.empty((sizeA, sizeB))
    i = 0
    for j in range(sizeB):
        a[0][j] = j * MismatchScore
    for i in range(sizeA):
        a[i][0] = i * MismatchScore

    for i in range(1, sizeA):
        for j in range(1, sizeB):
            if seqA[i] == seqB[j]:
                score = MatchScore
            else:
                score = MismatchScore

    choice1 = a[i-1][j-1] + score     # If characters are aligned    #Move diagonal - right & down
    choice2 = a[i-1][j] + GapPenalty       # Gap in seqB             #Move right
    choice3 = a[i][j-1] + GapPenalty       # Gap in seqA             #Move down
    a[i][j] = max(choice1, choice2, choice3)
    maxi = -1000
    for k in range(sizeA):
        for l in range(sizeB):
            if a[k][l] > maxi:
                maxi = a[k][l]

    return maxi

    # cost = 1

    # prev = [0 for i in range(len(seqB) + 1)]
    # curr = [0 for i in range(len(seqB) + 1)]

    # for i in range(0, len(seqA) + 1):
    #     for j in range(0, len(seqB) + 1):
    #         if i == 0 and j == 0:
    #             curr[0] = 0
    #         elif i == 0:
    #             curr[j] = curr[j-1] + cost
    #         elif j == 0:
    #             curr[j] = prev[j] + cost
    #         else:
    #             if seqA[i-1] == seqB[j-1]:
    #                 curr[j] = prev[j-1]
    #             else:
    #                 curr[j] = min(prev[j] + cost, curr[j-1] + cost, 1 + prev[j-1])
    #     prev = curr.copy()

    # return prev[len(seqB)]


def find(cluster_A, cluster_B, Z):
    '''
    Uses MIN linkage to find the distance b/w two clusters
    i.e the min distance between any 2 points,
    1 point from each cluster
    '''
    temp = []
    for i in range(len(cluster_A)):
        for j in range(len(cluster_B)):
            temp.append(Z[cluster_A[i]][cluster_B[j]])
    return min(temp)


def chooseMinValue(clusters, Z):
    '''
    Choosing the min value in Z
    and get the indices to get the
    points or clusters about to be merged
    '''
    minProx = 0
    indexOfminZvali = -1
    indexOfminZvalj = -1

    for i in range(len(clusters)):
        for j in range(len(clusters)):
            Z[i][j] = find(clusters[i], clusters[j], Z)
            if i != j and minProx > Z[i][j]:
                minProx = Z[i][j]
                indexOfminZvali = i
                indexOfminZvalj = j
    return Z, indexOfminZvali, indexOfminZvalj, minProx


def enterIntoX(temp, clusters, indexOfminZvali, indexOfminZvalj, minProx, X):
    '''
    For drawing dendrogram,
    storing info of clusters in X,
    in the following format
    Append to matrix Z
        #1 - Cluster A
        #2 - CLuster B
        #3 - Proximity value
        #4 - numOfpts in new cluster formed i.e Total pts in clusters A & B
    '''
    #Part 1
    X_temp = []
    for i in range(len(temp)):
        if temp[i] == clusters[indexOfminZvali]:
            X_temp.append(i)
            break
    #Part 2
    for i in range(len(temp)):
        if temp[i] == clusters[indexOfminZvalj]:
            X_temp.append(i)
            break
    #Part 3
    X_temp.append(minProx)
    #Part 4
    X_temp.append(len(clusters[indexOfminZvali]) + len(clusters[indexOfminZvalj]))

    X.append(X_temp)
    return X


def startClustering(keys, Z):
    '''
    Perform agglomerative clustering
    1. Choose the min value and hence the clusters to be merged,
        and update the proximity matrix
    2. Enter the values into X for plotting dendrogram
    3. Update cluster list, by removing 1 of the clusters,
        and merging the 2 clusters and naming it as 1 of them
    '''
    #numOfClusters = N
    clusters = [[i] for i in range(len(keys))]
    temp = [[i] for i in range(len(keys))]
    X = []

    while len(clusters) > 1:
        print(len(clusters))
        Z, indexOfminZvali, indexOfminZvalj, minProx = chooseMinValue(clusters, Z)
        X = enterIntoX(temp, clusters, indexOfminZvali, indexOfminZvalj, minProx, X)
        temp.append(clusters[indexOfminZvali] + clusters[indexOfminZvalj])
        clusters[indexOfminZvali] += clusters[indexOfminZvalj]
        clusters.remove(clusters[indexOfminZvalj])
    return X


def main():
    '''
    main - reads data, writes to pickle or loads from pickle
    Calls divisive and agglomerative clustering funcs
    '''
    data = {}
    keys = []
    data = ReadData()
    N = len(data)
    Z = np.empty((N, N))
    odata = collections.OrderedDict(sorted(data.items()))
    for key in odata.keys():
        keys.append(key)
    print("Coverted to ordered dict, created list keys")
    if first_time is True:
        for a, b in itertools.combinations(keys, 2):
            count1 = 0
            count2 = 0
            i = j = 0
            for k in odata.keys():
                if k == a:
                    i = count1
                    count2 = count2 + 1
                elif k == b:
                    j = count2
                    count1 = count1 + 1
                else:
                    count1 = count1 + 1
                    count2 = count2 + 1
            print("seqA:" + str(i) + " seqB:" + str(j))
            Z[i][j] = NeedlemanWunsch(data[a], data[b])

        #z = np.matrix(Z)
        print(Z)
        f = open('matrix2', 'wb')
        pickle.dump(Z, f)
        f.close()

    else:
        f = open('matrix2', 'rb')
        Z = pickle.load(f)
        f.close()

    # for i in range(len(data)):
    #     for j in range(len(data)):
    #         print(Z[i][j])
    #     print()
    for i in range(len(data)):
        for j in range(i, len(data)):
            Z[j][i] = Z[i][j]
    DivisiveClustering(keys, Z)
    X = startClustering(keys, Z)
    drawDendrogram(X)
    #print(DataFrame(Z))
    #Triangular matrix Z


def drawDendrogram(X):
    """
    Method to draw and save the dendrogram using linkage matrix.
    """
    #print("Here")
    X = np.array(X).astype(float)
    X = np.clip(X, 1, 1000)
    dendrogram(X, color_threshold=1, orientation='right')
    #clusters = fastcluster.fcluster(X)
    plt.show()
    #plt.savefig('plots/agglomerative.jpg')


def DivisiveClustering(keys, Z):
    '''
    Performs Divisive Clustering
    '''
    #print(Z)
    clusters = [[i for i in range(len(keys))]]
    X = []
    temp = []
    while len(clusters) != len(keys):
        index = chooseCluster(clusters, Z)
        origCluster = clusters[index]
        newClusterIndex = createNewCluster(origCluster, Z)
        newCluster = [origCluster[newClusterIndex]]
        origCluster.remove(origCluster[newClusterIndex])

        while True:
            newClusterIndex = splitPoints(newCluster, origCluster, Z)
            if newClusterIndex is None:
                break
            newCluster.append(origCluster[newClusterIndex])
            origCluster.remove(origCluster[newClusterIndex])

        # distance_bw_clusters = get_distance_bw_clusters(origCluster, newCluster)
        # length_new = len(origCluster) + len(newCluster)

        clusters[index] = origCluster.copy()
        clusters.append(newCluster.copy())

        origCluster.sort()
        newCluster.sort()
        temp.append([origCluster.copy(), newCluster.copy()])
        print("length : ", len(clusters))
        #index = get_cluster_to_split(clusters, Z)
    drawDendrogram2(temp, X, Z)


def chooseCluster(clusters, Z):
    '''
    Finds the index of the cluster to split
    i.e the cluster with the max dissimilarity
    between its points
    '''
    index = 0
    max_value = -inf
    for idx, cluster in enumerate(clusters):
        value = get_max_dissimilarity(cluster, Z)
        if len(clusters) > 1 and value > max_value:
            max_value = value
            index = idx
    return index


def get_max_dissimilarity(cluster, Z):
    """
    Computes maximum dissimilarity
    between any 2 points of a cluster
    """
    max_distance = 0
    for i in cluster:
        for j in cluster:
            #print(Z[i][j])
            max_distance = max(max_distance, Z[i][j])
    #print(max_distance)
    return max_distance


def createNewCluster(cluster, Z):
    """
    Gets the index of the element in the cluster which
    will be split from the cluster,
    to form a new cluster
    """
    maxSumDist = 0
    newClusterIndex = 0
    index = 0

    for i in cluster:
        sumDist = 0
        for j in cluster:
            sumDist += Z[i][j]
        if sumDist > maxSumDist:
            maxSumDist = sumDist
            newClusterIndex = index
        index += 1

    return newClusterIndex


def splitPoints(clusterA, clusterB, Z):
    """
    Splits points in the original cluster
    between the old and the new clusters
    clusterA - new
    clusterB - old
    """
    if len(clusterB) == 1:
        return None

    index = None
    maxDiff, idx = 0, 0

    for i in clusterB:
        sumdistToA, sumdistToB = 0, 0
        for _ in clusterA:
            sumdistToA += Z[i][_]
        for _ in clusterB:
            sumdistToB += Z[i][_]

        avdistToA = sumdistToA/len(clusterA)
        avdistToB = sumdistToB/(len(clusterB)-1)
        diff = avdistToB - avdistToA

        if maxDiff < diff:
            maxDiff = diff
            index = idx
        idx += 1
    return index


def distBetwClusters(clusterA, clusterB, Z):
    """
    Computes the distance between 2 clusters
    min linkage
    """
    dist = 1500
    for i in clusterA:
        for j in clusterB:

            if Z[i][j] != 0:
                #print("Z[i][j]:" + str(Z[i][j]))
                dist = min(dist, Z[i][j])
    return dist


def compute_linkage(temp, X, Z):
    """
    Computes the linkage matrix
    for plotting dendrogram
    """
    temp = list(reversed(temp))
    idx = 0

    for i in range(0, len(temp)):
        Xtemp = []
        Xtemp.append(idx)
        Xtemp.append(idx+1)
        Xtemp.append(distBetwClusters(temp[i][0].copy(), temp[i][1], Z))
        Xtemp.append(len(temp[i][0]) + len(temp[i][1]))

        X.append(Xtemp)
        idx += 2
    return X


def drawDendrogram2(temp, X, Z):
    """
    Plotting the dendrogram
    """
    X2 = compute_linkage(temp, X, Z)
    X2 = np.array(X).astype(float)
    #print(X2)
    dendrogram(X2, color_threshold=1, orientation='right')
    plt.savefig('divisiveClustering.png')
    plt.show()

if __name__ == "__main__":
    main()
