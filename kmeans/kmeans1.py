#1. Initialize the k centroids randomly to k genes as starting points.
#2. Assign each gene to the cluster that has the closest centroid.
#3. After all genes have been assigned to clusters, recalculate the centroids for each cluster (as averages
#of all genes in the cluster).
#4. Repeat the gene assig



from scipy.spatial import distance
from itertools import izip
import sys


def outputFile(pts, distor):
    outfile = open(sys.argv[3], "w")
    outfile.seek(0)
    outfile.truncate()

    for index, pt in enumerate(pts):
        minDist = float('inf')  # Initial min_distance (infinite)
        bcluster = None
        for p, centroid in enumerate(centroids):
            dist = distance.sqeuclidean(pt, centroid)
            if dist < minDist:
                minDist = dist
                bcluster = p
        outfile.write(", ".join(gene[index]) + " " + str(bcluster) + "\n")
    outfile.close()


def checkEqual(centroids, nCentroids):
    if len(centroids) != len(nCentroids):
        return False
    for centroid, newCentroid in izip(centroids, nCentroids):
        if centroid[0] != newCentroid[0] or centroid[1] != newCentroid[1]:
            return False
    return True


import numpy as np
def kcluster (pts, k, centroids):
    clusters = [[] for _ in range(k)]

    for point in pts:
        minDist = float('inf')   #minimum distance

        bcluster = None
        for j, centroid in enumerate(centroids):
            dist = distance.sqeuclidean(point, centroid)        # sqeuclidean Computes the squared Euclidean distance
            if dist < minDist:
                minDist = dist
                bcluster = j
        clusters[bcluster].append(point)
    centroidsNew = [np.mean(cluster, axis=0) for cluster in clusters]

    if not checkEqual(centroids, centroidsNew):                #checks centroid equal or not and assign the gene
                                                               # to the cluster that has the closest centroid
        clusters = kcluster(pts, k, centroidsNew )
    else:
        distortionDist = 0.0
        for i in range(len(point)):
            closest_dist = sys.maxsize
            for j in range(len(centroids)):
                dist = distance.sqeuclidean(point[i], centroids[j])
                if dist < closest_dist:
                    closest_dist = dist
            distortionDist += closest_dist * closest_dist
            distor = distortionDist / len(point)
        print "Squared Error Distortion--> ", distortionDist/len(point)  #print the Squared Error Distortion
        outputFile(pts, distor )                                                #printing the output to a file
    return clusters


def readData(inp, gene):
    f = map(str.rstrip, file(sys.argv[1], 'r').readlines())
    for line in f:
        columns = line.split('\t')  # make the lines to columns
        for n, i in enumerate(columns):
            if i == '':
                # empty value to 0
                columns[n] = '0'
        inp.append(columns)
    for x in inp:
        gene.append([x[0]])
        del x[0]
        for n, i in enumerate(x):
            x[n] = float(x[n])
    k = int(sys.argv[2])
    return k


import random
if __name__ == "__main__":
    inp = []
    gene = []
    k =readData(inp, gene)
    centroids = random.sample(inp, k)           # Randomly picking k points as centroids
    clusters = kcluster(inp, k, centroids)
