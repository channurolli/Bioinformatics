from random import randint

def Profile(DNA, s, k, dist):
    """
     Given a motif, scan each DNA sequence to find best match
     :param DNA:
     :param motif:
     :return:
     """
    totalDist = 0
    bestAlignment = []
    pb = 1.0
    for i, base in enumerate(DNA[s:s + k]):
        pb *= dist[i][base]
    # for seq in DNA:
    #     minHammingDist = k + 1
    #     for s in xrange(len(seq) - k + 1):
    #         HammingDist = sum([1 for i in xrange(k) if DNA[i] != seq[s + i]])
    #         if (HammingDist < minHammingDist):
    #             bestS = s
    #             minHammingDist = HammingDist
    #             bestAlignment.append(bestS)
    #     totalDist += minHammingDist
    return pb

def Consensus(s, DNA, k):
    """ compute the consensus k-Motif of an alignment given offsets into each DNA string.
            s = list of starting indices, DNA = list of nucleotide strings,
            k = Target Motif length """
    consensus = ''
    profile = dict(zip("acgt", ([], [], [], [])))
    for i in xrange(k):
        # loop over string positions
        cnt = dict(zip("acgt", (0, 0, 0, 0)))
        for j, sval in enumerate(s):
            # loop over DNA strands
            base = DNA[j][sval+i]
            cnt[base] += 1
        for base, c in cnt.iteritems():
            profile[base].append(round(float(c) / len(s), 2))

        consensus += max(cnt.iteritems(), key=lambda tup: tup[1])[0]

    for j, sval in enumerate(s):
        print "DNA[{}][{}:{}]\t{}\n".format(j, sval, sval+k, DNA[j][sval:sval+k])

    for base, p in profile.iteritems():
        print "{}\t{}\n".format(base, ",".join(map(str, p)))

    print "Consensus\t{}".format(consensus)

    return consensus


def GibbsProfileMotifSearch(seq, k):
    bestAlign = 0.0
    m = [randint(0,len(seq[x])-k+1) for x in xrange(len(seq))]
    notBest = 0
    while True:
        rm = randint(0,len(seq)-1)
        m[rm] = -1
        dist = [dict([(base, 0.1) for base in "acgt"]) for i in xrange(k)]
        for t in xrange(len(seq)):
            if (m[t] < 0):
                continue
            for i, base in enumerate(seq[t][m[t]:m[t] + k]):
                dist[i][base] += 1.0
        for i in xrange(k):
            total = sum(dist[i].values())
            for base in "acgt":
                dist[i][base] /= total
        mdist = dist
        score = 0.0
        for t in xrange(len(seq)):
            if (m[t] < 0):
                Scorem = 0.0
                for i in xrange(len(seq[rm])-k+1):
                    score = Profile(seq[rm], i, k, mdist)
                    if (score > Scorem):
                        StartM, Scorem = i, score
                score += Scorem
                m[t] = StartM
            else:
                score += Profile(seq[t], m[t], k, mdist)
        if (score > bestAlign):
            bestAlign = score
            notBest = 0
        else:
            notBest += 1
            if (notBest > len(seq)):
                break
    return score, m


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


with open('DNASamples.fa') as fp:
    seqApprox = []
    for name, seq in read_fasta(fp):
        seqApprox.append(seq)


def display(k):
    N = 100
    temp = 0.0
    for i in xrange(N):
        # print "Motif length -->  " ,ki
        scoreS, start = GibbsProfileMotifSearch(seqApprox, k)
        if scoreS > temp:
            motif = [s for s in start]
            temp = scoreS
    Consensus(motif, seqApprox, k)

if __name__ == "__main__":
    """
    main program
    """
    # k = [8,9,10,11]       # motif length
    for k in range(8,12):
        if k == k:
            print "Motif Length --> ", k
            display(k)