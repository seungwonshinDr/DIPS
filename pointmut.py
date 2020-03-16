import numpy as np

def pointmut(population, mutationProb):
    whereMutate = np.random.rand(np.shape(population)[0], np.shape(population)[1])
    for i in range(np.shape(whereMutate)[0]):
        for j in range(np.shape(whereMutate)[1]):
            if whereMutate[i, j] >= mutationProb:
                whereMutate[i, j] = 1
            else:
                whereMutate[i, j] = 0

    newPop = population * whereMutate
    #print(newPop)

    # population[np.where(whereMutate < self.mutationProb)] = 0

    for i in range(np.shape(newPop)[0]):
        for j in range(np.shape(newPop)[1]):
            if newPop[i, j] == 0:
                seqcode = np.random.choice(4) + 1
                newPop[i, j] = seqcode

    return newPop