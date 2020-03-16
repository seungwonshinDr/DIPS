import numpy as np

def shiftmut(newPop, shiftmutProb):

    for i in range(np.shape(newPop)[0]):
        SeqAddPoint = np.random.choice(np.shape(newPop)[1] + 1)
        direction = np.random.choice(2)
        chance = np.random.rand()
        #direction = np.random.choice(4)

        if chance <= shiftmutProb and direction == 0 and SeqAddPoint != 0:
            newPop[i, :SeqAddPoint - 1] = newPop[i, 1:(SeqAddPoint)]
            newPop[i, SeqAddPoint - 1] = np.random.choice(4) + 1
        elif chance <= shiftmutProb and direction == 1 and SeqAddPoint != np.shape(newPop)[1]:
            newPop[i, SeqAddPoint + 1:] = newPop[i, SeqAddPoint:-1]
            newPop[i, SeqAddPoint] = np.random.choice(4) + 1

    return newPop