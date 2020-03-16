import numpy as np
import math

class Keq:
    def __init__(self,gibbs,temp):
        self.gibbs = gibbs
        self.temp = temp

    def keq(self):
        keq = np.zeros((np.shape(self.gibbs)))
        for i in range(np.shape(self.gibbs)[0]):
            for j in range(np.shape(self.gibbs)[1]):
                keq[i,j] = math.exp(-self.gibbs[i,j]/0.001987/self.temp)/55.14

        return keq