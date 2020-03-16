import numpy as np

def cross(population, popSize):
    if np.mod(np.shape(population)[0], 2) == 0:
        newPopulation = np.zeros(np.shape(population))
        crossoverPoint = np.random.randint(0, np.shape(population)[1], popSize)
        for i in range(0, np.shape(population)[0], 2):
            # print(np.shape(population))
            newPopulation[i, :crossoverPoint[i]] = population[i, :crossoverPoint[i]]
            newPopulation[i + 1, :crossoverPoint[i]] = population[i + 1, :crossoverPoint[i]]
            newPopulation[i, crossoverPoint[i]:] = population[i + 1, crossoverPoint[i]:]
            newPopulation[i + 1, crossoverPoint[i]:] = population[i, crossoverPoint[i]:]
    else:
        newPopulation = np.zeros(np.shape(population))
        crossoverPoint = np.random.randint(0, np.shape(population)[1], popSize)
        for i in range(0, np.shape(population)[0] - 1, 2):
            # print(np.shape(population))
            newPopulation[i, :crossoverPoint[i]] = population[i, :crossoverPoint[i]]
            newPopulation[i + 1, :crossoverPoint[i]] = population[i + 1, :crossoverPoint[i]]
            newPopulation[i, crossoverPoint[i]:] = population[i + 1, crossoverPoint[i]:]
            newPopulation[i + 1, crossoverPoint[i]:] = population[i, crossoverPoint[i]:]
        newPopulation[-1, :] = population[-1, :]

    return newPopulation