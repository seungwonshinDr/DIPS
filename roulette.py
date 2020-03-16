import numpy as np

def roulette(elites, population, fitness):
    fitness_reduced = fitness - fitness.min()
    for i in range(np.shape(elites)[0]):
        order = np.argsort(fitness_reduced, axis= 0)
        population = np.delete(population, order[0,0], 0)
        fitness_reduced = np.delete(fitness_reduced, order[0,0], 0)

    probability = fitness_reduced / fitness_reduced.sum()
    select_range = np.zeros((np.shape(population)[0],2))
    for i in range(np.shape(select_range)[0]):
        if i == 0:
            select_range[i,0] = 0
            select_range[i,1] = probability[i,0]
        else:
            select_range[i,0] = select_range[i-1,1]
            select_range[i,1] = select_range[i,0] + probability[i,0]

    newPopulation = np.zeros((0,np.shape(population)[1]))
    for i in range(np.shape(population)[0]):
        ticket = np.random.rand()
        select_code = np.where(select_range < ticket)[0][-1]
        newPopulation = np.concatenate((newPopulation, np.array([population[select_code,:]])))

    return newPopulation