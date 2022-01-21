#from __future__ import absolute_import, division, print_function, unicode_literals
import array
import random
import json
# import operator  #  el objeto HallOfFame recibe un argumento llamdo "similar" que por defecto usa el metodo "operator.eq" de este modulo
import numpy as np
#import ctypes
from math import sqrt

#import networkx as nx
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import show

from deap import base
#from deap import benchmarks
from deap import creator
from deap import tools
from deap import algorithms
#from deap.benchmarks.tools import diversity, convergence, hypervolume

#from scipy.stats import pearsonr
#import os
import sys

import corr1 # 
import corr2 #
import corr3 #

sys.setrecursionlimit(100000) # ???

#global pcc_d      #
#global cont_pcc   #
#global cont_unif  #
global array_r    #
global train_char #
#global num_scores #
global index_n    #
global label      #
global hasha      #
global exp_scores #
global count_size #  Number of peptides for training  PCC and AUC

###########################
# Se crean los objetos para definir el problema y el tipo de indiviuos
# "FitnessMin"
#  base.Fitness
# weights=(1.0, 1.0)
creator.create("FitnessMin", base.Fitness, weights=(1.0, 1.0))
# "Individual"
#  array.array, typecode='d'
#  fitness=creator.FitnessMin
creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)

# Funcion que genera individuos aleatorios
def create_individual(low, up, size=None):
    try:
        return [random.uniform(a, b) for a, b in zip(low, up)]
    except TypeError:
        return [random.uniform(a, b) for a, b in zip([low] * size, [up] * size)]

# Este objeto (toolbox) permite registrar funciones que se utilizaran durante la operacion del GA
# El registro de funciones se realiza mediante el metodo "register" de la clase "base.Toolbox"                    
toolbox = base.Toolbox()

# defining the limits of the search for problem                                                                   
BOUND_LOW, BOUND_UP = -1.0, 1.0
NDIM = 9

# Generacion de individuos y poblacion inicial
toolbox.register("attr_float", create_individual, BOUND_LOW, BOUND_UP, NDIM)
toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_float)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

# Funcion de fitness: multiobjetivo 1) pcc 2) auc
def corr_pcc(individual):
    # CORR1
    # ????????
    # inputs:
    #   np.array(individual):
    #   label:
    #   train_char:
    #   hasha:
    #   array_r:
    #   exp_scores:
    #   index_n:
    #   count_size:
    #
    # outputs:
    #   pcc:
    #   auc:
    pcc, auc = corr1.corr_nsga2(np.array(individual), label, train_char, hasha, array_r, exp_scores, index_n, count_size)
    return pcc, auc

# Registro de operadores geneticos y funcion objetivo
# cruce
# tools.cxSimulatedBinaryBounded: 
toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0)
# mutacion
# tools.mutPolynomialBounded:
toolbox.register("mutate", tools.mutPolynomialBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0, indpb=1.0/NDIM)
# seleccion
# tools.selNSGA2:
toolbox.register("select", tools.selNSGA2)
# fitness
toolbox.register("evaluate", corr_pcc)

# Ejecucion del algoritmo multiobjetivo
# inputs:
#   NGEN: num_generations
#   MU:   num_individuals
#   seed: random_seed

def main(NGEN, MU, seed):
    print("*")
    print("*Starting evaluation")
    print("")
    random.seed(seed)

    CXPB = 0.7 # probabilidad de cruce (hiperparametro)  0.7 mas utilizado
    MUTPB = 0.3 # probabilidad de mutacion (hiperparametro)  0.3 mas utilizado
    MU = MU   # numero de individuos a seleccionar para cada generacion
    LAMBDA = MU  # numero de hijos a producir en cada generacion

    pop = toolbox.population(n=MU)
    # objeto tipo HallOfFame que almacena el mejor individuo encontrado a los largo de las generaciones del GA
    # es la estrategia que usa Deap para no perder al mejor individuos, debido a las operaciones geneticas que se aplican en la poblacion
    # HallOfFame se encuentra definida en submodulo tools
    # El objeto "hof" es actualizado mediante el metodo "update" en cada genracion del GA
    # se puede obtener el mejor indiviuo asi :  hof[0]
    # y su fitnes así: hof[0].fitness.values
    # recibe dos parametros: 
    #   maxsize: numero de indioviduos a almacenar
    #   similar: una funcion para comprar si dos individuos son igfuales. Por defecto utliza el metodo "operator.eq" del modulo operator
    hof = tools.HallOfFame(maxsize=1)##********************

    # Se define un objeto para generar las estadisticas de la poblacion a lo largo de las generaciones del algoritmo
    # Al crear el objeto se debe indicar sobre que atributo del os individuos se van a genrear las estadisticas
    # en este caso el "fitness"
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)

    # Este objeto permite almacenar todos los datos de evolucion del GA en un registro.
    # la informacion se almacena mediante diccionarios de Python
    # tiene metodos, y dos importantes son: "record" y "select"    
    logbook = tools.Logbook()


    # Este objeto ....
    pareto = tools.ParetoFront()

    # Ejecucion del algoritmo multiobjetivo
    #
    # algorithms.eaMuPlusLambda():
    #
    # inputs:
    #   pop:
    #   toolbox,
    #   mu:
    #   lambda_:
    #   cxpb:
    #   mutpb:
    #   ngen:
    #   stats:
    #   halloffame:
    #   verbose: default true, show stadistics for each iteration
    #   
    # outputs:
    #   pop: Poblacion final del algoritmo
    #   logbook: Registro de la evolucion
    pop, logbook = algorithms.eaMuPlusLambda(pop,
                                             toolbox,
                                             mu=MU,
                                             lambda_=LAMBDA,
                                             cxpb=CXPB,
                                             mutpb=MUTPB,
                                             ngen=NGEN,
                                             stats=stats,
                                             halloffame=pareto,
                                             verbose=True,
                                             )

#    history = tools.History()##********************
#    toolbox.decorate("mate", history.decorator)##********************
#    toolbox.decorate("mutate", history.decorator)##********************
#
#    # Evaluate the individuals with an invalid fitness
#    invalid_ind = [ind for ind in pop if not ind.fitness.valid]
#    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
#
#    for ind, fit in zip(invalid_ind, fitnesses):
#        ind.fitness.values = fit
#
#    # This is just to assign the crowding distance to the individuals
#    # no actual selection is done
#    pop = toolbox.select(pop, len(pop))
#
#    record = stats.compile(pop)
#    logbook.record(gen=0, evals=len(invalid_ind), **record)
#    print(logbook.stream)
#
#    print("      ")
#    print(" Individuals in generation 0: ")
#    print("      ")
#    j=0
#    for i in pop:
#        j+=1
#        print(j, ":", i, "PCC, AUC=", i.fitness)
#    print("      ")
#
#    # Begin the generational process
#    for gen in range(1, NGEN):
#          # Vary the population
#        offspring = tools.selTournamentDCD(pop, len(pop))
#        offspring = [toolbox.clone(ind) for ind in offspring]
#
#        for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
#            if random.random() <= CXPB:
#                toolbox.mate(ind1, ind2)
#
#            toolbox.mutate(ind1)
#            toolbox.mutate(ind2)
#            del ind1.fitness.values, ind2.fitness.values
#
#        # Evaluate the individuals with an invalid fitness
#        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
#        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
#        for ind, fit in zip(invalid_ind, fitnesses):
#              ind.fitness.values = fit
#
#        # Select the next generation population
#        pop = toolbox.select(pop + offspring, MU)
#        record = stats.compile(pop)
#        logbook.record(gen=gen, evals=len(invalid_ind), **record)
#        print(logbook.stream)
#
#        hof.update(pop)##********************
#
#        print("      ")
#        print(" Individuals in generation",gen,": ")
#        print("      ")
#        j=0
#        for i in pop:
#            j+=1
#            print(j, ":", i, "PCC, AUC=", i.fitness)
#        print("      ")
#
#    return pop, logbook, history, hof
    return pop, logbook, pareto

if __name__ == "__main__":
    # parsing command line arguments
    # revisar a futuro  https://docs.python.org/3/library/argparse.html
    # https://www.geeksforgeeks.org/command-line-arguments-in-python/
    num_generations=int(sys.argv[1])
    num_individuals=int(sys.argv[2])
    random_seed=int(sys.argv[3])

    train_txt = 'train1_DRB1_0101_e.txt'
    matrix_txt = 'matrix_DR1_label.txt'
    #cont_pcc = 0
    #cont_unif = 0
    if len(sys.argv)==5:
        train_txt = sys.argv[4]

    if len(sys.argv)==6:
        train_txt = sys.argv[4]
        matrix_txt = sys.argv[5]

    # PEPTIDE FILE SIZE
    count_size = len(open(train_txt).readlines())
    print("*")
    print("*Number of peptides (data to train): ", count_size)
    print("*")
    print("************Genetic Algorithm************")
    print("************Binding Predictor************")
    print("*")
    print("*Number of Generations: ", num_generations)
    print("*Number of Individuals: ", num_individuals)
    print("*Random seed: ", random_seed)

###########################################################################
    array_r = np.zeros((9,20), dtype='f')

    # CORR2 
    # ??????
    # inputs:
    #   count_size:
    #   train_txt:
    #   matrix_txt:
    #
    # outputs:
    #  array_r:
    #   label:
    #   index_n:
    array_r, label, index_n = corr2.store_nsga2(count_size, train_txt, matrix_txt)

    train_char = np.array((index_n,9), dtype='U')
    hasha = np.array((index_n), dtype='U')
    exp_scores = np.zeros((count_size), dtype='f')

    # CORR3
    train_char, hasha, exp_scores = corr3.store_nsga2(count_size, index_n, train_txt, matrix_txt)

#########################################################################
#    pop, stats, history, hof = main(num_generations, num_individuals, random_seed)
    pop, log, pareto = main(num_generations, num_individuals, random_seed)
#########################################################################
    #print("Resume")
    #print(stats)
    #print("--------")


    #print(toolbox.individual())
    #print(toolbox.individual())
    #print("Convergence: ", convergence(pop, optimal_front))
    #print("Diversity: ", diversity(pop, optimal_front[0], optimal_front[-1]))

    # import matplotlib.pyplot as plt
    # import numpy

    #from matplotlib.pyplot import imshow, show##********************
    #*************************
    #gen_best = history.getGenealogy(hof[0])##********************
    #print(gen_best)
    #print(hof[0])
    #graph = nx.DiGraph(gen_best).reverse()##********************networkx
    #nx.draw(graph)
    #show()
    #*************************

    #pos = nx.spring_layout(graph)

    #nx.draw_networkx_nodes(graph, pos, cmap=plt.get_cmap('jet'), node_size = 500)
                        #node_color = values, node_size = 500)
    #nx.draw_networkx_labels(graph, pos)
    #nx.draw_networkx_edges(graph, pos, edgelist=red_edges, edge_color='r', arrows=True)
    #nx.draw_networkx_edges(graph, pos, edgelist=black_edges, arrows=False)
    #nx.draw(graph)##********************
    #imshow(stats.avg[0], origin="lower")##********************
    #show()##********************

    #front = np.array([ind.fitness.values for ind in pop])
    #optimal_front = np.array(optimal_front)
    #plt.scatter(optimal_front[:,0], optimal_front[:,1], c="r")
    #plt.scatter(front[:,0], front[:,1], c="b")
    #plt.axis("tight")
    #plt.show()
