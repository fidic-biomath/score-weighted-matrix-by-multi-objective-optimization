import array
import random
import json
import numpy as np
from math import sqrt

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

import sys

import corr1 # 
import corr2 #
import corr3 #

sys.setrecursionlimit(100000) # ???

global array_r    #
global data_peptide_core_chars #
global index_n    #
global label      #
global hasha      #
global exp_scores #
global dataset_size #  Number of peptides for training  PCC and AUC

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
    #   data_peptide_core_chars:
    #   hasha:
    #   array_r:
    #   exp_scores:
    #   index_n:
    #   dataset_size:
    #
    # outputs:
    #   pcc:
    #   auc:
    pcc, auc = corr1.corr_nsga2(np.array(individual), label, data_peptide_core_chars, hasha, array_r, exp_scores, index_n, dataset_size)
    #pcc = np.person
    #auc = np.roc.auc
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

#######################################################################
# Funciones extras antes de la funcion main()

def plot_frente():
    """
    Representación del frente de Pareto que hemos obtenido
    """
    datos_pareto = np.loadtxt("fitness.log", delimiter=",")
    plt.scatter(datos_pareto[:, 0], datos_pareto[:, 1], s=30)

    plt.xlabel("PCC")
    plt.ylabel("AUC")
    plt.grid(True)
    plt.legend(["Pareto obtenido"], loc="upper right")
    #plt.savefig("ParetoBenchmark.eps", dpi=300, bbox_inches="tight")
    plt.show()

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
    return pop, logbook, pareto

if __name__ == "__main__":
    # parsing command line arguments
    # revisar a futuro  https://docs.python.org/3/library/argparse.html
    # https://www.geeksforgeeks.org/command-line-arguments-in-python/
    num_generations=int(sys.argv[1])
    num_individuals=int(sys.argv[2])
    random_seed=int(sys.argv[3])
    data_peptides = sys.argv[4]
    score_matrix = sys.argv[5]

    # PEPTIDE FILE SIZE
    dataset_size = len(open(data_peptides).readlines())
    print("*")
    print("*Number of peptides (data to train): ", dataset_size)
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
    #   dataset_size:
    #   data_peptides:
    #   score_matrix:
    #
    # outputs:
    #  array_r <class 'numpy.ndarray'>:  score_matrix
    #   label <class 'bytes'>:   columns labels of array_r
    #   index_n <class 'int'>:  number of 9-mer possible peptides
    array_r, label, index_n = corr2.store_nsga2(dataset_size, data_peptides, score_matrix)
    print( array_r )
    print( label )
    print( index_n )


    data_peptide_core_chars = np.array((index_n,9), dtype='U')
    hasha = np.array((index_n), dtype='U')
    exp_scores = np.zeros((dataset_size), dtype='f')

    # CORR3
    # ??????
    # inputs:
    #   dataset_size:
    #   index_n:
    #   data_peptides: 
    #   score_matrix:
    # outputs:
    #   data_peptide_core_chars <class 'numpy.ndarray'>:  9-mer possible peptides by chars
    #   hasha <class 'numpy.ndarray'>:  ????
    #   exp_scores <class 'numpy.ndarray'>:  experimental peptide binding scores 
    data_peptide_core_chars, hasha, exp_scores = corr3.store_nsga2(dataset_size, index_n, data_peptides, score_matrix)
    np.set_printoptions(threshold=sys.maxsize)
    print( type(data_peptide_core_chars) )
    print( data_peptide_core_chars.shape )
    print( data_peptide_core_chars )
    print( type(hasha) )
    print( hasha.shape )
    print( hasha )
    print( type(exp_scores) )
    print( exp_scores.shape )
    print( exp_scores )


#########################################################################
    pop, log, pareto = main(num_generations, num_individuals, random_seed)
#########################################################################
    res_individuals = open("individuals.log", "w")
    res_fitness = open("fitness.log", "w")
    for ind in pareto:
        res_individuals.write(str(ind))
        res_individuals.write("\n")
        res_fitness.write(str(ind.fitness.values[0]))
        res_fitness.write(",")
        res_fitness.write(str(ind.fitness.values[1]))
        res_fitness.write("\n")
    res_individuals.close()
    res_fitness.close()

