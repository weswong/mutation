import itertools
import collections
import  numpy
import numpy as np
import math



class Simulation:
    def __init__(self):
        self.sim_duration = 365   # days
        self.sim_tstep = 2        # 2 days/mitotic generation
        self.sim_N0 = 10
        self.population = Population.instantiate_population(self.sim_N0)
        self.output = 'simulation'
        self.day = 0
        
    
    
    def update(self):
        dt = self.sim_tstep
        self.day += dt
        
        self.population.advance_generation(self.day)
        print self.day
        numpy.save(self.output + '/' + str(self.day), self.population.get_mutation_freqs())
        
    def run(self):
        for t in range(self.sim_duration / self.sim_tstep):
            self.update()
            


class Genome:
    '''
    Our genome. We only care about the Genome ID and the mutations it carries.
    For simplicities sake, we are assuming infinite alleles model. Only way for two genomes
    to have the same mutation is through inheritance
    '''

    newid = itertools.count().next
    mutation_id = itertools.count().next
    
    def __init__(self, id=None, mutations=None, parent=None): 
        if not id and not mutations:
            self.id = Genome.newid()
            self.mutations=[]
            self.parent = 'origin'
        else:
            self.id = id
            self.mutations = mutations
            self.parent = parent
        
    
    @classmethod
    def add_mutation(cls, genome):
        # assume infinite sites model
        id = cls.newid()
        mutation = genome.mutations + [cls.mutation_id()]
        parent = genome
        return cls(id, mutation, parent)


def count_proportions(data):
    counts = collections.Counter(data)
    proportions = dict(zip(counts.keys(), counts.values() / np.sum(counts.values(), dtype=float)))
    return proportions


def asexual_demo_function(t):
    '''functions that determine the blood parasitemia within the host at any given time (post emergence). 
    Time is measured in days, after emergence of merozoites from liver
    Functions obtained from fitting equations to the average results of Johnston Fidock Model (100X)'''
    if t <= 13:
        coefficients = [-0.03159309,0.82910821, -1.30758949]
    elif t > 13 and t <=23:
        coefficients = [3.20993043e-03,-1.53006247e-01,2.27611666e+00,-7.01357758e+00]
    elif t >23 and t<=50:
        coefficients = [ 4.55947621e-06,  -8.20047138e-04,   5.80259918e-02,-2.01844752e+00,3.44679750e+01,  -2.27557373e+02]
    elif t> 50 and t <= 80:
        coefficients = [1.20307345e-04,  -2.36885079e-02,1.50310628e+00,-2.87256521e+01]
    elif t > 80:
        coefficients = [8.11449765e-11, -1.27763731e-07, 6.28393586e-05, -2.97606465e-02, 3.57669339e+00]
    log_N = numpy.poly1d(coefficients)(t)
    N = numpy.power(10, log_N) * 5000000
    #N = numpy.power(10, log_N)
    return math.floor(N)


# skeleton code for liver stage dynamics...if needed? ignoring this for now
'''def liver_stage(N0):
    #pre-baked demographic histories that take place in the liver
    
    # assuming exponential growth, total merozoite population given by Johnston-Fidock Model at time 0 ~2.5e5
    total_N = demo_function(0)
    # each sporozoite makes its own little subpopulation. Assume they are all equally likely to survive and population size of each is the same
    subpopulation_Nfinal = demo_function/N0
    # merozoite emergence happens in 
    r = np.log(subpopulation_Nfinal / N0) / 5.0'''
    



class Population:
    ''' our population'''
    def __init__(self, genomes):
        #genomes is a dictionary that shows the frequency of each genome object in the population
        self.genomes = genomes
        self.genome_stats = count_proportions(genomes)
        
        self.unique_genomes = self.genome_stats.keys()
        self.id = [genome.id for genome in self.genome_stats.keys()]
        self.props = self.genome_stats.values()

    def get_population_stats(self):
        return zip(self.id, self.props)
    
    def get_mutation_freqs(self):
        mutations = list(itertools.chain(*[genome.mutations for genome in self.genomes]))
        return np.log(count_proportions(mutations).values())    
    
    def update(self):
        self.unique_genomes = self.genome_stats.keys()
        self.id = [genome.id for genome in self.genome_stats.keys()]
        self.props = self.genome_stats.values()
        self.population_genomes = self.genomes
        
    @classmethod
    def merge_population(cls, *arg):
        genomes = []
        for population in args:
            genomes += population.genomes
        return cls(genomes)
    
    @classmethod
    def instantiate_population(cls, initialN):
        population = [Genome() for _ in range(initialN)]
        population_stats = count_proportions(population)
        return cls(population)
    
    def resample_population(self, N):
        # Can't go all the way up to 1e12, will determine optimum subsample size later
        if N > 1e6:
            N = 1e6
            
        sampling = np.random.choice(self.unique_genomes, N, p=self.props)
        return sampling
        #self.genome_stats = count_proportions(sampling)
        #self.genomes = sampling
        #self.update()
        #return cls(population_stats)
    
    def mutate_population(self, N):
        N=float(N)
        mu = 4.85e-9
        #mu = 0
        if N > 1e6:
            N = 1e6
        mutation_events = numpy.random.poisson(mu * 2.3e7 * N)
        if mutation_events ==0:
            mutants = []
        else:
            vectorized_func = np.vectorize(Genome.add_mutation)
            mutants = vectorized_func(np.random.choice(self.unique_genomes, mutation_events, p = self.props))
        self.genomes = np.concatenate([mutants, self.resample_population(N-mutation_events)])
                                                                                                                        
        '''for idx in np.random.randint(0, len(self.genomes), mutation_events):
            #print str(population.genomes[idx].id) + ' is mutating to ',
            self.genomes[idx] = Genome.add_mutation(self.genomes[idx])
            #print str(population.genomes[idx].id)'''
        self.genome_stats = count_proportions(self.genomes)
        self.update()
        
    def advance_generation(self, dt):
        #dt = the day we are currently on, post emergence from liver
        #self.resample_population(asexual_demo_function(dt))
        self.mutate_population(asexual_demo_function(dt))
        
        

simulation = Simulation()
simulation.run()


