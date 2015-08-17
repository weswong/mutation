import itertools
import collections
import  numpy
import numpy as np
import math

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

def count_proportions(data):
    counts = collections.Counter(data)
    proportions = dict(zip(counts.keys(), counts.values() / np.sum(counts.values(), dtype=float)))
    return proportions

# skeleton code for liver stage dynamics...if needed? ignoring this for now
'''def liver_stage(N0):
    #pre-baked demographic histories that take place in the liver
    
    # assuming exponential growth, total merozoite population given by Johnston-Fidock Model at time 0 ~2.5e5
    total_N = demo_function(0)
    # each sporozoite makes its own little subpopulation. Assume they are all equally likely to survive and population size of each is the same
    subpopulation_Nfinal = demo_function/N0
    # merozoite emergence happens in 
    r = np.log(subpopulation_Nfinal / N0) / 5.0'''
    

class Genome:
    '''
    Our genome. We only care about the Genome ID and the mutations it carries.
    For simplicities sake, we are assuming infinite alleles model. Only way for two genomes
    to have the same mutation is through inheritance
    '''

    newid = itertools.count().next
    mutation_id = itertools.count().next
    
    def __init__(self, freq,id=None, mutations=None, ancestry=None): 
        # add frequency here
        if not id and not mutations:
            self.id = Genome.newid()
            self.mutations=[]
        else:
            self.id = id
            self.mutations = mutations
        if ancestry:
            self.ancestry =[ancestry]
        else:
            self.ancestry = []
        
        self.freq = freq
    
    @classmethod
    def make_children(cls, genome, n_children):
        genome.n_children = n_children
        genome.n_mutant_children = numpy.random.poisson(4.85e-9*2.3e7*genome.n_children) #mu = 4.85e-9, genome_size = 2.3e7
        genome.n_nonmutant_children = genome.n_children - genome.n_mutant_children
            
    @classmethod
    def add_mutation(cls, genome, N_next, number=1):
        # assume infinite sites model
        id = cls.newid()
        mutation = genome.mutations + [cls.mutation_id() for _ in range(number)]
        ancestry = genome.id
        return cls(1./N_next, id, mutation, ancestry)
    
    @classmethod
    def create_mutant_pool(cls, genome, mutants, N):
        return [Genome.add_mutation(genome, N) for x in range(mutants)]

class Population:
    ''' Following the population dynamics of a single LINEAGE
    a population here is defined as the family tree of a single genetically distinct
    sporozoite
    I do not simulate '''
    
    def __init__(self, strains, N_current):
        self.strains = strains
        self.generations = 0
        #genomes is a dictionary that shows the frequency of each genome object in the population
        
        self.N_current = N_current
        
        self.strain_id = [strain.id for strain in self.strains]
        self.strain_freqs = [strain.freq for strain in self.strains]
    
    def __iter__(self):
       for strain in self.strains:
          yield strain
    
    def get_population_stats(self):
        return zip([strain.id for strain in self.strains], [strain.freq for strain in self.strains])
    
    def get_mutation_freqs(self):
        mutations = list(itertools.chain(*[strain.mutations for strain in self.strains]))
        mutation_id, counts = np.unique(mutations, return_counts=True)
        return mutation_id, counts / float(self.N_current)
    
    
    @classmethod
    def instantiate_population(cls, n_ihepatocytes):
        # assume all success infected hepatocytes survive and are at the same frequency, change later
        population = [Genome(1./n_ihepatocytes) for _ in range(n_ihepatocytes)]
        #for strain in population:
        #    strain.make_children(N_next)
        
        # start at emergence from liver
        return cls(population, n_ihepatocytes)
    
    @classmethod
    def advance_generation(cls, population, N_next):
        mutant_pool = []
        failed_genomes = []
        n_sampling_population = 0
        
        sampling = numpy.random.choice(population.strains, size=N_next, p=population.strain_freqs)
        sampled_strains, counts = numpy.unique(sampling, return_counts=True)
        
        for strain, count in zip(sampled_strains,counts):
            Genome.make_children(strain, count)
            mutants = Genome.create_mutant_pool(strain, strain.n_mutant_children, N_next)
            mutant_pool += mutants
            strain.freq = float(strain.n_nonmutant_children) / N_next

        #updating stats
        population.strains = list(sampled_strains) +mutant_pool    
        population.N_current = N_next
        population.generations += 1
        population.strain_id = [strain.id for strain in population.strains]
        population.strain_freqs = [strain.freq for strain in population.strains]

class Simulation:
    def __init__(self):
        self.sim_duration = 4   # days
        self.sim_tstep = 2        # 2 days/mitotic generation
        self.sim_N0 = 10
        self.output = 'simulation'
        self.day = 0
        self.population = Population.instantiate_population(10)
    
    
    def update(self):
        dt = self.sim_tstep
        self.day += dt
        
        #Population.advance_generation(self.population, asexual_demo_funcion(self.day))
        Population.advance_generation(self.population, 100)
        print self.day
        numpy.save(self.output + '/' + str(self.day) + '_' + str(self.population.N_current),self.population.get_mutation_freqs())
        
    def run(self):
        for t in range(self.sim_duration / self.sim_tstep):
            self.update()
            
s = Simulation()
s.run()