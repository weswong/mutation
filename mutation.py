import itertools
import collections
import  numpy
import numpy as np
import math
import json

#general functions --------------------------------------------
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
    #N = numpy.power(10, log_N) * 5000000
    N = numpy.power(10, log_N)
    return math.floor(N)


def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

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

#Classes ------------------------------------------------------------------

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
            self.ancestry =ancestry
        else:
            self.ancestry = [self.id]
        
        self.freq = freq
    
    
    
    @classmethod
    def add_mutation(cls, genome, N_next, number=1):
        # assume infinite sites model
        id = cls.newid()
        mutation = genome.mutations + [cls.mutation_id() for _ in range(number)]
        ancestry = genome.ancestry + [id]
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
        
        self.strain_ids = [strain.id for strain in self.strains]
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
    
    def get_mutation_classes(self):
        mutation_class_dict = {}
        for strain in self.strains:
            mutation_class = len(strain.mutations)
            if mutation_class in mutation_class_dict.keys():
                mutation_class_dict[mutation_class] += 1 * strain.freq * self.N_current
            else:
                mutation_class_dict[mutation_class] = 1 * strain.freq * self.N_current
        return mutation_class_dict
    
    @classmethod
    def instantiate_population(cls, n_ihepatocytes):
        # assume all success infected hepatocytes survive and are at the same frequency, change later
        population = [Genome(1./n_ihepatocytes) for _ in range(n_ihepatocytes)]
        #for strain in population:
        #    strain.make_children(N_next)
        
        # start at emergence from liver
        return cls(population, n_ihepatocytes)
    
    def sampling(self, n_sample):
        sampling = np.random.choice(self.strain_ids, size=n_sample, p=self.strain_freqs)
        sampled_strains, counts = np.unique(sampling, return_counts=True)
        return sampled_strains, counts
    
    
    @classmethod
    def advance_generation(cls, population, N_next):        
        pop_mutants = numpy.random.poisson(4.85e-9*2.3e7*N_next)
        if pop_mutants < N_next:
            pop_nonmutants = N_next - pop_mutants
            
        else:
            pop_nonmutants = 0.
            pop_mutants = N_next        
        

        nonmutant_sampled_strains, nonmutant_counts = population.sampling(pop_nonmutants)
        nommutant_sampling_strains = []
        for strain_id, count in zip(nonmutant_sampled_strains, nonmutant_counts):
            strain = population.strains[population.strain_ids.index(strain_id)]
            nommutant_sampling_strains.append(strain)
            strain.freq = float(count) / N_next
        
        mutant_sampled_strains, mutant_counts = population.sampling(pop_mutants)
        mutant_pool = []
        for mutant_id, count in zip(mutant_sampled_strains, mutant_counts): 
            mutant_strain = population.strains[population.strain_ids.index(mutant_id)]
            mutants = Genome.create_mutant_pool(mutant_strain, count, N_next)
            mutant_pool += mutants
        
        
        
        #updating stats
        population.strains = list(nommutant_sampling_strains) + mutant_pool   
        for strain in population.strains[0:10]:
            print strain.__dict__ 
        population.N_current = N_next
        population.generations += 1
        population.strain_ids = [strain.id for strain in population.strains]
        population.strain_freqs = [strain.freq for strain in population.strains]
	#print zip(population.strain_id, population.strain_freqs)
	
class Simulation:
    def __init__(self):
        self.sim_duration = 366   # days
        self.sim_tstep = 2        # 2 days/mitotic generation
        self.sim_N0 = 10
        self.output = 'simulation'
        self.day = 0
        self.population = Population.instantiate_population(10)
        self.capture_days = [2, 20, 46, 12, 25, 80]
    
    def update(self):
        dt = self.sim_tstep
        self.day += dt
        
        #Population.advance_generation(self.population, asexual_demo_funcion(self.day))
        if self.day > self.sim_duration:
            print 'finish'
        elif asexual_demo_function(self.day) < 1:
            print 'no more parasites'
        else:
            print 'starting day ', 
            print self.day, asexual_demo_function(self.day)
            Population.advance_generation(self.population,asexual_demo_function(self.day))
            #if self.day == self.sim_duration or self.day in self.capture_days:
            if self.day:
                numpy.save(self.output + '/' + str(self.day) + '_' + str(self.population.N_current),self.population.get_mutation_freqs())
                
                with open(self.output + '/' + 'mutation_classes_{day}.txt'.format(day = self.day), 'w') as outfile:
                    json.dump(self.population.get_mutation_classes(), outfile)
    
    def run(self):
        for t in range(self.sim_duration / self.sim_tstep):
            self.update()
            
if __name__ == "__main__":
    s = Simulation()
    s.run()
    
    
