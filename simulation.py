#!/usr/bin/python2.7
'''
Die Klasse erstellt ein Objekt 'Simulation' und benoetigt zu aller erst ein
zeitabhaengiges Netzwerk. Die Funktionen: sir(), sis(), si() und leslie() fuehren
die entsprechenden Simulationen durch und geben die aktuelle und akkumulierte
Pfaddichte zurueck.
'''
import numpy as np
from scipy import sparse as sp
from einlesen import readdata, load_sparse_csr

class Simulation:
    def __init__(self):
        self.size                    = 0
        self.runtime                 = 0
        self.simulation_time         = 0
        self.graph                   = ""
        self.isjob                   = False
        self.model                   = ""
        self.memory                  = 0
        self.mean_cumulative   = 0
        self.mean_prevalence = 0
        self.r0                      = False
        self.memory_efficient        = False
        
    def load(self, file):
        results = np.load(file)
        self.graph = results['graph']
        self.size = results['size']
        self.runtime = results['runtime']
        self.model = results['model']
        self.memory = results['memory']
        if 'mean_cumulative' in results:
            self.mean_cumulative = results['mean_cumulative']
        if 'mean_prevalence' in results:
            self.mean_prevalence = results['mean_prevalence']
        
    def save(self,path="/home/andreasko/reachability/simulationen/", fname = "simulation"):
        print         
        print "save to ..."
        self.path = str(path)
        self.fname = str(fname)
        print self.path + str(self.fname)
        np.savez(self.path + str(self.fname), graph=self.graph, size=self.size, runtime=self.runtime, model=self.model, memory=self.memory, mean_cumulative=self.mean_cumulative, mean_prevalence=self.mean_prevalence, r0=self.r0, simulation_time=self.simulation_time)
        
    def check(self,A):
        cumulative_SI = A[0].copy() + sp.eye(self.size, self.size, 0, format='csr')
        cumulative_SI.data = np.ones_like(cumulative_SI.data)
        
        incidence_class    = [sp.csr_matrix((self.size,self.size)) for ii in range(self.runtime)]
        incidence_class[0] = A[0]
        incidence_class[1] = sp.eye(self.size,self.size,format='csr')
        cumulative_SIR  = incidence_class[0] + incidence_class[1]
        prevalence  = incidence_class[0] + incidence_class[1]
        check       = np.empty(shape=(self.runtime,),dtype=bool)
        check[0]    = (cumulative_SI != cumulative_SIR).nnz > 0                        # check == 1 wenn wenigsten ein Element verschieden ist
        
        for ii in range(1,self.runtime):
            try:
                cumulative_SI      = cumulative_SI + A[ii]*cumulative_SI
                cumulative_SI.data = np.ones_like(cumulative_SI.data)
                
                incidence      = A[ii].dot(prevalence)
                incidence.data = np.ones_like(incidence.data)
                incidence      = incidence - incidence.multiply(cumulative_SIR)
                
                prevalence = prevalence - incidence_class[-1]
                prevalence = prevalence + incidence
                cumulative_SIR = cumulative_SIR     + incidence
                
                for jj in xrange(self.runtime-1,0,-1):
                    incidence_class[jj] = incidence_class[jj-1]
                incidence_class[0]      = incidence
                check[ii]    = (cumulative_SI != cumulative_SIR).nnz > 0
                print ii, check[ii]
            except:
                print 'Break at t = ', ii
                break
        return check
        
    def si(self,A):
        cumulative = A[0] + sp.eye(self.size, self.size, 0, format='csr')
        cumulative.data = np.ones_like(cumulative.data)
        mean_cumulative = np.zeros((1,self.runtime))[0]
        
        for ii in range(1,self.runtime):
            cumulative = (A[ii] + sp.eye(self.size, self.size, 0, format='csr'))*cumulative
            cumulative.data = np.ones_like(cumulative.data)
            mean_cumulative[ii] = cumulative.nnz
        self.mean_cumulative = mean_cumulative.astype(float)/self.size**2
        return
        
    def sir(self,A):
        print "Start the Simulation ... \n"
        #=======================================================================
        # sir() berechnet ein Ausbruchsszenario nach dem SIR-Modell.
        # Als Adjazenzmatrizen wird eine gegebene Serie unter 'graph.dat' geladen.
        # Bei einer infektioesen Periode von N werden N+1 Infektions-Klassen initialisiert.
        # Die Erste ist die Einheitsmatrix und alle anderen sind Nullmatrizen, wobei die
        # Letzte den R-Zustand darstellt. Kein Gedaechtnis (memory=0) entspricht einer I-
        # und einer R-Klasse.
        #=======================================================================
        #=======================================================================
        # Die zeitabhaengige Zugangsmatrix wird berechnet
        #=======================================================================
        if self.r0==False and self.memory_efficient==False:
            incidence_class    = [sp.csr_matrix((self.size,self.size)) for ii in range(self.memory+1)]
            incidence_class[0] = sp.eye(self.size,self.size,format='csr')
            cumulative      = sp.eye(self.size,self.size,format='csr')
            prevalence  = sp.eye(self.size,self.size,format='csr')
            
            mean_cumulative   = np.zeros((1,self.runtime))[0]
            mean_prevalence = np.zeros((1,self.runtime))[0]
            memory = self.memory
            for jj in range(0,self.runtime):
                
                #incidence = reduce(sparse_add,(A[jj].dot(cl) for cl in incidence_class[0:-1]))
                incidence      = A[jj].dot(prevalence)
                incidence.data = np.ones_like(incidence.data)
                incidence      = incidence - incidence.multiply(cumulative)
                
                prevalence = prevalence - incidence_class[-1]
                prevalence = prevalence + incidence
                cumulative     = cumulative     + incidence
                
                for ii in xrange(memory,0,-1):
                    incidence_class[ii] = incidence_class[ii-1]
                incidence_class[0]      = incidence
                
                mean_cumulative[jj]   = cumulative.nnz
                mean_prevalence[jj] = prevalence.nnz
                
            self.mean_cumulative   = mean_cumulative.astype(float)/self.size**2
            self.mean_prevalence = mean_prevalence.astype(float)/self.size**2
            
        elif self.r0==False and self.memory_efficient==True:
            mean_cumulative = np.zeros((self.runtime,))
            mean_prevalence = np.zeros((self.runtime,))
            memory = self.memory
            
            for ii in range(self.size):
                if ii%10000 == 0:
                    print 'node ', ii, ' of ', self.size
                incidence_class    = [sp.csr_matrix((self.size,1)) for kk in range(self.memory+1)]
                startvektor        = np.zeros((self.size,1), dtype=float)
                startvektor[ii]    = 1.
                incidence_class[0] = sp.csr_matrix(startvektor)
                cumulative         = sp.csr_matrix(startvektor)
                prevalence         = sp.csr_matrix(startvektor)
                
                for jj in range(0,self.simulation_time):
                    incidence      = A[jj].dot(prevalence)
                    incidence.data = np.ones_like(incidence.data)
                    incidence      = incidence - incidence.multiply(cumulative)
                    
                    prevalence = prevalence - incidence_class[-1]
                    prevalence = prevalence + incidence
                    cumulative = cumulative + incidence
                    
                    for ii in xrange(memory,0,-1):
                        incidence_class[ii] = incidence_class[ii-1]
                    incidence_class[0] = incidence
                    
                    mean_cumulative[jj] += cumulative.nnz
                    mean_prevalence[jj] += prevalence.nnz
            
            self.mean_cumulative = mean_cumulative.astype(float)/self.size**2
            self.mean_prevalence = mean_prevalence.astype(float)/self.size**2
        elif self.r0==True:
            incidence_class    = [sp.csr_matrix((self.size,self.size)) for ii in range(self.memory+1)]
            incidence_class[0] = sp.eye(self.size,self.size,format='csr')
            cumulative      = sp.eye(self.size,self.size,format='csr')
            prevalence  = sp.eye(self.size,self.size,format='csr')
            
            mean_cumulative   = np.zeros((1,self.runtime))[0]
            mean_prevalence = np.zeros((1,self.runtime))[0]
            
            r0   = np.zeros((self.runtime,))
            size = self.size**2
            memory = self.memory
            for jj in range(0,self.runtime):
                print jj
                #incidence = reduce(sparse_add,(A[jj].dot(cl) for cl in incidence_class[0:-1]))
                incidence      = A[jj].dot(prevalence)
                incidence.data = np.ones_like(incidence.data)
                incidence      = incidence - incidence.multiply(cumulative)
                
                prevalence = prevalence - incidence_class[-1]
                prevalence = prevalence + incidence
                cumulative     = cumulative     + incidence
                
                for ii in xrange(memory,0,-1):
                    incidence_class[ii] = incidence_class[ii-1]
                incidence_class[0]      = incidence
                
                susceptible = float(size - cumulative.nnz)
                infected    = float(prevalence.nnz)
                r0[jj]      = float(incidence.nnz) / susceptible / infected
            self.r0 = r0
        else:
            print self.r0, type(self.r0)
            print self.memory_efficient, type(self.memory_efficient)
        return
        
    def sis(self,A):
        
        incidence_class    = [sp.csr_matrix((self.size,self.size)) for ii in range(self.memory+1)]
        incidence_class[0] = sp.eye(self.size,self.size,format='csr')
        cumulative      = sp.eye(self.size,self.size,format='csr')
        prevalence  = sp.eye(self.size,self.size,format='csr')
        
        mean_cumulative   = np.zeros((1,self.runtime))[0]
        mean_prevalence = np.zeros((1,self.runtime))[0]
        
        #=======================================================================
        # Die zeitabhaengige Zugangsmatrix wird berechnet
        #=======================================================================
        memory = self.memory
        for jj in range(0,self.runtime):
            
            #incidence = reduce(sparse_add,(A[jj].dot(cl) for cl in incidence_class[0:-1]))
            incidence      = A[jj].dot(prevalence)
            incidence.data = np.ones_like(incidence.data)
            incidence      = incidence - incidence.multiply(prevalence)
            
            prevalence  = prevalence - incidence_class[-1]
            prevalence  = prevalence + incidence
            cumulative      = cumulative     + incidence
            cumulative.data = np.ones_like(cumulative.data)
            
            incidence_class[-1]     = incidence_class[-1] + incidence_class[-2]
            for ii in xrange(memory,0,-1):
                incidence_class[ii] = incidence_class[ii-1]
            incidence_class[0]      = incidence
            
            mean_cumulative[jj] = cumulative.nnz
            mean_prevalence[jj] = prevalence.nnz
            
        self.mean_cumulative   = mean_cumulative.astype(float)/self.size**2
        self.mean_prevalence = mean_prevalence.astype(float)/self.size**2
        return
        
    def leslie(self,A):
        
        incidence_class    = [sp.csr_matrix((self.size,self.size)) for ii in range(self.memory+1)]
        incidence_class[0] = sp.eye(self.size,self.size,format='csr')
        cumulative      = sp.eye(self.size,self.size,format='csr')
        prevalence  = sp.eye(self.size,self.size,format='csr')
        
        mean_cumulative   = np.zeros((1,self.runtime))[0]
        mean_prevalence = np.zeros((1,self.runtime))[0]
        #=======================================================================
        # Die zeitabhaengige Zugangsmatrix wird berechnet
        #=======================================================================
        sparse_add = sp.compressed._cs_matrix.__add__
        memory     = self.memory
        for jj in range(0,self.runtime):
            
            incidence = A[jj].dot(prevalence)
            for ii in xrange(memory,0,-1):
                incidence_class[ii] = incidence_class[ii-1]
            incidence_class[0] = incidence
            
            prevalence  = reduce(sparse_add, incidence_class)
            cumulative      = cumulative + incidence
            cumulative.data = np.ones_like(cumulative.data)
            
            mean_cumulative[jj]   = cumulative.nnz
            mean_prevalence[jj] = prevalence.nnz
            
        self.mean_cumulative   = mean_cumulative.astype(float)/self.size**2
        self.mean_prevalence = mean_prevalence.astype(float)/self.size**2
        return
        
    def aggregateNetwork(self,A):
        AA = A[0].astype(float)
        for ii in range(1,self.runtime):
            AA = A[ii].astype(float) + AA
            AA.data = np.ones_like(AA.data)
        
        A = []
        for ii in range(self.runtime):
            A.append(AA)
            
        return A
        
    def runSimulation(self,
                      graph="sexual_contacts.dat",
                      isjob=False,
                      model="SIR",
                      memory=1,
                      r0=False,
                      memory_efficient=True,
                      simulation_time=0):
        
        if isinstance(graph, basestring):
            if ((str(graph) == "sociopatterns_hypertext.dat") or (str(graph) == "sexual_contacts.dat") or (str(graph) == "testnetwork.txt")):
                A = readdata(str(graph), symmetrisch=True)
                self.graph = str(graph)
            elif str(graph)=="sociopattern_mlist.p" or str(graph)=="sociopattern_mlist_clean.p":
                A = load_sparse_csr(graph)
                self.graph = str(graph)
            elif str(graph)=="pig_trade_11-14_uvdw_from_v0_2.dat":
                A = readdata("HIT/"+str(graph), symmetrisch=False)
                self.graph = str(graph)
            elif str(graph)=="pig_trade_11-14_uvdw_from_v0_2.p":
                A = load_sparse_csr("HIT/"+str(graph))
                self.graph = str(graph)
        else:
            A = graph
            self.graph = "test_matrix"
        
        self.memory           = int(memory)
        self.size             = A[0].shape[0]
        self.runtime          = len(A)
        self.isjob            = bool(int(isjob))
        self.model            = str(model)
        self.r0               = bool(int(r0))
        self.memory_efficient = bool(int(memory_efficient))
        if int(simulation_time)  == 0:
            self.simulation_time = self.runtime
        else:
            self.simulation_time = int(simulation_time)
        
        if self.model == 'SI':
            self.memory = 0
            self.si(A)
        elif self.model == 'SIS':
            self.sis(A)
        elif self.model == 'SIR':
            self.sir(A)
        elif self.model == 'Leslie':
            self.leslie(A)
        elif self.model == 'check':
            check = self.check(A)
            if np.sum(check) == 0:
                print 'Die Unfolding cumulativeibility und das SI-Modell sind identisch'
            else:
                print 'Es traten ', np.sum(check), ' Abweichungen auf'

#===============================================================================
# Es folgt die Testumgebung:
#===============================================================================
#'load either sociopatterns_hypertext.dat, pig_trade_11-14_uvdw_from_v0_2.dat,  sexual_contacts.dat or nothing: []')
if __name__ == "__main__":
    A = []
    A.append(sp.csr_matrix(([1,1], [[1,0],[0,1]]),shape=(4,4)))
    A.append(sp.csr_matrix(([1,1], [[2,1],[1,2]]),shape=(4,4)))
    A.append(sp.csr_matrix(([1,1], [[3,2],[2,3]]),shape=(4,4)))
    A.append(sp.csr_matrix(([1,1], [[3,0],[0,3]]),shape=(4,4)))
    test = Simulation()
    test.runSimulation(graph           = "pig_trade_11-14_uvdw_from_v0_2.dat",
                       model           = "check",
                       memory          = 350,
                       memory_efficient= True,
                       simulation_time = 0)
    test.save("simulationen/")
#    test.load('/users/stud/koher/master/python/2014_02_26/results/m=3000_SIR.npz')
    #test.get_mean_cumulative()
    #test.get_mean_prevalence()
#    test.get_infection_time()
#    import matplotlib.pyplot as plt
#    Tbin = range(1,len(test.infection_time_distribution)+1)
#    plt.bar(Tbin, test.infection_time_distribution, width=1, color='black',)
#    plt.show()
