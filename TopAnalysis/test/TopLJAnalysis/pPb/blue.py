from numpy import *

class BLUE:
    """General implementation of BLUE for N measurements, by Andrea Giammanco (2016)"""
    """Caveat: here I am assuming that correlation is the same for any pair-wise combination of the N measurements for any given systematic."""
    """That's ok in many cases, but keep in mind: not always legit."""
 
    def __init__(self,dim=2):
        self.values = matrix(zeros((1,dim)))
        self.ematrix = matrix(zeros((dim,dim)))
        self.n = dim

    def AddMeasurement(self,v):
        self.values = v

    def AddUncertainty(self,v,corr):
        e = zeros(self.n)
        uncmatrix = matrix(zeros((self.n,self.n)))
        for i in range(self.n):
            e[i] = v[i]
            for j in range(self.n):
                e[j] = v[j]
                if i == j:
                    uncmatrix[i,j] = e[i]*e[j]
                else:
                    uncmatrix[i,j] = e[i]*e[j]*corr

        self.ematrix = self.ematrix + uncmatrix

    def Simple(self,scaleFactors=None):

        if (scaleFactors == None): print 'ematrix: \n', self.ematrix # reason: too verbose for the iterative case

        unitvector = matrix(ones((1,self.n))).T 

        # the following lines are motivated by the iterative case, otherwise there would be no need to define newmatrix, we could just use self.ematrix
        # note: probably numpy provides a more compact way to scale by column and by row?
        newmatrix = matrix(zeros((self.n,self.n)))
        for i in range(self.n):
            for j in range(self.n):
                if scaleFactors != None: 
                    newmatrix[i,j] = self.ematrix[i,j]*scaleFactors[i]*scaleFactors[j]
                else:
                    newmatrix[i,j] = self.ematrix[i,j]

        self.weights = newmatrix.I*unitvector/(unitvector.T*newmatrix.I*unitvector)

        self.weights = self.weights/sum(self.weights)
        print 'weights: ', self.weights.T#, ' (sum =', sum(self.weights), ')'
        
        self.average = self.values*self.weights
        print 'average: ', self.average

        self.error = sqrt(self.weights.T*newmatrix*self.weights)
        print 'error: ', self.error

        return [self.average,self.error]

    def Iterative(self, numiter): 
        """Read motivation here: https://arxiv.org/abs/1405.3425"""
        """Beware: this assumes that all uncertainties are proportional to the average; in real life, some are proportional and some are absolute"""
        """My personal rule of thumb is to run both the plain and the iterative version, and check how distant are the results. If not much, one can safely assume that these subtleties are unimportant."""
        
        scaleFactors = ones(self.n)
        for i in range(numiter):
            print "----- iteration #", i
            tempAverage = self.Simple(scaleFactors)[0]
            for j in range(self.n):
                scaleFactors[j] = tempAverage/(self.values[j])
            if i < numiter-1: # don't print it for the last iteration
                print "Scale factor on the errors:", scaleFactors, ", to be applied at iteration", i+1

        return [self.average,self.error]

class BLUE2x2:
    """Original class in collaboration with Steffen Roecker, hardcoded for 2 measurements"""
    """It may be useful for testing"""
 
    def __init__(self):
        self.values = matrix([0,0])
        self.ematrix = matrix([[0,0],[0,0]])

    def AddMeasurement(self,v1,v2):
        self.values = matrix([v1,v2])

    def AddUncertainty(self,e1,e2,corr):
        uncmatrix = matrix([[e1**2,e1*e2*corr],[e1*e2*corr,e2**2]])
        self.ematrix = self.ematrix + uncmatrix

    def Simple(self):
        print 'ematrix: \n', self.ematrix

        unitvector = matrix([1,1]).T
        self.weights = self.ematrix.I*unitvector/(unitvector.T*self.ematrix.I*unitvector)

        print 'weights: ', self.weights.T
        
        self.average = self.values*self.weights
        print 'average: ', self.average

        self.error = sqrt(self.weights.T*self.ematrix*self.weights)
        print 'error: ', self.error

        return [self.average,self.error]
