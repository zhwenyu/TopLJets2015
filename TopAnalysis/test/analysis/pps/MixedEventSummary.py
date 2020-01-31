class MixedEventSummary:

    """ 
    holds reconstructed proton information for an event and a set of pileup discriminant variables
    puDiscr           - a vector with floats
    {pos,neg}_protons - array with 3 columns of variable size
                        first column are protons from multiRP
                        second column are protons from far detectors (pixels in 2017)
                        third column are protons from near detectors (strips in 2017)
   
    """

    def __init__(self,puDiscr,pos_protons,neg_protons):
        self.puDiscr=puDiscr
        self.pos_protons=pos_protons
        self.neg_protons=neg_protons

    def getPuDiscriminants(self):

        """ getter for  pileup discriminants """
        
        return self.puDiscr

    def getProtons(self,pos=True):
        
        """ getter for protons """

        return self.pos_protons if pos else self.neg_protons
        

def main():
    print 'Defines MixedEventSummary class'

if __name__ == "__main__":
    sys.exit(main())
