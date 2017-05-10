from blue import BLUE

blue = BLUE(2)
blue.AddMeasurement([54.6,  43.3])
blue.AddUncertainty([15.3,  6.6],  0) #stat
blue.AddUncertainty([6.5,   5.2],  1) #lumi
blue.AddUncertainty([2.0,   1.6],  0) #lepton eff
blue.AddUncertainty([5.07,  1.95], 1) #JES
blue.AddUncertainty([0.11,  0.04], 1) #JER
blue.AddUncertainty([2.1,   1.6],  1) #theory


blue.Simple()

blue.Iterative(5)
