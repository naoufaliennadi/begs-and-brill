import numpy as np

gascr = .00685    # m^3/s
oilcr = .00368    # m^3/s
diameter = .0635  # m

rohg = 42   # kg/m^3
roho = 800  # kg/m^3

muo = .002     # Pa*s
mug = .0000131  # Pa*s

theta = 0  # 0: horizontal, -: downhill, +: uphill
sigma = .03  # N/m
epsillon = 0 


def get_lambdal_nfr(oil, gas, diameter, volumetric=False):
    '''
        Function to return lambdal and nfr
        Takes in oil and gas volumetric rate or rate 
        if volumetric needs diameter to finish calculation
    '''
    if not volumetric:
        return (oil / (oil + gas), (oil + gas) ** 2 / (9.81 * diameter))
    else:
        oil = oil / ((np.pi /4) * diameter ** 2) 
        gas = gas / ((np.pi /4) * diameter ** 2)
        return (oil / (oil + gas), oil + gas)


def det_flow_pattern(lambdal, nfr, theta):
    '''
        Function to return the flow pattern variables
        a,b,c,d,e,f,g from lambdal, vm and theta


        !!!!!there is probably a better way but here
        goes 50 elifs !!!!!!!!
    '''
    l1 = 316 * lambdal ** .302
    l2 = .0009252 * lambdal ** -2.4684
    l3 = .10 * lambdal ** -1.4516
    l4 = .5 * lambdal ** -6.738
    #  segregated flow
    if (lambdal < .01 and nfr < l1) or (lambdal >= .01 and nfr < l2):
        if theta >= 0:
            return 
        else:
            return 
    #  transition flow 
    elif lambdal >= .001 and l2 < nfr <= l3:
        if theta >= 0:
            return 
        else:
            return 
    #  interminent flow 
    elif not l4:
        return
    return
