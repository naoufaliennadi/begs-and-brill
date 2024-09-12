import numpy as np


class pgcalc():
    def __init__(self):
        self.G = 9.81  # m/s^2
        self.FP = {
            "segup" : (.98, .04846, .0868, .0110, -3.7680, 3.5390, -1.6140),
            "segdo" : (.98, .04846, .0868, 4.7000, -.3692, .1244, -.05056),
            "intup" : (.845, .5351, .0173, 2.9600, .3050, -.4473, .0978),
            "intdo" : (.845, .5351, .0173, 4.7000, -.3692, .1244, -.5056),
            "ditup" : (1.065, .05824, .0609, .00000001, 0, 0, 0),  # should make it negative so C=0
            "ditdo" : (1.065, .5824, .0609, 4.7000, -.3692, .1244, -.5056)
        }

    def switch_mode(self):
        self.G = 32.2  # ft/s^2


    def get_dimensionless(self, oil, gas, diameter, sigma, roho):
        '''
            Function to return lambdal, nfr and nvl
            Takes in oil/gas rate, diameter, sigma and oil density 
        '''
        return (oil / (oil + gas), (oil + gas) ** 2 / (self.G * diameter), oil * (roho / (self.G * sigma) ** 1 / 4))


    def get_mix_properties(self, lambdal, roho, rohg, muo, mug):
        '''
            Function to return mixture density and viscosity
            from gas and oil properties and lambdal
        '''
        rohn = roho * lambdal + rohg * (1 - lambdal) 
        mun = muo * lambdal + mug * (1 - lambdal)
        return rohn, mun


    def det_flow_pattern(self, lambdal, nfr, theta):
        '''
            Function to return the flow pattern
            from lambdal, nfr and theta
        '''
        l1 = 316 * lambdal ** .302
        l2 = .0009252 * lambdal ** -2.4684
        l3 = .10 * lambdal ** -1.4516
        l4 = .5 * lambdal ** -6.738
        #  segregated flow
        if (lambdal < .01 and nfr < l1) or (lambdal >= .01 and nfr < l2):
            if theta >= 0:
                return "segup"
            else:
                return "segdo"
        #  transition flow 
        elif lambdal >= .001 and l2 < nfr <= l3:
            if theta >= 0:
                return "tranup"
            else:
                return "trando"
        #  interminent flow 
        elif (.01 <= lambdal < .4 and l3 < nfr <= l1) or (lambdal > .4 and l3 < nfr <= l4):
            if theta >= 0:
                return "intup"
            else:
                return "intdo"
        #  distributed flow 
        elif (lambdal < .4 and nfr >= l1) or (lambdal >= .4 and nfr > l4):
            if theta >= 0:
                return "ditup"
            else:
                return "ditdo"
        else:
            raise "DIDNT FIND A FLOW PATTERN FROM DATA"


    def get_liquid_holdup(self, lambdal, nfr, nvl, theta, flop):
        """
            Function to get the liquid holdup from
            dimensionless variables and flow pattern

            !!!probably should find a better way to do
            transition flow!!!
        """
        tran = False
        if flop == "tranup":
            a,b,c,e,f,g,h = self.FP["segup"]
            a2,b2,c2,e2,f2,g2,h2 = self.FP["intup"]
            tran = True
        elif flop == "trando":
            a,b,c,e,f,g,h = self.FP["segdo"]
            a2,b2,c2,e2,f2,g2,h2 = self.FP["intdo"]
            tran = True
        else:
            a,b,c,e,f,g,h = self.FP[flop]

        if tran:
            l2 = .0009252 * lambdal ** -2.4684
            l3 = .10 * lambdal ** -1.4516
            trana = (l3 - nfr) / (l3 - l2)
            hlo01 = (a * lambdal ** b) / (nfr ** c) 
            hlo02 = (a2 * lambdal ** b2) / (nfr ** c2)

            if hlo01 < lambdal:
                hlo01 = lambdal
            if hlo02 < lambdal:
                hlo02 = lambdal

            cor1 = (1.0 - lambdal) * np.log(e * lambdal ** f * nvl ** g * nfr ** h) 
            cor2 = (1.0 - lambdal) * np.log(e2 * lambdal ** f2 * nvl ** g2 * nfr ** h2)

            if cor1 < 0:
                cor1 = 0
            if cor2 < 0:
                cor2 = 0

            tri1 = 1.0 + cor1 * (np.sin(1.8 * theta) - .333 * np.sin(1.8 * theta) ** 3)
            tri2 = 1.0 + cor2 * (np.sin(1.8 * theta) - .333 * np.sin(1.8 * theta) ** 3)
            hlot1 = hlo01 * tri1
            hlot2 = hlo02 * tri2
            if theta >= 0:
                hlot1 *= .924
                hlot2 *= .924
            else:
                hlot1 *= .685
                hlot2 *= .685
            hlot = trana * hlot1 + (1 - trana) * hlot2
        else:
            hlo0 = (a * lambdal ** b) / (nfr ** c)

            if hlo0 < lambdal:
                hlo0 = lambdal

            cor = (1.0 - lambdal) * np.log(e * lambdal ** f * nvl ** g * nfr ** h) 

            if cor < 0:
                cor = 0

            tri = 1.0 + cor * (np.sin(1.8 * theta) - .333 * np.sin(1.8 * theta) ** 3)
            hlot = hlo0 * tri
            if theta >= 0:
                hlot *= .924
            else:
                hlot *= .685
        return hlot


    def get_reynold_num(self, lambdal, rohn, mun, diameter, mr):
        '''
            Function to return reynold's number
        '''
        return (rohn * mr * diameter) / mun


    def get_friction_factor(self, ren, diameter, epsilon):
        '''
            Function to get the friction factor
            uses the c(Re)^-n method and the haaland 
            correlation
        '''
        if ren < 2000:
            return 64 / ren
        elif ren > 2000 and epsilon == 0:
            return .184 * ren ** -.2
        else:
            return (-1.8 * np.log10(((epsilon / diameter) / 3.7) ** 1.11 + 6.9 / ren)) ** -2


    def correct_friction_f(self, lambdal, hlo, ff):
        '''
            Function to correct the friction factor
            using the liquid holdup
        '''
        y = lambdal / (hlo ** 2)
        if 1 < y < 1.2:
            s = np.log(2.2 * y - 1.2)
        else:
            s = np.log(y) / (-.0523 + 3.182 * np.log(y) + -.8725 * np.log(y) ** 2 + .01853 * np.log(y) ** 4)
        return ff * np.exp(s)


    def get_rohs(self, hlo, roho, rohg):
        ''' Just for rohs '''
        return roho * hlo + rohg * (1 - hlo)


    def get_pressure_gradient(self, ff, rohn, mr, diameter, rohs, theta):
        ''' Just for pressure gradient '''
        return -(ff * rohn * mr ** 2) / (2 * diameter) - rohs * self.G * np.sin(theta)
