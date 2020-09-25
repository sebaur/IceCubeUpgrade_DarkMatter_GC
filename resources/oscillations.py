import numpy as np

def PMNS_factory(t12, t13, t23, d):
    s12 = np.sin(t12)
    c12 = np.cos(t12)
    s23 = np.sin(t23)
    c23 = np.cos(t23)
    s13 = np.sin(t13)
    c13 = np.cos(t13)
    cp  = np.exp(1j*d)
    return np.array([[ c12*c13, s12*c13, s13*np.conj(cp) ],
                  [-s12*c23 - c12*s23*s13*cp, c12*c23 - s12*s23*s13*cp, s23*c13],
                  [ s12*s23 - c12*s23*s13*cp,-c12*s23 - s12*c23*s13*cp, c23*c13]])

class infOsz():

    ##Probability of flavor to change when L->inf
    def prob(self,a, b, U):
        s = 0
        for i in range(3):
                s += (np.conj(U[a,i])*U[b,i]*U[a,i]*np.conj(U[b,i])).real
        return s

    #Define Oscillation Matrix
    def prob_matrix(self,U):
        return np.array([[self.prob(0, 0, U), self.prob(0, 1, U), self.prob(0,2,U)],
                     [self.prob(1, 0, U), self.prob(1, 1, U), self.prob(1,2,U)],
                     [self.prob(2, 0, U), self.prob(2, 1, U), self.prob(2,2,U)]])


    def __init__(self):
        self.t12 = 0.5934
        self.t13 = 0.1514
        self.t23 = 0.785398

        self.matrix = self.prob_matrix(PMNS_factory(self.t12, self.t13, self.t23, 0))


    def oscillate(self,dNdE_e,dNdE_ebar,dNdE_mu,dNdE_mubar,dNdE_tau,dNdE_taubar):
        dNdE_nu_e_osc = []
        dNdE_nu_mu_osc = []
        dNdE_nu_tau_osc = []
        for i in range(len(dNdE_e)):
            dNdE_nu_e_osc.append(np.dot(self.matrix,np.array([dNdE_e[i] + dNdE_ebar[i],dNdE_mu[i] + dNdE_mubar[i],dNdE_tau[i] + dNdE_taubar[i]]))[0])
            dNdE_nu_mu_osc.append(np.dot(self.matrix,np.array([dNdE_e[i] + dNdE_ebar[i],dNdE_mu[i] + dNdE_mubar[i],dNdE_tau[i] + dNdE_taubar[i]]))[1])
            dNdE_nu_tau_osc.append(np.dot(self.matrix,np.array([dNdE_e[i] + dNdE_ebar[i],dNdE_mu[i] + dNdE_mubar[i],dNdE_tau[i] + dNdE_taubar[i]]))[2])

        dNdE_nu_e_osc = np.array(dNdE_nu_e_osc)
        dNdE_nu_mu_osc = np.array(dNdE_nu_mu_osc)
        dNdE_nu_tau_osc = np.array(dNdE_nu_tau_osc)

        return dNdE_nu_e_osc, dNdE_nu_mu_osc, dNdE_nu_tau_osc



class threeFlavorVacuumOsc():


    def __init__(self):
        self.t12 = 0.5934
        self.t13 = 0.1514
        self.t23 = 0.785398

        self.dm2 = [ 7.58E-05, 2.27E-03, 2.35E-03]

        self.U = PMNS_factory(self.t12, self.t13, self.t23, 0)


    def posc(self, a, b, L, E):
        """
        Gives the oscillation probability for nu(a) -> nu(b)
        for PMNS matrix U, and L in km and E in GeV
        """
        s = 0
        for i in range(2):
            for j in range(i+1, 3):
                arg = 5.068*self.dm2[i+j-1]*L/E
                mxe = np.conj(self.U[a,i])*self.U[b,i]*self.U[a,j]*np.conj(self.U[b,j])
                s += -4*mxe.real*np.sin(0.25*arg)**2 + 2*mxe.imag*np.sin(0.50*arg)**2
        if a == b: s += 1.0
        return s
