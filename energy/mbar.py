import numpy as np
nlambdas = 31
niterations = 100


def run(epsilon_r):
        
        b = np.loadtxt('u_kln.txt')
        u_kln = b.reshape((nlambdas,nlambdas,niterations))
        nstates = len(u_kln)
        print("u_kln shape: ",u_kln.shape)

        # Estimate free energy of Lennard-Jones particle insertion
        from pymbar import MBAR, timeseries

        ## Subsample data to extract uncorrelated equilibrium timeseries
        print('The timeseries part')
        N_k = np.zeros([nstates], np.int32) # number of uncorrelated samples
        for k in range(nstates):
                [nequil, g, Neff_max] = timeseries.detect_equilibration(u_kln[k,k,:])
                indices = timeseries.subsample_correlated_data(u_kln[k,k,:], g=g)
                N_k[k] = len(indices)
                u_kln[k,:,0:N_k[k]] = u_kln[k,:,indices].T

        # Compute free energy differences
        print('The mbar part')
        mbar = MBAR(u_kln, N_k)

        # dont compute uncertainties here, if you do it may fail with an error for
        # pymbar versions > 3.0.3. See this issue: https://github.com/choderalab/pymbar/issues/419
        dict = mbar.compute_free_energy_differences(compute_uncertainty=True)
        DeltaF_ij = dict['Delta_f']
        Error = dict['dDelta_f']

        print("Free energy change to insert a particle = ",DeltaF_ij[nstates-1][0])
        print("Error = ",Error[nstates-1][0])


run(15)
