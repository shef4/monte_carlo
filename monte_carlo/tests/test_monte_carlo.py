"""
Unit and regression test for the monte_carlo package.
"""

# Import package, test suite, and other packages as needed
import monte_carlo
import pytest
import sys
import random
import numpy as np
import copy as cp
import matplotlib as mpl
import monte_carlo
import networkx as nx

def test_monte_carlo_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "monte_carlo" in sys.modules

def build_1d_graph(N, Jval):
    """
    Build a 1D graph with a single J value (Jval)
    """
    G = nx.Graph()
    G.add_nodes_from([i for i in range(N)])
    G.add_edges_from([(i,(i+1)% G.number_of_nodes() ) for i in range(N)])
    # G.add_edge(2,5)
    # G.add_edge(4,8)
    # G.add_edge(4,0)
    for e in G.edges:
        G.edges[e]['weight'] = Jval
    return G


def get_IsingHamiltonian(G, mus=None):
    if mus == None:
        mus = np.zeros(len(G.nodes()))

    if len(G.nodes()) != len(mus):
        print("error: Dimension Mismatch")
        return None

    if len(G.nodes()) != len(mus):
        print("error: Dimension Mismatch")
        return None
    
    J = [[] for i in G.nodes()]
    for e in G.edges:
        J[e[0]].append((e[1], G.edges[e]['weight']))
        J[e[1]].append((e[0], G.edges[e]['weight']))
    return monte_carlo.IsingHamiltonian(J,mus)

def test_ising():
    random.seed(2)
    #testcase 1
    N = 8
    Jval = 1
    G = build_1d_graph(N, Jval)
    ham = get_IsingHamiltonian(G, mus=[.1 for i in range(N)])
    conf = monte_carlo.BitString(N=N)
    conf.set_config([0, 0, 0, 0, 0, 0, 1, 1])

    Ei = ham.energy(conf)
    assert(abs(Ei-3.6) < 1e-12)

    conf.set_int_config(106)
    Ei = ham.energy(conf)
    assert(abs(Ei+4.0) < 1e-12)

    #testcase 2
    N = 6
    conf = monte_carlo.BitString(N=N)
    G = build_1d_graph(N, 2)
    ham = get_IsingHamiltonian(G, mus=[1.1 for i in range(N)])
    E, M, HC, MS = ham.compute_average_values(conf, 1)
    assert(np.isclose(E,  -11.90432015))
    assert(np.isclose(M,  -0.02660820))
    assert(np.isclose(HC, 0.59026994))
    assert(np.isclose(MS, 0.05404295))

    #testcase 3
    N = 8
    e_list = []
    e2_list = []
    m_list = []
    m2_list = []
    T_list = []
    conf = monte_carlo.BitString(N=N)
    G = build_1d_graph(N, 1)
    ham = get_IsingHamiltonian(G, mus=[.1 for i in range(N)])

    for Ti in range(1,100):
        T = .1*Ti
        
        E, M, HC, MS = ham.compute_average_values(conf, T)
        
        e_list.append(E)
        m_list.append(M)
        e2_list.append(HC)
        m2_list.append(MS)
        T_list.append(T)

    Tc_ind = np.argmax(m2_list)
    Tc2 = T_list[np.argmax(e2_list)]
    assert(np.isclose(T_list[Tc_ind],  2.00000000))
    assert(np.isclose(e_list[Tc_ind],  -3.73231850))
    assert(np.isclose(m_list[Tc_ind],  -0.14658168))
    assert(np.isclose(e2_list[Tc_ind], 1.64589165))
    assert(np.isclose(m2_list[Tc_ind], 1.46663062))
    assert(np.isclose(Tc2,  1.00000000))

    #testcase 4
    conf = monte_carlo.BitString(N=N)
    conf.initialize(M=4)
    T = 2
    E, M, EE, MM = monte_carlo.metropolis_monte_carlo(ham, conf, T=T, nsweep=80000, nburn=1000)
    HC = (EE[-1] - E[-1]*E[-1])/T/T
    MS = (MM[-1] - M[-1]*M[-1])/T
    assert(np.isclose(E[-1], -3.70852500))
    assert(np.isclose(M[-1], -0.14625000))
    assert(np.isclose(HC, 1.63372658))
    assert(np.isclose(MS, 1.48850546))

def test_average_values():
    N=10
    conf = monte_carlo.BitString(N=N)
    #conf.initialize(M=20)

    T = 2.0

    
    # now test the general ising hamiltonian
    Jval = 1.0
    mu = [.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = monte_carlo.IsingHamiltonian(J=J, mu=mu)
    E, M, HC, MS = ham2.compute_average_values(conf, T) 
 
    assert(np.isclose(E, -4.6378514858094695))
    assert(np.isclose(M, -0.1838233606011354 ))
    assert(np.isclose(HC, 1.9883833749653714 ))
    assert(np.isclose(MS, 1.8391722085614428))


def test_metropolis():
    random.seed(2)
    N=20
    conf = monte_carlo.BitString(N=N)
    T = 2

    J = []
    Jval = 1.0
    mu = [.1 for i in range(N)]
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham = monte_carlo.IsingHamiltonian(J=J, mu=mu)

    conf = monte_carlo.BitString(N=N)
    E, M, EE, MM = monte_carlo.metropolis_monte_carlo(ham, conf, T=T, nsweep=8000, nburn=1000)
    HC = (EE[-1] - E[-1]*E[-1])/T/T
    MS = (MM[-1] - M[-1]*M[-1])/T
    print
    print("     E:  %12.8f" %(E[-1]))
    print("     M:  %12.8f" %(M[-1]))
    print("     HC: %12.8f" %(HC))
    print("     MS: %12.8f" %(MS))

    assert(np.isclose(-9.3096250, E[-1]))
     

def test_classes():
    random.seed(2)
    conf = monte_carlo.BitString(N=10)
    conf.initialize(M=5)
    assert(all(conf.config == [1, 1, 1, 0, 0, 0, 0, 1, 1, 0]))

    N = 10
    J = []
    Jval = -1.0
    mu = [-.001 for i in range(N)]
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = monte_carlo.IsingHamiltonian(J=J, mu=mu)

    e = ham2.energy(conf)
    print(" Energy = ", e)
    assert(np.isclose(e,-2))

    conf.flip_site(3)
    print(conf.config)
    print(" Energy = ", e)
    assert(np.isclose(ham2.energy(conf),-2.002))
    
    # now flip back
    conf.flip_site(3)
    print(conf.config)
    e = ham2.energy(conf)
    print(" Energy = ", e)
    assert(np.isclose(ham2.energy(conf),-2.00))


    conf_old = cp.deepcopy(conf)
   
    ham2.mu = np.array([-1.1 for i in range(N)])
    ham2.metropolis_sweep(conf, T=.9)
    print(conf_old, " --> ", conf)
    print("Energy: %12.8f --> %12.8f" %(e, ham2.energy(conf)))  
    assert(all(conf.config == np.ones(10)))
    

    random.seed(2)
    conf.set_int_config(44)
    conf_old = cp.deepcopy(conf)

    
    Jval = 1.0
    mu = [-.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = monte_carlo.IsingHamiltonian(J=J, mu=mu)
    random.seed(2)
    conf.set_int_config(44)
    conf_old = cp.deepcopy(conf)
    ham2.metropolis_sweep(conf, T=.9)
    print(conf_old, " --> ", conf)
    print(conf_old)
    print(ham2.energy(conf_old))  
    print("Energy: %12.8f --> %12.8f" %(e, ham2.energy(conf)))  
    assert(all(conf.config == [1,1,1,0,1,0,0,1,1,0]))
   
def test_delta_e():
    random.seed(2)
    N = 20
    conf = monte_carlo.BitString(N=N)
    conf.initialize(M=10)

    Jval = -1.0

    mu = [-.001 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])


    ham = monte_carlo.IsingHamiltonian(J=J, mu=mu)

    e1 = ham.energy(conf)
    print(conf)
    print(" Energy = ", e1)

    delta_e1 = ham.delta_e_for_flip(3, conf)

    conf.flip_site(3)
    e2 = ham.energy(conf)
    print(conf)
    print(" Energy = ", e2)

    print(" delta E: %12.8f" %(e2-e1))
    print(" delta E: %12.8f" %(delta_e1))

    assert(np.isclose(e2-e1, delta_e1))
    
if __name__== "__main__":
    test_monte_carlo_imported()
    test_ising()
    test_classes()
    test_average_values()
    test_metropolis()
    test_delta_e()

