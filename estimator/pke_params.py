from dataclasses import dataclass
import sys
sys.path.append('./lattice-estimator')

from estimator import LWE, ND
from estimator.nd import sigmaf,stddevf
from sage.all import var, log, ceil, floor, sqrt, Infinity, euler_phi, round, pi, n
import os

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# bytes pretty-printing
UNITS_MAPPING = [
    # (1<<53, ' PB'),
    # (1<<43, ' TB'),
    # (1<<33, ' GB'),
    # (1<<23, ' MB'),
    (1<<13, ' KB'),
    (1<<3, ' B'),
]

def pretty_size(bits, units=UNITS_MAPPING):
    """Get human-readable file sizes.
    simplified version of https://pypi.python.org/pypi/hurry.filesize/
    """
    for factor, suffix in units:
        if bits >= factor:
            break
    # amount = int(bits / factor)
    amount = n(bits/factor, digits=4)

    if isinstance(suffix, tuple):
        singular, multiple = suffix
        if amount == 1:
            suffix = singular
        else:
            suffix = multiple
    return f'{amount:5.1f}' + suffix

@dataclass
class PKEParams:
    secpar: int
    f: int # conductor
    n: int # module rank, height of matrix A 
    q: int # modulus
    ell: int # message dimension
    p: int # message modulus

    def __repr__(self):
        return f"PKE( secpar: {self.secpar:3d}, f: {self.f:3d}, phi: {self.phi():3d}, n: {self.n:1d}, q: 2^{ceil(log(self.q,2)):2d}, ell: {self.ell: 1d}, |ct|: {self.ct():6s} )"

    def phi(self):
        return euler_phi(self.f)

    def ct(self):
        return pretty_size((self.n + self.ell)* self.phi() * log(self.q,2))

def gen_pke_params(secpar, f, base_n = 1, ell = 1, p = 2, b = 2):
    chi = ND.CenteredBinomial(b)
    # chi = ND.UniformMod(2*b+1)
    phi = euler_phi(f)
    for n in range(base_n,5*base_n,1):
        q = 8 * phi * n + 5
        
        lwe_param = LWE.Parameters(
            n=phi*n, m = phi*(n + ell), q=q, Xs=chi, Xe=chi
        )
        
        lwe_param_2 = LWE.Parameters(
            n=phi*n, m = phi*n, q=q, Xs=chi, Xe=chi
        )
        
        lwe_params = [ # comment out LWE instance if estimation is not needed
            lwe_param,
            # lwe_param_2
        ]
        
        with HiddenPrints():
            lwe_costs = []
            for lwe in lwe_params:
                costs = LWE.estimate(lwe, deny_list=("arora-gb", "bkw", "bdd_hyrbid", "bdd_mitm_hybrid"))
                lwe_costs += [costs]
        
        min_lwe_costs = []
        for costs in lwe_costs:
            min_cost = log(min(cost["rop"] for cost in costs.values()),2)
            min_cost = floor(min_cost) if min_cost < Infinity else Infinity
            min_lwe_costs += [min_cost]
        
        min_lwe_cost = min(min_lwe_costs)    
                   
        print(f'secpar: {min_lwe_costs}, n: {n}, q: 2^{ceil(log(q,2)):2d}')
                   
        if (
            min_lwe_cost >= secpar 
        ):
            pke = PKEParams(min_lwe_cost, f, n, q, ell, p)
            print(f"{pke}")
            return pke
