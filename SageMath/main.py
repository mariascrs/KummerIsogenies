from fp_arith import reset_counter, get_cost, update_p
from hash import KummerHash
from random import randint
from kummer_arithmetic import normalise


# ###########################################################################
#     Main file, used for benchmarking 
# ###########################################################################


def BenchmarkHash(secp, num_samples):
    """ Function to benchmark the hash function without optimal strategies using our cost metric """
    costs = []

    for i in range(num_samples):
        print(f"Running sample number {i+1}")
        HASH = KummerHash(secp, optimal=False)
        update_p(HASH.p)
        # We set the message to be the scalars defining the kernel (for now); need to fix the top bit of the scalars
        m1 = 2**HASH.l + randint(1,3**(HASH.e-1))
        m2 = 2**HASH.l + randint(1,3**(HASH.e-1))
        m3 = 2**HASH.l + randint(1,3**(HASH.e-1))
        msg = [m1,m2,m3]
        print(f'Message is {msg}.')
        
        reset_counter()
        h = HASH.hash(msg)
        print(f"Hash is {h}.\n")
        
        cost_run = get_cost()
        print(f'Total cost: {cost_run}')
        costs.append(cost_run)
        print("")
    

    cost_avg = sum(costs)/len(costs)
    print(f'Averge cost across {num_samples} runs is {cost_avg} using our cost metric.')

def BenchmarkHashOptimal(secp, num_samples):
    """ Function to benchmark the hash function with optimal strategies using our cost metric """
    
    costs = []
    print("Running hash with optimal strategies.")
    for i in range(num_samples):
        print(f"Running sample number {i+1}")
        HASH = KummerHash(secp)
        update_p(HASH.p)
        # We set the message to be the scalars defining the kernel (for now); need to fix the top bit of the scalars to be 1
        m1 = 2**HASH.l + randint(1,3**(HASH.e-1))
        m2 = 2**HASH.l + randint(1,3**(HASH.e-1))
        m3 = 2**HASH.l + randint(1,3**(HASH.e-1))
        msg = [m1,m2,m3]
        print(f'Message is {msg}.')
        
        reset_counter()
        h = HASH.hash_optimal(msg)
        print(f"Hash is {h}.\n")

        cost_run = get_cost()
        print(f'Total cost: {cost_run}')
        costs.append(cost_run)
        print("")


    cost_avg = sum(costs)/len(costs)
    print(f'Averge cost across runs is {cost_avg} using our cost metric.')
 
            
if __name__ == "__main__":
    secp = '128'
    # secp = '192'
    # secp = '256'
    num_samples = 100
    print(f"Security parameter is {secp}")
    # BenchmarkHash(secp, num_samples)
    BenchmarkHashOptimal(secp, num_samples)
    