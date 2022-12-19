from bayes_opt import BayesianOptimization, UtilityFunction
import subprocess
import multiprocess as mp

pool = mp.Pool(mp.cpu_count() - 1)

alpha = 1.0
beta = 1.0
rho = 0.05
phi = 0.05
tau0 = 1
nAnts = 15
minPheromone = 10e-6

problems_full = ["./problems/ch130.tsp", "./problems/d198.tsp", "./problems/eil76.tsp", "./problems/fl1577.tsp", "./problems/kroA100.tsp",
                "./problems/lin318.tsp", "./problems/pcb442.tsp", "./problems/pr439.tsp", "./problems/rat783.tsp", "./problems/u1060.tsp"]

problems_reduced = ["./problems/ch130.tsp", "./problems/d198.tsp", "./problems/eil76.tsp", "./problems/fl1577.tsp", "./problems/pcb442.tsp", "./problems/rat783.tsp", "./problems/u1060.tsp"]

problems = problems_full

def process_problem(data):
    #print(f"Processing problem {data[0]}")
    problem_name = data["problem"].split("/")[-1]
    program_arg_results_file = "python_result_{}.tsp".format(problem_name) # The program will modify it...
    results_file = "python_result_{}.opt.tour".format(problem_name) # So it will be renamed to this.
    with open("out_{}".format(problem_name), "w") as f:
        p = subprocess.Popen(["./implementations/tsp_run", data["problem"], program_arg_results_file, str(
            data["alpha"]), str(data["beta"]), str(data["rho"]), str(data["phi"]), str(data["tau0"]), str(data["nAnts"]), str(data["minPheromone"])], stdout=f)
        p.wait()
    #print(f"Finished processing {data[0]}")
    return_val = 0.0
    with open(results_file, "r") as f:
        # Scan the file until we get to a line of the form "TOUR_LENGTH : <value>".
        for line in f:
            if line.startswith("ERROR"):
                # Return the value.
                return_val = float(line.split(":")[1])
                break
    subprocess.call(["rm", results_file])
    return return_val


def program(alpha, tau0, minPheromone):
    # Create an empty file to store the result. If it exists, it will be overwritten.
    open("python_result.tmp", "w").close()
    #p = subprocess.Popen(["parallel", "./implementations/tsp_run", "{}", "python_result.tmp", str(alpha), str(beta), ":::", ], stdout = subprocess.PIPE)
    #problems = ["./problems/eil76.tsp"]
    #problems = ["./problems/u1060.tsp"]
    #processes = [subprocess.Popen(["./implementations/tsp_run", p, "python_result.tmp", str(alpha), str(beta)]) for p in problems]
    #[p.wait() for p in processes]
    values = pool.map(process_problem, [{'problem': p, 'alpha': alpha, 'beta': beta,
             'rho': rho, 'phi': phi, 'tau0': tau0, 'nAnts': round(nAnts), 'minPheromone': minPheromone} for p in problems])
    #print("Finished all problems.")
    average = sum(values) / float(len(values))
    return 1.0 / average  # Convert to a maximization problem.


# Set range of C to optimize for.
# bayes_opt requires this to be a dictionary.
pbounds = {
    "nAnts": [5, 50],
    "beta": [10e-10, 5],
    "rho": [0.01, 0.2],
    "phi": [0.01, 0.2],
}
pbounds = {
    "alpha": [0.8, 2],
    "tau0": [1, 3],
    "minPheromone": [10e-7, 10e-4],
}

# Create a BayesianOptimization optimizer,
# and optimize the given black_box_function.
optimizer = BayesianOptimization(f=program,
                                 pbounds=pbounds, verbose=2,
                                 random_state=1)

optimizer.maximize(init_points=10, n_iter=9000)


print("Best result: {}; f(x) = {}.".format(
    optimizer.max["params"], optimizer.max["target"]))
