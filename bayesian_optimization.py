from bayes_opt import BayesianOptimization, UtilityFunction
import subprocess
import multiprocess as mp

pool = mp.Pool(mp.cpu_count() - 1)

alpha = 8.95666582
beta = 8.95666582
rho = 1.202
phi = 0.8452
tau0 = 0.3449


def process_problem(data):
    #print(f"Processing problem {data[0]}")
    problem_name = data["problem"].split("/")[-1]
    program_arg_results_file = "python_result_{}.tsp".format(problem_name) # The program will modify it...
    results_file = "python_result_{}.opt.tour".format(problem_name) # So it will be renamed to this.
    with open("out_{}".format(problem_name), "w") as f:
        p = subprocess.Popen(["./implementations/tsp_run", data["problem"], program_arg_results_file, str(
            data["alpha"]), str(data["beta"]), str(data["rho"]), str(data["phi"]), str(data["tau0"])], stdout=f)
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


def program(alpha, beta, rho, phi, tau0):
    # Create an empty file to store the result. If it exists, it will be overwritten.
    open("python_result.tmp", "w").close()
    #p = subprocess.Popen(["parallel", "./implementations/tsp_run", "{}", "python_result.tmp", str(alpha), str(beta), ":::", ], stdout = subprocess.PIPE)
    problems = ["./problems/ch130.tsp", "./problems/d198.tsp", "./problems/eil76.tsp", "./problems/fl1577.tsp", "./problems/kroA100.tsp",
                "./problems/lin318.tsp", "./problems/pcb442.tsp", "./problems/pr439.tsp", "./problems/rat783.tsp", "./problems/u1060.tsp"]
    #problems = ["./problems/eil76.tsp"]
    #problems = ["./problems/u1060.tsp"]
    #processes = [subprocess.Popen(["./implementations/tsp_run", p, "python_result.tmp", str(alpha), str(beta)]) for p in problems]
    #[p.wait() for p in processes]
    values = pool.map(process_problem, [{'problem': p, 'alpha': alpha, 'beta': beta,
             'rho': rho, 'phi': phi, 'tau0': tau0} for p in problems])
    #print("Finished all problems.")
    average = sum(values) / float(len(values))
    return 1.0 / average  # Convert to a maximization problem.


# Set range of C to optimize for.
# bayes_opt requires this to be a dictionary.
pbounds = {
    "alpha": [0, 10],
    "beta": [0, 10],
#}
#pbounds = {
    "rho": [1e-10, 0.1],
    "phi": [1e-10, 0.1],
    "tau0": [0, 1.5],
}

# Create a BayesianOptimization optimizer,
# and optimize the given black_box_function.
optimizer = BayesianOptimization(f=program,
                                 pbounds=pbounds, verbose=2,
                                 random_state=1)

optimizer.maximize(init_points=7, n_iter=90)


print("Best result: {}; f(x) = {}.".format(
    optimizer.max["params"], optimizer.max["target"]))
