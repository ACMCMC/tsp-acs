from bayes_opt import BayesianOptimization, UtilityFunction
import subprocess
import multiprocess as mp

pool = mp.Pool(mp.cpu_count() - 1)

def process_problem(data):
  #print(f"Processing problem {data[0]}")
  with open("out_{}".format(data[0].split("/")[-1]), "w") as f:
    p = subprocess.Popen(["./implementations/tsp_run", data[0], "python_result.tmp", str(data[1]), str(data[2])], stdout = f)
    p.wait()
  #print(f"Finished processing {data[0]}")

def program(alpha, beta):
  # Create an empty file to store the result. If it exists, it will be overwritten.
  open("python_result.tmp", "w").close()
  #p = subprocess.Popen(["parallel", "./implementations/tsp_run", "{}", "python_result.tmp", str(alpha), str(beta), ":::", ], stdout = subprocess.PIPE)
  problems = ["./problems/ch130.tsp", "./problems/eil76.tsp", "./problems/fl1577.tsp", "./problems/pr439.tsp", "./problems/rat783.tsp", "./problems/u1060.tsp", "./problems/pcb442.tsp"]
  #problems = ["./problems/eil76.tsp"]
  #problems = ["./problems/u1060.tsp"]
  #processes = [subprocess.Popen(["./implementations/tsp_run", p, "python_result.tmp", str(alpha), str(beta)]) for p in problems]
  #[p.wait() for p in processes]
  pool.map(process_problem, [(p, alpha, beta) for p in problems])
  #print("Finished all problems.")
  with open("python_result.tmp", "r") as f:
    # The file contains several CSV lines with the return scores. Sum them up and get the average.
    values = [float(line.split(",")[2]) for line in f.readlines()]
  average = sum(values) / float(len(values))
  return 1.0 / average # Convert to a maximization problem.

# Set range of C to optimize for.
# bayes_opt requires this to be a dictionary.
pbounds = {
  "alpha": [0.1, 10],
  "beta": [0.1, 10],
}

# Create a BayesianOptimization optimizer,
# and optimize the given black_box_function.
optimizer = BayesianOptimization(f = program,
                                 pbounds = pbounds, verbose = 2,
                                 random_state = 1)

optimizer.maximize(init_points = 7, n_iter = 45)

# Delete the temporary file.
subprocess.call(["rm", "python_result.tmp"])

print("Best result: {}; f(x) = {}.".format(optimizer.max["params"], optimizer.max["target"]))