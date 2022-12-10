from bayes_opt import BayesianOptimization, UtilityFunction
import subprocess
import os

def program(alpha, beta):
  # Create an empty file to store the result. If it exists, it will be overwritten.
  open("python_result.tmp", "w").close()
  p = subprocess.Popen(["implementations/tsp_run", "./problems/lin318.tsp", "python_result.tmp", str(alpha), str(beta)], stdout = subprocess.PIPE)
  p.wait()
  with open("python_result.tmp", "r") as f:
    # The file contains a CSV line with the return score.
    return - float(f.readline().split(",")[1]) # Convert to a maximization problem.

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
                                 random_state = 4)

optimizer.maximize(init_points = 5, n_iter = 20)

# Delete the temporary file.
subprocess.call(["rm", "python_result.tmp"])

print("Best result: {}; f(x) = {}.".format(optimizer.max["params"], optimizer.max["target"]))