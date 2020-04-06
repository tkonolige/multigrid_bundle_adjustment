using bamg
using BALUtils

ba = readbal(ARGS[1])
_, jacobian, _ = generate_problem(ba)


