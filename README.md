# BLDE
Differential Evolution for Bilevel Programming (BLDE)\
When using this code, please refer to 
"J. S. Angelo, E. Krempser and H. J. C. Barbosa, "Differential evolution for bilevel programming," 2013 IEEE Congress on Evolutionary Computation, Cancun, Mexico, 2013, pp. 470-477, doi: 10.1109/CEC.2013.6557606."

# Compile and create executable file
g++ -O2 blde.cpp funcoes* -o blde.exe

# An example to run the code
./blde.exe -genL 30 -genF 100 -popL 20 -popF 30 -CR 0.7 -F 0.9 -Var 1 -seed 2 -func 1

-- popL and popF: leader and follower population size\
-- genL and genF: leader and follower number of generations\
-- CR and F: DE's parameters\
-- seed: seed\
-- Var: DE's variant (**1**. DE/rand/1/bin, **2**. DE/best/1/bin, **3**.DE/target-to-rand/1/bin\
-- func: test functions (1001 to 1006: unconstrained SMD's test functions), ( 1 to 19: "simple" test problems)\
