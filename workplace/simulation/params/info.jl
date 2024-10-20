"""
["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]

"kaczmarz"        - kaczmarz method (the default)
"cgnr"            - CGNR
"direct"          - A direct solver using the backslash operator
"daxkaczmarz"     - Dax algorithm (with Kaczmarz) for unconstrained problems
"daxconstrained"  - Dax algorithm for constrained problems
"pseudoinverse"   - approximates a solution using the More-Penrose pseudo inverse
"fusedlasso"      - solver for the Fused-Lasso problem
"fista"           - Fast Iterative Shrinkage Thresholding Algorithm
"optista"         - "Optimal" ISTA
"pogm"            - Proximal Optimal Gradient Method
"admm"            - Alternating Direcion of Multipliers Method
"splitBregman"    - Split Bregman method for constrained & regularized inverse problems
"primaldualsolver"- First order primal dual method
"""

"""
admm         ["L2", "L1", "L21", "TV", "Positive", "Proj"]
cgnr         L2
fista        ["L2", "L1", "L21", "TV", "Positive", "Proj"]
optista      ["L2", "L1", "L21", "TV", "Positive", "Proj"] 
pogm         ["L2", "L1", "L21", "TV", "Positive", "Proj"] 
splitBregman ["L2", "L1", "L21", "TV", "Positive", "Proj"] 
"""