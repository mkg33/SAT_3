# SAT_3

Project 3

Deadline: 15.03.2021

UPDATE (24.01)

As of today, our solver can solve almost all test cases well within the 1m CPU time limit. There's only one timeout (hole9.cnf). All tests performed with VSIDS.
There's also one borderline case that takes ca 57-58s but doesn't exceed the limit.
I'll post the full results soon. It's promising that VSIDS actually outperforms the well-worn FIRST strategy. Preprocessing (tautology, subsumption, duplicates) probably looks ugly but does what it says on the tin and doesn't hog the resources.

TODO:

* clause deletion strategy (P)
* VSIDS / VMTF branching heuristic (M)
* restarts + phase saving + different restart policies (P)
* preprocessing techniques (M)
* proof logging in DRUP (or DRAT) format
* tests
