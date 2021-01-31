# SAT_3

Project 3

Deadline: 15.03.2021 (2:00 PM)

**UPDATE (31.01)**

ALl tests pass well within the time limit (current version). Turns out that the VSIDS heuristic works just fine (provided there are restarts with the right parameters).

EVSIDS also works but the conflictsNumber counter needs to be reset with every
restart. I thought the counter was supposed to be incremented irrespective of
restarts. Anyway, I'm glad this version works. Now we can consolidate our updates
and get even better results. I suppose from now on it's just stress-free work :-).

**UPDATE (24.01)**

As of today, our solver can solve almost all test cases well within the 1m CPU time limit. There's only one timeout (hole9.cnf). All tests performed with VSIDS.
<del>There's also one borderline case that takes ca 57-58s but doesn't exceed the limit.
I'll post the full results soon.</del> It's promising that VSIDS actually outperforms the well-worn FIRST strategy. Preprocessing (<del>tautology</del>, subsumption, duplicates) probably looks ugly but does what it says on the tin and doesn't hog the resources.

TODO:

* clause deletion strategy (P)
* VSIDS / VMTF branching heuristic (M) :white_check_mark:
* phase saving (P)
* restarts + different restart policies (M) :white_check_mark:
* preprocessing techniques (M) :white_check_mark: [might add more]
* proof logging in DRUP (or DRAT) format (M) :white_check_mark:
* tests + cactus plots (M)
