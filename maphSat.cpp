// The solver is largely based on the following paper: http://poincare.matf.bg.ac.rs/~filip//phd/sat-tutorial.pdf
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

#include "maphSat.hpp"

// Elimiate pure literals.
void MaphSAT::pureLiteral() {
    std::vector<int> trackLiterals; // keep track of literals occurring with unique polarity
    std::vector<int> erasedLiterals; // keep track of the erased literals
    trackLiterals.reserve(numberVariables*2);
    erasedLiterals.reserve(numberVariables*2);

    for (const auto & clause : formula) {
        for (int literal : clause) {

            if (std::find_if(trail.begin(), trail.end(), [&](const auto & lit) {
                return lit.first == literal || lit.first == -literal;
            }) == trail.end()) {

                auto it = std::find_if(trackLiterals.begin(), trackLiterals.end(), [&](const auto & lit) {
                    return lit == -literal;
                });
                auto erased = std::find_if(erasedLiterals.begin(), erasedLiterals.end(), [&](const auto & lit) {
                    return lit == literal || lit == -literal;
                });
                if (it == trackLiterals.end() && erased == erasedLiterals.end()) {
                    trackLiterals.push_back(literal);
                    std::sort(trackLiterals.begin(), trackLiterals.end());
                    auto pairs = unique(trackLiterals.begin(), trackLiterals.end()); // this seems to be faster than erase(unique ...); tested w.r.t. CPU time
                    trackLiterals.erase(pairs, trackLiterals.end());
                }
                else if (it != trackLiterals.end()) {
                    erasedLiterals.push_back(literal);
                    trackLiterals.erase(it);
                }
            }
        }
    }
    if (!trackLiterals.empty()) {
        for (const auto & lit : trackLiterals) {
            assertLiteral(lit, true);
            #ifdef DEBUG
            std::cout << "Pure literal: "<< lit << '\n';
            #endif
        }
        trackLiterals.clear();
        erasedLiterals.clear();
    }
}

// EVSIDS branching heuristic
int MaphSAT::selectEVSIDS() {
    int maxLit = 0;

    for (size_t i = 0; i < VSIDSvec.size(); ++i) {
        auto it = std::find_if(trail.begin(), trail.end(), [&](const auto & lit) { //check that the literal is not in the trail
            return lit.first == VSIDSvec[i].first || lit.first == -VSIDSvec[i].first;
        });
        if (it == trail.end()) { //pick the literal with the highest score
            maxLit = VSIDSvec[i].first;
            break;
        }
    }

    #ifdef DEBUG
    std::cout << "\nmaxLit: " << maxLit << '\n';
    #endif

    return maxLit;
}

// VSIDS branching heuristic
int MaphSAT::selectVSIDS() {
    int maxLit = 0;

    for (size_t i = 0; i < VSIDSvec.size(); ++i) {
        auto it = std::find_if(trail.begin(), trail.end(), [&](const auto & lit) { //check that the literal is not in the trail
            return lit.first == VSIDSvec[i].first || lit.first == -VSIDSvec[i].first;
        });
        if (it == trail.end()) { //pick the literal with the highest score
            maxLit = VSIDSvec[i].first;
        }
        if (maxLit != 0 && VSIDSinterval != 256) {
            break;
        }
        else if (VSIDSinterval == 256) { // decay by 20% (initially)
            if (VSIDScounter >= 2000 && decay > 0) { // update the decay factor periodically; similar technique implemented in Glucose
                VSIDScounter = 0;
                decay -= 0.15; // decrement the decay factor
            }
           VSIDSvec[i].second *= decay;
       }
    }

    if (VSIDSinterval == 256) { // interval used in GRASP
        VSIDSinterval = 0;
    }

    #ifdef DEBUG
    std::cout << "\nmaxLit: " << maxLit << '\n';
    #endif

    return maxLit;
}

// Remove tautologies from the formula. Right now, it leads to contradiction (in rare cases).
// We should come back to it once we have proper clause deletion strategies.
void MaphSAT::removeTautologies() {
    for (const auto & clause : formula) {
        for (int literal : clause) {
            auto it = std::find_if(clause.begin(), clause.end(), [&](const auto & lit) { //check that the literal is not in the trail
                return lit == -literal;
            });
            if (it != clause.end()) {
                auto trailIt = std::find_if(trail.begin(), trail.end(), [&](const auto & lit){
                    return lit.first == -literal;
                });
                if (trailIt == trail.end()) {
                    assertLiteral(literal, true);
                    std::cout << "Asserted " << literal << '\n';
                    #ifdef DEBUG
                    std::cout << "Found tautology in clause: ";
                    for (int lit : clause) {
                        std::cout << lit << ' ';
                    }
                    std::cout << '\n';
                    #endif
                }
            }
        }
    }
}

// Assert a literal as a decision literal or as a non-decision literal.
void MaphSAT::assertLiteral(int literal, bool decision) {
    trail.emplace_back(literal, decision);
    notifyWatches(-literal);
}

// Select a literal that is not yet asserted and assert it as a decision literal.
void MaphSAT::applyDecide() {
    int literal = 0;

    switch (heuristic) {
    case MaphSAT::Heuristic::VSIDS:
        literal = selectVSIDS();
        break;
    case MaphSAT::Heuristic::EVSIDS:
        literal = selectEVSIDS();
        break;
    case MaphSAT::Heuristic::WEVSIDS:
        // this variant will be available within EVSIDS()
        break;
    }

    if (literal == 0)
        return;

    assertLiteral(literal, true);
    ++numberDecisions;
}

// If there are any unit literals due to the current partial evaluation, assert
// them as non-decision literals. Repeat until there are no more unit literals.
void MaphSAT::applyUnitPropagate() {
    while (!unitQueue.empty() && !conflict) {
        const auto literal = unitQueue.back();
        unitQueue.pop_back();
        assertLiteral(literal, false);
    }
}

// Returns the number of decision literals in the trail that precede the first
// occurrence of 'literal', including 'literal' itself if it is a decision literal.
std::size_t MaphSAT::level(int literal) const {
    std::size_t decisions = 0;
    for (const auto & lit : trail) {
        if (lit.second)
            ++decisions;
        if (lit.first == literal)
            break;
    }
    return decisions;
}

// Returns a literal from 'clause' that is in the trail such that no other
// literal from 'clause' comes after it in the trail.
int MaphSAT::lastAssertedLiteral() const {
    auto last = std::find_first_of(trail.rbegin(), trail.rend(), backjumpClause.begin(), backjumpClause.end(),
        [](const auto & lit1, const auto & lit2) { return lit1.first == -lit2; });
    if (last == trail.rend())
        return 0;
    return last->first;
}

// Check if the backjump clause satisfies the first UIP condition, which is the
// case if the backjump clause contains exactly one literal of the current decision level.
bool MaphSAT::isUIP() const {
    const int literal = lastAssertedLiteral();
    for (int lit : backjumpClause) {
        if (-lit != literal && level(-lit) == level(literal))
            return false;
    }
    return true;
}

// Perform a single resolution step between the backjump clause and a clause
// that is the reason for the propagation of -'literal'.
void MaphSAT::applyExplain(int literal) {
    // The index of the clause that forced the propagation of 'literal'.
    auto it = reason.find(literal);
    if (it == reason.end())
        return;

    // Add the remaining literals from the reason clause to the backjump clauses.
    for (int lit : formula[it->second]) {
        if (lit != literal)
            backjumpClause.push_back(lit);
    }

    backjumpClause.erase(std::remove(backjumpClause.begin(), backjumpClause.end(), literal), backjumpClause.end());
    backjumpClause.erase(std::remove(backjumpClause.begin(), backjumpClause.end(), -literal), backjumpClause.end());
    std::sort(backjumpClause.begin(), backjumpClause.end());
    auto pairs = unique(backjumpClause.begin(), backjumpClause.end());
    backjumpClause.erase(pairs, backjumpClause.end());
}

// Construct the backjump clause by repeatedly explaining a literal that lead to a
// conflict until the backjump clause satisfies the first UIP condition.
void MaphSAT::applyExplainUIP() {
    while (!isUIP())
        applyExplain(lastAssertedLiteral());
}

// Add a learned clause to the formula to prevent the same conflict from happening again.
void MaphSAT::applyLearn() {
    formula.push_back(backjumpClause);

    // Add the clause to the watch list.
    watchList[backjumpClause[0]].push_back(formula.size() - 1);
    watchList[backjumpClause[1]].push_back(formula.size() - 1);

    if (proofLogging) { // Add the clause to the proof log.
        proofClauses.push_back(backjumpClause);
    }
}

// Return an iterator to the first literal in the trail that has a decision level greater than 'level'.
std::vector<std::pair<int, bool> >::iterator MaphSAT::firstLiteralPast(int level) {
    int decisions = 0;
    for (auto it = trail.begin(); it != trail.end(); ++it) {
        if (it->second)
            ++decisions;
        if (decisions > level)
            return it;
    }
    return trail.end();
}

// Remove any literals from the trail that have a decision level greater than 'level'.
void MaphSAT::removePast(int level) {
    auto first = firstLiteralPast(level);
    trail.erase(first, trail.end());
    numberDecisions = level;
}

// Return the greatest decision level of the backjump clause exluding 'literal'.
int MaphSAT::getBackjumpLevel(int literal) {
    if ((literal != 0 && backjumpClause.size() > 1) ||
        (literal == 0 && backjumpClause.size() > 0)) {
        int maxLvl = 0;
        for (int lit : backjumpClause) {
            if (-lit != literal) {
                const int lvl = level(-lit);
                if (lvl > maxLvl)
                    maxLvl = lvl;
            }
        }
        return maxLvl;
    } else
        return 0;
}

// Backtrack literals from the trail until the backjump clause becomes a unit
// clause and then assert the unit literal.
void MaphSAT::applyBackjump() {
    const int literal = lastAssertedLiteral();
    const int level = getBackjumpLevel(literal);
    removePast(level);

    conflict = false;
    unitQueue.clear();
    unitQueue.push_front(-literal);
    reason[-literal] = formula.size() - 1;
    VSIDScounter++;
}

 // Notify clauses that a literal has been asserted.
void MaphSAT::notifyWatches(int literal) {
    // Are there any clauses that watch 'literal'?
    if (watchList.find(literal) == watchList.end())
        return;

    std::vector<std::size_t> newWL;
    newWL.reserve(watchList[literal].size());

    const auto list = watchList[literal];

    for (std::size_t clauseIndex : list) {

        auto & clause = formula[clauseIndex];
        // Swap the watched literals if the first watched literal was falsified.
        if (clause[0] == literal)
            std::swap(clause[0], clause[1]);

        // Is the clause already satisfied? Only check the first watched literal.
        if (std::find_if(trail.begin(), trail.end(), [&clause](const auto & p) { return p.first == clause[0]; }) != trail.end()) {
            newWL.push_back(clauseIndex);
            continue;
        }

        // Are there any other unfalsified literals in the clause?
        std::vector<int>::iterator other = clause.end();
        for (auto it = clause.begin() + 2; it != clause.end(); ++it) {
            if (std::find_if(trail.begin(), trail.end(), [it](const auto & p) { return p.first == -*it; }) == trail.end()) {
                other = it;
                break;
            }
        }
        // If there is, swap the unfalsified literal with the second watched literal.
        if (other != clause.end()) {
            std::iter_swap(clause.begin() + 1, other);
            watchList[clause[1]].push_back(clauseIndex);
            continue;
        }

        // If there is no other unfalsified literal and the first watched literal is
        // also falsified, then there is a conflict.
        newWL.push_back(clauseIndex);
        if (std::find_if(trail.begin(), trail.end(), [&clause](const auto & p) { return p.first == -clause[0]; }) != trail.end()) {
            conflict = true;
            ++numberConflicts;
            ++VSIDSinterval;
            backjumpClause.clear();
            for (const int literal : clause) {
                backjumpClause.push_back(literal);
            }
            for (size_t i = 0; i < backjumpClause.size(); ++i) {
                for (size_t j = 0; j < VSIDSvec.size(); ++j) {
                    if (VSIDSvec[j].first == backjumpClause[i]) {
                        VSIDSvec[j].second++; // we need a conditional here for VSIDS vs EVSIDS (if we want to keep both - I think it's good for benchmarking)

                        /*if (VSIDSvec[i].first == literal) {
                            VSIDSvec[i].second += pow(1/0.2, numberConflicts); // EVSIDS
                            break;
                        }*/

                        /*if (i == backjumpClause.size() - 1) { // weighted EVSIDS variant
                            VSIDSvec[j].second  += 0.5*pow(1/0.2, numberConflicts);
                        }
                        else {
                            VSIDSvec[j].second += (1-((i+1)/10000))*pow(1/0.2, numberConflicts);
                        }*/
                        break;
                    }
                }
            }
            std::sort(VSIDSvec.begin(), VSIDSvec.end());
            auto pairs = unique(VSIDSvec.begin(), VSIDSvec.end());
            VSIDSvec.erase(pairs, VSIDSvec.end());
            std::sort(VSIDSvec.begin(), VSIDSvec.end(), [](auto &pair1, auto &pair2) {
                return pair1.second > pair2.second; //sort by score in descending order
            });
        } else if (std::find(unitQueue.begin(), unitQueue.end(), clause[0]) == unitQueue.end()) {
            // If the first watched literal is not falsified, it is a unit literal.
            unitQueue.push_front(clause[0]);
            reason[clause[0]] = clauseIndex;
        }
    }

    watchList[literal] = newWL;
}

// A geometric policy used in MiniSat v1.14.
// Initial restart interval: 100 conflicts.
// Increase: factor of 1.5 after each restart.
// Set restartLimit to 100 in .hpp.
void MaphSAT::restartPolicy100() {
    if (numberConflicts >= restartLimit) {
        numberConflicts = 0;
        restartLimit *= 1.5;
        removePast(0);
    }
}

// Fixed-restart interval policy.
void MaphSAT::fixedRestart(int restartInterval) {
    if (numberConflicts >= restartInterval) {
        numberConflicts = 0;
        removePast(0);
    }
}

// From https://baldur.iti.kit.edu/sat/files/2018/l06.pdf.
// My own implementation is below (LubySequence).
// Need to measure which is faster.
unsigned int MaphSAT::Luby(unsigned int index) {
    for (unsigned int k = 1; k < 32; k++) {
        if (index == (1 << k) - 1) {
            return 1 << (k-1);
        }
    }
    for (unsigned int k = 1;;k++) {
        if ((1 << (k-1)) <= index && index < (1 << k) - 1) {
            return Luby(index-(1 << (k-1)) + 1);
        }
    }
}

/* Computes Luby's sequence: 1,1,2,1,1,2,4,1,1,2,1,1,2,4,8,...
   Recursive definition:
   if i=2ˆk - 1, then t_i=2ˆ(k-1),
   if 2ˆ(k-1) <= i < 2ˆk - 1, then t_i = t_{i-2ˆ(k-1)+1},
   where i is a positive integer.
   Takes an index (i).
   Returns the ith element of the sequence (t_i). */
int MaphSAT::LubySequence(int index) {
    int result = 0;
    double parameter = log2(index+1);
    /* Because:
       i = 2ˆparameter - 1
       2ˆparameter = i + 1
       parameter = log_2(i+1) */

    if (index == 1) {
        return 1;
    }
    if (floor(parameter) == ceil(parameter)) { // check if it's an integer
        result = pow(2, parameter-1);
    }
    else if (pow(2, floor(parameter)-1) == index) { // check if 2ˆ(floor(parameter)-1)) is an integer
            result = LubySequence(index - pow(2, floor(parameter)-1) + 1);
    }
    else { // it holds that 2ˆ(ceil(parameter-1)) is an integer
        result = LubySequence(index - pow(2, ceil(parameter)-1) + 1);
    }

    return result;
}

// Restart policy based on Luby's sequence.
// A unit interval is defined as 32 conflicts.
// The restart intervals are thus: 32, 32, 64, 32, 32, 64, 128...
// The choice of unit intervals based on:
// Huang, J. (2007). The effect of restarts on the efficiency of clause learning.
// Note: unit interval of 512 recommended in https://baldur.iti.kit.edu/sat/files/2018/l06.pdf.
// MiniSAT uses unit intervals of 100.
// Set restartLimit=32 in .hpp.
void MaphSAT::restartLuby() {
    if (numberConflicts >= restartLimit) {
        numberConflicts = 0;
        removePast(0);
        restartLimit = 8*LubySequence(++numberRestarts);
        numberConflicts = 0;
    }
}

// Parse a CNF formula and throw invalid_argument() if unsuccessful.
MaphSAT::MaphSAT(std::istream & stream, MaphSAT::Heuristic heuristic, bool proofLogging, std::string proofName) :
    heuristic(heuristic), state(MaphSAT::State::UNDEF), numberVariables(0),
    numberClauses(0), numberDecisions(0), conflict(false), proofLogging(proofLogging), proofName(proofName) {
    // Skip optional comments and the mandatory 'p cnf' appearing at the top of the CNF formula.
    char c;
    while (stream >> c) {
        if (c == 'c')
            stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        else if (c == 'p') {
            stream.ignore(4, '\0');
            break;
        }
    }

    stream >> std::skipws;
    // Parse the number of variables and the number of clauses.
    if (!(stream >> numberVariables))
        throw std::invalid_argument("Error parsing DIMACS.");
    if (!(stream >> numberClauses))
        throw std::invalid_argument("Error parsing DIMACS.");

    // Reserve memory for the clauses and the trail.
    formula.reserve(numberClauses);
    trail.reserve(numberVariables);
    VSIDSvec.reserve(numberVariables*2); //don't know if it makes sense to reserve more space than most likely needed (?)

    // Parse all clauses.
    int literal;
    std::vector<int> clause;
    for (std::size_t i = 0; i < numberClauses; ++i) {
        while (stream >> literal) {
            if (literal == 0 && clause.size() > 1) {
                std::sort(clause.begin(), clause.end());
                #ifdef DEBUG
                bool subsumed = false;
                #endif
                // experimental: I tried to remove duplicate clauses and (backward)-subsumed clauses at the same time
                auto it = std::find_if(formula.begin(), formula.end(), [&](const auto & formulaClause) {
                    if (clause == formulaClause) { // check for duplicates
                        return true;
                    }

                    for (int lit : clause) { // check for subsumption
                        if (std::find(formulaClause.begin(), formulaClause.end(), lit) == formulaClause.end()) {
                            return false;
                        }
                    }
                    #ifdef DEBUG
                    subsumed = true;
                    #endif

                    return true;
                });
                if (it == formula.end()) {
                    formula.push_back(clause);
                    watchList[clause[0]].push_back(formula.size() - 1);
                    watchList[clause[1]].push_back(formula.size() - 1);
                    for (int literal : clause) {
                        VSIDSvec.push_back(std::make_pair(literal, 0)); //initialize the VSIDS vector
                    }
                    clause.clear();
                    break;
                }
                else {
                    clause.clear(); // don't insert the clause if duplicate/subsumption case
                    break;

                    #ifdef DEBUG
                    if (subsumed) {
                        std::cout << "Found subsumed clause: ";
                    }
                    else {
                    std::cout << "Found duplicate: ";
                    }
                    for (int lit : clause) {
                        std::cout << lit << ' ';
                    }
                    std::cout << '\n';
                    #endif
                }
            } else if (literal == 0 && clause.size() == 1) {
                unitQueue.push_front(clause[0]);
                clause.clear();
                break;
            } else if (literal != 0 && std::find(clause.begin(), clause.end(), literal) == clause.end()) {
                clause.push_back(literal);
            }
        }
        if (stream.fail())
            throw std::invalid_argument("Error parsing DIMACS.");
    }
}

// Solve the CNF formula.
bool MaphSAT::solve() {

    // Are there any conflicts with the unit literals?
    for (std::size_t i = 0; i < unitQueue.size(); ++i) {
        const int unitLiteral = unitQueue[i];
        if (std::any_of(unitQueue.begin() + i, unitQueue.end(), [unitLiteral](int literal) { return literal == -unitLiteral; })) {
            state = MaphSAT::State::UNSAT;
        return false;
        }
    }

    //removeTautologies(); // only in the preprocessing stage; doesn't work for now

    applyUnitPropagate(); // before we can apply pureLiteral(), we have to call applyUnitPropagate() [cf. the case with unit.cnf if this line is uncommented]
    pureLiteral(); // eliminate pure literals only in the preprocessing stage [considerable speedup]

    // Until the formula is satisfiable or unsatisfiable, the state of the solver is undefined.
    while (state == MaphSAT::State::UNDEF) {

        // Assert any unit literals.
        applyUnitPropagate();
        // Do the current assignments lead to a conflict?
        if (conflict) {
            // Can we backtrack to resolve the conflict?
            if (numberDecisions == 0)
                state = MaphSAT::State::UNSAT;
            else {
                applyExplainUIP();
                applyLearn();
                applyBackjump();
            }
        } else {
            // Does every variable have an assignment? If that is the case, we are done.
            // Otherwise assign a value to a variable that has no assignment yet.
            if (trail.size() == numberVariables) {
                state = MaphSAT::State::SAT;

                #ifdef DEBUG
                if (verify()) {
                    std::cout << "\nVerifier: OK\n";
                }
                else {
                    std::cout << "\nVerifier: FAIL\n";
                }
                #endif
            }
            else {
                restartLuby();
                applyDecide();
            }
        }
    }

    if (proofLogging) { // write proof to file
        std::ofstream file;
        file.open(proofName + "_proofLog.txt");
        for (const auto & clause : proofClauses) {
            for (int lit : clause) {
                file << lit << ' ';
            }
            file << "0\n";
        }
        file.close();
    }

    // If the formula is satisfiable, the trail represents the satisfying assignment.
    // Sort the trail before it gets printed.
    if (state == MaphSAT::State::SAT) {
        std::sort(trail.begin(), trail.end(), [](const auto & l1, const auto & l2) {
            const int abs1 = std::abs(l1.first);
            const int abs2 = std::abs(l2.first);
            return abs1 < abs2;
        });
        return true;
    }

    return false;
}

// Verifier used to check satisfiable cases.
// Checks that every clause contains a literal that is true in the assignment trail.
// Returns true if the variable assignment is correct, false otherwise.
bool MaphSAT::verify() {
    bool clauseValid; // flag for checking the validity of every clause

    for (auto const & clause : formula) {
        clauseValid = false; // reset the flag for the next clause
        for (auto const literal : clause) {
            auto it = std::find_if(trail.begin(), trail.end(), [literal] (const auto & lit) {
                return literal == lit.first;
            });
            if (it != trail.end()) {
                clauseValid = true;
                break;
            }
        }
        if (!clauseValid){
            return false;
        }
    }
    return true;
}

// Print the current state of the SAT solver.
std::ostream & operator<<(std::ostream & out, const MaphSAT & maph) {
    switch (maph.state) {
    case MaphSAT::State::UNDEF:
        out << "UNDEF\n";
        break;
    case MaphSAT::State::SAT:
        //out << "s SATISFIABLE\n";
        out << "S";
        // to use the runTests.py script, change the output to "SAT"
        break;
    case MaphSAT::State::UNSAT:
        //out << "s UNSATISFIABLE\n";
        out << "U";
        // to use the runTests.py script, change the output to "UNSAT"
        break;
    }

    // to use the runTests.py script, uncomment the following output
    /*if (maph.state == MaphSAT::State::SAT) {
        out << "v ";
        for (const auto & literal : maph.trail)
            out << literal.first << ' ';
    }*/

    return out;
}
