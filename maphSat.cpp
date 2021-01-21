// The solver is largely based on the following paper: http://poincare.matf.bg.ac.rs/~filip//phd/sat-tutorial.pdf
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

#include "maphSat.hpp"

// Helper for random selection heuristics.
// Takes a lower and upper bound and returns a random index within the bounds.
// It's mainly used for vectors, therefore the upper bound is decremented.
int getRandomIndex(int lowerBound, int upperBound) {
    std::random_device randDevice;
    std::mt19937 rng(randDevice());
    std::uniform_int_distribution<int> result(lowerBound,upperBound-1);
    int randIndex = result(rng);
    return randIndex;
}

// Helper for computing the combined sum of occurrences (both polarities).
// If the first parameter is true, it computes the combined sum for the MOMS heuristic
// and uses the int to determine the 'cutoffLength' of a clause.
// If false, it computes the usual combined sum for DLCS.
void MaphSAT::combinedSum(std::vector<std::pair<int, int> > & pos, std::vector<std::pair<int, int> > & neg, bool constraint, std::size_t cutoffLength) const {
    int counterPos = 0;
    int counterNeg = 0;
    for (const auto & clause : formula){
        for (int literal : clause) {
            if (constraint && clause.size() > cutoffLength)
                break;
            if (std::find_if(trail.begin(), trail.end(), [literal](const auto & p) { return p.first == literal || p.first == -literal; }) == trail.end()) {
                auto itPos = std::find_if(pos.begin(), pos.end(), [&counterPos, literal](const auto & p) {
                    counterPos = p.second;
                    return p.first == literal || p.first == -literal;
                });
                auto itNeg = std::find_if(neg.begin(), neg.end(), [&counterNeg, literal](const auto & p) {
                    counterNeg = p.second;
                    return p.first == literal || p.first == -literal;
                });

                if (literal > 0) {
                    if (itPos == pos.end() && literal > 0)
                        pos.emplace_back(std::make_pair(literal, 1));
                    else {
                        pos.erase(itPos);
                        const int newPosCounter = ++counterPos;
                        pos.emplace_back(literal, newPosCounter);
                    }
                } else if (literal < 0) {
                    if (itNeg == neg.end())
                        neg.emplace_back(literal, 1);
                    else {
                        neg.erase(itNeg);
                        const int newNegCounter = ++counterNeg;
                        neg.emplace_back(literal, newNegCounter);
                    }
                }
            }
        }
    }
}

// Select the first literal that is not yet asserted.
int MaphSAT::selectFirst() const {
    for (const auto & clause : formula) {
        for (int literal : clause) {
            if (std::find_if(trail.begin(), trail.end(), [literal](const auto & p) { return p.first == literal || p.first == -literal; }) == trail.end())
                return literal;
        }
    }
    return 0;
}

// Selection heuristic: Pick a random literal.
int MaphSAT::selectRandom() const {
    int maxLit = 0;
    std::vector<int> randCandidates;

    for (const auto & clause : formula) {
        for (int literal : clause) {
            if (std::find_if(trail.begin(), trail.end(), [literal](const auto & p) { return p.first == literal || p.first == -literal; }) == trail.end())
                randCandidates.push_back(literal);
        }
    }
    if (!randCandidates.empty()) {
        const int randIndex = getRandomIndex(0, randCandidates.size());
        maxLit = randCandidates[randIndex];
    }
    if (maxLit > 0)
        return maxLit;
    else
        return -maxLit;
}

// Selection heuristic: Dynamic Largest Individual Sum.
// Picks the literal with the highest number of occurrences in the unsatisfied clauses.
// Sets value to true if the literal is positive.
// If the literal is negative, sets the value of its negation to true.
// If randomized true, it runs the randomized DLIS variant.
int MaphSAT::selectDLIS(bool random) const {
    int counter = 0;
    int maxNumber = 0;
    int maxLit = 0;
    std::vector<int> randCandidates;
    std::vector<std::pair<int, int> > vCount;

    for (const auto & clause : formula) {
        for (int literal : clause) {
            if (std::find_if(trail.begin(), trail.end(),
                [literal](const auto & p) { return p.first == literal || p.first == -literal; }) == trail.end()) {
                auto it = std::find_if(vCount.begin(), vCount.end(), [&counter, literal](const auto & p) {
                    counter = p.second;
                    return p.first == literal || p.first == -literal;
                });
                if (it == vCount.end()) {
                    vCount.emplace_back(literal, 1);
                    if (maxNumber > 1) {
                        maxNumber = 1;
                        maxLit = literal;
                    } else if (random && maxNumber == 1)
                        randCandidates.push_back(literal);
                } else {
                    vCount.erase(it);
                    const int newValue = ++counter;
                    vCount.emplace_back(literal, newValue);
                    if (newValue > maxNumber) {
                        maxNumber = newValue;
                        maxLit = literal;
                    } else if (random && newValue == maxNumber)
                        randCandidates.push_back(literal);
                }
            }
        }
    }

    if (random && !randCandidates.empty()) {
        const int randIndex = getRandomIndex(0, randCandidates.size());
        maxLit = randCandidates[randIndex];
    }
    if (maxLit > 0)
        return maxLit;
    else
        return -maxLit;
}

// Selection heuristic: Dynamic Largest Combined Sum.
// Picks the variable with the highest number of occurrences of its positive and negative literals (combined).
// If randomized true, it runs the randomized DLCS variant.
int MaphSAT::selectDLCS(bool random) const {
    int counterNeg = 0;
    int maxScore = 0;
    int maxLit = 0;
    std::vector<int> randCandidates;
    std::vector<std::pair<int, int> > posVariableCount;
    std::vector<std::pair<int, int> > negVariableCount;

    combinedSum(posVariableCount, negVariableCount, false, 0);

    for (auto & posLit : posVariableCount) {
        auto itNeg = std::find_if(negVariableCount.begin(), negVariableCount.end(), [&counterNeg, posLit](const auto & negLit) {
            counterNeg = negLit.second;
            return negLit.first == -posLit.first;
        });
        if (itNeg != negVariableCount.end()) {
            if (counterNeg + posLit.second > maxScore) {
                maxScore = counterNeg + posLit.second;
                maxLit = posLit.first;
                negVariableCount.erase(itNeg);
            } else if (random && (counterNeg + posLit.second) == maxScore) {
                randCandidates.push_back(posLit.first);
                negVariableCount.erase(itNeg);
            }
        }
    }

    if (maxLit == 0) {
        if (!posVariableCount.empty())
            maxLit = posVariableCount.front().first;
        else if (!negVariableCount.empty()) {
            maxLit = -negVariableCount.front().first;
        }
    }

    if (random && !randCandidates.empty()) {
        const int randIndex = getRandomIndex(0, randCandidates.size());
        maxLit = randCandidates[randIndex];
    }

    return maxLit;
}

// Selection heuristic: the Jeroslow-Wang method.
// If randomized true, it runs the randomized J-W variant.
int MaphSAT::selectJW(bool random) const {
    double score = 0;
    double maxScore = std::numeric_limits<double>::min();
    int maxLit = 0;
    std::vector<int> randCandidates;
    std::vector<std::pair<int, double> > JWcount;

    for (const auto & clause : formula) {
        for (int literal : clause) {
            if (std::find_if(trail.begin(), trail.end(), [literal](const auto & p) { return p.first == literal || p.first == -literal; }) == trail.end()) {
                auto it = std::find_if(JWcount.begin(), JWcount.end(), [&score, literal](const auto & p) {
                    score = p.second;
                    return p.first == literal || p.first == -literal;
                });
                if (it == JWcount.end()) {
                    double initialScore = std::pow(2.0, -clause.size());
                    JWcount.emplace_back(literal, initialScore);
                    if (initialScore > maxScore) {
                        maxScore = initialScore;
                        maxLit = literal;
                    } else if (random && initialScore == maxScore)
                        randCandidates.push_back(literal);
                } else {
                    JWcount.erase(it);
                    score = score + pow(2.0, -clause.size());
                    JWcount.emplace_back(literal, score);
                    if (score > maxScore) {
                        maxScore = score;
                        maxLit = literal;
                    } else if (random && score == maxScore)
                        randCandidates.push_back(literal);
                }
            }
        }
    }

    if (maxLit != 0) {
        if (random && !randCandidates.empty()) {
            const int randIndex = getRandomIndex(0, randCandidates.size());
            maxLit = randCandidates[randIndex];
        }
        if (maxLit > 0)
            return maxLit;
        else
            return -maxLit;
    }

    return 0;
}

// Selection heuristic: Maximum [number of] Occurrences in Minimum [length] Clauses.
// If randomized true, it runs the randomized MOMS variant.
int MaphSAT::selectMOMS(bool random) const {
    int maxLit = 0;
    int counterNeg = 0;
    int maxScore = 0;
    int parameter = 10; // as suggested in: J. Freeman, “Improvements to propositional satisfiability search algorithms” , PhD thesis, The University of Pennsylvania, 1995.
    int cutoffLength = 0;
    std::size_t totalClauseLength = 0;
    std::vector<int> randCandidates;
    std::vector<std::pair<int, int> > posVariableCount;
    std::vector<std::pair<int, int> > negVariableCount;

    for (const auto & clause : formula)
        totalClauseLength += clause.size();

    cutoffLength = totalClauseLength / formula.size();
    if (cutoffLength >= 2)
        --cutoffLength;

    combinedSum(posVariableCount, negVariableCount, true, cutoffLength);

    for (auto & posLit : posVariableCount) {
        auto itNeg = std::find_if(negVariableCount.begin(), negVariableCount.end(), [&counterNeg, posLit](const auto & negLit) {
            counterNeg = negLit.second;
            return negLit.first == -posLit.first;
        });
        if (itNeg != negVariableCount.end()) {
            const int tempScore = (counterNeg + posLit.second) * std::pow(2, parameter) + (counterNeg * posLit.second);
            if (tempScore > maxScore) {
                maxScore = tempScore;
                maxLit = posLit.first;
                negVariableCount.erase(itNeg);
            } else if (random && tempScore == maxScore) {
                randCandidates.push_back(posLit.first);
            }
        }
    }

    if (maxLit == 0)
        maxLit = selectFirst();
    if (random && !randCandidates.empty()) {
        const int randIndex = getRandomIndex(0, randCandidates.size());
        maxLit = randCandidates[randIndex];
    }

    return maxLit;
}

// Elimiate pure literals.
void MaphSAT::pureLiteral() {
    std::vector<int> trackLiterals; // keep track of literals occurring with unique polarity
    std::vector<int> erasedLiterals; // keep track of the erased literals
    for (const auto & clause : formula) {
        for (int literal : clause) {

            if (std::find_if(trail.begin(), trail.end(), [&](const auto & lit) {
                return lit.first == literal || lit.first == -literal;
            }) == trail.end()) {

                auto it = std::find_if(trackLiterals.begin(), trackLiterals.end(), [&](const auto & lit) {
                    return lit == -literal;
                });
                auto duplicate = std::find_if(trackLiterals.begin(), trackLiterals.end(), [&](const auto & lit) {
                    return lit == literal;
                });
                auto erased = std::find_if(erasedLiterals.begin(), erasedLiterals.end(), [&](const auto & lit) {
                    return lit == literal || lit == -literal;
                });
                if (it == trackLiterals.end() && duplicate == trackLiterals.end() && erased == erasedLiterals.end()) {
                    trackLiterals.push_back(literal);
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
        }
        trackLiterals.clear();
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
    case MaphSAT::Heuristic::FIRST:
        literal = selectFirst();
        break;
    case MaphSAT::Heuristic::RANDOM:
        literal = selectRandom();
        break;
    case MaphSAT::Heuristic::DLIS:
        literal = selectDLIS(false);
        break;
    case MaphSAT::Heuristic::RDLIS:
        literal = selectDLIS(true);
        break;
    case MaphSAT::Heuristic::DLCS:
        literal = selectDLCS(false);
        break;
    case MaphSAT::Heuristic::RDLCS:
        literal = selectDLCS(true);
        break;
    case MaphSAT::Heuristic::JW:
        literal = selectJW(false);
        break;
    case MaphSAT::Heuristic::RJW:
        literal = selectJW(true);
        break;
    case MaphSAT::Heuristic::MOMS:
        literal = selectMOMS(false);
        break;
    case MaphSAT::Heuristic::RMOMS:
        literal = selectMOMS(true);
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
    backjumpClause.erase(std::unique(backjumpClause.begin(), backjumpClause.end()), backjumpClause.end());
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
            backjumpClause.clear();
            for (int literal : clause)
                backjumpClause.push_back(literal);
        } else if (std::find(unitQueue.begin(), unitQueue.end(), clause[0]) == unitQueue.end()) {
            // If the first watched literal is not falsified, it is a unit literal.
            unitQueue.push_front(clause[0]);
            reason[clause[0]] = clauseIndex;
        }
    }

    watchList[literal] = newWL;
}

// Parse a CNF formula and throw invalid_argument() if unsuccessful.
MaphSAT::MaphSAT(std::istream & stream, MaphSAT::Heuristic heuristic) :
    heuristic(heuristic), state(MaphSAT::State::UNDEF), numberVariables(0),
    numberClauses(0), numberDecisions(0), conflict(false) {
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

    // Parse all clauses.
    int literal;
    std::vector<int> clause;
    for (std::size_t i = 0; i < numberClauses; ++i) {
        while (stream >> literal) {
            if (literal == 0 && clause.size() > 1) {
                // Add the clause to the formula.
                formula.push_back(clause);
                // Add the clause to the watch list.
                watchList[clause[0]].push_back(formula.size() - 1);
                watchList[clause[1]].push_back(formula.size() - 1);
                clause.clear();
                break;
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

    // Until the formula is satisfiable or unsatisfiable, the state of the solver is undefined.
    while (state == MaphSAT::State::UNDEF) {
        // Assert any unit literals.
        applyUnitPropagate();
        // Eliminate pure literals. It slowed out solver down so we uncommented it.
        // pureLiteral();
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
            if (trail.size() == numberVariables)
                state = MaphSAT::State::SAT;
            else
                applyDecide();
        }
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

// Print the current state of the SAT solver.
std::ostream & operator<<(std::ostream & out, const MaphSAT & maph) {
    switch (maph.state) {
    case MaphSAT::State::UNDEF:
        out << "UNDEF\n";
        break;
    case MaphSAT::State::SAT:
        out << "s SATISFIABLE\n";
        // to use the runTests.py script, change the output to "SAT"
        break;
    case MaphSAT::State::UNSAT:
        out << "s UNSATISFIABLE\n";
        // to use the runTests.py script, change the output to "UNSAT"
        break;
    }

    // to use the runTests.py script, uncomment the following output
    if (maph.state == MaphSAT::State::SAT) {
        out << "v ";
        for (const auto & literal : maph.trail)
            out << literal.first << ' ';
    }

    return out;
}
