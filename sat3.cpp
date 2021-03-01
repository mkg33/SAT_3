#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

using Literal =  int32_t;   // Literal data type.
using   CSize = uint32_t;   // Small size data type.
using   Size  = uint64_t;   // Large size data type.

// Comparator based on activity level of a literal.
struct ActivityCompare {
    bool operator()(const std::pair<Literal, Size> & p1, const std::pair<Literal, Size> & p2) const {
        return p1.second > p2.second || (p1.second == p2.second && p1.first < p2.first);
    }
};

struct Clause {
    CSize   activity;       // Activity level of the clause.
    CSize       size;       // Size of the clause.
    Size       index;       // Start of the literals of the clause.
    Literal cache[4];       // First four literals of the clause.

    template <typename... T>
    Clause(CSize activity, CSize size, Size index, T... cache) :
        activity(activity), size(size), index(index), cache{cache...} {}
    
    Size begin() const { return index; }
    Size   end() const { return index + size; }
};

struct Variable {
    Size      level;     // Decision level of the variable. Meaningless if 'trail' is false.
    bool      trail;     // Is the variable in the trail?
    bool       sign;     // Is the variable positive or negative? Meaningless if 'trail' is false.
    bool   decision;     // Is the variable a decision variable? Meaningless if 'trail' is false.

    Variable(Size level, bool trail, bool sign, bool decision) :
        level(level), trail(trail), sign(sign), decision(decision) {}
};

// Return the index of the variable of the literal.
Size getIndex(Literal literal) {
    return std::abs(literal) - 1;
}

// Place literals over the same variable close to each other.
bool litAbsLessThan(Literal a, Literal b) {
    Literal s = std::abs(a);
    Literal t = std::abs(b);
    return s < t || (s == t && a < b);
}

class SAT {
public:

    SAT(std::istream &);    // Constructor parses a CNF file.
    bool solve();           // Solve the SAT problem.

private:

    enum class State {
        UNDEF,              // No solution yet.
          SAT,              // Formula is satisfiable.
        UNSAT               // Formula is unsatisfiable.
    } state;

    struct Statistics {
        Size   varNum;   // Number of variables.
        Size   clsFrm;   // Number of clauses of the formula.
        Size   clsLrn;   // Number of learned clauses.
        Size   clsTrs;   // Threshold on the length of learned clauses.
        Size   decNum;   // Number of decisions.
        Size   conNum;   // Number of conflicts.
        Size   conRes;   // Number of conflicts since the last restart.
        Size   resNum;   // Number of restarts.
        Size   resTrs;   // Do not restart if 'conRes' is below this threshold.
        Size   litHgh;   // Number of literals of the backjump clause of the highest decision level.
        double lrnTrs;   // The limit of learned clauses is a factor of the size of the formula.
        double lrnFac;   // 'lrnTrs' is multiplied with this factor each restart.

        Statistics(Size varNum, Size clsFrm, Size clsLrn, Size clsTrs, Size decNum, Size conNum,
            Size conRes, Size resNum, Size resTrs, Size litHgh, double lrnTrs, double lrnFac) :
            varNum(varNum), clsFrm(clsFrm), clsLrn(clsLrn), clsTrs(clsTrs), decNum(decNum),
            conNum(conNum), conRes(conRes), resNum(resNum), resTrs(resTrs), litHgh(litHgh),
            lrnTrs(lrnTrs), lrnFac(lrnFac) {}
    } statistics;

    bool    conflict;    // Conflict flag.
    Literal lastAss;     // Last asserted literal of the backjump clause.

    std::vector<Literal>                              literals;   // All literals.
    std::vector<Variable>                            variables;   // All variables.
    std::vector<Clause>                                formula;   // All clauses.
    std::set<std::pair<Literal, Size>, ActivityCompare> scores;   // Sorted list of literals based on activity level.
    std::vector<Literal>                                 trail;   // Current partial assignment.
    std::vector<Literal>                              backjump;   // Backjump clause.
    std::vector<Literal>                           backjumpLow;   // Literals of the backjump clause that are not of the current decision level.
    std::unordered_map<Literal, bool>              backjumpMap;   // Is a literal in the backjump clause?
    std::deque<Literal>                              unitQueue;   // Queue of literals to be propagated.
    std::unordered_map<Literal, std::vector<Size>>       watch;   // Watch lists.
    std::unordered_map<Literal, Size>                   reason;   // What is the clause that propagated a literal?

    void addClause(std::vector<Literal> &);                       // Add a clause to the formula.
    void assertLiteral(Literal, bool);                            // Add a literal as a decision or non-decision literal to the trail.
    void assertPureLiterals();                                    // Check if there are any pure literals and if so, assert them.
    void backjumpApply();                                         // Analyze the conflict and execute the backjump.
    void backjumpInit();                                          // Initialize the conflict analysis.
    void decide();                                                // Decide on the next decision literal.
    void findLastAsserted();                                      // Find the last asserted literal of the backjump clause.
    void forget();                                                // Forget unused learned clauses.
    void notify(Literal);                                         // Notify watch lists that a literal was asserted.
    void prefixTo(Size);                                          // Remove any literals from the trail past a certain level.
    void resolve(Literal, const Clause &);                        // Explain a literal.
    void restart();                                               // Restart using Luby's sequence.
    bool satisfied(Literal);                                      // Check if a literal is satisfied by the current partial assignment.
    bool shouldKeep(std::vector<Literal> &);                      // Check if a clause should be kept in the formula.
    void subsumeApply();                                          // Check if a clause subsumes the backjump clause.
    void unitPropagate();                                         // Propagate unit literals.

    friend std::ostream & operator<<(std::ostream &, SAT &);      // Print the solver state.
};

// OK
SAT::SAT(std::istream & stream) :
    state(State::UNDEF),
    statistics(0, 0, 0, 12, 0, 0, 0, 0, 1024, 0, 0.3, 1.1),
    conflict(false) {
    // Skip optional comments and 'p cnf' appearing at the top of the file.
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
    Size clauses;
    if (!(stream >> statistics.varNum) || !(stream >> clauses))
        throw std::invalid_argument("Error while parsing CNF formula.");

    // Reserve some memory.
    variables.assign(statistics.varNum, {0, false, false, false});
    trail.reserve(statistics.varNum);
    formula.reserve(clauses);

    // Parse all clauses.
    std::vector<Literal> clause;
    Literal literal;
    for (Size i = 0; i < clauses; ++i) {
        while (stream >> literal) {
            if (literal == 0) break;
            clause.push_back(literal);
        }
        addClause(clause);
        clause.clear();
    }

    if (stream.fail())
        throw std::invalid_argument("Error while parsing CNF formula.");
}

// OK
bool SAT::solve() {
    assertPureLiterals();
    while (state == State::UNDEF) {
        unitPropagate();
        if (conflict) {
            if (statistics.decNum == 0)
                state = State::UNSAT;
            else {
                backjumpInit();
                backjumpApply();
            }
        } else {
            if (trail.size() == statistics.varNum)
                state = State::SAT;
            else {
                if (statistics.conRes > statistics.resTrs) restart();
                decide();
            }
        }
    }
    return state == State::SAT;
}

// OK
bool SAT::shouldKeep(std::vector<Literal> & clause) {
    // Sort the clause to make preprocessing easier.
    std::sort(clause.begin(), clause.end());

    // Remove duplicate literals.
    clause.erase(std::unique(clause.begin(), clause.end()), clause.end());

    // Remove falsified literals.
    clause.erase(std::remove_if(clause.begin(), clause.end(),
        [this](Literal literal){ return satisfied(-literal); }), clause.end());

    // If the clause contains a true literal we do not keep it.
    if (std::any_of(clause.begin(), clause.end(),
        [this](Literal literal){ return satisfied(literal); })) return false;

    // If the clause is now of size 0 or 1 we do not keep it.
    if (clause.size() < 2)
        return false;

    // Sort the clause again using absolute values to make the tautological check easier.
    std::sort(clause.begin(), clause.end(), litAbsLessThan);

    // If the clause if tautological we do not keep it.
    for (Size i = 0; i < clause.size() - 1; ++i) {
        if (clause[i] == -clause[i + 1])
            return false;
    }

    return true;
}

// OK
void SAT::addClause(std::vector<Literal> & clause) {
    static Size nextLit = 0;   // Index of the next literal.
    static Size nextCls = 0;   // Index of the next clause.
    const bool     keep = shouldKeep(clause);
    const CSize    size = clause.size();

    if (keep) {
        // Add the literals.
        for (Literal literal : clause) {
            literals.push_back(literal);
            scores.emplace(literal, 0);
        }
        // Add the clause.
        if      (size >= 4) formula.emplace_back(0, size, nextLit, clause[0], clause[1], clause[2], clause[3]);
        else if (size == 3) formula.emplace_back(0, size, nextLit, clause[0], clause[1], clause[2]);
        else                formula.emplace_back(0, size, nextLit, clause[0], clause[1]);
        // Update the watch list.
        watch[clause[0]].push_back(nextCls);
        watch[clause[1]].push_back(nextCls);
        // Update indexes.
        nextLit += size;
        ++nextCls;
        ++statistics.clsFrm;
    } else if (size == 1 && !satisfied(clause[0])) {
        // Assert a pure literal.
        assertLiteral(clause[0], false);
        unitPropagate();
    } else if (size == 0) {
        // A conflict occured at decision level 0.
        state = State::UNSAT;
    }
}

// OK
void SAT::assertLiteral(Literal literal, bool decision) {
    Variable & variable = variables[getIndex(literal)];
    variable.level      = statistics.decNum;
    variable.trail      = true;
    variable.sign       = (literal < 0);
    variable.decision   = decision;
    trail.push_back(literal);
    notify(-literal);
}

// OK
void SAT::assertPureLiterals() {
    // We note for each variable whether it occurs in the positive or negative.
    // If it occurs in only one form, then it can be propagated as a pure literal.
    std::vector<std::pair<bool, bool> > vars(statistics.varNum, {false, false});

    // Go through the formula and add the literals of each clause to the variable list.
    for (const auto & clause : formula) {
        for (Size i = clause.begin(); i < clause.end(); ++i) {
            const Literal literal  = literals[i];
            const Literal absolute = std::abs(literal);
            if (literal > 0) vars[absolute - 1].first  = true;
            else             vars[absolute - 1].second = true;
        }
    }

    // Assert pure literals and check if any variable was deleted during the process.
    // If a variable was deleted, we can add it to the trail in any form.
    for (Size i = 0; i < vars.size(); ++i) {
        const auto & variable = vars[i];
        const Literal literal = i + 1;
        if (!satisfied(literal) && !satisfied(-literal)) {
            if (variable.first && !variable.second) {
                assertLiteral(literal, false);
                unitPropagate();
            } else if (!variable.first && variable.second) {
                assertLiteral(-literal, false);
                unitPropagate();
            } else if (!variable.first && !variable.second) {
                assertLiteral(literal, false);
            }
        }
    }
}

// OK
void SAT::findLastAsserted() {
    for (auto it = trail.rbegin(); it != trail.rend(); ++it) {
        if (backjumpMap.find(-*it) != backjumpMap.end() && backjumpMap[-*it]) {
            lastAss = *it;
            break;
        }
    }
}

// OK
void SAT::backjumpApply() {
    // Explain literals of the backjump clause until the first UIP condition is satisfied.
    while (statistics.litHgh != 1) {
        resolve(lastAss, formula[reason[lastAss]]);
        findLastAsserted();
    }

    // Build the initial backjump clause.
    backjump = backjumpLow;
    backjump.push_back(-lastAss);

    subsumeApply();

    // Build the final backjump clause.
    backjump = backjumpLow;
    backjump.push_back(-lastAss);

    const CSize size = backjump.size();
    Size       level = 0;

    // Do not learn clauses of size 1.
    if (size > 1) {
        // Set the watched literals of the clause.
        std::swap(backjump[0], backjump[size - 1]);
        for (auto it = trail.rbegin(); it != trail.rend(); ++it) {
            auto last = std::find(backjump.begin() + 1, backjump.end(), -*it);
            if (last != backjump.end()) {
                level = variables[getIndex(*it)].level;
                std::swap(*(backjump.begin() + 1), *last);
                break;
            }
        }
        watch[backjump[0]].push_back(formula.size());
        watch[backjump[1]].push_back(formula.size());

        reason[-lastAss] = formula.size();

        // Learn the clause.
        const CSize activity = (size <= statistics.clsTrs) ? size : size + statistics.clsTrs;
        if      (size >= 4) formula.emplace_back(activity, size, literals.size(), backjump[0], backjump[1], backjump[2], backjump[3]);
        else if (size == 3) formula.emplace_back(activity, size, literals.size(), backjump[0], backjump[1], backjump[2]);
        else                formula.emplace_back(activity, size, literals.size(), backjump[0], backjump[1]);
        for (Literal literal : backjump) literals.push_back(literal);
        ++statistics.clsLrn;
    }

    // Execute the backjump.
    conflict  = false;
    unitQueue = {-lastAss};
    prefixTo(level);
}

// OK
void SAT::backjumpInit() {
    // Clear the contents of the previous conflict analysis.
    backjumpMap.clear();
    backjumpLow.clear();
    statistics.litHgh = 0;
    ++statistics.conNum;
    ++statistics.conRes;

    // Set up the backjump containers. Do not explain literals of decision level 0.
    for (Literal literal : backjump) {
        const Size level = variables[getIndex(-literal)].level;
        if (level == statistics.decNum) {
            backjumpMap[literal] = true;
            ++statistics.litHgh;
        } else if (level != 0) {
            backjumpMap[literal] = true;
            backjumpLow.push_back(literal);
        }
    }

    findLastAsserted();
}

// OK
void SAT::decide() {
    for (auto it = scores.begin(); it != scores.end(); ++it) {
        const Literal literal = it->first;
        if (!satisfied(literal) && !satisfied(-literal)) {
            ++statistics.decNum;
            assertLiteral(literal, true);
            return;
        }
    }
}

// OK
void SAT::forget() {
    // Delete 50% of the clauses with the highest activity level.
    double remaining = statistics.clsLrn / 2;

    // Sort the learned clauses by activity level.
    std::set<std::pair<Size, CSize>, ActivityCompare> learned;
    for (Size i = statistics.clsFrm; i < formula.size(); ++i)
        learned.emplace(i, formula[i].activity);

    // For now, only remove the clauses from the watch list.
    for (auto it = learned.begin(); it != learned.end() && remaining > 0; ++it, --remaining) {
        const auto & clause = formula[it->first];
        if (clause.size == 2) {
            --remaining;
            continue;
        }

        --statistics.clsLrn;
        auto & list1 = watch[clause.cache[0]];
        auto & list2 = watch[clause.cache[1]];
        list1.erase(std::remove(list1.begin(), list1.end(), it->first), list1.end());
        list2.erase(std::remove(list2.begin(), list2.end(), it->first), list2.end());
    }

    // Increase the allowed amount of learned clauses.
    statistics.lrnTrs *= statistics.lrnFac;
}

// OK
void SAT::notify(Literal literal) {
    // Check if there is any clause that is watching the literal.
    if (watch.find(literal) == watch.end()) return;

    // New updated watch list for the literal.
    std::vector<Size> newWL;

    for (Size index : watch[literal]) {

        // Clause that is watching the literal.
        Clause & clause = formula[index];

        // If the literal is the first watched literal, swap with the second watched literal.
        if (clause.cache[0] == literal) {
            std::swap(clause.cache[0], clause.cache[1]);
            std::swap(literals[clause.index], literals[clause.index + 1]);
        }

        // Check if the clause is satisfied by looking at the cached literals.
        if (satisfied(clause.cache[0]) || satisfied(clause.cache[1]) ||
            satisfied(clause.cache[2]) || satisfied(clause.cache[3])) {
            newWL.push_back(index);
            continue;
        }

        const Literal watch1 = clause.cache[0];

        // Is there any other unfalsified literal?
        for (Size i = clause.begin() + 2; i < clause.end(); ++i) {
            if (!satisfied(-literals[i])) {
                std::swap(literals[clause.index + 1], literals[i]);
                for (Size j = 0; j < (clause.size > 4 ? 4 : clause.size); ++j)
                    clause.cache[j] = literals[clause.index + j];
                watch[clause.cache[1]].push_back(index);
                goto next;
            }
        }

        newWL.push_back(index);

        // All literals of the clause are falsified if the first watched literal is falsified.
        if (satisfied(-watch1)) {
            // Update the scores.
            for (Size i = clause.begin(); i < clause.end(); ++i) {
                for (auto it = scores.begin(); it != scores.end(); ++it) {
                    if (it->first == literals[i]) {
                        const Literal literal = it->first;
                        const Size   activity = it->second + 1;
                        scores.erase(it);
                        scores.emplace(literal, activity);
                        break;
                    }
                }
            }
            // Set up the conflict.
            if (!conflict) {
                conflict = true;
                backjump.clear();
                for (Size i = clause.begin(); i < clause.end(); ++i)
                    backjump.push_back(literals[i]);
            }
            continue;
        }

        // The first watched literal is unit.
        if (std::find(unitQueue.begin(), unitQueue.end(), watch1) == unitQueue.end()) {
            unitQueue.push_back(watch1);
            reason[watch1] = index;
            // If the clause is a learned clause update its activity.
            // TODO: Check if statistics.clsFrm needs updating after preprocessing.
            if (index >= statistics.clsFrm && (statistics.clsTrs + statistics.decNum) < clause.activity)
                clause.activity = statistics.clsTrs + statistics.decNum;
        }

        next:;
    }

    // Update the watch list for the literal.
    watch[literal] = newWL;
}

// OK
void SAT::prefixTo(Size level) {
    // Update variables that are going to be deleted from the trail.
    Size i = trail.size();
    while (i-- > 0) {
        Variable & variable = variables[getIndex(trail[i])];
        if (variable.level == level) break;
        variable.level    = 0;
        variable.trail    = false;
        variable.sign     = false;
        variable.decision = false;
    }

    // Delete the variables from the trail.
    trail.erase(trail.begin() + i + 1, trail.end());
    statistics.decNum = level;
}

// OK
void SAT::resolve(Literal literal, const Clause & clause) {
    // Delete the literal from the backjump clause.
    backjumpMap[-literal] = false;
    if (variables[getIndex(literal)].level == statistics.decNum) --statistics.litHgh;
    else backjumpLow.erase(std::remove(backjumpLow.begin(), backjumpLow.end(), -literal), backjumpLow.end());

    // Add every literal of the clause, except 'literal', to the backjump clause.
    for (Size i = clause.begin(); i < clause.end(); ++i) {
        const Literal lit = literals[i];
        if (lit != literal && !backjumpMap[lit]) {
            const Size level = variables[getIndex(-lit)].level;
            if (level == statistics.decNum) {
                backjumpMap[lit] = true;
                ++statistics.litHgh;
            } else if (level != 0) {
                backjumpMap[lit] = true;
                backjumpLow.push_back(lit);
            }
        }
    }
}

Size lubySeq(Size rest) {
    if (rest == 1)
        return 1;

    const double param  = std::log2(rest + 1);  

    if (std::floor(param) == param)
        return std::pow(2, param - 1);

    const double powf = std::pow(2, std::floor(param) - 1);
    if (powf == rest)
        return lubySeq(rest + powf - 1);
    
    return lubySeq(rest - std::pow(2, std::ceil(param) - 1) + 1);
}

// OK
void SAT::restart() {
    statistics.conRes = 0;
    statistics.resTrs = 64 * lubySeq(++statistics.resNum);
    prefixTo(0);
    reason.clear();
    scores.clear();
    for (Size i = 1; i <= statistics.varNum; ++i) {
        scores.emplace(i, 0);
        scores.emplace(-i, 0);
    }
    if (statistics.clsLrn > statistics.clsFrm * statistics.lrnTrs) forget();
}

// OK
bool SAT::satisfied(Literal literal) {
    if (literal == 0) return false;
    const Variable & variable = variables[getIndex(literal)];
    return variable.trail && variable.sign == (literal < 0);
}

// OK
void SAT::subsumeApply() {
    for (Size i = 0; i < backjump.size(); ++i) {
        // Skip decision literals as they have no reason clause.
        const Literal literal = -backjump[i];
        if (!satisfied(literal) || variables[getIndex(literal)].decision) continue;

        // Does the reason clause of the literal subsume the backjump clause?
        const Clause & clause = formula[reason[literal]];
        if (clause.size - 1 > backjump.size()) continue;
        if (std::all_of(literals.begin() + clause.begin(), literals.begin() + clause.end(),
            [this, literal] (Literal l) {
            return l == literal || std::find(backjump.begin(), backjump.end(), l) != backjump.end();
        })) {
            // The reason clause subsumes the backjump clause. Explain the literal.
            backjumpMap[-literal]     = false;
            if (variables[getIndex(literal)].level == statistics.decNum) --statistics.litHgh;
            else backjumpLow.erase(std::remove(backjumpLow.begin(), backjumpLow.end(), -literal), backjumpLow.end());

            // Add every literal of the clause, except 'literal', to the backjump clause.
            for (Size j = clause.begin(); j < clause.end(); ++j) {
                const Literal l = literals[j];
                if (l != literal && !backjumpMap[l]) {
                    const Size level = variables[getIndex(-l)].level;
                    if (level == statistics.decNum) {
                        backjumpMap[l] = true;
                        ++statistics.litHgh;
                    } else if (level != 0) {
                        backjumpMap[l] = true;
                        backjumpLow.push_back(l);
                    }
                }
            }

            findLastAsserted();
        }
    }
}

// OK
void SAT::unitPropagate() {
    while (!unitQueue.empty() && !conflict) {
        assertLiteral(unitQueue.front(), false);
        unitQueue.pop_front();
    }
}

// OK
std::ostream & operator<<(std::ostream & out, SAT & sat) {
    switch (sat.state) {
    case SAT::State::UNDEF:
        out << "UNDEF\n";
        break;
    case SAT::State::SAT:
        out << "s SATISFIABLE\n";
        break;
    case SAT::State::UNSAT:
        out << "s UNSATISFIABLE\n";
        break;
    }
    if (sat.state == SAT::State::SAT) {
        std::sort(sat.trail.begin(), sat.trail.end(), [] (Literal l1, Literal l2) { return std::abs(l1) < std::abs(l2); });
        out << "v ";
        for (Literal literal : sat.trail) out << literal << ' '; out << "0";
    }
    return out;
}

// OK
int main(int argc, char ** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <CNF file>\n";
        return 1;
    }

    std::ifstream stream(argv[1]);
    if (stream.fail()) {
        std::cerr << "Usage: " << argv[0] << " <CNF file>\n";
        return 1;
    }

    SAT solver(stream);
    if (solver.solve())
        std::cout << "SAT" << std::endl;
    else
        std::cout << "UNSAT" << std::endl;
    //std::cout << solver << std::endl;
}