#ifndef __MAPHSAT_HPP__
#define __MAPHSAT_HPP__

#include <cstdlib>
#include <deque>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

class MaphSAT {

public:

    enum class Heuristic {
        FIRST,
        RANDOM,
        DLIS,
        RDLIS,
        DLCS,
        RDLCS,
        JW,
        RJW,
        MOMS,
        RMOMS
    };
    Heuristic heuristic;

private:

    enum class State {
        UNDEF,  // The formula has not yet been found to be satisfiable or unsatisfiable.
        SAT,    // The formula is satisfiable.
        UNSAT   // The formula is unsatisfiable.
    };
    State state;

    std::size_t numberVariables;
    std::size_t numberClauses;
    std::size_t numberDecisions;

    bool conflict;

    // The formula in CNF format. Each inner vector represents a clause.
    std::vector<std::vector<int> > formula;

    // The trail represents the current partial evaluation, with each pair being
    // a literal and a boolean denoting whether it is a decision literal or not.
    std::vector<std::pair<int, bool> > trail;

    // The variable assignment that lead to a conflict and its opposite.
    std::vector<int> backjumpClause;

    // Maps a propagated literal to the clause that forced its propagation.
    std::unordered_map<int, std::size_t> reason;

    // Literals than can be unit propagated and the clause that forced the propagation.
    std::deque<int> unitQueue;

    // Maps a literal to the clauses that are watching the literal.
    std::unordered_map<int, std::vector<std::size_t> > watchList;

    void combinedSum(std::vector<std::pair<int, int> > &, std::vector<std::pair<int, int> > &, bool, std::size_t) const;
    int selectFirst() const;
    int selectRandom() const;
    int selectDLIS(bool) const;
    int selectDLCS(bool) const;
    int selectJW(bool) const;
    int selectMOMS(bool) const;

    // Elimiate pure literals.
    void pureLiteral();

    // Assert a literal as a decision literal or as a non-decision literal.
    void assertLiteral(int, bool);

    // Select a literal that is not yet asserted and assert it as a decision literal.
    void applyDecide();

    // If there are any unit literals due to the current partial evaluation, assert
    // them as non-decision literals. Repeat until there are no more unit literals.
    void applyUnitPropagate();

    // START CONFLICT RESOLUTION

    // Once a conflict occurs, the conflict resolution process starts. A conflict occurs
    // when unit propagation forces the propagation of a literal l, even though literal
    // -l is already asserted. When such a conflict occurs, we need to find the variable
    // assignments that are incompatible and lead to the conflict. Let's say literal 1
    // and literal 3 are incompatible and lead to a conflict. Then the clause [-1, -3]
    // represents this incompatability and we can add this clause to our formula to prevent
    // the same conflict to happen repeatedly. To find this clause, we can imagine a
    // directed graph. Each node of the graph is an asserted literal and an edge from
    // node a to node b means that node a was a reason for the propagation of node b. Let's
    // say we have the clause [1, 2, 3] and -1 and -2 is asserted. Then the graph
    // contains the nodes -1, -2 and 3 with an edge from -1 to 3 and an edge from -2
    // to 3 due to -1 and -2 forcing the propagation of 3. When we detect a conflict,
    // meaning there are literals l and -l as nodes in the graph, we split the graph
    // into a reason side and a conflict side. The conflict side contains literals l
    // and -l. We add all nodes to the conflict side, that we reach before we get
    // to the first occurance of a UIP (first UIP) while going backwards to the node of the current
    // decision literal. A UIP (unique implication point) is any node of the current decision
    // level such that any path from the current decision literal to the conflict node
    // must pass through it. Then the conflict clause consists of all nodes, belonging
    // to the reason side, that have an edge into the conflict side.

    // Returns the number of decision literals in the trail that precede the first
    // occurance of 'literal', including 'literal' itself if it is a decision literal.
    std::size_t level(int literal) const;

    // Returns a literal from 'clause' that is in the trail such that no other
    // literal from 'clause' comes after it in the trail.
    int lastAssertedLiteral() const;

    //int lastAssertedLiteralNonDecision(const std::vector<int> &) const;

    // Check if the backjump clause satisfies the first UIP condition, which is the
    // case if the backjump clause contains exactly one literal of the current decision level.
    bool isUIP() const;

    // Perform a single resolution step between the backjump clause and a clause
    // that is the reason for the propagation of -'literal'.
    void applyExplain(int);

    // Construct the backjump clause by repeatedly explaining a literal that lead to a
    // conflict until the backjump clause satisfies the first UIP condition.
    void applyExplainUIP();

    // Add a learned clause to the formula to prevent the same conflict from happening again.
    void applyLearn();

    // Return an iterator to the first literal in the trail that has a decision level greater than 'level'.
    std::vector<std::pair<int, bool> >::iterator firstLiteralPast(int);

    // Remove any literals from the trail that have a decision level greater than 'level'.
    void removePast(int);

    // Return the greatest decision level of the backjump clause exluding 'literal'.
    int getBackjumpLevel(int);

    // Backtrack literals from the trail until the backjump clause becomes a unit
    // clause and then assert the unit literal.
    void applyBackjump();

    // Notify clauses that a literal has been asserted.
    void notifyWatches(int);

    //bool pureLiteral();

public:

    // Parse a CNF formula and throw invalid_argument() if unsuccessful.
    MaphSAT(std::istream &, Heuristic);

    // Solve the CNF formula.
    bool solve();

    // Print the current state of the SAT solver.
    friend std::ostream & operator<<(std::ostream &, const MaphSAT &);

public:

};

#endif
