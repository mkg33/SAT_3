#include <cstdlib>
#include <fstream>
#include <iostream>

#include "maphSat.hpp"

void printError(char * prog) {
    std::cerr << R"(
███╗   ███╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██╗     ██╗   ██╗███████╗██████╗
████╗ ████║██╔══██╗██╔══██╗██║  ██║██╔════╝██╔═══██╗██║     ██║   ██║██╔════╝██╔══██╗
██╔████╔██║███████║██████╔╝███████║███████╗██║   ██║██║     ██║   ██║█████╗  ██████╔╝
██║╚██╔╝██║██╔══██║██╔═══╝ ██╔══██║╚════██║██║   ██║██║     ╚██╗ ██╔╝██╔══╝  ██╔══██╗
██║ ╚═╝ ██║██║  ██║██║     ██║  ██║███████║╚██████╔╝███████╗ ╚████╔╝ ███████╗██║  ██║
╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚══════╝  ╚═══╝  ╚══════╝╚═╝  ╚═╝
)"  << "\nUsage: " << prog << " <DIMACS file>" << " selection heuristic:\n <FIRST=0 | RANDOM=1 | DLIS=2 | RDLIS=3 | DLCS=4 | RDLCS=5 | JW=6 | RJW=7 | MOMS=8 | RMOMS=9 | VSIDS=10>"
    << " optional DRAT proof logging: <Y|N>\n\n"
    << "Available selection heuristics: \n" << "- FIRST: select the first available literal\n"
    << "- RANDOM: select a random literal\n" << "- DLIS: Dynamic Largest Individual Sum\n" << "- RDLIS: randomized Dynamic Largest Individual Sum\n"
    << "- DLCS: Dynamic Largest Combined Sum\n" << "- RDLCS: randomized Dynamic Largest Combined Sum\n" << "- JW: Jeroslow-Wang heuristic\n"
    << "- RJW: randomized Jeroslow-Wang heuristic\n" << "- MOMS: Maximum [number of] Occurrences in Minimum [length] Clauses\n"
    << "- RMOMS: randomized Maximum [number of] Occurrences in Minimum [length] Clauses\n";
}

int main(int argc, char ** argv) {
    bool proofLogging = false;
    if (argc < 3) {
        printError(argv[0]);
        return 1;
    }
    else if (argc == 4) {
        if (argv[3][0] == 'Y' || argv[3][0] == 'y') {
            proofLogging = true;
        }
        else if (argv[3][0] == 'N' || argv[3][0] == 'n') {
            proofLogging = false;
        }
        else {
            printError(argv[0]);
        }
    }

    std::ifstream stream(argv[1]);
    const int heuristic = atoi(argv[2]);

    if (stream.fail() || heuristic < 0 || heuristic > 10) {
        printError(argv[0]);
        return 1;
    }

    MaphSAT solver(stream, static_cast<MaphSAT::Heuristic>(heuristic), proofLogging);

    solver.solve();
    std::cout << solver;

    if (std::cout.bad()) {
        std::cerr << "Error while printing.\n";
        return 1;
    }
}
