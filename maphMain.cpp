#include <cstdlib>
#include <fstream>
#include <iostream>

#include "maphSat.hpp"

namespace fs = std::__fs::filesystem;

void printError(char * prog) {
    std::cerr << R"(
███╗   ███╗ █████╗ ██████╗ ██╗  ██╗███████╗ ██████╗ ██╗     ██╗   ██╗███████╗██████╗
████╗ ████║██╔══██╗██╔══██╗██║  ██║██╔════╝██╔═══██╗██║     ██║   ██║██╔════╝██╔══██╗
██╔████╔██║███████║██████╔╝███████║███████╗██║   ██║██║     ██║   ██║█████╗  ██████╔╝
██║╚██╔╝██║██╔══██║██╔═══╝ ██╔══██║╚════██║██║   ██║██║     ╚██╗ ██╔╝██╔══╝  ██╔══██╗
██║ ╚═╝ ██║██║  ██║██║     ██║  ██║███████║╚██████╔╝███████╗ ╚████╔╝ ███████╗██║  ██║
╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚══════╝  ╚═══╝  ╚══════╝╚═╝  ╚═╝
)"  << "\nUsage: " << prog << " <DIMACS file>" << " selection heuristic:\n <VSIDS=0 | EVSIDS=1>"
    << " optional DRAT proof logging: <Y>\n\n"
    << "Available selection heuristics:\n" << "- VSIDS: Variable State Independent Decaying Sum\n"
    << "- EVSIDS: Exponential Variable State Independent Decaying Sum\n"
    << "- WEVSIDS: Weighted Exponential Variable State Independent Decaying Sum\n";
    // I suggest we git rid of the choice of heuristics (switching between cases is too slow when updating scores).
    // I also think we should leave them in code and allow the user to comment them out to check that they really work.
    // We should also use all three heuristics for overlaying cactus plots.
}

int main(int argc, char ** argv) {
    bool proofLogging = false;
    std::string proofName = "";

    if (argc < 3) {
        printError(argv[0]);
        return 1;
    }
    else if (argc == 4) {
        if (argv[3][0] == 'Y' || argv[3][0] == 'y') {
            proofLogging = true;
            proofName = fs::path(argv[1]).filename();
        }
        else {
            printError(argv[0]);
        }
    }

    std::ifstream stream(argv[1]);
    const int heuristic = atoi(argv[2]);

    if (stream.fail() || heuristic < 0 || heuristic > 1) {
        printError(argv[0]);
        return 1;
    }

    MaphSAT solver(stream, static_cast<MaphSAT::Heuristic>(heuristic), proofLogging, proofName);

    solver.solve();
    std::cout << solver;

    if (std::cout.bad()) {
        std::cerr << "Error while printing.\n";
        return 1;
    }
}
