import os
import subprocess
import sys

class bcolors:
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    TIMEOUT = '\033[93m'
    BOLD = '\033[1m'
    ENDC = '\033[0m'

def runForAll(exe, dir, res):
    out = []
    heuristics = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    files = sorted(os.listdir(dir))
    for heuristic in heuristics:
        print("HEURISTIC: ", heuristic)
        for filename in files:
            if filename.endswith(".cnf"):
                fOut    = []
                try:
                    result = subprocess.run(["time", "./" + exe, dir + "/" + filename, str(heuristic)], capture_output=True, text=True, timeout = 600)
                    execTime = result.stderr.strip().split()[2]
                    output   = result.stdout.strip().split()[0]
                    if output != res:
                        fOut.append("[" + bcolors.FAIL + "FAIL" + bcolors.ENDC + "]")
                    else:
                        fOut.append("[" + bcolors.OKGREEN + "OK" + bcolors.ENDC + "]")
                    fOut.append(output)
                    fOut.append(execTime)
                    fOut.append(filename)

                except subprocess.TimeoutExpired:
                    fOut.append("[" + bcolors.TIMEOUT + "TIMEOUT" + bcolors.ENDC + "]")
                    fOut.append(filename)
                print(" ".join(fOut))

def main(argv):
    if len(argv) < 4:
        print(argv[0], " <executable> <directory of cnf files> <expected result>")
        sys.exit(2)
    else:
        runForAll(argv[1], argv[2], argv[3])

if __name__ == "__main__":
    main(sys.argv)
