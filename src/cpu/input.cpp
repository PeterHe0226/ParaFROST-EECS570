/***********************************************************************[input.cpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Technische Universiteit Eindhoven (TU/e).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************************/

#include "sort.h"
#include "input.h"
#include "control.h"
#include "version.h"

namespace ParaFROST {

    Vec<ARG*, int> options; // container for all options available

    void printUsage(int argc, char** argv, bool verbose)
    {
        PFNAME("ParaFROST (Parallel Formal ReasOning about Satisfiability)", version());
        PFAUTHORS("Muhammad Osama Mahmoud");
        PFLOG0("");
        uint64 sysmem = 0;
        getCPUInfo(sysmem);
        getBuildInfo();
        PFLRULER('-', RULELEN);
        PFLOG0("");
        PFLOG1(" Usage: %sparafrost%s [ <cnf> ][ <option> ... ]", CSOLVER, CNORMAL);
        PFLOG0("");
        Sort(options.data(), options.size(), ARG::ARG_CMP());
        arg_t prev_type = NULL;
        PFLOG0("");
        PFLOG1(" Options (simplification + solve):");
        for (int i = 0; i < options.size(); i++) {
            if (options[i]->type != prev_type) {
                PRINT(PREFIX);
                PUTCH('\n');
            }
            options[i]->help(verbose);
            prev_type = options[i]->type;
        }
        PFLOG0("");
        PFLOG1("  %s-h or --help  print available options.%s", CHELP, CNORMAL);
        PFLOG1("  %s--helpmore    print available options with verbose message.%s", CHELP, CNORMAL);
        PFLOG0("");
        PFLRULER('-', RULELEN);
        exit(EXIT_SUCCESS);
    }

    bool parseArguments(int& argc, char** argv)
    {
        if (argc <= 1) return false;
        const char* arg = argv[1];
        int dashes = (arg[0] == '-') + (arg[1] == '-');
        if ((dashes & 1) && arg[1] == 'h')
            printUsage(argc, argv);
        else if ((dashes & 2) && hasstr(arg, "help")) {
            if (hasstr(arg, "more"))
                printUsage(argc, argv, true);
            else
                printUsage(argc, argv);
        }
        struct stat st;
        bool ispath = canAccess(arg, st);
        bool ret = false;
        for (int i = 1 + ispath; i < argc; i++) {
            const size_t arglen = strlen(argv[i]);
            if (arglen == 1) 
                PFLOGE("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
            else if (arglen > 1) {
                const char* arg = argv[i];
                int dashes = (arg[0] == '-') + (arg[1] == '-');
                if (!dashes) 
                    PFLOGE("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
                else if ((dashes & 1) && arg[1] == 'h')
                    printUsage(argc, argv);
                else if ((dashes & 2) && hasstr(arg, "help")) {
                    if (hasstr(arg, "more")) 
                        printUsage(argc, argv, true);
                    else 
                        printUsage(argc, argv);
                }
                else {
                    int k = 0;
                    bool parsed = false;
                    while (k < options.size() && !(parsed = options[k++]->parse(argv[i])));
                    if (!parsed)  PFLOGE("unknown input \"%s\". Use '-h or --help' for help.", argv[i]);
                    if (!ret && parsed) ret = true;
                }
            }
        }
        return ret;
    }

    void ARG::insert(ARG* opt) { options.push(this); }

}