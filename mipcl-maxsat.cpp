// g++ -Wall -g -DMIP_API= -I/usr/local/mipcl-1.1.2/mipcl/headers mipcl-maxsat.cpp -L/usr/local/mipcl-1.1.2/lib -lmipcl -o mipcl-maxsat
/*
Copyright 2016 Masahiro Sakai. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
   3. The name of the author may not be used to endorse or promote
      products derived from this software without specific prior
      written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#define __STDC_FORMAT_MACROS

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <vector>
#include <map>
#include <iostream>
#include <cmip.h>
#include <except.h>

class CMaxSat : public CMIP {

public:
  
  CMaxSat(const char* name) : CMIP(name) {
  }
  
  virtual void infoMessage(const char* msg, int level=0) {
      std::cout << "c " << msg << std::endl;
  }

  /*
  virtual void lpInfo(const char *method, const char *time, int itNum, int degItNum, double objVal) {
  }

  virtual void mipInfo(char *timeElapsed, int nodeNum, int leafNum,
               double bestObjVal, double objBound, double gap, int solsFound, bool sense, bool header) {
  }

  virtual void cutInfo(__LONG time, int round, double objVal, int fracNum, int cutNum) {
  }

  virtual void probingInfo(char* timeStr, int round,
                   int probeVarNum,int varFixed,int ctrTightened,int varBdAdded, int implications) {
  }
  */
  /*
  virtual void solStatistics(std::ostream &out, const char* MIPCLver,
                             const char* solTime, bool timeLimit, int nodeNum,
                             bool feasible, bool hasSolution, double objVal,
                             bool opt, double gap, bool gapLimit, double bound,
                             int difficultNodes)
  {
      if (out == std::cout)
        out << "c ";
      CMIP::solStatistics(out,MIPCLver,solTime,timeLimit,nodeNum,feasible,hasSolution,objVal,opt,gap,gapLimit,bound,difficultNodes);
  }
  */
};

static
int read_wcnf(CMaxSat &prob, const char *filename)
{
    FILE *file = fopen(filename, "r");
    char line[1024*128];
    int nv, nc;
    int64_t top = -1;
    bool isWCNF = 0;

    while (1) {
        fgets(line, sizeof(line), file);
        if (line[0] == 'c')
            continue;
        else if (line[0] == 'p') {
            int ret = sscanf(line, "p cnf %d %d", &nv, &nc);
            if (ret == 2) goto BODY;

            ret = sscanf(line, "p wcnf %d %d %"PRId64, &nv, &nc, &top);
            if (ret >= 2) {
                isWCNF = 1;
                goto BODY;
            }
        }

        fprintf(stderr, "unexpected line: %s\n", line);
        exit(1);
    }

BODY:
    std::vector<std::map<int,double> > rows;
    std::vector<int> lbs;
    int nz = 0;
    std::vector<int64_t> obj(1+nv);

    for (int i = 1; i <= nc; i++) {
        int64_t cost = 1;
        if (isWCNF) fscanf(file, " %"PRId64, &cost);

        std::vector<int> lits;
        while (1) {
            int lit;
            fscanf(file, "%d", &lit);
            if (lit == 0) break;
            lits.push_back(lit);
        }

        if (cost != top && lits.size() == 1) {
            int lit = lits[0];
            if (lit > 0) {
                int v = lit;
                // obj += cost*(1 - v)
                obj[0] += cost;
                obj[v] -= cost;
            } else {
                int v = - lit;
                // obj += cost*v
                obj[v] += cost;
            }
            continue;
        }

        std::map<int,double> coeffs;

        if (cost != top) {
	    int r = obj.size();
            obj.push_back(cost);
            coeffs[r] = 1;
        }

        // we need to sum up coefficients of a same variable.
        int lb = 1;
        for (std::vector<int>::iterator j = lits.begin(); j != lits.end(); j++) {
            int lit = *j;
            if (lit > 0) {
                coeffs[lit] += 1;
            } else {
                coeffs[-lit] -= 1;
                lb--;
            }
        }

        nz += coeffs.size();
        rows.push_back(coeffs);
        lbs.push_back(lb);
    }

    fclose(file);

    printf("openMatrix(%d, %d, %d)\n", (int)rows.size(), (int)obj.size(), (int)nz);
    prob.openMatrix(rows.size(), obj.size(), nz);
    prob.setObjSense(false);
    prob.addVar(0, CMIP::VAR_INT, 0, 1.0, 1.0); // constant 1.0
    for (int i = 1; i < obj.size(); i++) {
        printf("addVer(%d)\n",i);
        prob.addVar(i, CMIP::VAR_BIN, obj[i], 0.0, 1.0);
    }
    for (int i = 0; i < rows.size(); i++) {
        printf("addCtr(%d,0,%d,CLP::INF)\n",i,lbs[i]);
	prob.addCtr(i, 0, lbs[i], CLP::INF);
        std::map<int,double> coeffs = rows[i];
        for (std::map<int,double>::iterator j = coeffs.begin(); j != coeffs.end(); j++) {
            printf("addEntry(%f, %d, %d)\n",j->second, i, j->first);
	    prob.addEntry(j->second, i, j->first);
	}
    }
    prob.closeMatrix();

    return nv;
}

static
int output_comment(void *info, const char *s)
{
    printf("c %s", s);
    fflush(stdout);
    return 1;
}

/*
static
void print_model(glp_prob *prob, int nv)
{
    for (int i = 1; i <= nv; i++) {
        if (i % 10 == 1) {
            if (i != 1) puts(""); // new line
            fputs("v", stdout);
        }
        if (glp_mip_col_val(prob, i) >= 0.5)
            printf(" %d", i);
        else
            printf(" -%d", i);
    }
    puts(""); // new line
}

static
void intopt_callback(glp_tree *tree, void *info)
{
    switch (glp_ios_reason(tree)) {
    case GLP_IBINGO: // Better integer solution found
        {
            glp_prob *prob = glp_ios_get_prob(tree);
            int64_t val = (int64_t)llround(glp_mip_obj_val(prob));
            printf("o %"PRId64 "\n", val);
            fflush(stdout);
        }
        break;
    default:
        break;
    }
    return;
}
*/

int main(int argc, char **argv)
{
    if (1 >= argc) {
        fprintf(stderr, "USAGE: mipcl-maxsat [file.cnf|file.wcnf]");
        exit(1);
    }
    char *filename = argv[1];

    try {
        puts("B");
        CMaxSat prob(filename); // 1
	puts("C");
        int nv = read_wcnf(prob, filename);      
	puts("D");
        prob.optimize();
	puts("E");

	if (prob.isSolution()) {
	    double *dpX;
	    int *ipHd;
	    int n = prob.getSolution(dpX, ipHd);
	    // dpX,ipHd	dpX[j] is value of variable whose handle is ipHd[j], j=1,...,n, where n is return value.
            //prob.printSolution("primer.sol"); // (10)
	} else {
	    puts("s UNKNOWN");
	    exit(1);
	}
    }
    catch(CException* pe) {
        std::cerr << pe->getErrorMessage() << std::endl;
        delete pe;
        return 1;
    }
    return 0;
}
