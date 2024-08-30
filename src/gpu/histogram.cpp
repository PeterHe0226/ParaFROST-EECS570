/***********************************************************************[histogram.cpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

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

#include "solver.hpp" 

using namespace ParaFROST;

inline bool Solver::isBinary(const C_REF& r, uint32& first, uint32& second)
{
	assert(!DL());
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	first = 0, second = 0;
	forall_clause(c, k) {
		const uint32 lit = *k;
		CHECKLIT(lit);
		const LIT_ST val = sp->value[lit];
		if (UNASSIGNED(val)) {
			if (second) return false; // not binary
			if (first) second = lit;
			else first = lit;
		}
		else if (val) {
			// satisfied
			removeClause(c, r);
			return false;
		}
	}
	if (!second) return false; // all falsified except 'first'
	return true;
}

inline void	Solver::resetoccurs()
{
	forall_vector(OCCUR, occurs, o) {
		o->ps = o->ns = 0;
	}
}

inline void	Solver::countoccurs(CLAUSE& c)
{
	assert(!c.deleted());
	forall_clause(c, k) {
		const uint32 lit = *k;
		CHECKLIT(lit);
		if (SIGN(lit)) occurs[ABS(lit)].ns++;
		else occurs[ABS(lit)].ps++;
	}
}

void Solver::histClause(CLAUSE& c)
{
	assert(!c.deleted());
	forall_clause(c, k) {
		CHECKLIT(*k);
		vhist[*k]++;
	}
}

void Solver::histBins(BCNF& cnf)
{
	forall_cnf(cnf, i) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		uint32 a, b;
		if (isBinary(r, a, b)) {
			CHECKLIT(a), CHECKLIT(b);
			vhist[a]++;
			vhist[b]++;
		}
	}
}

void Solver::histCNF(BCNF& cnf, const bool& reset) {
	if (cnf.empty()) return;
	if (reset) resetoccurs();
	forall_cnf(cnf, i) {
		const C_REF ref = *i;
		if (cm.deleted(ref)) continue;
		countoccurs(cm[*i]);
	}
}

void Solver::histCNFFlat(BCNF& cnf, const bool& reset) {
	if (cnf.empty()) return;
	if (reset) {
		forall_vector(uint32, vhist, h) { 
			*h = 0;
		}
	}
	forall_cnf(cnf, i) {
		const C_REF ref = *i;
		if (cm.deleted(ref)) continue;
		histClause(cm[*i]);
	}
}