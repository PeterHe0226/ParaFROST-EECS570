/***********************************************************************[function.cuh]
Copyright(c) 2021, Muhammad Osama - Anton Wijs,
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

#ifndef __FUN_
#define __FUN_
#include "elimination.cuh"

namespace ParaFROST {
	constexpr uint32  TABLE_LENGTH = 64;
	typedef uint64 TruthTable[TABLE_LENGTH];

	__constant__ uint64 BITMASKS[6] = {
    	0x5555555555555555ULL, 0x3333333333333333ULL, 0x0F0F0F0F0F0F0F0FULL,
    	0x00FF00FF00FF00FFULL, 0x0000FFFF0000FFFFULL, 0x00000000FFFFFFFFULL
	};

	constexpr int MAX_VARIABLE_COUNT = 12;
	constexpr uint64 FULL_MASK = ~0ULL;

	_PFROST_D_ void clause2fun(const int& variable, const bool& polarity, TruthTable table)
	{
		if (variable >= 6) {
			uint64 alternatingMask = polarity ? FULL_MASK : 0ULL;
			int segmentCounter = 0;
			int segmentSize = 1 << (variable - 6);

			#pragma unroll
			for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
				table[i] |= alternatingMask;
				if (++segmentCounter >= segmentSize) {
					alternatingMask = ~alternatingMask;
					segmentCounter = 0;
				} else {
					continue;
				}
			}
    	} else {
			uint64 mask = BITMASKS[variable];
			if (polarity) mask = ~mask;
			#pragma unroll
			for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
				table[i] |= mask;
			}
    	}
	}

	_PFROST_D_ void orfun(TruthTable dest, const TruthTable source)
	{
		#pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			dest[i] |= source[i];
		}
	}

	_PFROST_D_ void andfun(TruthTable dest, const TruthTable source)
	{
		#pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			dest[i] &= source[i];
		}
	}

	_PFROST_D_ void copyfun(TruthTable dest, const TruthTable source)
	{
		#pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			dest[i] = source[i];
		}
	}

	_PFROST_D_ void falsefun(TruthTable table)
	{
		#pragma unroll
    	for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
        	table[i] = 0ULL;
    	}
	}

	_PFROST_D_ void truefun(TruthTable table)
	{
		 #pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			table[i] = FULL_MASK;
		}
	}

	_PFROST_D_ uint64 collapsefun(const TruthTable firstTable, const TruthTable  secondTable)
	{
		uint64 combinedResult = 0;
		#pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			combinedResult |= (firstTable[i] & secondTable[i]);
		}
		return combinedResult;
	}

	

	_PFROST_D_ bool isfalsefun(const TruthTable table)
	{
		#pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			if (table[i]!= 0ULL) {
				return false;
			}
		}
		return true;
	}

	_PFROST_D_ bool istruefun(const TruthTable table)
	{
		#pragma unroll
		for (uint32 i = 0; i < TABLE_LENGTH; ++i) {
			if (table[i] != FULL_MASK) {
				return false;
			}
		}
		return true;
	}


	_PFROST_D_ bool buildfuntab(const uint32& literal, const uint32* variableMap, CNF& formula, OT& occurrenceTracker, TruthTable sharedTable, TruthTable localTable)
	{
		truefun(localTable);
		OL& occurrences = occurrenceTracker[literal];
		forall_occurs(occurrences, i) {
			SCLAUSE& clause = formula[*i];
        	if (clause.learnt()) {
				continue;
			}
			falsefun(sharedTable);
			forall_clause(clause, j) {
				const uint32 otherLiteral = *j;
				if (otherLiteral == literal) continue;
				const uint32 mappedVariable = variableMap[ABS(otherLiteral)];
				if (mappedVariable >= MAX_VARIABLE_COUNT) {
					return false;
				} else {
					clause2fun(mappedVariable, SIGN(otherLiteral), sharedTable);
				}
			}
			andfun(localTable, sharedTable);
		}
		return true;
	}

	_PFROST_D_ void buildfuntab(const uint32& literal,const uint32* variableMapping,const int& limit,const OL& occurrences,CNF& formula,TruthTable sharedMemTable,TruthTable localTable,bool& isCore)
	{
		 for (int index = 0; index < limit; ++index) {
			SCLAUSE& clause = formula[occurrences[index]];
			if (clause.learnt()) {
				continue;
			}
			falsefun(sharedMemTable);
			forall_clause(clause, iter) {
				const uint32 currentLiteral = *iter;
				if (currentLiteral == literal) {
					continue;
				}
				const uint32 mappedVariable = variableMapping[ABS(currentLiteral)];
				clause2fun(mappedVariable, SIGN(currentLiteral), sharedMemTable);
			}
			andfun(localTable, sharedMemTable);
    	}

		if (!isfalsefun(localTable)) {
			return;
    	} else{
			isCore = true;
        	formula[occurrences[limit]].melt();
		}
	}


	_PFROST_D_ bool countCoreSubstituted(const uint32& variable,const uint32& initialClauseCount,CNF& formula,OL& positiveOccurrences,OL& negativeOccurrences,uint32& coreCount,uint32& addedClauses,uint32& addedLiterals)
	{
		const int clauseLimit = kOpts->ve_clause_max;
		if (kOpts->proof_en) {
			uint32 proofSize = 0;
			forall_occurs(positiveOccurrences, posIndex) {
				SCLAUSE& positiveClause = formula[*posIndex];
				if (positiveClause.learnt()) {
					continue;
				}
				const bool isMolten = positiveClause.molten();
				forall_occurs(negativeOccurrences, negIndex) {
					SCLAUSE& negativeClause = formula[*negIndex];
					if (negativeClause.original() && (!isMolten || !negativeClause.molten())) {
						const int proofLength = mergeProof(variable, positiveClause, negativeClause, proofSize);
						if (proofLength == 1)
							coreCount++;
						else if (proofLength) {
							if (++addedClauses > initialClauseCount || (clauseLimit && proofLength > clauseLimit))
								return true;
							addedLiterals += proofLength;
						}
					}
				}
        	}
			if (coreCount > ADDEDCLS_MAX  || proofSize > ADDEDPROOF_MAX){
				return true;
			} 

        	coreCount = ENCODEPROOFINFO(coreCount, proofSize);

		} else {
			forall_occurs(positiveOccurrences, posIndex) {
				SCLAUSE& positiveClause = formula[*posIndex];
				if (positiveClause.learnt()) {
					continue;
				}
				const bool isMolten = positiveClause.molten();
				forall_occurs(negativeOccurrences, negIndex) {
					SCLAUSE& negativeClause = formula[*negIndex];
					if (negativeClause.original() && (!isMolten || !negativeClause.molten())) {
						const int clauseSize = merge(variable, positiveClause, negativeClause);
						if (clauseSize == 1)
							coreCount++;
						else if (clauseSize) {
							if (++addedClauses > initialClauseCount || (clauseLimit && clauseSize > clauseLimit)){
								return true;
							}
							addedLiterals += clauseSize;
						}
					}
				}
        	}
		}

		if (addedClauses > ADDEDCLS_MAX || addedLiterals > ADDEDLITS_MAX) {
			return true;
		}else if (kOpts->ve_lbound_en) {
			uint32 initialLiteralCount = 0;
			countLitsBefore(formula, positiveOccurrences, initialLiteralCount);
			countLitsBefore(formula, negativeOccurrences, initialLiteralCount);
			if (addedLiterals > initialLiteralCount) {
				return true;
			} else {
				return false;
			}
    	}

    	return false;
	}

	_PFROST_D_ bool find_fun_gate(const uint32& positiveLiteral,const uint32& negativeLiteral,const uint32& originalClauseCount,const uint32* variableMapping,CNF& formula,OT& occurrenceTracker,uint32* sharedMemory,uint32& coreElementCount,uint32& newClauses,uint32& newLiterals)
	{
		TruthTable positiveTable, negativeTable;
		uint64* sharedTable = (uint64*)sharedMemory;

		if (buildfuntab(positiveLiteral, variableMapping, formula, occurrenceTracker, sharedTable, positiveTable) &&
        buildfuntab(negativeLiteral, variableMapping, formula, occurrenceTracker, sharedTable, negativeTable)) {
			if (!collapsefun(positiveTable, negativeTable)) {
				TruthTable& gateTable = positiveTable;
				OL& positiveOccurrenceList = occurrenceTracker[positiveLiteral];
				bool isCore = false;

				for (int i = positiveOccurrenceList.size() - 1; i >= 0; i--) {
					copyfun(gateTable, negativeTable);
					if (formula[positiveOccurrenceList[i]].original()){
						buildfuntab(positiveLiteral, variableMapping, i, positiveOccurrenceList, formula, sharedTable, gateTable, isCore);
					}
				}
				OL& negativeOccurrenceList = occurrenceTracker[negativeLiteral];
				for (int i = negativeOccurrenceList.size() - 1; i >= 0; i--) {
					truefun(gateTable);
					if (formula[negativeOccurrenceList[i]].original()){
						buildfuntab(negativeLiteral, variableMapping, i, negativeOccurrenceList, formula, sharedTable, gateTable, isCore);
					}
				}

				coreElementCount = 0;
				newClauses = 0;
				newLiterals = 0;
				if (countCoreSubstituted(ABS(positiveLiteral), originalClauseCount, formula, positiveOccurrenceList, negativeOccurrenceList, coreElementCount, newClauses, newLiterals)) {
					if (isCore) {
						freezeClauses(formula, positiveOccurrenceList, negativeOccurrenceList);
					}
					return false;
				}

				return true;
			}
    	}
    	return false;
	}
} 
#endif