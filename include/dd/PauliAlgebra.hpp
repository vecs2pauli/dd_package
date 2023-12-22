//
// Created by lieuwe on 13/12/2021.
//

#ifndef DDPACKAGE_PAULIALGEBRA_HPP
#define DDPACKAGE_PAULIALGEBRA_HPP

#include "ComplexNumbers.hpp"
#include "LimTable.hpp"
#include "Log.hpp"
#include "PauliUtilities.hpp"

#include <algorithm>
#include <array>
#include <iostream>

namespace dd {

// TODO write a test
    // TODO refactor to use recoverElement
    // precondition: group is in column echelon form
    template<std::size_t NUM_QUBITS>
    inline phase_t recoverPhase(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* a, const Qubit nQubits = NUM_QUBITS - 1) {
        if (a == LimEntry<NUM_QUBITS>::noLIM) {
            throw std::runtime_error("[recoverPhase] a is noLIM.\n");
        }
        // startProfile(recoverPhase)
        //std::cout << "in recoverphase\n";
        LimEntry<NUM_QUBITS> A(a);
        LimEntry<NUM_QUBITS> B;
        for (std::size_t g = 0; g < G.size(); g++) {
            //std::cout << "in forloop recoverphase\n";
            //std::cout << LimEntry<>::to_string(G[g], nQubits) << ", ";
            auto const pivot = G[g]->pivotPosition(nQubits);
            //std::cout << "got pivotposition: " << pivot << std::endl;
            if (A.paulis.test(pivot)) {
                //std::cout << "pivot getest\n";
                A.multiplyBy(*G[g], nQubits);
                B.multiplyBy(*G[g], nQubits);
            }
            //std::cout << "dan niet\n";
        }
        // endProfile(recoverPhase)
        //std::cout << "klaar met recoverphase\n";
        return B.getPhase();
    }

    // precondition: group is in column echelon form
    template<std::size_t NUM_QUBITS>
    inline phase_t recoverPhase(const std::vector<LimEntry<NUM_QUBITS>>& G, const LimEntry<NUM_QUBITS>* a, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(recoverPhase)
        LimEntry<NUM_QUBITS> A(a);
        LimEntry<NUM_QUBITS> B;
        for (std::size_t g = 0; g < G.size(); g++) {
            auto const pivot = G[g].pivotPosition();
            if (A.paulis.test(pivot)) {
                A.multiplyBy(G[g], nQubits);
                B.multiplyBy(G[g], nQubits);
            }
        }
        // endProfile(recoverPhase)
        return B.getPhase();
    }

    // precondition: group is in column echelon form
    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> recoverElement(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* a, const Qubit nQubits = NUM_QUBITS - 1) {
        if (a == LimEntry<NUM_QUBITS>::noLIM) {
            throw std::runtime_error("[recoverPhase] a is noLIM.\n");
        }
        // startProfile(recoverPhase)
        LimEntry<NUM_QUBITS> A(a);
        LimEntry<NUM_QUBITS> B;
        for (std::size_t g = 0; g < G.size(); g++) {
            auto const pivot = G[g]->pivotPosition();
            if (A.paulis.test(pivot)) {
                A.multiplyBy(G[g], nQubits);
                B.multiplyBy(G[g], nQubits);
            }
        }
        // endProfile(recoverPhase)
        return B;
    }


    // Performs Gaussian elimination on G
    // We assume that G is not stored in the LimTable.
    // In more detail: the elements of G are modified in place
    // Therefore, the LimEntry objects should NOT be stored in the LimTable;
    template<std::size_t NUM_QUBITS>
    inline void GaussianElimination(std::vector<LimEntry<NUM_QUBITS>*>& G) {
        // startProfile(gaussianElimination)
        if (G.size() <= 1) {
            return;
        }
        //        Log::log << "[Gaussian Elimination] start. |G| = " << G.size() << ".\n"; Log::log.flush();
        constexpr std::size_t pauli_height = 2 * NUM_QUBITS; // length of the columns as far as they contain Pauli operators
        for (std::size_t h = 0; h < pauli_height; h++) {
            // Step 1: Find a column with a '1' at position h
            std::size_t reducingColId = std::numeric_limits< decltype(reducingColId)>::max();
            for (std::size_t i = 0; i < G.size(); i++) {
                if (G[i]->pivotPosition() == h) {
                    reducingColId = i;
                    break;
                }
            }
            if (reducingColId == std::numeric_limits< decltype(reducingColId)>::max()) {
                continue;
            }
            // Step 2: Reduce all other columns via G[reducingColId]
            for (std::size_t reduceColId = 0; reduceColId < G.size(); reduceColId++) {
                if (reduceColId == reducingColId) {
                    continue;
                }
                if (G[reduceColId]->paulis.test(h)) {
                    //                    Log::log << "[Gaussian Elimination] Multiplying col " << reduceColId << " with col " << reducingColId << ".\n"; Log::log.flush();
                    G[reduceColId] = LimEntry<NUM_QUBITS>::multiply(G[reduceColId], G[reducingColId]);
                    //G[reduceColId]->multiplyBy(G[reducingColId]); // TODO this would solve the memory leak; but it leads to failed tests
                }
            }
        }
        // endProfile(gaussianElimination)
    }

    //    // precondition: G is sorted
    //    template<std::size_t NUM_QUBITS>
    //    inline void GaussianEliminationSortedFast(std::vector<LimEntry<NUM_QUBITS>*>& G) {
    //        if (G.size() <= 1) return;
    //        unsigned int pivot;
    //
    //        for (unsigned int g = 0; g + 1 < G.size(); g++) {
    //            pivot = G[g]->pivotPosition();
    //            if (pivot >= 2 * NUM_QUBITS) continue; // In this case G[g] is an all-zero column (i.e., is the identity)
    //            for (unsigned int h = g + 1; h < G.size(); h++) {
    //                if (G[h]->paulis.test(pivot)) {
    //                    G[h] = LimEntry<NUM_QUBITS>::multiply(G[h], G[g]);
    //                }
    //            }
    //        }
    //    }

    template<std::size_t NUM_QUBITS>
    inline void GaussianEliminationSortedFast(std::vector<LimEntry<NUM_QUBITS>>& G, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(gaussianElimination)
        if (G.size() <= 1) {
            return;
        }

        for (std::size_t g = 0; g + 1 < G.size(); g++) {
            auto const pivot = G[g].pivotPosition(nQubits);
            if (pivot >= 2 * NUM_QUBITS) {
                continue; // In this case G[g] is an all-zero column (i.e., is the identity)
            }
            for (std::size_t h = g + 1; h < G.size(); h++) {
                if (G[h].paulis.test(pivot)) {
                    G[h].multiplyBy(G[g], nQubits);
                }
            }
        }
        // endProfile(gaussianElimination)
    }

    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline void GaussianEliminationModuloPhaseSortedFast(std::vector<LimBitset<NUM_QUBITS, NUM_BITS>*>& G) {
        // startProfile(gaussianElimination)
        if (G.size() <= 1) {
            return;
        }

        for (std::size_t g = 0; g + 1 < G.size(); g++) {
            auto const pivot = G[g]->lim.pivotPosition();
            if (pivot >= 2 * NUM_QUBITS) {
                continue;
            }
            for (std::size_t h = g + 1; h < G.size(); h++) {
                if (G[h]->lim.paulis.test(pivot)) {
                    G[h] = LimBitset<NUM_QUBITS, NUM_BITS>::multiply(G[h], G[g]);
                }
            }
        }
        // endProfile(gaussianElimination)
    }

    // TODO make faster by telling the multiply() routine how many qubits there are (!)
    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline void GaussianEliminationModuloPhaseSortedFast(std::vector<LimBitset<NUM_QUBITS, NUM_BITS>>& G, Qubit nQubits = NUM_QUBITS-1) {
        // startProfile(gaussianElimination)
        if (G.size() <= 1) {
            return;
        }

        for (std::size_t g = 0; g + 1 < G.size(); g++) {
            auto const pivot = G[g].lim.pivotPosition(nQubits);
            if (pivot >= 2 * NUM_QUBITS) {
                continue;
            }
            for (std::size_t h = g + 1; h < G.size(); h++) {
                if (G[h].lim.paulis.test(pivot)) {
                    //G[h] = LimBitset<NUM_QUBITS, NUM_BITS>::multiply(G[h], G[g]); // TODO use multiplyBy (this avoids a copy constructor)
                    G[h].multiplyBy(G[g], nQubits);
                }
            }
        }
        // endProfile(gaussianElimination)
    }

    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline bool isAbelian(const std::vector<LimBitset<NUM_QUBITS, NUM_BITS>*>& G) {
        for (std::size_t i = 0; i < G.size(); i++) {
            for (std::size_t j = i + 1; j < G.size(); j++) {
                if (!G[i]->lim.commutesWith(G[j]->lim)) {
                    return false;
                }
            }
        }
        return true;
    }

    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline bool isAbelian(const std::vector<LimBitset<NUM_QUBITS, NUM_BITS>>& G) {
        for (std::size_t i = 0; i < G.size(); i++) {
            for (std::size_t j = i + 1; j < G.size(); j++) {
                if (!G[i].lim.commutesWith(G[j].lim)) {
                    return false;
                }
            }
        }
        return true;
    }

    // Performs Gaussian Elimination on the group G, ignoring the phase of the LIMs involved
    // todo it is possible to write a faster procedure, if we are allowed to assume that G is sorted
    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline void GaussianEliminationModuloPhase(std::vector<LimBitset<NUM_QUBITS, NUM_BITS>*>& G) {
        if (G.size() <= 1) {
            return;
        }
        //        Log::log << "[Gaussian Elimination Bitset] start. |G| = " << G.size() << ".\n"; Log::log.flush();
        constexpr std::size_t pauli_height = 2 * NUM_QUBITS; // length of the columns as far as they contain Pauli operators
        for (std::size_t row = 0; row < pauli_height; row++) {
            // Step 1: Find a column with a '1' at position row
            std::size_t reducingColId = std::numeric_limits<decltype(reducingColId)>::max();

            for (std::size_t i = 0; i < G.size(); i++) {
                if (G[i]->lim.pivotPosition() == row) {
                    reducingColId = i;
                    break;
                }
            }
            if (reducingColId == std::numeric_limits<decltype(reducingColId)>::max()) {
                continue;
            }
            // Step 2: Reduce all other columns via G[reducingColId]
            for (std::size_t reduceColId = 0; reduceColId < G.size(); reduceColId++) {
                if (reduceColId == reducingColId) {
                    continue;
                }
                if (G[reduceColId]->lim.paulis.test(row)) {
                    //                    Log::log << "[Gaussian Elimination Bitset] Multiplying col " << reduceColId << " with col " << reducingColId << ".\n"; Log::log.flush();
                    G[reduceColId] = LimBitset<NUM_QUBITS, NUM_BITS>::multiply(G[reduceColId], G[reducingColId]);
                }
            }
        }
    }

    //    // Puts the stabilizer group in column echelon form; specifically:
    //    //   1. performs Gaussian elimination on G
    //    //   2. prunes the all-zero columns
    //    //   3. sorts the columns lexicographically, i.e., so that 'pivots' appear in the matrix
    //    inline void toColumnEchelonForm(StabilizerGroup& G) {
    //        std::sort(G.begin(), G.end(), LimEntry<>::geq);
    //        GaussianEliminationSortedFast(G);
    //        //        Log::log << "[toColumnEchelonForm] After Gaussian Elimination, before pruning zero col's, group is:\n";Log::log.flush();
    //        //        printStabilizerGroup(G);
    //        pruneZeroColumns(G);
    //        // To obtain a lower triangular form, we now sort the vectors descending lexicographically, descending
    //        std::sort(G.begin(), G.end(), LimEntry<>::geq);
    //        //        Log::log << "[toColumnEchelonForm] After CEF, group is:\n"; Log::log.flush();
    //        //        printStabilizerGroup(G);
    //    }

    inline void toColumnEchelonForm(StabilizerGroupValue& G, const Qubit nQubits = NUM_QUBITS - 1) {
        std::sort(G.begin(), G.end(), LimEntry<>::greaterValue);
        GaussianEliminationSortedFast(G, nQubits);
        pruneZeroColumns(G);
        std::sort(G.begin(), G.end(), LimEntry<>::greaterValue);
    }

    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline void toColumnEchelonFormModuloPhase(std::vector<LimBitset<NUM_QUBITS, NUM_BITS>>& G, Qubit nQubits = NUM_QUBITS - 1) {
        std::sort(G.begin(), G.end(), LimBitset<NUM_QUBITS, NUM_BITS>::greaterValue);
        GaussianEliminationModuloPhaseSortedFast(G, nQubits);
        pruneZeroColumnsModuloPhase(G);
        std::sort(G.begin(), G.end(), LimBitset<NUM_QUBITS, NUM_BITS>::greaterValue);
    }

    // Reduces a vector 'x' via a group 'G' via the Gram-Schmidt procedure
    // Returns the reduced vector
    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> GramSchmidt(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x) {
        // startProfile(gramSchmidt)
        //        Log::log << "[GramSchmidt] |G|=" << G.size() << "  x = " << LimEntry<>::to_string(x) << "\n"; Log::log.flush();
        LimEntry<NUM_QUBITS> y(x); // = new LimEntry<NUM_QUBITS>(x);
        if (G.empty()) {
            return y;
        }
        constexpr std::size_t height = 2 * NUM_QUBITS;
        for (std::size_t h = 0; h < height; h++) {
            if (y.paulis[h]) {
                //                Log::log << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (std::size_t v = 0; v < G.size(); v++) {
                    if (G[v]->pivotPosition() == h) {
                        //                        Log::log << "[GramSchmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << LimEntry<>::to_string(G[v]) << "\n";
                        y.multiplyBy(*G[v]);
                    }
                }
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    // Reduces a vector 'x' via a group 'G' via the Gram-Schmidt procedure
    // Returns the reduced vector
    // Precondition: the group G is in column echelon form
    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> GramSchmidtFastSorted(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x) {
        // startProfile(gramSchmidt)
        LimEntry<NUM_QUBITS> y(x);
        for (std::size_t g = 0; g < G.size(); g++) {
            auto const pivot = G[g]->pivotPosition();
            if (pivot == std::numeric_limits<decltype(pivot)>::max()) continue;
            if (y.paulis.test(pivot)) {
                y.multiplyBy(*G[g]);
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> GramSchmidtFastSorted(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x, Qubit nQubits) {
        // startProfile(gramSchmidtLazy)
        LimEntry<NUM_QUBITS> y(x);
        for (std::size_t g = 0; g < G.size(); g++) {
            auto const pivot = G[g]->pivotPosition();
            if (pivot == std::numeric_limits<decltype(pivot)>::max()) continue;
            if (y.paulis.test(pivot)) {
                y.multiplyBy(*G[g], nQubits);
            }
        }
        // endProfile(gramSchmidtLazy)
        return y;
    }

    // Reduces a vector 'x' via a group 'G' via the Gram-Schmidt procedure
    // Returns the reduced vector
    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> GramSchmidt(const std::vector<LimEntry<NUM_QUBITS>>& G, const LimEntry<NUM_QUBITS>* x) {
        // startProfile(gramSchmidt)
        //        Log::log << "[GramSchmidt] |G|=" << G.size() << "  x = " << LimEntry<>::to_string(x) << "\n"; Log::log.flush();
        LimEntry<NUM_QUBITS> y(x); // = new LimEntry<NUM_QUBITS>(x);
        if (G.size() == 0) {
            return y;
        }
        constexpr std::size_t height = 2 * NUM_QUBITS;
        for (std::size_t h = 0; h < height; h++) {
            if (y.paulis[h]) {
                //                Log::log << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (std::size_t v = 0; v < G.size(); v++) {
                    if (G[v].pivotPosition() == h) {
                        //                        Log::log << "[GramSchmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << LimEntry<>::to_string(G[v]) << "\n";
                        y.multiplyBy(G[v]);
                    }
                }
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> GramSchmidtFastSorted(const std::vector<LimEntry<NUM_QUBITS>>& G, const LimEntry<NUM_QUBITS>* x, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(gramSchmidt)
        LimEntry<NUM_QUBITS> y(x);
        for (std::size_t g = 0; g < G.size(); g++) {
            auto const pivot = G[g].pivotPosition(nQubits);
            if (pivot == std::numeric_limits<decltype(pivot)>::max()) {
                continue;
            }
            if (y.paulis.test(pivot)) {
                y.multiplyBy(G[g], nQubits);
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline LimBitset<NUM_QUBITS, NUM_BITS> GramSchmidt(const std::vector<LimBitset<NUM_QUBITS, NUM_BITS>*>& G, const LimBitset<NUM_QUBITS, NUM_BITS>* x) {
        // startProfile(gramSchmidt)
        LimBitset<NUM_QUBITS, NUM_BITS> y(x);
        if (G.size() == 0) return y;
        constexpr std::size_t height = 2 * NUM_QUBITS;
        //        Log::log << "[Gram Schmidt] start y = " << y << "\n";
        for (std::size_t h = 0; h < height; h++) {
            if (y.lim.paulis[h]) {
                //                Log::log << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (std::size_t v = 0; v < G.size(); v++) {
                    if (G[v]->lim.pivotPosition() == h) {
                        //                        Log::log << "[Gram Schmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << *G[v] << "\n";
                        y.multiplyBy(*G[v]);
                        //                        Log::log << "[Gram Schmidt] after multiplication, y = " << y << "\n";
                    }
                }
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline LimBitset<NUM_QUBITS, NUM_BITS> GramSchmidt(const std::vector<LimBitset<NUM_QUBITS, NUM_BITS>>& G, const LimBitset<NUM_QUBITS, NUM_BITS>& x) {
        // startProfile(gramSchmidt)
        LimBitset<NUM_QUBITS, NUM_BITS> y(x);
        if (G.size() == 0) return y;
        constexpr std::size_t height = 2 * NUM_QUBITS;
        //        Log::log << "[Gram Schmidt] start y = " << y << "\n";
        for (std::size_t h = 0; h < height; h++) {
            if (y.lim.paulis[h]) {
                //                Log::log << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (std::size_t v = 0; v < G.size(); v++) {
                    if (G[v].lim.pivotPosition() == h) {
                        //                        Log::log << "[Gram Schmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << *G[v] << "\n";
                        y.multiplyBy(G[v]);
                        //                        Log::log << "[Gram Schmidt] after multiplication, y = " << y << "\n";
                    }
                }
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    // Precondition: G is in column echelon form
    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline LimBitset<NUM_QUBITS, NUM_BITS> GramSchmidtFastSorted(const std::vector<LimBitset<NUM_QUBITS, NUM_BITS>>& G, const LimBitset<NUM_QUBITS, NUM_BITS>& x, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(gramSchmidt)
        LimBitset<NUM_QUBITS, NUM_BITS> y(x);
        for (std::size_t g = 0; g < G.size(); g++) {
            auto const pivot = G[g].lim.pivotPosition(nQubits);
            if (y.lim.paulis.test(pivot)) {
                y.multiplyBy(G[g], nQubits);
            }
        }
        // endProfile(gramSchmidt)
        return y;
    }

    // Performs the GramSchmidt algorithm,, i.e.,
    //   given a group G and a vector x,
    //   reduces the vector x via G, and returns this reduced vector
    //   The decomposition that is found, is recorded in the bitset 'indicator'
    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline void GramSchmidt(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x, std::bitset<NUM_BITS>& indicator) {
        // startProfile(gramSchmidt)
        //        Log::log << "[GramSchmidt] |G|=" << G.size() << "  x = " << LimEntry<>::to_string(x) << "\n";
        //        LimEntry<NUM_QUBITS>* y = new LimEntry<NUM_QUBITS>(x);
        LimEntry<NUM_QUBITS> y(x);
        constexpr std::size_t          height = 2 * NUM_QUBITS;
        for (std::size_t i = 0; i < height && i < NUM_QUBITS; i++) {
            indicator.set(i, 0);
        }
        if (G.size() == 0) {
            return;
        }
        for (std::size_t h = 0; h < height; h++) {
            if (y.paulis[h]) {
                //                Log::log << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (std::size_t v = 0; v < G.size(); v++) {
                    if (G[v]->pivotPosition() == h) {
                        //                        Log::log << "[GramSchmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << LimEntry<>::to_string(G[v]) << "\n";
                        //                        y = LimEntry<NUM_QUBITS>::multiply(*G[v], *y);
                        y.multiplyBy(*G[v]);
                        //                        Log::log << "[GramSchmidt] after multiplication, y = " << y << Log::endl;
                        indicator.set(v, 1);
                    }
                }
            }
        }
        // endProfile(gramSchmidt)
    }

    // Given a group G and a 0/1 indicator vector,
    //   returns the product of the indicated elements of G
    //   e.g., with G={ZIZ, IZZ, IXY} and indicator = '110', we return ZZI
    template<std::size_t NUM_QUBITS, std::size_t NUM_BITS>
    inline LimEntry<NUM_QUBITS> getProductOfElements(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::bitset<NUM_BITS>& indicator, const Qubit nQubits = NUM_QUBITS - 1) {
        LimEntry<NUM_QUBITS> g = LimEntry<NUM_QUBITS>();
        assert(G.size() <= NUM_BITS);
        for (std::size_t i = 0; i < G.size(); i++) {
            if (indicator.test(i)) {
                g.multiplyBy(G[i], nQubits);
            }
        }
        return g;
    }

    template<std::size_t NUM_QUBITS>
    inline std::vector<std::bitset<NUM_QUBITS>> getKernelZ([[maybe_unused]] const std::vector<LimEntry<NUM_QUBITS>*>& G) {
        std::vector<std::bitset<NUM_QUBITS>> kernel;
        return kernel;
        //        std::vector<LimBitset<NUM_QUBITS>*> G_Id = appendIdentityMatrixBitset(G);
        //        GaussianElimination(G_Id);
        //        //        Log::log << "[getKernelZ] after Gaussian elimination, G_Id:\n";
        //        //        printStabilizerGroup(G_Id);
        //        for (unsigned int i = 0; i < G_Id.size(); i++) {
        //            if (G_Id[i]->lim.isIdentityOperator()) {
        //                kernel.push_back(G_Id[i]->bits);
        //            }
        //        }
        //        //        Log::log << "[getKernelZ] found kernel:\n";
        //        //        printKernel(kernel);
        //        // free / deallocate G_Id and its elements
        //        for (unsigned int i = 0; i < G_Id.size(); i++) {
        //            delete G_Id[i];
        //        }
        //        return kernel;
    }

    // Returns the kernel of the group G modulo phase, as a vector<bitset>
    // we assume the width of G is at most 2*NUM_QUBITS
    template<std::size_t NUM_QUBITS>
    inline std::vector<std::bitset<2 * NUM_QUBITS>> getKernelModuloPhase(const std::vector<LimEntry<NUM_QUBITS>>& G, const Qubit nQubits = NUM_QUBITS - 1) {
        std::vector<LimBitset<NUM_QUBITS, 2 * NUM_QUBITS>> G_Id = appendIdentityMatrixBitsetBig(G);

        std::sort(G_Id.begin(), G_Id.end(), LimBitset<NUM_QUBITS, 2 * NUM_QUBITS>::greaterValue);
        GaussianEliminationModuloPhaseSortedFast(G_Id, nQubits);
        std::vector<std::bitset<2 * NUM_QUBITS>> kernel;
        for (std::size_t i = 0; i < G_Id.size(); i++) {
            if (G_Id[i].lim.isIdentityModuloPhase()) {
                kernel.push_back(G_Id[i].bits);
            }
        }
        return kernel;
    }

    // Assume: G_Id has been sorted, and then Gaussian elimination was performed. (Zero columns have NOT been pruned)
    template<std::size_t NUM_QUBITS>
    inline std::vector<std::bitset<2 * NUM_QUBITS>> getKernelModuloPhase2(const std::vector<LimEntry<NUM_QUBITS>>& G_Id) {
        //std::vector<LimBitset<NUM_QUBITS, 2 * NUM_QUBITS>> G_Id = appendIdentityMatrixBitsetBig(G);

        //std::sort(G_Id.begin(), G_Id.end(), LimBitset<NUM_QUBITS, 2 * NUM_QUBITS>::greaterValue);
        //GaussianEliminationModuloPhaseSortedFast(G_Id);
        std::vector<std::bitset<2 * NUM_QUBITS>> kernel;
        for (std::size_t i = 0; i < G_Id.size(); i++) {
            if (G_Id[i].lim.isIdentityModuloPhase()) {
                kernel.push_back(G_Id[i].bits);
            }
        }
        return kernel;
    }

    // Given two groups G, H, computes the intersection, <G> intersect <H>
    // TODO refactor with NUM_QUBITS template parameter
    //    inline StabilizerGroup intersectGroupsZ(const StabilizerGroup& G, const StabilizerGroup& H) {
    //        StabilizerGroup                          intersection;
    //        StabilizerGroup                          GH     = groupConcatenate(G, H);
    //        std::vector<std::bitset<dd::NUM_QUBITS>> kernel = getKernelZ(GH);
    //        LimEntry<>                               g;
    //        for (unsigned int i = 0; i < kernel.size(); i++) {
    //            g = getProductOfElements(G, kernel[i]);
    //            intersection.push_back(new LimEntry<NUM_QUBITS>(g));
    //        }
    //        //        Log::log << "[intersectGroupsZ] found intersection:\n";
    //        //        printStabilizerGroup(intersection);
    //        return intersection;
    //    }

    // Returns a generating set J for the intersection of G and H, so <J>= <G> intersect <H>
    //    J is not necessarily in Column Echelon Form
    //    J may contain elements that are equal up to phase
    inline StabilizerGroupValue intersectGroupsModuloPhase(const StabilizerGroup& G, const StabilizerGroupValue& H, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(groupIntersect)
        StabilizerGroupValue                         intersection;
        StabilizerGroupValue                         concat = groupConcatenate(G, H);
        std::vector<std::bitset<2 * dd::NUM_QUBITS>> kernel = getKernelModuloPhase(concat, nQubits);
        //        Log::log << "[intersectGroups mod phase] |kernel| = " << kernel.size() << "\n";
        for (std::size_t i = 0; i < kernel.size(); i++) {
            auto g = getProductOfElements(G, kernel[i], nQubits);
            intersection.push_back(g);
        }
        // endProfile(groupIntersect)
        return intersection;
    }

    inline StabilizerGroupValue intersectGroupsModuloPhaseValue(const StabilizerGroup& G, const StabilizerGroup& H, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(groupIntersect)
        StabilizerGroupValue                         intersection; // TODO reserve some storage? or use std::array instead of std::vector?
        StabilizerGroupValue                         concat = groupConcatenateValue(G, H);
        std::vector<std::bitset<2 * dd::NUM_QUBITS>> kernel = getKernelModuloPhase(concat, nQubits);
        for (unsigned int i = 0; i < kernel.size(); i++) {
            auto g = getProductOfElements(G, kernel[i], nQubits);
            intersection.push_back(g);
        }
        // endProfile(groupIntersect)
        return intersection;
    }

    // TODO (low priority) prevent the use of the oppositePhaseGenerators vector
    //   instead, use some bookkeeping variables to add the appropriately constructed objects to the intersection vector
    //   the purpose is to have less dynamically allocated memory. In this case the use of the vector oppositePhaseGenerators's DAM is prevented
    inline StabilizerGroupValue intersectGroupsPauli(const StabilizerGroup& G, const StabilizerGroupValue& H, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(groupIntersect)
        StabilizerGroupValue intersection = intersectGroupsModuloPhase(G, H, nQubits);
        StabilizerGroupValue oppositePhaseGenerators;
        toColumnEchelonForm(intersection);
        //Log::log << "[intersect groups Pauli] intersection mod phase = " << groupToString(intersection, n) << '\n';
        // Remove all elements from intersection where the G-phase is not equal to the H-phase
        std::size_t g = 0;
        std::size_t intersectionSize = intersection.size();
        for (std::size_t i = 0; i < intersectionSize; i++) {
            auto const phaseG = recoverPhase(G, &intersection[g], nQubits);
            auto const phaseH = recoverPhase(H, &intersection[g], nQubits);
            if (phaseG == phaseH) {
                intersection[g].setPhase(phaseG);
                //Log::log << "[intersect groups Pauli] Adding " << LimEntry<>::to_string(&intersection[g], 3) << " to intersection.\n";
                g++;
            } else {
                // add it to the list of opposite phases
                //Log::log << "[intersect groups Pauli] Element " << LimEntry<>::to_string(&intersection[g], 3) << " has phase(G)=" << phaseToString(phaseG) << " and phaseH=" << phaseToString(phaseH) << ".\n";
                oppositePhaseGenerators.push_back(intersection[g]);
                // remove this from the intersection
                intersection[g] = intersection[intersection.size() - 1];
                intersection.pop_back();
            }
        }
        for (std::size_t i = 1; i < oppositePhaseGenerators.size(); i++) {
            auto a      = LimEntry<>::multiply(oppositePhaseGenerators[0], oppositePhaseGenerators[i]);
            auto const phaseG = recoverPhase(G, &a, nQubits);
            a.setPhase(phaseG);
            intersection.push_back(a);
        }
        //Log::log << "[intersect groups Pauli] intersection = " << groupToString(intersection, n) << '\n';
        // endProfile(groupIntersect)
        return intersection;
    }

    // Precondition: G and H are in column echelon form
    // TODO refactor so that oppositePhaseGenerators are not allocated dynamically but are integrated in the main for-loop
    //    inline StabilizerGroup intersectGroupsPauli(const StabilizerGroup& G, const StabilizerGroup& H) {
    //        //        Qubit n = 5;
    //        //Log::log << "[intersect groups Pauli] G = " << groupToString(G, n) << " H = " << groupToString(H, n) << '\n';
    //        StabilizerGroup intersection = intersectGroupsModuloPhase(G, H);
    //        StabilizerGroup oppositePhaseGenerators;
    //        toColumnEchelonForm(intersection);
    //        //Log::log << "[intersect groups Pauli] intersection mod phase = " << groupToString(intersection, n) << '\n';
    //        // Remove all elements from intersection where the G-phase is not equal to the H-phase
    //        unsigned int s = intersection.size();
    //        phase_t      phaseG, phaseH;
    //        unsigned int g = 0;
    //        for (unsigned int i = 0; i < s; i++) {
    //            phaseG = recoverPhase(G, intersection[g]);
    //            phaseH = recoverPhase(H, intersection[g]);
    //            if (phaseG == phaseH) {
    //                intersection[g]->setPhase(phaseG);
    //                //Log::log << "[intersect groups Pauli] Adding " << LimEntry<>::to_string(intersection[g], 3) << " to intersection.\n";
    //                g++;
    //            } else {
    //                // add it to the wrong phase stuff
    //                //Log::log << "[intersect groups Pauli] Element " << LimEntry<>::to_string(intersection[g], 3) << " has phase(G)=" << phaseToString(phaseG) << " and phaseH=" << phaseToString(phaseH) << ".\n";
    //                oppositePhaseGenerators.push_back(intersection[g]);
    //                // remove this from the intersection
    //                intersection[g] = intersection[intersection.size() - 1];
    //                intersection.pop_back();
    //            }
    //        }
    //        LimEntry<>* a;
    //        for (unsigned int i = 1; i < oppositePhaseGenerators.size(); i++) {
    //            a      = LimEntry<>::multiply(oppositePhaseGenerators[0], oppositePhaseGenerators[i]);
    //            phaseG = recoverPhase(G, a);
    //            a->setPhase(phaseG);
    //            intersection.push_back(a);
    //        }
    //        //Log::log << "[intersect groups Pauli] intersection = " << groupToString(intersection, n) << '\n';
    //        return intersection;
    //    }

    //    inline StabilizerGroup conjugateGroup(const StabilizerGroup& G, const LimEntry<>* a) {
    //        StabilizerGroup H;
    //        for (unsigned int i = 0; i < G.size(); i++) {
    //            H.push_back(new LimEntry<>(G[i]));
    //            if (!H[i]->commutesWith(a)) {
    //                H[i]->setPhase(multiplyPhases(H[i]->getPhase(), phase_t::phase_minus_one));
    //            }
    //        }
    //        return H;
    //    }

    inline StabilizerGroupValue conjugateGroupValue(const StabilizerGroup& G, const LimEntry<>* a) {
        StabilizerGroupValue H;
        for (unsigned int i = 0; i < G.size(); i++) {
            H.push_back(*G[i]);
            if (!H[i].commutesWith(a)) {
                // TODO implement and use a new function in LimEntry which directly multiplies the phase by -1
                H[i].setPhase(multiplyPhases(H[i].getPhase(), phase_t::phase_minus_one));
            }
        }
        return H;
    }

    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> getCosetIntersectionElementModuloPhase(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::vector<LimEntry<NUM_QUBITS>*>& H, const LimEntry<NUM_QUBITS>& a, bool& foundElement, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(cosetIntersectModP)
        std::vector<LimBitset<NUM_QUBITS, 2 * NUM_QUBITS>> GH_Id = concatenateAndAppendIdentityMatrix(G, H);
        toColumnEchelonFormModuloPhase(GH_Id, nQubits);

        std::bitset<NUM_QUBITS>               decomposition; // decomposition of 'a'
        LimBitset<NUM_QUBITS, 2 * NUM_QUBITS> a_bitset(a);
        a_bitset = GramSchmidtFastSorted(GH_Id, a_bitset, nQubits);
        std::bitset<NUM_QUBITS> decomposition_G, decomposition_H; // these bitsets are initialized to 00...0, according to the C++ reference
        bitsetCopySegment(decomposition_G, a_bitset.bits, 0, 0, G.size());
        bitsetCopySegment(decomposition_H, a_bitset.bits, 0, G.size(), G.size() + H.size());
        LimEntry<NUM_QUBITS> a_G = getProductOfElements(G, decomposition_G, nQubits);
        //        Log::log << "[coset intersection P] got first product. Computing second product.\n"; Log::log.flush();
        LimEntry<NUM_QUBITS> a_H     = getProductOfElements(H, decomposition_H, nQubits);
        LimEntry<NUM_QUBITS> a_prime = LimEntry<NUM_QUBITS>::multiply(a_G, a_H, nQubits);
        if (!LimEntry<NUM_QUBITS>::EqualModuloPhase(a, a_prime)) {
            foundElement = false;
        } else {
            foundElement = true;
        }
        // endProfile(cosetIntersectModP)
        return a_G;
    }

    // Given Pauli groups G,H, and Pauli strings a,b, and a phase lambda,
    // Finds an element in the set G intersect lambda a H b,
    // or returns LimEntry::noLIM, if this set is empty
    // TODO refactor to allocate less dynamic memory
    template<std::size_t NUM_QUBITS>
    std::pair<LimEntry<NUM_QUBITS>, bool> getCosetIntersectionElementPauli(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::vector<LimEntry<NUM_QUBITS>*>& H, const LimEntry<NUM_QUBITS>* a, const LimEntry<NUM_QUBITS>* b, phase_t lambda, [[maybe_unused]] const Qubit nQubits = NUM_QUBITS - 1) {
        //std::cout << "helooo\n";
        if (lambda == phase_t::no_phase) {
            return {LimEntry<NUM_QUBITS>(), false};
        }
        // startProfile(cosetIntersectPauli)
        // find an element in G intersect abH modulo phase
        //std::cout << "in getCosetIntersectionElementPauli\n";
        LimEntry<NUM_QUBITS> ab = LimEntry<NUM_QUBITS>::multiply(*a, *b, nQubits);
        bool                 foundCIEMP;
        //std::cout << "in getCosetIntersectionElementPauli\n";
        LimEntry<NUM_QUBITS> c = getCosetIntersectionElementModuloPhase(G, H, ab, foundCIEMP, nQubits);
        //std::cout << "na getCosetIntersectionElementModuloPhase\n";
        if (!foundCIEMP) {
            //            std::cout << "[get coset intersection] Even modulo phase there is no element.\n";
            //            std::cout << "[coset intersection] a = " << LimEntry<>::to_string(a, nQubits) << " b = " << LimEntry<>::to_string(b, nQubits) << " c = " << LimEntry<>::to_string(c, nQubits) << " ab = " << LimEntry<>::to_string(ab, nQubits) << " lambda = " << phaseToString(lambda) << '\n';
            //            std::cout << "[coset intersection] G = " << groupToString(G, nQubits) << "  H = " << groupToString(H, nQubits) << "\n";
            // endProfile(cosetIntersectPauli)
            return {LimEntry<NUM_QUBITS>(), false};
        }
        //std::cout << "coset1 voor recoverphase: (";
        // for (int i = 0; i < int(G.size()); i++){
        //     std::cout << LimEntry<>::to_string(G[i], nQubits) << ", ";
        // }
        c.setPhase(recoverPhase(G, &c));
        //std::cout << "na setPhase\n";
        LimEntry<NUM_QUBITS> acb = LimEntry<NUM_QUBITS>::multiply(*a, c);
        acb                      = LimEntry<NUM_QUBITS>::multiply(acb, *b);
        phase_t alpha            = multiplyPhases(acb.getPhase(), getPhaseInverse(lambda));
        //std::cout << "na multiplyphases\n";
        // Retrieve the phase of acb in H
        // std::cout << "coset2 voor recoverphase: (";
        // for (int i = 0; i < int(H.size()); i++){
        //     std::cout << LimEntry<>::to_string(H[i], nQubits) << ", ";
        // }
        phase_t tau = recoverPhase(H, &acb);
        //std::cout << "na recoverphase\n";
        //Log::log << "[coset intersection] a = " << LimEntry<>::to_string(a, nQubits) << " b = " << LimEntry<>::to_string(b, nQubits) << " c = " << LimEntry<>::to_string(c, nQubits) << " ab = " << LimEntry<>::to_string(ab, nQubits) << " abc = " << LimEntry<>::to_string(acb, nQubits) << " lambda = " << phaseToString(lambda) << " alpha = " << phaseToString(alpha) << " tau = " << phaseToString(tau) << '\n';
        //Log::log << "[coset intersection] G = " << groupToString(G, nQubits) << "  H = " << groupToString(H, nQubits) << "\n";
        if (alpha == tau) {
            // endProfile(cosetIntersectPauli)
            return {c, true};
        }
        // TODO we should just be able to say 'else', because ALWAYS alpha == -tau in this case.
        //    Check if this conjecture is true.
        else if (alpha == multiplyPhases(tau, phase_t::phase_minus_one)) {
            // See if some element of J has xy = -1
            //std::cout << "na multiplyphases2\n";
            std::vector<LimEntry<NUM_QUBITS>> GintersectH = intersectGroupsModuloPhaseValue(G, H);
            //std::cout << "na intersectgroupsmodulophasevalue\n";
            for (std::size_t i = 0; i < GintersectH.size(); i++) {
                if ((!GintersectH[i].commutesWith(b)) ^ (recoverPhase(G, &GintersectH[i]) != recoverPhase(H, &GintersectH[i]))) {
                    // endProfile(cosetIntersectPauli)
                    return {LimEntry<NUM_QUBITS>::multiply(c, recoverElement(G, &GintersectH[i])), true};
                }
            }
        }
        // endProfile(cosetIntersectPauli)
        return {c, false}; // dummy element
    }

    // Given Pauli groups G,H, and Pauli strings a,b, and a phase lambda,
    // Finds an element in the set G intersect lambda a H b,
    // GH_Id is the concatenation of G and H, to which an identity matrix is appended, and on which Gaussian elimination modulo phase has been performed
    // TODO use the data in struct MemoizedData, instead of taking parameters
    template <std::size_t NUM_QUBITS>
    std::pair<LimEntry<NUM_QUBITS>, bool> getCosetIntersectionElementPauli2(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::vector<LimEntry<NUM_QUBITS>*>& H, const LimEntry<NUM_QUBITS>* a, const LimEntry<NUM_QUBITS>* b, phase_t lambda, const std::vector<LimBitset<NUM_QUBITS, 2*NUM_QUBITS> >& GH_Id_CEF, std::vector<LimEntry<NUM_QUBITS>>& GintersectH, bool& memoizedGintersectH, CachingStrategy cachingStrategy, [[maybe_unused]] const  Qubit nQubits = NUM_QUBITS - 1) {
        if (lambda == phase_t::no_phase) {
            return {LimEntry<NUM_QUBITS>(), false};
        }
        // startProfile(cosetIntersectPauli)
        //Log::log << "[cosetIntersection] G = ";
        //printStabilizerGroup(G, nQubits);
        //Log::log << " H = ";
        //printStabilizerGroup(H, nQubits);
        //Log::log << " a=" << a->to_string(nQubits) << " b=" << LimEntry<>::to_string(b, nQubits) << " lambda = " << phaseToString(lambda) << ". Looking for element in G intersect lambda a H b.\n";
        // find an element in G intersect abH modulo phase
        LimEntry<NUM_QUBITS> ab = LimEntry<NUM_QUBITS>::multiply(*a, *b, nQubits);
        bool                 foundCIEMP;
        LimEntry<NUM_QUBITS> c = getCosetIntersectionElementModuloPhase2(G, H, ab, foundCIEMP, GH_Id_CEF, nQubits);
        if (!foundCIEMP) {
            //            std::cout << "[get coset intersection] Even modulo phase there is no element.\n";
            //            std::cout << "[coset intersection] a = " << LimEntry<>::to_string(a, nQubits) << " b = " << LimEntry<>::to_string(b, nQubits) << " c = " << LimEntry<>::to_string(c, nQubits) << " ab = " << LimEntry<>::to_string(ab, nQubits) << " lambda = " << phaseToString(lambda) << '\n';
            //            std::cout << "[coset intersection] G = " << groupToString(G, nQubits) << "  H = " << groupToString(H, nQubits) << "\n";
            // endProfile(cosetIntersectPauli)
            return {LimEntry<NUM_QUBITS>(), false};
        }
        c.setPhase(recoverPhase(G, &c, nQubits));
        LimEntry<NUM_QUBITS> acb = LimEntry<NUM_QUBITS>::multiply(*a, c, nQubits);
        acb                      = LimEntry<NUM_QUBITS>::multiply(acb, *b, nQubits);
        phase_t alpha            = multiplyPhases(acb.getPhase(), getPhaseInverse(lambda));
        // Retrieve the phase of acb in H
        phase_t tau = recoverPhase(H, &acb, nQubits);
        //Log::log << "[coset intersection] a = " << LimEntry<>::to_string(a, nQubits) << " b = " << LimEntry<>::to_string(b, nQubits) << " c = " << LimEntry<>::to_string(c, nQubits) << " ab = " << LimEntry<>::to_string(ab, nQubits) << " abc = " << LimEntry<>::to_string(acb, nQubits) << " lambda = " << phaseToString(lambda) << " alpha = " << phaseToString(alpha) << " tau = " << phaseToString(tau) << '\n';
        //Log::log << "[coset intersection] G = " << groupToString(G, nQubits) << "  H = " << groupToString(H, nQubits) << "\n";
        if (alpha == tau) {
            // endProfile(cosetIntersectPauli)
            //Log::log << "[cosetIntersection] found " << c.to_string(nQubits) << " in G intersect lambda a H b\n";
            return {c, true};
        }
            // TODO we should just be able to say 'else', because ALWAYS alpha == -tau in this case.
            //    Check if this conjecture is true.
        else if (alpha == multiplyPhases(tau, phase_t::phase_minus_one)) {
            if ((!memoizedGintersectH && usingLazyMemoizationGroupIntersect(cachingStrategy)) || !usingLazyMemoizationGroupIntersect(cachingStrategy)) {
                GintersectH = intersectGroupsModuloPhaseValue(G, H, nQubits);
                memoizedGintersectH = true;
            } else {
                //intersectionMemoizationHits++;
            }
            // See if some element of J has xy = -1
            for (std::size_t i = 0; i < GintersectH.size(); i++) {
                if ((!GintersectH[i].commutesWith(b)) ^ (recoverPhase(G, &GintersectH[i], nQubits) != recoverPhase(H, &GintersectH[i], nQubits))) {
                    // endProfile(cosetIntersectPauli)
                    //Log::log << "[cosetIntersection] found " << LimEntry<NUM_QUBITS>::multiply(c, recoverElement(G, &GintersectH[i])).to_string(nQubits) << "\n";
                    return {LimEntry<NUM_QUBITS>::multiply(c, recoverElement(G, &GintersectH[i], nQubits)), true};
                }
            }
        }
        // endProfile(cosetIntersectPauli)
        return {c, false}; // dummy element
    }

    template<std::size_t NUM_QUBITS>
    inline LimEntry<NUM_QUBITS> getCosetIntersectionElementModuloPhase2(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::vector<LimEntry<NUM_QUBITS>*>& H, const LimEntry<NUM_QUBITS>& a, bool& foundElement, const std::vector<LimBitset<NUM_QUBITS, 2*NUM_QUBITS>>& GH_Id_CEF, const Qubit nQubits = NUM_QUBITS - 1) {
        // startProfile(cosetIntersectModP)
        //std::vector<LimBitset<NUM_QUBITS, 2 * NUM_QUBITS>> GH_Id = concatenateAndAppendIdentityMatrix(G, H);
        //toColumnEchelonFormModuloPhase(GH_Id);
        //pruneZeroColumnsModuloPhase(GH_Id);
        //std::sort(GH_Id.begin(), GH_Id.end(), LimBitset<NUM_QUBITS, 2*NUM_QUBITS>::greaterValue);

        std::bitset<NUM_QUBITS>               decomposition; // decomposition of 'a'
        LimBitset<NUM_QUBITS, 2 * NUM_QUBITS> a_bitset(a);
        a_bitset = GramSchmidtFastSorted(GH_Id_CEF, a_bitset, nQubits);
        std::bitset<NUM_QUBITS> decomposition_G, decomposition_H; // these bitsets are initialized to 00...0, according to the C++ reference
        bitsetCopySegment(decomposition_G, a_bitset.bits, 0, 0, G.size());
        bitsetCopySegment(decomposition_H, a_bitset.bits, 0, G.size(), G.size() + H.size());
        LimEntry<NUM_QUBITS> a_G = getProductOfElements(G, decomposition_G, nQubits);
        //        Log::log << "[coset intersection P] got first product. Computing second product.\n"; Log::log.flush();
        LimEntry<NUM_QUBITS> a_H     = getProductOfElements(H, decomposition_H, nQubits);
        LimEntry<NUM_QUBITS> a_prime = LimEntry<NUM_QUBITS>::multiply(a_G, a_H, nQubits);
        if (!LimEntry<NUM_QUBITS>::EqualModuloPhase(a, a_prime)) {
            foundElement = false;
        } else {
            foundElement = true;
        }
        // endProfile(cosetIntersectModP)
        return a_G;
    }



} // namespace dd
#endif //DDPACKAGE_PAULIALGEBRA_HPP
