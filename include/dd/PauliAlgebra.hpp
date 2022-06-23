//
// Created by lieuwe on 13/12/2021.
//

#ifndef DDPACKAGE_PAULIALGEBRA_HPP
#define DDPACKAGE_PAULIALGEBRA_HPP

#include "Edge.hpp"
#include "Nodes.hpp"
#include "LimTable.hpp"
#include "ComplexNumbers.hpp"
#include <iostream>

// note: my package won't compile unless I put my functions in a class
// for now, I've called this class Pauli
// - Lieuwe

namespace dd {

typedef std::vector<LimEntry<>*> StabilizerGroup;

enum LIMDD_group {
	QMDD_group,
	Z_group,
	Pauli_group
};

class Pauli {
public:

	static void superFlush() {
		for (unsigned int i=0; i<1000; i++) {
			std::cout.flush();
		}
	}

    // todo find an appropriate place for this utility function
    template <std::size_t N, std::size_t M>
    static void bitsetCopySegment(std::bitset<N>& x, const std::bitset<M> y, unsigned int begin_x, unsigned int begin_y, unsigned int end_y) {
        assert(end_y <= N && end_y <= M);
        for (unsigned int i=begin_y; i<end_y; i++) {
            x.set(i+begin_x-begin_y, y.test(i));
        }
    }

    // Returns whether the two groups have the same vector of generators
    // Note that, mathematically speaking, two generator sets G and H can still generate the same groups,
    // even if they have different orders
    static bool stabilizerGroupsEqual(const StabilizerGroup& G, const StabilizerGroup& H) {
//        std::cout << "[stabilizerGroupsEqual] start.\n";  std::cout.flush();
        if (G.size() != H.size()) return false;
        for (unsigned int i=0; i<G.size(); i++) {
            if (G[i]->operator!=(*H[i])) return false;
        }
//        std::cout << "[stabilizerGroupsEqual] groups are equal.\n"; std::cout.flush();
        return true;
    }

    // Returns whether the group G is sorted descending,
    // i.e., whether all the zeros are in the top right quadrant
    static bool stabilizerGroupIsSorted(const StabilizerGroup& G) {
        for (unsigned int i=0; i+1<G.size(); i++) {
            if (LimEntry<>::leneq(G[i], G[i+1]))
                return false;
        }
        return true;
    }

    static void printStabilizerGroup(const StabilizerGroup& G) {
        std::cout << "Stabilizer group (" << G.size() << " elements)\n";  std::cout.flush();
        for (unsigned int i=0; i<G.size(); i++) {
            std::cout << LimEntry<>::to_string(G[i]) << std::endl;  std::cout.flush();
        }
        std::cout.flush();
    }

    template <std::size_t NUM_QUBITS>
    static void printStabilizerGroup(const std::vector<LimBitset<NUM_QUBITS>*>& G) {
        if (G.size() == 0) {
            std::cout << "  (empty group)\n";
        }
        for (unsigned int i=0; i<G.size(); i++) {
            std::cout << *G[i] << std::endl; std::cout.flush();
        }
    }

    template <std::size_t NUM_QUBITS>
    static void printKernel(const std::vector<std::bitset<NUM_QUBITS> >& kernel) {
        for (unsigned int i=0; i<kernel.size(); i++) {
            for (unsigned int j=0; j<NUM_QUBITS; j++) {
                std::cout << kernel[i].test(j);
            }
            std::cout << '\n';
        }
    }

    static StabilizerGroup groupConcatenate(const StabilizerGroup& G, const StabilizerGroup& H) {
        StabilizerGroup concat = G;
        for (unsigned int i=0; i<H.size(); i++) {
            concat.push_back(H[i]);
        }
        return concat;
    }

//    static StabilizerGroup deepcopy(const StabilizerGroup& G) {
//        StabilizerGroup copy;
//        for (unsigned int i=0; i<G.size(); i++) {
//            copy.push_back(new LimEntry<>(G[i]));
//        }
//        return copy;
//    }

    template <std::size_t NUM_QUBITS>
    static std::vector<LimBitset<NUM_QUBITS>*> appendIdentityMatrixBitset(const std::vector<LimEntry<NUM_QUBITS>*>& G) {
        std::vector<LimBitset<NUM_QUBITS>*> GI;
        LimBitset<NUM_QUBITS>* col;
        for (unsigned int i=0; i<G.size(); i++) {
            col = new LimBitset<NUM_QUBITS>();
            col->lim = *G[i];
            col->bits.set(i, 1);
            GI.push_back(col);
        }
        return GI;
    }

    // Performs Gaussian elimination on G
    // We assume that G is not stored in the LimTable.
    // In more detail: the elements of G are modified in place
    // Therefore, the LimEntry objects should NOT be stored in the LimTable;
    // todo this algorithm allocates many new LIMs; by modifying in place, less memory can be allocated, and we solve a memory leak
    template <std::size_t NUM_QUBITS>
    static void GaussianElimination(std::vector<LimEntry<NUM_QUBITS>*>& G) {
        if (G.size() <= 1) return;
//        std::cout << "[Gaussian Elimination] start. |G| = " << G.size() << ".\n"; std::cout.flush();
        unsigned int pauli_height = 2*NUM_QUBITS; // length of the columns as far as they contain Pauli operators
        unsigned int reducingColId;
        for (unsigned int h=0; h<pauli_height; h++) {
            // Step 1: Find a column with a '1' at position h
            reducingColId = -1;
            for (unsigned int i=0; i<G.size(); i++) {
                if (G[i]->pivotPosition() == h) {
                    reducingColId = i;
                    break;
                }
            }
            if (reducingColId == (unsigned int)-1) continue;
            // Step 2: Reduce all other columns via G[reducingColId]
            for (unsigned int reduceColId =0; reduceColId < G.size(); reduceColId++) {
                if (reduceColId == reducingColId) continue;
                if (G[reduceColId]->paulis.test(h)) {
//                    std::cout << "[Gaussian Elimination] Multiplying col " << reduceColId << " with col " << reducingColId << ".\n"; std::cout.flush();
                    G[reduceColId] = LimEntry<NUM_QUBITS>::multiply(G[reduceColId], G[reducingColId]);
                }
            }
        }
    }

    template <std::size_t NUM_QUBITS>
    static bool isAbelian(const std::vector<LimBitset<NUM_QUBITS>*>& G) {
    	for (unsigned int i=0; i<G.size(); i++) {
    		for (unsigned int j=i+1; j<G.size(); j++) {
    			if (!G[i]->lim.commutesWith(&(G[j]->lim))) {
    				return false;
    			}
    		}
    	}
    	return true;
    }

    // does not initialize 'decomposition' to an identity matrix
    // TODO don't reallocate so much memory
    // TODO there is a faster version if the group is sorted
    //    (calls to GaussianElimination from getKernel do not sort their group before calling)
    template <std::size_t NUM_QUBITS>
    static void GaussianElimination(std::vector<LimBitset<NUM_QUBITS>*>& G) {
        if (G.size() <= 1) return;
//        std::cout << "[Gaussian Elimination Bitset] start. |G| = " << G.size() << ".\n"; std::cout.flush();
        unsigned int pauli_height = 2*NUM_QUBITS; // length of the columns as far as they contain Pauli operators
        unsigned int reducingColId;
        if (!isAbelian(G)) {
        	// Add -I
        	G.push_back(new LimBitset<NUM_QUBITS>(LimEntry<NUM_QUBITS>::getMinusIdentityOperator()));
        }
        for (unsigned int h=0; h<pauli_height; h++) {
            // Step 1: Find a column whose first '1' is at position h
            reducingColId = -1;
            for (unsigned int i=0; i<G.size(); i++) {
                if (G[i]->lim.pivotPosition() == h) {
                    reducingColId = i;
                    break;
                }
            }
            if (reducingColId == (unsigned int)-1) continue;
            // Step 2: Reduce all other columns via G[reducingColId]
            for (unsigned int reduceColId =0; reduceColId < G.size(); reduceColId++) {
                if (reduceColId == reducingColId) continue;
                if (G[reduceColId]->lim.paulis.test(h)) {
//                    std::cout << "[Gaussian Elimination Bitset] Multiplying col " << reduceColId << " with col " << reducingColId << ".\n"; std::cout.flush();
                    G[reduceColId] = LimBitset<NUM_QUBITS>::multiply(G[reduceColId], G[reducingColId]);
                }
            }
        }
        // the pauli's are done; handle the phase
        reducingColId = -1;
        for (unsigned int i=0; i<G.size(); i++) {
            if (G[i]->lim.isIdentityModuloPhase() && G[i]->lim.getPhase() == phase_t::phase_minus_one) {
                reducingColId = i;
                break;
            }
        }
        if (reducingColId == (unsigned int) -1) return;
        for (unsigned int reduceColId = 0; reduceColId < G.size(); reduceColId++) {
            if (reduceColId == reducingColId) continue;
            if (G[reduceColId]->lim.getPhase() == phase_t::phase_minus_one || G[reduceColId]->lim.getPhase() == phase_t::phase_minus_i) {
                G[reduceColId] = LimBitset<NUM_QUBITS>::multiply(G[reduceColId], G[reducingColId]);
            }
        }
    }

    // Performs Gaussian Elimination on the group G, ignoring the phase of the LIMs involved
    // todo it is possible to write a faster procedure, if we are allowed to assume that G is sorted
    template <std::size_t NUM_QUBITS>
    static void GaussianEliminationModuloPhase(std::vector<LimBitset<NUM_QUBITS>*>& G) {
        if (G.size() <= 1) return;
//        std::cout << "[Gaussian Elimination Bitset] start. |G| = " << G.size() << ".\n"; std::cout.flush();
        unsigned int pauli_height = 2*NUM_QUBITS; // length of the columns as far as they contain Pauli operators
        unsigned int reducingColId;
        for (unsigned int row=0; row<pauli_height; row++) {
            // Step 1: Find a column with a '1' at position row
            reducingColId = -1;
            for (unsigned int i=0; i<G.size(); i++) {
                if (G[i]->lim.pivotPosition() == row) {
                    reducingColId = i;
                    break;
                }
            }
            if (reducingColId == (unsigned int)-1) continue;
            // Step 2: Reduce all other columns via G[reducingColId]
            for (unsigned int reduceColId =0; reduceColId < G.size(); reduceColId++) {
                if (reduceColId == reducingColId) continue;
                if (G[reduceColId]->lim.paulis.test(row)) {
//                    std::cout << "[Gaussian Elimination Bitset] Multiplying col " << reduceColId << " with col " << reducingColId << ".\n"; std::cout.flush();
                    G[reduceColId] = LimBitset<NUM_QUBITS>::multiply(G[reduceColId], G[reducingColId]);
                }
            }
        }
    }

    // TODO this procedure should reduce the refcount of the LimEntry objects it removes from G
    //   or should deallocate these objects in some other way
    template <std::size_t NUM_QUBITS>
    static void pruneZeroColumns(std::vector<LimEntry<NUM_QUBITS>*>& G) {
        unsigned int i=0;
        while (i<G.size()) {
            if (G[i]->isAllZeroVector()) {
                // Remove this all-zero vector from the matrix
                // Step 1: replace it with the last vector in the matrix
                G[i] = G[G.size()-1];
                G.pop_back();
            } else {
                i++;
            }
        }
    }

    template <std::size_t NUM_QUBITS>
    static void pruneZeroColumns(std::vector<LimBitset<NUM_QUBITS>*>& G) {
        unsigned int i=0;
        while (i<G.size()) {
            if (G[i]->lim.isAllZeroVector()) {
                // Remove this all-zero vector from the matrix
                // Step 1: replace it with the last vector in the matrix
                G[i] = G[G.size()-1];
                G.pop_back();
            } else {
                i++;
            }
        }
    }

    template <std::size_t NUM_QUBITS>
    static void pruneZeroColumnsModuloPhase(std::vector<LimBitset<NUM_QUBITS>*>& G) {
        unsigned int i=0;
        while (i<G.size()) {
            if (G[i]->lim.isIdentityModuloPhase()) {
                // Remove this all-zero vector from the matrix
                // Step 1: replace it with the last vector in the matrix
                G[i] = G[G.size()-1];
                G.pop_back();
            } else {
                i++;
            }
        }
    }

    // Puts the stabilizer group in column echelon form; specifically:
    //   1. performs Gaussian elimination on G
    //   2. prunes the all-zero columns
    //   3. sorts the columns lexicographically, i.e., so that 'pivots' appear in the matrix
    static void toColumnEchelonForm(StabilizerGroup& G) {
//        std::cout << "[toColumnEchelonForm] Before CEF, group is:\n"; std::cout.flush();
//        printStabilizerGroup(G);
        GaussianElimination(G);
//        std::cout << "[toColumnEchelonForm] After Gaussian Elimination, before pruning zero col's, group is:\n";std::cout.flush();
//        printStabilizerGroup(G);
        pruneZeroColumns(G);
        // To obtain a lower triangular form, we now sort the vectors descending lexicographically, descending
        std::sort(G.begin(), G.end(), LimEntry<>::geq);
//        std::cout << "[toColumnEchelonForm] After CEF, group is:\n"; std::cout.flush();
//        printStabilizerGroup(G);
    }

    template <std::size_t NUM_QUBITS>
    static void toColumnEchelonForm(std::vector<LimBitset<NUM_QUBITS>*>& G) {
//        std::cout << "[toColumnEchelonForm] start, |G| = " << G.size() << "\n";
        GaussianElimination(G);
        pruneZeroColumns(G);
        std::sort(G.begin(), G.end(), LimBitset<NUM_QUBITS>::geq);
    }

    template <std::size_t NUM_QUBITS>
    static void toColumnEchelonFormModuloPhase(std::vector<LimBitset<NUM_QUBITS>* >& G) {
        GaussianEliminationModuloPhase(G);
        pruneZeroColumnsModuloPhase(G);
        std::sort(G.begin(), G.end(), LimBitset<NUM_QUBITS>::geq);
    }

    // Reduces a vector 'x' via a group 'G' via the Gram-Schmidt procedure
    // Returns the reduced vector
    template<std::size_t NUM_QUBITS>
    static LimEntry<NUM_QUBITS>* GramSchmidt(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x) {
//        std::cout << "[GramSchmidt] |G|=" << G.size() << "  x = " << LimEntry<>::to_string(x) << "\n"; std::cout.flush();
        LimEntry<NUM_QUBITS> y(x); // = new LimEntry<NUM_QUBITS>(x);
        if (G.size() == 0) return new LimEntry<NUM_QUBITS>(y);
        std::size_t height = 2*NUM_QUBITS;
        for (unsigned int h=0; h<height; h++) {
            if (y.paulis[h]) {
//                std::cout << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (unsigned int v=0; v<G.size(); v++) {
                    if (G[v]->pivotPosition() == h) {
//                        std::cout << "[GramSchmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << LimEntry<>::to_string(G[v]) << "\n";
//                        y = LimEntry<NUM_QUBITS>::multiply(*y, *G[v]);
                        y.multiplyBy(*G[v]);
                    }
                }
            }
        }
        return new LimEntry<NUM_QUBITS>(y);
    }

    // todo this algorithm can be sped up if we are allowed to assume that the group G is sorted
    // todo this version uses right multiplication; refactor to left multiplication
    template <std::size_t NUM_QUBITS>
    static LimBitset<NUM_QUBITS> GramSchmidt(const std::vector<LimBitset<NUM_QUBITS>*>& G, const LimBitset<NUM_QUBITS>* x) {
//        LimBitset<NUM_QUBITS>* y = new LimBitset<NUM_QUBITS>(x);
        LimBitset<NUM_QUBITS> y(x);
        if (G.size() == 0) return y;
        std::size_t height = 2*NUM_QUBITS;
//        std::cout << "[Gram Schmidt] start y = " << y << "\n";
        for (unsigned int h=0; h<height; h++) {
            if (y.lim.paulis[h]) {
//                std::cout << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (unsigned int v=0; v<G.size(); v++) {
                    if (G[v]->lim.pivotPosition() == h) {
//                        std::cout << "[Gram Schmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << *G[v] << "\n";
                        y.multiplyBy(*G[v]);
//                        std::cout << "[Gram Schmidt] after multiplication, y = " << y << "\n";
                    }
                }
            }
        }
        return y;
    }

    // Performs the GramSchmidt algorithm,, i.e.,
    //   given a group G and a vector x,
    //   reduces the vector x via G, and returns this reduced vector
    //   The decomposition that is found, is recorded in the bitset 'indicator'
    // todo this algorithm can be sped up if the group G is sorted
    template <std::size_t NUM_QUBITS>
    static void GramSchmidt(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x, std::bitset<NUM_QUBITS>& indicator) {
//        std::cout << "[GramSchmidt] |G|=" << G.size() << "  x = " << LimEntry<>::to_string(x) << "\n";
//        LimEntry<NUM_QUBITS>* y = new LimEntry<NUM_QUBITS>(x);
        LimEntry<NUM_QUBITS> y(x);
        std::size_t height = 2*NUM_QUBITS;
        for (unsigned int i=0; i<height && i<NUM_QUBITS; i++) {
            indicator.set(i, 0);
        }
        if (G.size() == 0) return;
        for (unsigned int h=0; h<height; h++) {
            if (y.paulis[h]) {
//                std::cout << "[GramSchmidt] h=" << h << ".\n";
                // Look for a vector whose first '1' entry is at position h
                for (unsigned int v=0; v<G.size(); v++) {
                    if (G[v]->pivotPosition() == h) {
//                        std::cout << "[GramSchmidt] found '1' in G[" << v << "][" << h << "]; multiplying by " << LimEntry<>::to_string(G[v]) << "\n";
//                        y = LimEntry<NUM_QUBITS>::multiply(*G[v], *y);
                        y.multiplyBy(*G[v]);
//                        std::cout << "[GramSchmidt] after multiplication, y = " << y << std::endl;
                        indicator.set(v, 1);
                    }
                }
            }
        }
//        return y;
    }


    // Given a group G and a 0/1 indicator vector,
    //   returns the product of the indicated elements of G
    //   e.g., with G={ZIZ, IZZ, IXY} and indicator = '110', we return ZZI
    template<std::size_t NUM_QUBITS>
    static LimEntry<NUM_QUBITS>* getProductOfElements(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::bitset<NUM_QUBITS>& indicator) {
        LimEntry<NUM_QUBITS>* g = LimEntry<NUM_QUBITS>::getIdentityOperator();
        assert(G.size() <= NUM_QUBITS);
        for (unsigned int i=0; i<G.size(); i++) {
            if (indicator.test(i)) {
                g->multiplyBy(G[i]);
            }
        }
        return g;
    }

    // TODO free / deallocate G_Id and its elements
    template <std::size_t NUM_QUBITS>
    static std::vector<std::bitset<NUM_QUBITS> > getKernelZ(const std::vector<LimEntry<NUM_QUBITS>* >& G) {
        std::vector<LimBitset<NUM_QUBITS>* > G_Id = appendIdentityMatrixBitset(G);
        GaussianElimination(G_Id);
//        std::cout << "[getKernelZ] after Gaussian elimination, G_Id:\n";
//        printStabilizerGroup(G_Id);
        std::vector<std::bitset<NUM_QUBITS> > kernel;
        for (unsigned int i=0; i<G_Id.size(); i++) {
            if (G_Id[i]->lim.isIdentityOperator()) {
                kernel.push_back(G_Id[i]->bits);
            }
        }
        // TODO free / deallocate G_Id and its elements
//        std::cout << "[getKernelZ] found kernel:\n";
//        printKernel(kernel);
        return kernel;
    }

    // Returns the kernel of the group G modulo phase, as a vector<bitset>
    // TODO free / deallocate G_Id and its elements
    template <std::size_t NUM_QUBITS>
    static std::vector<std::bitset<NUM_QUBITS> > getKernelModuloPhase(const std::vector<LimEntry<NUM_QUBITS>* >& G) {
        std::vector<LimBitset<NUM_QUBITS>* > G_Id = appendIdentityMatrixBitset(G);
        GaussianEliminationModuloPhase(G_Id);
        std::vector<std::bitset<NUM_QUBITS> > kernel;
        for (unsigned i=0; i<G_Id.size(); i++) {
            if (G_Id[i]->lim.isIdentityModuloPhase()) {
                kernel.push_back(G_Id[i]->bits);
            }
        }
        // TODO free / deallocate G_Id and its elements
        return kernel;
    }

    // Given two groups G, H, computes the intersection, <G> intersect <H>
    static StabilizerGroup intersectGroupsZ(const StabilizerGroup& G, const StabilizerGroup& H) {
        StabilizerGroup intersection;
        StabilizerGroup GH = groupConcatenate(G, H);
        std::vector<std::bitset<32> > kernel = getKernelZ(GH);
        LimEntry<>* g;
        for (unsigned int i=0; i<kernel.size(); i++) {
            g = getProductOfElements(G, kernel[i]);
//            g = LimEntry<>::getIdentityOperator();
//            for (unsigned int j=0; j<G.size(); j++) {
//                if (kernel[i].test(j)) {
//                    g->multiplyBy(*G[j]);
//                }
//            }
            intersection.push_back(g);
        }
//        std::cout << "[intersectGroupsZ] found intersection:\n";
//        printStabilizerGroup(intersection);
        return intersection;
    }

    // TODO Not yet implemented
    static StabilizerGroup intersectGroupsPauli(const StabilizerGroup& G, const StabilizerGroup& H) {
        return intersectGroupsZ(G, H);
    }

    // Returns a generating set J for the intersection of G and H, so <J>= <G> intersect <H>
    //    J is not necessarily in Column Echelon Form
    //    J may contain elements that are equal up to phase
    // todo verify: does this method work for Pauli groups? it certainly works for Z groups
    // TODO refactor loop body to getProductOfElements
    static StabilizerGroup intersectGroupsModuloPhase(const StabilizerGroup& G, const StabilizerGroup& H) {
//        std::cout << "[intersectGroups mod phase] start  |G|=" << G.size() << "  |H| = " << H.size() << std::endl;
//        std::cout << "[intersectGroups mod phase] Group G:\n";
//        printStabilizerGroup(G);
//        std::cout << "[intersectGroups mod phase] Group H:\n";
//        printStabilizerGroup(H);
        StabilizerGroup intersection;
        StabilizerGroup concat = groupConcatenate(G, H);
        std::vector<std::bitset<32> > kernel = getKernelModuloPhase(concat);
//        std::cout << "[intersectGroups mod phase] |kernel| = " << kernel.size() << "\n";
        LimEntry<>* g;
        for (unsigned int i=0; i<kernel.size(); i++) {
        	g = getProductOfElements(G, kernel[i]);
//            g = LimEntry<>::getIdentityOperator();
//            for (unsigned int j=0; j<G.size(); j++) {
//                if (kernel[i].test(j)) {
//                    g->multiplyBy(*G[j]);
//                }
//            }
            intersection.push_back(g);
        }
//        std::cout << "[intersectGroups mod phase] Found intersection: \n";
//        printStabilizerGroup(intersection);

        return intersection;
    }

    // Given sets G, H which generate Pauli groups <G> and <H>, respectively, and a Pauli string a,
    // Searches for an element in the set '<G> intersect (<H>+a)'
    // Note that '<G>' is a group, and '<H>+a' is a coset
    // If the intersection of these two sets is empty, 'LimEntry<..>::noLIM' is returned
    template <std::size_t NUM_QUBITS>
    static LimEntry<NUM_QUBITS>* getCosetIntersectionElementPauli(const std::vector<LimEntry<NUM_QUBITS>*>& G, const std::vector<LimEntry<NUM_QUBITS>*>& H, const LimEntry<NUM_QUBITS>* a) {
//        std::cout << "[coset intersection P] a = " << *a << "\n"; std::cout.flush();
        std::vector<LimEntry<NUM_QUBITS>*> GH = groupConcatenate(G, H);
        assert(GH.size() <= NUM_QUBITS); // todo fix me; Concatenating two StabilizerGroups of size NUM_QUBITS can result in a StabilizerGroups > NUM_QUBITS. This causes an error in line appendIdentityMatrixBitset
        std::vector<LimBitset<NUM_QUBITS>*> GH_Id = appendIdentityMatrixBitset(GH);
        toColumnEchelonFormModuloPhase(GH_Id);
//        std::cout << "[coset intersection P] Found GH +Id to column echelon mod phase:\n";
//        printStabilizerGroup(GH_Id);
        std::bitset<NUM_QUBITS> decomposition;   // decomposition of 'a'
//        std::cout << "[coset intersection P] Doing Gram-Schmidt with GH, a.\n"; std::cout.flush();
        LimBitset<NUM_QUBITS> a_bitset(a);
        // todo refactor this to the GramSchmidt(Group, LimEntry, std::bitset) version instead of the GramSchmidt(Group, LimBitset) version
        a_bitset = GramSchmidt(GH_Id, &a_bitset);
//        std::cout << "[coset intersection P] after Gram-Schmidt, a_bitset = " << a_bitset << "\n";
//        std::cout << "[coset intersection P] a.bits[0] = " << a_bitset.bits.test(0) << ", a.bits[1] = " << a_bitset.bits.test(1) << "\n";
//        std::cout << "[coset intersection P] computing product of elements (GH, decomposition).\n"; std::cout.flush();
//        LimEntry<NUM_QUBITS>* c = getProductOfElements(GH, decomposition);
//        std::cout << "[coset intersection P] Computed product of elements.\n"; std::cout.flush();
        std::bitset<NUM_QUBITS> decomposition_G, decomposition_H;  // these bitsets are initialized to 00...0, according to the C++ reference
//        std::cout << "[coset intersection P] copying bit segment 1...\n"; std::cout.flush();
        bitsetCopySegment(decomposition_G, a_bitset.bits, 0, 0, G.size());
        bitsetCopySegment(decomposition_H, a_bitset.bits, 0, G.size(), G.size()+H.size());
//        bitsetCopySegment<NUM_QUBITS, 2*NUM_QUBITS+2>(decomposition_G, c->paulis, 0, 2*nqubits, 2*nqubits+G.size());
//        std::cout << "[coset intersection P] copying bit segment 2...\n"; std::cout.flush();
//        bitsetCopySegment(decomposition_H, c->paulis, 0, 2*nqubits + G.size(), 2*nqubits+G.size() + H.size());
//        std::cout << "[coset intersection P] copy successful! Getting product of elements.\n"; std::cout.flush();
        LimEntry<NUM_QUBITS>* a_G = getProductOfElements(G, decomposition_G);
//        std::cout << "[coset intersection P] got first product. Computing second product.\n"; std::cout.flush();
        LimEntry<NUM_QUBITS>* a_H = getProductOfElements(H, decomposition_H);
        LimEntry<NUM_QUBITS>* a_prime = LimEntry<NUM_QUBITS>::multiply(a_G, a_H);
        phase_t phase_diff = (phase_t) ((a->getPhase() + a_prime->getPhase()) & 0x3);
//        std::cout << "[coset intersection P] a_G = " << *a_G << "\n";
//        std::cout << "[coset intersection P] a_H = " << *a_H << "\n";
//        std::cout << "[coset intersection P] decomposition G = " << decomposition_G << "\n";
//        std::cout << "[coset intersection P] decomposition H = " << decomposition_H << "\n";
//        std::cout << "[coset intersection P] aprime = " << *a_prime << "\n"; std::cout.flush();
        // todo replace this check with the more 'elegant' checking whether the vector returned by GramSchmidt, above, is identity modulo phase
        //   that has the advantage that the check can be done earlier
        if (!LimEntry<NUM_QUBITS>::EqualModuloPhase(a, a_prime)) {
            return LimEntry<NUM_QUBITS>::noLIM;
        }
        if (phase_diff == phase_t::phase_one) {
            return a_G;
        }
        else if (phase_diff == phase_t::phase_i || phase_diff == phase_t::phase_minus_i) {
            return LimEntry<NUM_QUBITS>::noLIM;
        }
        // phase_diff must be -1  (i.e., must be phase_t::phase_minus_one)
//        std::cout << "[coset intersection P] Phase difference is -1. Computing intersection...\n"; std::cout.flush();
        const std::vector<LimEntry<NUM_QUBITS>*> intersection = intersectGroupsModuloPhase(G, H);
        LimEntry<NUM_QUBITS>* k_G;
        LimEntry<NUM_QUBITS>* k_H;
        for (unsigned int i=0; i<intersection.size(); i++) {
            GramSchmidt(G, intersection[i], decomposition_G);
            GramSchmidt(H, intersection[i], decomposition_H);
            k_G = getProductOfElements(G, decomposition_G);
            k_H = getProductOfElements(H, decomposition_H);
            if (((k_G->getPhase() + k_H->getPhase()) & 0x3) == phase_t::phase_minus_one) {
                a_G->multiplyBy(k_G);
                return a_G;
            }
        }
//        std::cout << "[coset intersection P] no element in coset; returning noLIM.\n"; std::cout.flush();
        return LimEntry<NUM_QUBITS>::noLIM;
    }

    // We assume that only vNodes are passed
    // todo deallocate minus, m
    static StabilizerGroup constructStabilizerGeneratorSetZ(const vNode node) {
        //empty
        Edge<vNode> low, high;
        low  = node.e[0];
        high = node.e[1];
        unsigned int n = node.v;

        StabilizerGroup stabgenset;
        // Case 0: Check if this node is the terminal node (aka the Leaf)
        if (n == (unsigned int)-1) {
            // Return the trivial group.
            // This group is generated by the empty set; therefore, we just return the empty stabgenset
            return stabgenset;
        }
        // Case 1: right child is zero
        else if (high.isZeroTerminal()) {
//            std::cout << "[stab genZ] |0> knife case  n = " << n << ". Low stabilizer group is:\n";
            stabgenset = low.p->limVector; // copies the stabilizer group of the left child
//            printStabilizerGroup(stabgenset);
            LimEntry<>* idZ = LimEntry<>::getIdentityOperator();
            idZ->setOperator(n, 'Z');
            stabgenset.push_back(idZ);
//            std::cout << "[stab genZ] Added Z. Now stab gen set is:\n";
//            printStabilizerGroup(stabgenset);
            // the matrix set is already in column echelon form,
            // so we do not need to perform that step here
        }
        // Case 2: left child is zero
        else if (low.isZeroTerminal()) {
//            std::cout << "[stab genZ] |1> knife case. n = " << n << ". High stabilizer group is:\n";
            stabgenset = high.p->limVector;
//            printStabilizerGroup(stabgenset);
            LimEntry<>* minusIdZ = LimEntry<>::getMinusIdentityOperator();
            minusIdZ->setOperator(n, 'Z');
            stabgenset.push_back(minusIdZ);
//            std::cout << "[stab genZ] Added -Z. now stab gen set is:\n";
//            printStabilizerGroup(stabgenset);
        }
        // Case 3: the node is a 'fork': both its children are nonzero
        else {
            // Gather the stabilizer groups of the two children
//            std::cout << "[constructStabilizerGeneratorSet] Case fork n = " << n << ".\n";
            StabilizerGroup* stabLow  = &(low. p->limVector);
            StabilizerGroup* stabHigh = &(high.p->limVector);
            // Step 1: Compute the intersection
            stabgenset = intersectGroupsZ(*stabLow, *stabHigh);

            // Step 2: if some element v is in the set <G> intersect (<H> * -I),
            //   then add Z tensor v to the stabgenset
            LimEntry<>* minus = LimEntry<>::getMinusIdentityOperator();
            LimEntry<>* m = getCosetIntersectionElementPauli(*stabLow, *stabHigh, minus);
            if (m != LimEntry<>::noLIM) {
                m->setOperator(n, 'Z');
                stabgenset.push_back(m);
            }
            toColumnEchelonForm(stabgenset);
            // todo deallocate minus, m
        }

        return stabgenset;
    }

    // Construct the stabilizer generator set of 'node' in the Pauli group
	// TODO limdd store stab in LimTable
    static StabilizerGroup constructStabilizerGeneratorSetPauli(const vNode node) {
        Edge<vNode> low, high;
        low  = node.e[0];
        high = node.e[1];
        unsigned int n = node.v;

        StabilizerGroup stabgenset;
        // Case 0: Check if this node is the terminal node (aka the Leaf)
        if (n == (unsigned int)-1) { // TODO replace with a direct check whether 'node' is a terminal node
            // Return the trivial group.
            // This group is generated by the empty set; therefore, we just return the empty stabgenset
            return stabgenset;
        }
        // Case 1: right child is zero
        else if (high.isZeroTerminal()) {
            std::cout << "[stab genPauli] |0> knife case  n = " << n << ". Low stabilizer group is:\n";
            stabgenset = low.p->limVector; // copies the stabilizer group of the left child
            printStabilizerGroup(stabgenset);
            LimEntry<>* idZ = LimEntry<>::getIdentityOperator();
            idZ->setOperator(n, 'Z');
            stabgenset.push_back(idZ);
            std::cout << "[stab genPauli] Added Z. Now stab gen set is:\n";
            printStabilizerGroup(stabgenset);
            // the matrix set is already in column echelon form,
            // so we do not need to perform that step here
        }
        // Case 2: left child is zero
        else if (low.isZeroTerminal()) {
            std::cout << "[stab genPauli] |1> knife case. n = " << n << ". High stabilizer group is:\n";
            stabgenset = high.p->limVector; // copy the stabilizer of the right child
            printStabilizerGroup(stabgenset);
            LimEntry<>* minusIdZ = LimEntry<>::getMinusIdentityOperator();
            minusIdZ->setOperator(n, 'Z');
            stabgenset.push_back(minusIdZ);
            std::cout << "[stab genPauli] Added -Z. now stab gen set is:\n";
            printStabilizerGroup(stabgenset);
        }
        // Case 3: the node is a 'fork': both its children are nonzero
        else {
            // Gather the stabilizer groups of the two children
            std::cout << "[constructStabilizerGeneratorSet] Case fork n = " << n << ".\n";
			StabilizerGroup* stabLow  = &(low. p->limVector);
			StabilizerGroup* stabHigh = &(high.p->limVector);
			// Step 1: Compute the intersection
			stabgenset = intersectGroupsZ(*stabLow, *stabHigh);
			// Step 2: find out whether an element P*P' should be added, where P acts on qubit 'n'
			LimEntry<>* stab = LimEntry<>::noLIM;
            if (low.p == high.p) {
            	// Step 2.1: Find out if we should add X*P' or Y*P'
            	Complex highWeight = high.w;
            	if (highWeight.approximatelyEquals(Complex::one)) {
            		stab = new LimEntry<>(high.l);
            		stab->setOperator(n, 'X');
            	}
            	else if (highWeight.approximatelyEquals(Complex::minus_one)) {
            		stab = new LimEntry<>(high.l);
            		stab->setOperator(n, 'X');
            		stab->setPhase(phase_t::phase_minus_one);
            	}
            	else if (highWeight.approximatelyEquals(Complex::complex_i)) {
            		stab = new LimEntry<>(high.l);
            		stab->setOperator(n, 'Y');
            	}
            	else if (highWeight.approximatelyEquals(Complex::minus_i)) {
            		stab = new LimEntry<>(high.l);
            		stab->setOperator(n, 'Y');
            		stab->setPhase(phase_t::phase_minus_one);
            	}
            } else {
				// Step 2.2: if some element v is in the set <G> intersect (<H> * -I),
				//   then add Z tensor v to the stabgenset
				LimEntry<>* minus = LimEntry<>::getMinusIdentityOperator();
				stab = getCosetIntersectionElementPauli(*stabLow, *stabHigh, minus);
				if (stab != LimEntry<>::noLIM) {
					stab->setOperator(n, 'Z');
				}
				delete minus;
            }
			if (stab != LimEntry<>::noLIM) {
				stabgenset.push_back(stab);
			}
			toColumnEchelonForm(stabgenset);
        }

        return stabgenset;
    }

    // Returns an isomorphism between u and v,
    // or LimEntry<>::noLim if u and v are not isomorphic
    // Assumes that the low edges of u and v have an Identity LIM
    // TODO should we add assertions that u and v do not represent zero vectors?
    // TODO this function does not take into account the different phases... but maybe it doesn't need to...
    static LimEntry<>* getIsomorphismZ(const vNode* u, const vNode* v) {
        assert( u != nullptr );
        assert( v != nullptr );
//        std::cout << "[getIsomorphismZ] Start.\n";
        LimEntry<>* iso;
//         TODO add assertion that the nodes are on the same number of qubits u->v == v->v
//        assert (u->v == v->v);
        Edge<vNode> uLow  = u->e[0];
        Edge<vNode> uHigh = u->e[1];
        Edge<vNode> vLow  = v->e[0];
        Edge<vNode> vHigh = v->e[1];
        assert (!(uLow.isZeroTerminal() && uHigh.isZeroTerminal()));
        assert (!(vLow.isZeroTerminal() && vHigh.isZeroTerminal()));
        // TODO this assertion is not necessarily true; in the normalizeLIMDD function, we hve vLow.l != nullptr
        assert (uLow.l == nullptr && vLow.l == nullptr);
        // Case 0.1: the nodes are equal
        if (u == v) {
//            std::cout << "[getIsomorphismZ] case u == v.\n"; std::cout.flush();
            // In this case, we return the Identity operator, which is represented by a null pointer
            return nullptr;
        }
        // Case 0.2: The leaf case.
        // TODO this case should already be covered by case 0.1, since in this case v is also the terminal node
        //   Do we need this extra check?
        else if (vNode::isTerminal(u)) {
//            std::cout << "[getIsomorphismZ] Case u is terminal.\n"; std::cout.flush();
            // Return the identity operator, which is represented by a null pointer
            return nullptr;
        }
        // Case 1 ("Left knife"): Left child is nonzero, right child is zero
        else if (uHigh.isZeroTerminal()) {
//            std::cout << "[getIsomorphismZ] Case uHigh is terminal\n";
            if (!vHigh.isZeroTerminal()) return LimEntry<>::noLIM;
            if (uHigh.p != vHigh.p) return LimEntry<>::noLIM;
            return LimEntry<>::multiply(*uHigh.l, *vHigh.l);
        }
        // Case 2 ("Right knife"): Left child is zero, right child is nonzero
        else if (uLow.isZeroTerminal()) {
//            std::cout << "[getIsomorphismZ] case uLow is terminal.\n";
            if (!vLow.isZeroTerminal()) return LimEntry<>::noLIM; // not isomorphic
            if (uLow.p != vLow.p) return LimEntry<>::noLIM;
            return nullptr;  // return the Identity isomorphism
        }
        // Case 3 ("Fork"): Both children are nonzero
        else {
//            std::cout << "[getIsomorphismZ] case Fork.\n"; std::cout.flush();
//            std::cout << "[getIsomorphismZ] ulw " << uLow.w << " uhw " << uHigh.w << " vlw " << vLow.w << " vhw " << vHigh.w << std::endl;
            // Step 1.2: check if the amplitudes satisfy uHigh = -1 * vHigh
            bool amplitudeOppositeSign = uHigh.w.approximatelyEqualsMinus(vHigh.w);
            // Step 1.1:  check if the amplitudes are equal, up to a sign
            if (!uLow.w.approximatelyEquals(vLow.w) || (!uHigh.w.approximatelyEquals(vHigh.w) && !amplitudeOppositeSign)) return LimEntry<>::noLIM;
//            std::cout << "[getIsomorphismZ] edge weights are approximately equal.\n"; std::cout.flush();
            // Step 2: Check if nodes u and v have the same children
            if (uLow.p != vLow.p || uHigh.p != vHigh.p) return LimEntry<>::noLIM;
//            std::cout << "[getIsomorphismZ] children of u and v are the same nodes.\n"; std::cout.flush();
            // Step 3: (optional) check if the automorphism groups are equal
//            if (!stabilizerGroupsEqual(u->limVector, v->limVector)) {
//                return LimEntry<>::noLIM;
//            }
//            std::cout << "[getIsomorphismZ] the stabilizer Groups of u and v are equal.\n"; std::cout.flush();
            // Step 4: If G intersect (H+isoHigh) contains an element P, then Id tensor P is an isomorphism
            LimEntry<>* isoHigh = LimEntry<>::multiply(uHigh.l, vHigh.l);
//            std::cout << "[getIsomorphismZ] multiplied high isomorphisms:" << LimEntry<>::to_string(isoHigh) << ".\n"; std::cout.flush();
            if (amplitudeOppositeSign) {
//                std::cout << "[getIsomorphismZ] multiplying Phase by -1 because high amplitudes had opposite signs\n"; std::cout.flush();
                isoHigh->multiplyPhaseBy(phase_t::phase_minus_one); // multiply by -1
            }
            iso = getCosetIntersectionElementPauli(uLow.p->limVector, uHigh.p->limVector, isoHigh);
//            std::cout << "[getIsomorphismZ] completed coset intersection element.\n"; std::cout.flush();
            if (iso != LimEntry<>::noLIM) {
//                std::cout << "[getIsomorphismZ] The coset was non-empty; returning element.\n"; std::cout.flush();
                return iso;
            }
            // Step 5: If G intersect (H-isomorphism) contains an element P, then Z tensor P is an isomorphism
//            std::cout << "[getIsomorphismZ] multiplying phase by -1.\n"; std::cout.flush();
            isoHigh->multiplyPhaseBy(phase_t::phase_minus_one);
            iso = getCosetIntersectionElementPauli(uLow.p->limVector, uHigh.p->limVector, isoHigh);
//            std::cout << "[getIsomorphismZ] found coset intersection element.\n"; std::cout.flush();
            if (iso != LimEntry<>::noLIM) {
//                std::cout << "[getIsomorphismZ] Coset was not empty; returning result.\n"; std::cout.flush();
                iso->setOperator(u->v, pauli_op::pauli_z); // TODO should we do this? write a test
                return iso;
            }
            else {
//                std::cout << "[getIsomorphismZ] Coset was empty; returning -1.\n"; std::cout.flush();
                return LimEntry<>::noLIM;
            }
        }
    }

    // Assumes that u and v are semi-reduced:
    // - low edge label is identity
    // TODO take the edge weights into account:
    //    in case 3.1
    //    in knife cases
    //    check if uhigh.w = 1 / vhigh.w
    static LimWeight<>* getIsomorphismPauli(const vNode* u, const vNode* v, ComplexNumbers& cn) {
        assert( u != nullptr );
        assert( v != nullptr );
        std::cout << "[getIsomorphismPauli] Start. states have " << (int) u->v+1 << " qubits.\n";
        assert (u->v == v->v);  // Assert u and v have the same nubmer of qubits
        Edge<vNode> uLow  = u->e[0];
        Edge<vNode> uHigh = u->e[1];
        Edge<vNode> vLow  = v->e[0];
        Edge<vNode> vHigh = v->e[1];
        // Assert that neither u nor v is the Zero vector
        assert (!(uLow.isZeroTerminal() && uHigh.isZeroTerminal()));
        assert (!(vLow.isZeroTerminal() && vHigh.isZeroTerminal()));
        std::cout << "[getIsomorphismPauli] uLow .l = " << LimEntry<>::to_string(uLow.l) << " vLow.l = " << LimEntry<>::to_string(vLow.l) << std::endl;
        std::cout << "[getIsomorphismPauli] uHigh.l = " << LimEntry<>::to_string(uHigh.l) << " vHigh.l = " << LimEntry<>::to_string(vHigh.l) << std::endl;
        assert (LimEntry<>::isIdentityOperator(uLow.l));
        assert (LimEntry<>::isIdentityOperator(vLow.l));
        LimWeight<>* iso = new LimWeight<>();
        // Case 0: the nodes are equal
        if (u == v) {
            std::cout << "[getIsomorphismPauli] case u == v.\n"; std::cout.flush();
            // In this case, we return the Identity operator, which is represented by a null pointer
            iso = new LimWeight<>((LimEntry<>*)nullptr);
        }
        // Case 1 ("Left knife"): Left child is nonzero, right child is zero
        else if (uHigh.isZeroTerminal()) {
            std::cout << "[getIsomorphismPauli] Case |u> = |0>|u'>, since uHigh is zero\n";
            if (vHigh.isZeroTerminal()) {
            	if (uLow.p == vLow.p) iso = new LimWeight<>((LimEntry<>*)nullptr);
            	else iso = LimWeight<>::noLIM;
            }
            else if (vLow.isZeroTerminal()) {
            	if (uLow.p == vHigh.p) {
					// TODO limdd inspect weight on high edge
            		iso->lim = new LimEntry<>(vHigh.l);
            		iso->lim->setOperator(u->v, 'X');
            	}
            	else iso = LimWeight<>::noLIM;
            }
            else {
            	iso = LimWeight<>::noLIM;
            }
        }
        // Case 2 ("Right knife"): Left child is zero, right child is nonzero
        else if (uLow.isZeroTerminal()) {
            std::cout << "[getIsomorphismPauli] case uLow is zero, so |u> = |1>|u'>.\n";
        	if (vLow.isZeroTerminal()) {
        		// TODO limdd inspect weights
        		if (uHigh.p == vHigh.p) return new LimWeight<>(LimEntry<>::multiply(uHigh.l, vHigh.l));
        	}
        	else if (vHigh.isZeroTerminal()) {
        		// TODO limdd inspect weights
        		if (uHigh.p == vLow.p) {
					iso->lim = new LimEntry<>(uHigh.l);
					iso->lim->setOperator(u->v, 'X');
        		}
        	}
        	else iso = LimWeight<>::noLIM;
        }
        // Case 3 ("Fork"): Both children are nonzero
        else {
        	// Case 3.1: uLow == vHigh, uHigh == vLow but uLow != uHigh, i.e., the isomorphism's first Pauli operator is an X or Y
        	// TODO by handling case 3.1 more efficiently, we can prevent unnecessary copying of u->limVector
			if (uLow.p == vHigh.p && uHigh.p == vLow.p && uLow.p != uHigh.p) {
				// Return lambda^-1 * R * (X tensor P), where
				//    P is the uHigh's edge label
				//    lambda is uHigh's weight
				//    R is an isomorphism between uPrime and v
				std::cout << "[getIsomorphismPauli] case 3.1: children of nodes are opposite pair. Qubits: " << (int)(u->v) << "\n";
				// TODO refactor this piece of code which swaps two edges
				vNode uPrime;
				uPrime.v = u->v;
				uPrime.limVector = u->limVector;
				uPrime.e[0]   = u->e[1];
				uPrime.e[0].l = nullptr;
				uPrime.e[0].w = u->e[1].w;
				uPrime.e[1]   = u->e[0];
				uPrime.e[1].l = u->e[1].l;
				uPrime.e[1].w = u->e[0].w;
				LimWeight<>* R = getIsomorphismPauli(&uPrime, v, cn);
				if (R == LimWeight<>::noLIM) return LimWeight<>::noLIM;
				LimEntry<> P = *(u->e[1].l);
				P.setOperator(u->v, pauli_op::pauli_x);
				R->multiplyBy(P);
				return R;
			}
        	// Case 3.2: uLow == vLow and uHigh == vHigh
            std::cout << "[getIsomorphismPauli] case Fork with weights ulw " << uLow.w << " uhw " << uHigh.w << " vlw " << vLow.w << " vhw " << vHigh.w << std::endl; superFlush();
            // Step 1.1: Check if uLow == vLow and uHigh == vHigh, i.e., check if nodes u and v have the same children
            if (uLow.p != vLow.p || uHigh.p != vHigh.p) return LimWeight<>::noLIM;
            std::cout << "[getIsomorphismPauli] children of u and v are the same nodes.\n"; std::cout.flush(); superFlush();
			// TODO should we refactor this last part and just call getIsomorphismZ?
			//      we could refactor ONLY this last part, and thereby make both this and the getIsomorphismZ functions more readable
            // Step 1.2: check if the weights satisfy uHigh = -1 * vHigh
            if (uLow.p == uHigh.p) {
            	// Then the high weights can be uHigh = 1 / vHigh or uHigh = - 1 / vHigh
            }

//            Complex div0 = cn.getTemporary();  // Eventually returned to cache
//            ComplexNumbers::div(div0, u->e[1].w, u->e[0].w);
//            Complex div1 = cn.getTemporary();  // Eventually returned to cache
//            ComplexNumbers::div(div1, v->e[1].w, v->e[0].w);
//            iso->weight  = cn.getTemporary();  // We ask the calling function to return this to cache
//            ComplexNumbers::div(iso->weight, u->e[0].w, v->e[1].w);

//            bool amplitudeOppositeSign = div0.approximatelyEqualsMinus(div1);
            bool edgeWeightsOppositeSign = u->e[1].w.approximatelyEqualsMinus(v->e[1].w);
//            bool edgeWeightsPlusMinus  = div0.approximatelyEqualsPlusMinus(div1);
            bool edgeWeightsPlusMinus  = u->e[1].w.approximatelyEquals(v->e[1].w) || edgeWeightsOppositeSign;
//            cn.returnToCache(div0);
//            cn.returnToCache(div1);
            // Step 1.3:  check if the edge weights are equal, up to a sign
            // TODO limdd check if uhigh.w = 1 / vhigh.w
            if (!edgeWeightsPlusMinus) return LimWeight<>::noLIM;
            std::cout << "[getIsomorphismPauli] edge weights are approximately equal.\n"; std::cout.flush(); superFlush();
            // Step 2: If G intersect (H+isoHigh) contains an element P, then Id tensor P is an isomorphism
            LimEntry<>* isoHigh = LimEntry<>::multiply(uHigh.l, vHigh.l);
            std::cout << "[getIsomorphismPauli] multiplied high isomorphisms:" << LimEntry<>::to_string(isoHigh) << ".\n"; std::cout.flush(); superFlush();
            if (edgeWeightsOppositeSign) {
                std::cout << "[getIsomorphismPauli] multiplying Phase by -1 because high weights had opposite signs\n"; std::cout.flush(); superFlush();
                isoHigh->multiplyPhaseBy(phase_t::phase_minus_one); // multiply by -1
            }
            iso->lim = getCosetIntersectionElementPauli(uLow.p->limVector, uHigh.p->limVector, isoHigh);
            std::cout << "[getIsomorphismPauli] completed coset intersection element.\n"; std::cout.flush();
            if (iso->lim != LimEntry<>::noLIM) {
                std::cout << "[getIsomorphismPauli] The coset was non-empty; returning element.\n"; std::cout.flush();
                return iso;
            }
            // Step 3: If G intersect (H-isomorphism) contains an element P, then Z tensor P is an isomorphism
            std::cout << "[getIsomorphismPauli] multiplying phase by -1.\n"; std::cout.flush();
            isoHigh->multiplyPhaseBy(phase_t::phase_minus_one);
            iso->lim = getCosetIntersectionElementPauli(uLow.p->limVector, uHigh.p->limVector, isoHigh);
            if (iso->lim != LimEntry<>::noLIM) {
                iso->lim->setOperator(u->v, pauli_op::pauli_z);
                std::cout << "[getIsomorphismPauli] Coset was not empty; returning result.\n"; std::cout.flush();
            }
            else {
                std::cout << "[getIsomorphismPauli] Coset was empty; returning -1.\n"; std::cout.flush();
//                cn.returnToCache(iso->weight);
                iso = LimWeight<>::noLIM;  // TODO limdd I added this, is this right? -LV
            }
        }
        return iso;
    }

    // Choose the label on the High edge, in the Z group
    static LimEntry<>* highLabelZ(const vNode* u, const vNode* v, LimEntry<>* vLabel, Complex& weight, bool& s) {
        // We assert that the LIM has phase +1  (we expect normalizeLIMDD to guarantee this)
        assert(LimEntry<>::getPhase(vLabel) == phase_t::phase_one);
//        std::cout << "[highLabelZWeight] Start; |Gu| = " << u->limVector.size() << " |Gv| = " << v->limVector.size() << ".\n"; std::cout.flush();
        StabilizerGroup GH = groupConcatenate(u->limVector, v->limVector);
//        std::cout << "[highLabelZWeight] Concatenated; |GH| = " << GH.size() << std::endl;
        toColumnEchelonForm(GH);
//        std::cout << "[highLabelZWeight] to CEF'ed; now |GH| = " << GH.size() << std::endl;
        LimEntry<>* newHighLabel = GramSchmidt(GH, vLabel);
        // Set the new phase to +1
        newHighLabel->setPhase(phase_t::phase_one);
        s = false;
        if (weight.lexLargerThanxMinusOne()) {
//            std::cout << "[highLabelZWeight] Multiplying weight by -1, since weight = " << weight << ".\n";
            weight.multiplyByMinusOne(false);
            s = true;
        }
//        std::cout << "[highLabelZWeight] end.\n"; std::cout.flush();
        return newHighLabel;
    }

    // Finds a high label for a node with low child u and high child v, with current high edge label vLabel, and current high LIM vLabel
    // Sets the
    // Here we demand that 'weight' and 'weightInv' are retrieved with ComplexTable.getTemporary(..),
    // since they will be assigned values but will not be looked up in the ComplexTable
    // TODO limdd:
    //   1. make NUM_QUBITS a template parameter
    //   2. Communicate that we need a new reduction rule for high edge weights. This procedure now implements a candidate for that.
    static LimEntry<>* highLabelPauli(const vNode* u, const vNode* v, LimEntry<>* vLabel, Complex& weight, bool& s, bool& x) {
    	LimEntry<>* newHighLabel;
    	if (u == v) {
    		newHighLabel = GramSchmidt(u->limVector, vLabel);

    		if (CTEntry::val(weight.r) < 0) {
    			weight.multiplyByMinusOne(true);
    			s = true;
    			std::cout << "[highLabelPauli] the high edge weight is flipped, so setting s:=true. New weight is " << weight << ".\n";
    		}
    		fp norm = ComplexNumbers::mag2(weight);
    		if (norm > 1) {
    			ComplexNumbers::div(weight, Complex::one, weight);
    			x = true;
    		}

//    		if (weight.lexSmallerThanxMinusOne()) {
//    			weight.multiplyByMinusOne(true);
//    			s = true;
//    		}
//    		ComplexNumbers::div(weightInv, Complex::one, weight); // temp := 1/weight
//    		if (weightInv.lexSmallerThanxMinusOne()) {
//    			weightInv.multiplyByMinusOne(true);
//    			sInv = true;
//    		}
//    		if (weightInv.lexSmallerThan(weight)) {
//    			weight.setVal(weightInv);
//    			x = true;
//    			s = sInv;
//    		}
    	}
    	else {
    		StabilizerGroup GH = groupConcatenate(u->limVector, v->limVector);
    		toColumnEchelonForm(GH);
    		newHighLabel = GramSchmidt(GH, vLabel);
    		if (weight.lexSmallerThanxMinusOne()) {
    			weight.multiplyByMinusOne(true);
    			s = true;
    		}
    		x = false;
    	}
    	newHighLabel->setPhase(phase_t::phase_one);
    	return newHighLabel;
    }

    static LimEntry<>* highLabelZ(const mNode* u, const mNode* v, LimEntry<>* vLabel) {
    	std::cout << "[highLabelZ] Called with matrix node. Throwing exception.\n" << u << v << vLabel << std::endl;
        throw std::exception();
    }

    static LimEntry<>* highLabelPauli(const mNode* u, const mNode* v, LimEntry<>* vLabel) {
        std::cout << "[highLabelPauli] Called with matrix node. Throwing exception.\n" << u << v << vLabel << std::endl;
        throw std::exception();
    }

    static LimEntry<>* getIsomorphismZ(const mNode* u, const mNode* v) {
        std::cout << "[getIsomorphismZ] called with matrix nodes. Throwing exception.\n" << u << v << std::endl;
        throw std::exception();
    }

    static LimEntry<>* getIsomorphismPauli(const mNode* u, const mNode* v) {
        std::cout << "[getIsomorphismPauli] called with matrix nodes. Throwing exception.\n" << u << v << std::endl;
        throw std::exception();
    }

    // Returns the lexicographically smallest LIM R such that R * |v> == lim * |v>
    // This is useful when a canonical edge is needed for a cache entry
    // TODO in Pauli LIMDD, we need to right-multiply the LIM here; whereas in other applications we need a left-multiplication
    //    make sure the left and right-handed multiplications go well
    static LimEntry<>* getRootLabel(const vNode* v, const LimEntry<>* lim) {
    	return GramSchmidt(v->limVector, lim);
    }

};



} // namespace dd
#endif //DDPACKAGE_PAULIALGEBRA_HPP
