//
// Created by lieuwe on 13/12/2021.
//

#ifndef DDPACKAGE_PAULIALGEBRA_HPP
#define DDPACKAGE_PAULIALGEBRA_HPP

#include "Edge.hpp"
#include "Nodes.hpp"
#include "LimTable.hpp"

#include <iostream>

// note: my package won't compile unless I put my functions in a class
// for now, I've called this class Pauli
// - Lieuwe

namespace dd {

typedef std::vector<LimEntry<>*> StabilizerGroup;

class Pauli {
public:

    // Returns whether the two groups have the same vector of generators
    // Note that, mathematically speaking, two generator sets G and H can still generate the same groups,
    // even if they have different orders
    static bool stabilizerGroupsEqual(const StabilizerGroup& G, const StabilizerGroup& H) {
        if (G.size() != H.size()) return false;
        for (unsigned int i=0; i<G.size(); i++) {
            if (G[i]->operator!=(*H[i])) return false;
        }
        return true;
    }

    // Returns whether the group G is sorted descending,
    // i.e., whether all the zeros are in the top right quadrant
    static bool stabilizerGroupIsSorted(const StabilizerGroup& G) {
        for (unsigned int i=0; i<G.size()-1; i++) {
            if (LimEntry<>::leneq(G[i], G[i+1]))
                return false;
        }
        return true;
    }

    static void printStabilizerGroup(const StabilizerGroup& G) {
        std::cout << "Stabilizer group (" << G.size() << " elements)\n";
        for (unsigned int i=0; i<G.size(); i++) {
            std::cout << LimEntry<>::to_string(G[i]) << std::endl;
        }
        std::cout.flush();
    }

    static StabilizerGroup groupConcatenate(const StabilizerGroup& G, const StabilizerGroup& H) {
        StabilizerGroup concat = G;
        for (unsigned int i=0; i<H.size(); i++) {
            concat.push_back(H[i]);
        }
        return concat;
    }

    // This function does nothing; it just 'catches' the mNode case
//    static LimEntry<>* highLabelZ(const mNode u, const mNode v, LimEntry<>* vLabel) {
//        return vlabel;
//    }

    // Returns the LIM a * b
    static LimEntry<>* multiply(const LimEntry<>& a, const LimEntry<>& b) {
        LimEntry<>* c = LimEntry<>::getIdentityOperator();
        c->multiplyBy(a);
        c->multiplyBy(b);
        return c;
    }

    static StabilizerGroup deepcopy(const StabilizerGroup& G) {
        StabilizerGroup copy;
        for (unsigned int i=0; i<G.size(); i++) {
            copy.push_back(new LimEntry<>(G[i]));
        }
        return copy;
    }

    // Appends an Identity matrix below the group G, i.e., returns the matrix
    //   [ G  ]
    //   [ Id ]
    static StabilizerGroup appendIdentityMatrix(const StabilizerGroup& G, unsigned int height) {
        StabilizerGroup GI = deepcopy(G);
        for (unsigned int i=0; i<G.size(); i++) {
            // Make column G[i][height ... height + width-1] into the column [0 ...0 1 0 ... 0]
            for (unsigned int h=height; h<height + i; h++) {
                GI[i]->paulis.set(h, 0);
            }
            GI[i]->paulis.set(height+i, 1);
            for (unsigned int h=height+i+1; h<G.size(); h++) {
                GI[i]->paulis.set(h, 0);
            }
        }
        return GI;
    }

    // Performs Gaussian elimination on G
    // We assume that G is not stored in the LimTable.
    // In more detail: the elements of G are modified in place
    // Therefore, the LimEntry objects should NOT be stored in the LimTable;
    template <std::size_t NUM_QUBITS>
    static void GaussianElimination(const std::vector<LimEntry<NUM_QUBITS>*>& G) {
        if (G.size() <= 1) return;
        unsigned int pauli_height = 2*NUM_QUBITS; // length of the columns as far as they contain Pauli operators
        unsigned int reducingColId;
        for (unsigned int h=0; h<pauli_height; h++) {
            // Step 1: Find a column with a '1' at position h
            reducingColId = -1;
            for (unsigned int i=0; i<G.size(); i++) {
                if (G[i]->paulis.test(h)) {
                    reducingColId = i;
                    break;
                }
            }
            if (reducingColId == (unsigned int)-1) continue;
            // Step 2: Reduce all other columns via G[reducingColId]
            for (unsigned int reduceColId =0; reduceColId < G.size(); reduceColId++) {
                if (reduceColId == reducingColId) continue;
                if (G[reduceColId]->paulis.test(h)) {
                    G[reduceColId]->multiplyBy(*G[reducingColId]);
                }
            }
        }
    }

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

    // Puts the stabilizer group in column echelon form; specifically:
    //   1. performs Gaussian elimination on G
    //   2. prunes the all-zero columns
    //   3. sorts the columns lexicographically, i.e., so that 'pivots' appear in the matrix
    static void toColumnEchelonForm(StabilizerGroup& G) {
        GaussianElimination(G);
        pruneZeroColumns(G);
        // To obtain a lower triangular form, we now sort the vectors descending lexicographically, descending
        std::sort(G.begin(), G.end(), LimEntry<>::leq);
    }

    // Reduces a vector 'x' via a group 'G' via the Gram-Schmidt procedure
    // Returns the reduced vector
    template<std::size_t NUM_QUBITS>
    static LimEntry<NUM_QUBITS>* GramSchmidt(const std::vector<LimEntry<NUM_QUBITS>*>& G, const LimEntry<NUM_QUBITS>* x) {
        LimEntry<NUM_QUBITS>* y = new LimEntry<NUM_QUBITS>(x);
        if (G.size() == 0) return y;
        std::size_t height = 2*NUM_QUBITS;
        for (unsigned int h=0; h<height; h++) {
            // Look for a vector with a '1' in place h
            for (unsigned int v=0; v<G.size(); v++) {
                if (G[v]->paulis.test(h)) {
                    y = Pauli::multiply(*y, *G[v]);
                }
            }
        }
        return y;
    }

    // Interprets G as a 0/1 matrix, where each operator (i.e., LimEntry object) forms a column
    // returns the kernel of G in column echelon form
    // TODO this method requires LimEntry objects of length 4*NUM_QUBITS,
    //    but the LimEntry objects only give it vectors of length 2*NUM_QUBITS
    //    so right now, we can 'only' do simulations of up to 15 qubits
    static StabilizerGroup getKernelZ(const StabilizerGroup& G, unsigned int nqubits) {
        StabilizerGroup kernel;
        unsigned int width = G.size();  // minor note: size of G may change during toColumnEchelonForm(G)
        StabilizerGroup GI = appendIdentityMatrix(G, nqubits);
        toColumnEchelonForm(GI);
        for (unsigned int i=0; i<GI.size(); i++) {
            // Copy bits [2*nqubits ... 2*nqubits + width -1] to kernel[i][0 ... width-1]
            kernel.push_back(new LimEntry<>());
            for (unsigned int b=0; b<width; b++) {
                kernel[i]->paulis.set(b, GI[i]->paulis.test(2*nqubits+b));
            }
        }
        toColumnEchelonForm(kernel);
        return kernel;
    }

    static StabilizerGroup intersectGroupsZ(const StabilizerGroup& G, const StabilizerGroup& H, unsigned int nqubits) {
        StabilizerGroup intersection;
        StabilizerGroup concat = groupConcatenate(G, H);
        StabilizerGroup kernel = getKernelZ(concat, nqubits);
        LimEntry<>* g;
        for (unsigned int i=0; i<kernel.size(); i++) {
            g = LimEntry<>::getIdentityOperator();
            for (unsigned int j=0; j<G.size(); j++) {
                if (kernel[i]->paulis.test(j)) {
                    g->multiplyBy(*G[j]);
                }
            }
            intersection.push_back(g);
        }

        return intersection;
    }

    // Returns an element in the coset <G> intersect (b + <H>)
    // Here <G> is the group generated by G; likewise for H
    // Note that b + <H> is therefore a coset.
    // If the coset <G> intersect (b + <H>) is empty, nullptr is returned
    // TODO implement this function
    static LimEntry<>* getCosetIntersectionElementZ(const StabilizerGroup& G, const StabilizerGroup& H, const LimEntry<>* b) {
        if (!stabilizerGroupIsSorted(G)) return nullptr;
        if (!stabilizerGroupIsSorted(H)) return nullptr;

        LimEntry<>* b2;
        LimEntry<>* bOrth;
        b2 = new LimEntry<>(b);
        b2 = GramSchmidt(H, b2);
        bOrth = new LimEntry<>(b2);
        bOrth = GramSchmidt(H, bOrth);
        if (bOrth->isIdentityOperator()) {
            bOrth->multiplyBy(*b2);
        }
        return nullptr;
//        if (!G.isSorted()) return nullptr; // TODO throw exception
//        if (!H.isSorted()) return nullptr; // TODO replace by assertions so that in Release, this doesn't happen?
//
    }

    // For now, we assume that only vNodes are passed
//    template <class Node>
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
            stabgenset = low.p->limVector; // copies the stabilizer group of the left child
            LimEntry<>* idZ = LimEntry<>::getIdentityOperator();
            idZ->setOperator(n, 'Z');
            stabgenset.push_back(idZ);
            // the matrix set is already in column echelon form,
            // so we do not need to perform that step here
        }
        // Case 2: left child is zero
        else if (low.isZeroTerminal()) {
            stabgenset = high.p->limVector;
            LimEntry<>* minusIdZ = LimEntry<>::getMinusIdentityOperator();
            minusIdZ->setOperator(n, 'Z');
            stabgenset.push_back(minusIdZ);
        }
        // Case 3: the node is a 'fork': both its children are nonzero
        else {
            // Gather the stabilizer groups of the two children
            StabilizerGroup* stabLow  = &(low. p->limVector);
            StabilizerGroup* stabHigh = &(high.p->limVector);
            // Step 1: Compute the intersection
            stabgenset = intersectGroupsZ(*stabLow, *stabHigh, n);

            // Step 2:
            LimEntry<>* minus = LimEntry<>::getMinusIdentityOperator();
            LimEntry<>* m = getCosetIntersectionElementZ(*stabLow, *stabHigh, minus);
            if (m != nullptr) {
                m->setOperator(n, 'Z');
                stabgenset.push_back(m);
            }
            // The matrix is now guaranteed to be in column echelon form
        }

        return stabgenset;
    }

    // Returns an isomorphism between u and v,
    // or -1 if u and v are not isomorphic
    // Assumes that the low edges of u and v have an Identity LIM
    // TODO should we add assertions that u and v do not represent zero vectors?
    // TODO this function does not take into account the different phases... but maybe it doesn't need to...
    static LimEntry<>* getIsomorphismZ(const vNode* u, const vNode* v) {
        LimEntry<>* iso = nullptr;
//         TODO add assertion that the nodes are on the same number of qubits u->v == v->v
//        assert (u->v == v->v);
        Edge<vNode> uLow  = u->e[0];
        Edge<vNode> uHigh = u->e[1];
        Edge<vNode> vLow  = v->e[0];
        Edge<vNode> vHigh = v->e[1];
        assert (!(uLow.isZeroTerminal() && uHigh.isZeroTerminal()));
        assert (!(vLow.isZeroTerminal() && vHigh.isZeroTerminal()));
        assert (uLow.l == nullptr && vLow.l == nullptr);
        // Case 0.1: the nodes are equal
        if (u == v) {
            // In this case, we return the Identity operator, which is represented by a null pointer
            return nullptr;
        }
        // Case 0.2: The leaf case.
        // TODO this case should already be covered by case 0.1, since in this case v is also the terminal node
        //   Do we need this extra check?
        else if (vNode::isTerminal(u)) {
            // Return the identity operator, which is represented by a null pointer
            return nullptr;
        }
        // Case 1 ("Left knife"): Left child is nonzero, right child is zero
        else if (uHigh.isZeroTerminal()) {
            if (!vHigh.isZeroTerminal()) return (LimEntry<>*)-1;
            if (uHigh.p != vHigh.p) return (LimEntry<>*) -1;
            return multiply(*uHigh.l, *vHigh.l);
        }
        // Case 2 ("Right knife"): Left child is zero, right child is nonzero
        else if (uLow.isZeroTerminal()) {
            if (!vLow.isZeroTerminal()) return (LimEntry<>*)-1; // not isomorphic
            if (uLow.p != vLow.p) return (LimEntry<>*) -1;
            return nullptr;  // return the Identity isomorphism
        }
        // Case 3 ("Fork"): Both children are nonzero
        else {
            // Step 1: check if the amplitudes are equal, up to a sign
            if (!uLow.w.approximatelyEquals(vLow.w) || !uHigh.w.approximatelyEquals(vHigh.w)) return (LimEntry<>*) -1;
            // TODO check if the amplitudes satisfy uHigh = -1 * vHigh
            bool amplitudeOppositeSign = false;
            // Step 2: Check if nodes u and v have the same children
            if (uLow.p != vLow.p || uHigh.p != vHigh.p) return (LimEntry<>*) -1;
            // Step 3: check if the automorphism groups are equal
            if (!stabilizerGroupsEqual(u->limVector, v->limVector)) {
                return (LimEntry<>*) -1;
            }
            // Step 4: If G intersect (H+isoHigh) contains an element P, then Id tensor P is an isomorphism
            LimEntry<>* isoHigh = multiply(*uHigh.l, *vHigh.l);
            if (amplitudeOppositeSign) {
                isoHigh->multiplyPhaseBy(2); // multiply by -1
            }
            iso = getCosetIntersectionElementZ(u->limVector, v->limVector, isoHigh);
            if (iso != nullptr) {
                return iso;
            }
            // Step 5: If G intersect (H-isomorphism) contains an element P, then Z tensor P is an isomorphism
            isoHigh->multiplyPhaseBy(2);
            iso = getCosetIntersectionElementZ(u->limVector, v->limVector, isoHigh);
            if (iso != nullptr) {
                return iso;
            }
            else
                return (LimEntry<>*) -1;
        }

        return iso;
    }

    // Choose the label on the High edge, in the Z group
//    template <class Node>
    static LimEntry<>* highLabelZ(const vNode* u, const vNode* v, LimEntry<>* vLabel) {
        StabilizerGroup GH = groupConcatenate(u->limVector, v->limVector);
        toColumnEchelonForm(GH);
        LimEntry<>* newHighLabel = GramSchmidt(GH, vLabel);
        // Set the new label's "-1 bit" to 'true'
        // This is the vector's last bit
        newHighLabel->paulis.set(LimEntry<>::NUM_BITSETBITS-1, 1);
        return newHighLabel;
    }

    static LimEntry<>* highLabelZ(const mNode* u, const mNode* v, LimEntry<>* vLabel) {
        throw std::exception();
    }

    static LimEntry<>* getIsomorphismZ(const mNode* u, const mNode* v) {
        throw std::exception();
    }

};



} // namespace dd
#endif //DDPACKAGE_PAULIALGEBRA_HPP
