/*
 * This file is part of the JKQ DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_COMPLEX_HPP
#define DD_PACKAGE_COMPLEX_HPP

#include "ComplexTable.hpp"
#include "ComplexValue.hpp"

#include <cstddef>
#include <iostream>
#include <utility>

namespace dd {
    using CTEntry = ComplexTable<>::Entry;

    struct Complex {
        CTEntry* r;
        CTEntry* i;

        static Complex zero;
        static Complex one;

        void setVal(const Complex& c) const {
            r->value = CTEntry::val(c.r);
            i->value = CTEntry::val(c.i);
        }

        [[nodiscard]] inline bool approximatelyEquals(const Complex& c) const {
            return CTEntry::approximatelyEquals(r, c.r) && CTEntry::approximatelyEquals(i, c.i);
        };

        [[nodiscard]] inline bool approximatelyZero() const {
            return CTEntry::approximatelyZero(r) && CTEntry::approximatelyZero(i);
        }

        [[nodiscard]] inline bool approximatelyOne() const {
            return CTEntry::approximatelyOne(r) && CTEntry::approximatelyZero(i);
        }

        static inline bool approximatelyEqual(const std::complex<fp>& x, const std::complex<fp>& y) {
            return CTEntry::approximatelyEquals(x.real(), y.real()) && CTEntry::approximatelyEquals(x.imag(), y.imag());
        }

        static inline bool approximatelyZero(const std::complex<fp> x) {
            return CTEntry::approximatelyZero(x.real()) && CTEntry::approximatelyZero(x.imag());
        }

        inline bool operator==(const Complex& other) const {
            return r == other.r && i == other.i;
        }

        inline bool operator!=(const Complex& other) const {
            return !operator==(other);
        }

        // TODO this function is only used in debugging; we could move it to a more suitable place
        static Complex minusOne() {
            Complex min = one;
            min.r = ComplexTable<>::Entry::flipPointerSign(min.r);
            return min;
        }

        // Sets this complex value z to '-z'
        void multiplyByMinusOne() {
            r = CTEntry::flipPointerSign(r);
            i = CTEntry::flipPointerSign(i);
        }

        [[nodiscard]] std::string toString(bool formatted = true, int precision = -1) const {
            return ComplexValue::toString(CTEntry::val(r), CTEntry::val(i), formatted, precision);
        }

        // todo this function is required once we implement Pauli functionality; now we only have <Z> functionality
//        bool lexSmallerThan(const Complex& other) const {
//        }

        // Returns whether z=a+bi is lexicographically smaller than -z
        //   which is true iff
        //   b > 0, or b=0 and a>0
        //   i.e., iff z != 0 and z=r e^(it) with 0 <= t < pi
        bool lexSmallerThanxMinusOne() const {
            if (!ComplexTable<>::Entry::approximatelyZero(i)) {
                std::cout << "[lexSmallerThanxMinusOne] imag is not zero in " << toString() << "\n";
//                return !ComplexTable<>::Entry::isNegativePointer(i);
                return ComplexTable<>::Entry::val(i) > 0;
            }
            if (!ComplexTable<>::Entry::approximatelyZero(r)) {
                std::cout << "[lexSmallerThanxMinusOne] real is not zero in " << toString() << "\n";
//                return !ComplexTable<>::Entry::isNegativePointer(r);
                return ComplexTable<>::Entry::val(r) > 0;
            }
            return false;
        }

        void writeBinary(std::ostream& os) const {
            CTEntry::writeBinary(r, os);
            CTEntry::writeBinary(i, os);
        }
    };

    inline std::ostream& operator<<(std::ostream& os, const Complex& c) {
        return os << c.toString();
    }

    inline Complex Complex::zero{&ComplexTable<>::zero, &ComplexTable<>::zero};
    inline Complex Complex::one{&ComplexTable<>::one, &ComplexTable<>::zero};
} // namespace dd

namespace std {
    template<>
    struct hash<dd::Complex> {
        std::size_t operator()(dd::Complex const& c) const noexcept {
            auto h1 = dd::murmur64(reinterpret_cast<std::size_t>(c.r));
            auto h2 = dd::murmur64(reinterpret_cast<std::size_t>(c.i));
            return dd::combineHash(h1, h2);
        }
    };
} // namespace std

#endif //DD_PACKAGE_COMPLEX_HPP
