#include <cmath>

namespace Details {
int gcdExpanded(int first, int second, int& coefficientFirst, int& coefficientSecond);
int gcd(int first, int second);
unsigned int EulerFunction(unsigned int number);

template <unsigned long long N, unsigned long long left_bound, unsigned long long diapason>
struct SquareRootInterval {
    static const unsigned int value = SquareRootInterval<N,
            (left_bound + diapason / 2) * (left_bound + diapason / 2) <= N
                    ? left_bound + diapason / 2
                    : left_bound,
            (left_bound + diapason / 2) * (left_bound + diapason / 2) <= N
                    ? diapason - diapason / 2
                    : diapason / 2>::value;
};

template <unsigned long long N, unsigned long long left_bound>
struct SquareRootInterval<N, left_bound, 1> {
    static const unsigned int value = left_bound;
};

template <unsigned int N>
struct SquareRoot {
    static const unsigned int value = SquareRootInterval<N, 1, N>::value;
};

template <unsigned int N, unsigned int prime>
struct IsPrimeDegreeHelper {
    static const bool value = IsPrimeDegreeHelper<
            N % prime == 0 ? N / prime : 0, prime>::value;
};

template <unsigned int prime>
struct IsPrimeDegreeHelper<1, prime> {
    static const bool value = true;
};

template <unsigned int prime>
struct IsPrimeDegreeHelper<0, prime> {
    static const bool value = false;
};

template <unsigned int N, unsigned int D>
struct IsPrimeHelper {
    static const bool value = N % D == 0 ? false : IsPrimeHelper<N, D - 1>::value;
    static const unsigned int firstDivisor =
            value ? N
                  : (IsPrimeHelper<N, D - 1>::value
                        ? D
                        : IsPrimeHelper<N, D - 1>::firstDivisor);
};

template <unsigned int N>
struct IsPrimeHelper<N, 1> {
    static const bool value = true;
    static const unsigned int firstDivisor = N;
};

template <bool isCompile>
struct CheckBoolCompile {
    static void check() {}
};

template <>
struct CheckBoolCompile<false> {
    static void check() = delete;
};
} //end of namespace Details

template <unsigned int N>
class Residue {
private:
    long long remainder = 0;
    static const unsigned int _phi;
public:
    explicit Residue(int number);
    explicit operator int() const;
    explicit operator bool() const;

    Residue<N>& operator +=(const Residue<N>& another);
    Residue<N>& operator *=(const Residue<N>& another);
    Residue<N>& operator -=(const Residue<N>& another);
    Residue<N>& operator /=(const Residue<N>& another);

    Residue<N> operator +(const Residue<N>& another) const;
    Residue<N> operator -(const Residue<N>& another) const;
    Residue<N> operator *(const Residue<N>& another) const;
    Residue<N> operator /(const Residue<N>& another) const;

    bool operator ==(const Residue<N>& another) const;
    bool operator !=(const Residue<N>& another) const;

    Residue<N> pow(unsigned int power) const;
    Residue<N> pow(signed int power) const = delete;

    Residue<N> getInverse() const;
    unsigned int order() const;
    static Residue<N> getPrimitiveRoot();
};

template <unsigned int N>
struct is_prime {
    static const bool value = Details::
            IsPrimeHelper<N, Details::SquareRoot<N>::value>::value;
    static const unsigned int firstDivisor = Details::
            IsPrimeHelper<N, Details::SquareRoot<N>::value>::firstDivisor;
};

template <>
struct is_prime<1> {
    static const bool value = false;
    static const unsigned int firstDivisor = 1;
};


template <unsigned int N>
struct IsPrimeDegree {
    static const bool value = Details::
            IsPrimeDegreeHelper<N, is_prime<N>::firstDivisor>::value;
};

template <>
struct IsPrimeDegree<0> {
    static const bool value = false;
};

template <unsigned int N>
struct has_primitive_root {
    static const bool value =
            IsPrimeDegree<N % 4 == 0 ? 0 : (N % 2 == 0 ? N / 2 : N)>::value;
};

template <>
struct has_primitive_root<2> {
    static const bool value = true;
};

template <>
struct has_primitive_root<4> {
    static const bool value = true;
};

template <unsigned int N>
const bool is_prime_v = is_prime<N>::value;

template <unsigned int N>
const bool has_primitive_root_v = has_primitive_root<N>::value;

template <unsigned int N>
const unsigned int Residue<N>::_phi = Details::EulerFunction(N);

template <unsigned int N>
Residue<N>::Residue(int number):
        remainder(number >= 0 ? number % N : (number % int(N)) + N) {}

template <unsigned int N>
Residue<N>::operator int() const {
    return static_cast<int>(remainder);
}

template <unsigned int N>
Residue<N>::operator bool() const {
    return remainder != 0;
}


template <unsigned int N>
Residue<N>& Residue<N>::operator +=(const Residue<N>& another) {
    remainder += another.remainder;
    remainder %= N;
    return *this;
}

template <unsigned int N>
Residue<N>& Residue<N>::operator *=(const Residue<N>& another) {
    remainder *= another.remainder;
    remainder %= N;
    return *this;
}

template <unsigned int N>
Residue<N>& Residue<N>::operator -=(const Residue<N>& another) {
    remainder -= another.remainder;
    remainder += N;
    remainder %= N;
    return *this;
}

template <unsigned int N>
Residue<N>& Residue<N>::operator /=(const Residue<N>& another) {
    *this *= another.getInverse();
    return *this;
}

template <unsigned int N>
Residue<N> Residue<N>::operator +(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary += another;
    return temporary;
}

template <unsigned int N>
Residue<N> Residue<N>::operator -(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary -= another;
    return temporary;
}

template <unsigned int N>
Residue<N> Residue<N>::operator *(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary *= another;
    return temporary;
}

template <unsigned int N>
Residue<N> Residue<N>::operator /(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary *= another.getInverse();
    return temporary;
}

template <unsigned int N>
bool Residue<N>::operator ==(const Residue<N>& another) const {
    return remainder == another.remainder;
}

template <unsigned int N>
bool Residue<N>::operator !=(const Residue<N>& another) const {
    return !(*this == another);
}

template <unsigned int N>
Residue<N> Residue<N>::pow(unsigned int power) const {
    if (power == 0) return Residue<N>(1);
    Residue<N> answer(1);
    if (power % 2 != 0) {
        answer *= *this;
    }
    Residue<N> half_power = pow(power / 2);
    answer *= half_power;
    answer *= half_power;
    return answer;
}

template <unsigned int N>
Residue<N> Residue<N>::getInverse() const {
    Details::CheckBoolCompile<is_prime_v<N>>::check();
    int x, y;
    int g = Details::gcdExpanded(remainder, N, x, y);
    if (g == 1) {
        return Residue<N>(x);
    }
    return Residue<N>(0);
}

template <unsigned int N>
unsigned int Residue<N>::order() const {
    if (remainder == 1) return 1;
    unsigned int last_divider = 2;
    for (unsigned int divider = 2; divider * divider <= _phi; ++divider) {
        if (pow(divider).remainder == 1) {
            return divider;
        }
        last_divider = divider;
    }
    for (unsigned int divider = last_divider; divider >= 1; --divider) {
        if (pow(_phi / divider).remainder == 1) {
            return _phi / divider;
        }
    }
    return _phi;
}

template <unsigned int N>
Residue<N> Residue<N>::getPrimitiveRoot() {
    Details::CheckBoolCompile<has_primitive_root_v<N>>::check();
    Residue<N> one(1);
    bool isPrimitiveRoot;
    for (unsigned int candidateRemainder = 2; candidateRemainder <= N; ++candidateRemainder) {
        isPrimitiveRoot = true;
        Residue<N> candidate(candidateRemainder);
        if (Details::gcd(candidateRemainder, N) == 1 && candidate.pow(_phi / 2) != one) {
            for (unsigned int power = 2; power <= sqrt(_phi); ++power) {
                if (candidate.pow(_phi / power) == one) {
                    isPrimitiveRoot = false;
                    break;
                }
            }
            if (isPrimitiveRoot) {
                return candidate;
            }
        }
    }
    return Residue<N>(1);
}

namespace Details {
int gcdExpanded(int first, int second, int& coefficientFirst, int& coefficientSecond) {
    if (first == 0) {
        coefficientFirst = 0;
        coefficientSecond = 1;
        return second;
    }
    int coefficientFirstNew;
    int coefficientSecondNew;
    int gcd = gcdExpanded(second % first, first, coefficientFirstNew, coefficientSecondNew);
    coefficientFirst = coefficientSecondNew - (second / first) * coefficientFirstNew;
    coefficientSecond = coefficientFirstNew;
    return gcd;
}

int gcd(int first, int second) {
    int x, y;
    return gcdExpanded(first, second, x, y);
}

unsigned int EulerFunction(unsigned int number) {
    unsigned int _phi = number;
    for (int divider = 2; divider * divider <= number; ++divider) {
        if (number % divider == 0) {
            while (number % divider == 0) {
                number /= divider;
            }
            _phi -= _phi / divider;
        }
    }
    // maximum divider with order 1
    if (number > 1) {
        _phi -= _phi / number;
    }
    return _phi;
}
} // end of namespace Details
