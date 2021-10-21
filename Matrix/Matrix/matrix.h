#include <vector>
#include <string>
#include <iostream>
#include <initializer_list>
#include <cmath>

class BigInteger;
class Rational;
template <size_t N>
class Residue;
template <size_t N, size_t M, typename Field = Rational>
class Matrix;

class BigInteger {
    friend bool operator<(const BigInteger&, const BigInteger&);
    friend std::ostream& operator<<(std::ostream& out, const BigInteger& number);
    friend std::istream& operator>>(std::istream& in, BigInteger& number);
    
private:
    static const long long OSN = 1e9;
    std::vector<long long> blocks;
    bool is_positive;

    void normalize_zeros();
    void normalize();
    BigInteger divide_by_base(size_t order) const;
    BigInteger multiplicate_by_base(size_t order) const;

public:
    BigInteger();
    BigInteger(long long number);
    
    explicit operator bool() const;
    std::string toString() const;

    BigInteger& operator+=(const BigInteger& number);
    BigInteger& operator-=(const BigInteger& number);
    BigInteger& operator*=(const BigInteger& number);
    BigInteger& operator/=(const BigInteger& number);
    BigInteger& operator%=(const BigInteger& number);
    bool isEven();

    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);

    BigInteger operator-() const;

    void halving();
    void doubling();
};

bool operator<(const BigInteger& first, const BigInteger& second);
bool operator>(const BigInteger& first, const BigInteger& second);
bool operator<=(const BigInteger& first, const BigInteger& second);
bool operator>=(const BigInteger& first, const BigInteger& second);
bool operator==(const BigInteger& first, const BigInteger& second);
bool operator!=(const BigInteger& first, const BigInteger& second);

BigInteger operator+(const BigInteger& first, const BigInteger& second);
BigInteger operator-(const BigInteger& first, const BigInteger& second);
BigInteger operator*(const BigInteger& first, const BigInteger& second);
BigInteger operator/(const BigInteger& first, const BigInteger& second);
BigInteger operator%(const BigInteger& first, const BigInteger& second);

std::ostream& operator<<(std::ostream& out, const BigInteger& number);
std::istream& operator>>(std::istream& in, BigInteger& number);

BigInteger gcd(const BigInteger& first, const BigInteger& second);
BigInteger abs(const BigInteger& number);
void swap(BigInteger& first, BigInteger& second);

class Rational {
    friend bool operator<(const Rational& first, const Rational& second);
    friend Rational operator-(const Rational& number);
private:
    BigInteger nominator;
    BigInteger denominator;

    void reduction();
public:
    Rational();
    Rational(BigInteger number);
    Rational(int number);
    explicit operator bool() const;
    explicit operator double() const;

    Rational& operator+=(const Rational&);
    Rational& operator-=(const Rational&);
    Rational& operator*=(const Rational&);
    Rational& operator/=(const Rational&);

    std::string toString() const;
    std::string asDecimal(size_t precision = 6) const;
};

bool operator<(const Rational& first, const Rational& second);
bool operator>(const Rational& first, const Rational& second);
bool operator<=(const Rational& first, const Rational& second);
bool operator>=(const Rational& first, const Rational& second);
bool operator==(const Rational& first, const Rational& second);
bool operator!=(const Rational& first, const Rational& second);

Rational operator-(const Rational& number);

Rational operator+(const Rational& first, const Rational& second);
Rational operator-(const Rational& first, const Rational& second);
Rational operator*(const Rational& first, const Rational& second);
Rational operator/(const Rational& first, const Rational& second);

std::ostream& operator<<(std::ostream& out, const Rational& number);
std::istream& operator>>(std::istream& in, Rational& number);

template <size_t N>
class Residue {
private:
    long long remainder;
    static const size_t _phi;
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

    Residue<N> pow(size_t power) const;
    Residue<N> pow(signed int power) const = delete;

    Residue<N> getInverse() const;
    size_t order() const;
    static Residue<N> getPrimitiveRoot();
};

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& number);

template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& number);

template <size_t N, size_t M, typename Field>
class Matrix {
private:
    std::vector<std::vector<Field>> table;
public:
    Matrix();
    explicit Matrix(const std::vector<std::vector<Field>>& matrix);
    template<typename AnotherField>
    Matrix(std::initializer_list<std::initializer_list<AnotherField>> matrix);
    template<size_t AnotherN, size_t AnotherM>
    explicit Matrix(const Matrix<AnotherN, AnotherM, Field>& another): Matrix(another.table) {}

    Matrix<N, M, Field>& operator+=(const Matrix<N, M, Field>& another);
    Matrix<N, M, Field>& operator-=(const Matrix<N, M, Field>& another);
    Matrix<N, M, Field>& operator*=(const Field& scalar);
    Matrix<N, N, Field>& operator*=(const Matrix<N, M, Field>& another);

    Matrix<N, M, Field> operator+(const Matrix<N, M, Field>& second) const;
    Matrix<N, M, Field> operator-(const Matrix<N, M, Field>& second) const;

    bool operator==(const Matrix<N, M, Field>& another) const;
    bool operator!=(const Matrix<N, M, Field>& another) const;

    std::vector<Field>& operator[](size_t index);
    const std::vector<Field>& operator[](size_t index) const;

    Field det() const;
    Matrix<M, N, Field> transposed() const;
    size_t rank() const;
    Field trace() const;
    void invert();
    Matrix<N, N, Field> inverted() const;
    std::vector<Field> getRow(size_t index) const;
    std::vector<Field> getColumn(size_t index) const;
    const std::vector<std::vector<Field>>& getTable() const;
};

template <size_t N, size_t M, typename Field>
std::ostream& operator<<(std::ostream& out,
                         const Matrix<N, M, Field>& matrix);

template <size_t N, size_t M, typename Field>
std::istream& operator>>(std::istream& in,
                         Matrix<N, M, Field>& matrix);

template <size_t N, size_t M, size_t K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field>& first,
                              const Matrix<M, K, Field>& second);

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Field& scalar,
                              const Matrix<N, M, Field>& second);

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& first,
                              const Field& scalar);

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> gaussExpanded(const Matrix<N, M, Field>& matrix,
                                  Matrix<N, M, Field>& addMatrix,
                                  bool returnAddMatrix = false);

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> gauss(const Matrix<N, M, Field>& matrix);

template <typename Field>
std::vector<std::vector<Field>> StrassenProduct(const std::vector<std::vector<Field>> first,
                                                const std::vector<std::vector<Field>> second);

template <typename Field>
std::vector<std::vector<Field>> tableSumAssign(std::vector<std::vector<Field>>& first,
                                               const std::vector<std::vector<Field>>& second);
template <typename Field>
std::vector<std::vector<Field>> tableSum(const std::vector<std::vector<Field>>& first,
                                         const std::vector<std::vector<Field>>& second);
template <typename Field>
std::vector<std::vector<Field>> tableDifAssign(std::vector<std::vector<Field>>& first,
                                               const std::vector<std::vector<Field>>& second);
template <typename Field>
std::vector<std::vector<Field>> tableDif(const std::vector<std::vector<Field>>& first,
                                         const std::vector<std::vector<Field>>& second);
template <typename Field>
std::vector<std::vector<Field>> tableScalarAssign(std::vector<std::vector<Field>>& first,
                                                  const Field& scalar);
template <typename Field>
std::vector<std::vector<Field>> tableProduct(const std::vector<std::vector<Field>>& first,
                                             const std::vector<std::vector<Field>>& second);
template <typename Field>
std::vector<std::vector<Field>> resizeTable(const std::vector<std::vector<Field>>& table,
                                            size_t N, size_t M);



template <size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;




//------------------------------- Definitions -------------------------------




namespace Details {

int gcdExpanded(int first, int second, int& coefficientFirst, int& coefficientSecond);
int gcd(int first, int second);
size_t EulerFunction(size_t number);
size_t upperPowerofTwo(size_t n);
std::string block_to_string(long long number, long long osn);

template <unsigned long long N, unsigned long long left_bound, unsigned long long diapason>
struct SquareRootInterval {
    static const size_t value = SquareRootInterval<N,
            (left_bound + diapason / 2) * (left_bound + diapason / 2) <= N
                    ? left_bound + diapason / 2
                    : left_bound,
            (left_bound + diapason / 2) * (left_bound + diapason / 2) <= N
                    ? diapason - diapason / 2
                    : diapason / 2>::value;
};

template <unsigned long long N, unsigned long long left_bound>
struct SquareRootInterval<N, left_bound, 1> {
    static const size_t value = left_bound;
};

template <size_t N>
struct SquareRoot {
    static const size_t value = SquareRootInterval<N, 1, N>::value;
};

template <size_t N, size_t prime>
struct IsPrimeDegreeHelper {
    static const bool value = IsPrimeDegreeHelper<
            N % prime == 0 ? N / prime : 0, prime>::value;
};

template <size_t prime>
struct IsPrimeDegreeHelper<1, prime> {
    static const bool value = true;
};

template <size_t prime>
struct IsPrimeDegreeHelper<0, prime> {
    static const bool value = false;
};

template <size_t N, size_t D>
struct IsPrimeHelper {
    static const bool value = N % D == 0 ? false : IsPrimeHelper<N, D - 1>::value;
    static const size_t firstDivisor =
            value ? N
                  : (IsPrimeHelper<N, D - 1>::value
                        ? D
                        : IsPrimeHelper<N, D - 1>::firstDivisor);
};

template <size_t N>
struct IsPrimeHelper<N, 1> {
    static const bool value = true;
    static const size_t firstDivisor = N;
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


bool operator<(const BigInteger& first, const BigInteger& second) {
    if (first.is_positive && !second.is_positive) return false;
    if (!first.is_positive && second.is_positive) return true;
    bool delta_differ = false;

    if (!first.is_positive && !second.is_positive) {
        delta_differ = true;
    }

    if (first.blocks.size() != second.blocks.size()) {
        return (first.blocks.size() < second.blocks.size()) != delta_differ;
    }
    size_t i = first.blocks.size();
    while (i > 0) {
        auto delta = first.blocks[i - 1] - second.blocks[i - 1];
        if (delta || !(i - 1)) return (delta < 0) != delta_differ;
        --i;
    }
    return false;
}

bool operator>(const BigInteger& first, const BigInteger& second) {
    return second < first;
}

bool operator==(const BigInteger& first, const BigInteger& second) {
    return !(first < second) && !(second < first);
}

bool operator!=(const BigInteger& first, const BigInteger& second) {
    return (first < second) || (second < first);
}

bool operator>=(const BigInteger& first, const BigInteger& second) {
    return !(first < second);
}

bool operator<=(const BigInteger& first, const BigInteger& second) {
    return !(second < first);
}

BigInteger BigInteger::operator-() const {
    BigInteger result = *this;
    if (*this != 0) {
        result.is_positive = !is_positive;
    }
    return result;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result += second;
    return result;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
    return first + (-second);
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result *= second;
    return result;
}

BigInteger operator/(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result /= second;
    return result;
}

BigInteger operator%(const BigInteger& first, const BigInteger& second) {
    BigInteger result = first;
    result %= second;
    return result;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& number) {
    char symbol;
    long long block = 0;
    long long order = 1;
    std::vector<long long> signQueue;
    number.blocks.clear();
    signQueue.clear();
    number.is_positive = true;
    bool isStarted = false;

    while ((symbol = in.get()) && (!iswspace(symbol) || !isStarted)) {
        if (symbol == '-') {
            number.is_positive = false;
            continue;
        }
        if (!iswspace(symbol)) {
            signQueue.push_back(symbol - '0');
            isStarted = true;
        }
    }

    for (size_t i = signQueue.size(); i > 0; --i) {
        block += order * signQueue[i - 1];
        order *= 10;
        if (order >= number.OSN) {
            number.blocks.push_back(block);
            order = 1;
            block = 0;
        }
    }
    if (order != 1) {
        number.blocks.push_back(block);
    }
    number.normalize_zeros();
    return in;
}

void BigInteger::normalize_zeros() {
    while (blocks.size() > 1 && blocks.back() == 0) {
        blocks.pop_back();
    }
    if (blocks.size() == 1 && blocks.back() == 0) {
        is_positive = true;
    }
}

void BigInteger::normalize() {
    long long next_digit;
    for (size_t i = 0; i < blocks.size(); ++i) {
        if (blocks[i] >= 0) {
            next_digit = blocks[i] / OSN;
            blocks[i] = blocks[i] % OSN;
            if (i == blocks.size() - 1) {
                if (next_digit != 0) {
                    blocks.push_back(next_digit);
                }
                continue;
            }
            blocks[i + 1] += next_digit;
        }
        else{
            if (i == blocks.size() - 1) {
                is_positive = !is_positive;
                for (size_t j = 0; j < blocks.size(); ++j) {
                    blocks[j] *= -1;
                }
                normalize();
                return;
            }

            next_digit = (-blocks[i]) / OSN;
            blocks[i] = blocks[i] % OSN;
            if (blocks[i] < 0) {
                blocks[i] += OSN;
                ++next_digit;
            }
            blocks[i + 1] -= next_digit;
        }
    }
}

BigInteger BigInteger::divide_by_base(size_t order) const {
    BigInteger res = 0;
    if (order < blocks.size()) res.blocks.clear();
    for (size_t i = order; i < blocks.size(); ++i) {
        res.blocks.push_back(blocks[i]);
    }
    return res;
}

BigInteger BigInteger::multiplicate_by_base(size_t order) const {
    BigInteger res = 0;
    res.blocks.clear();
    for (size_t i = 0; i < order; ++i) {
        res.blocks.push_back(0);
    }
    for (size_t i = 0; i < blocks.size(); ++i) {
        res.blocks.push_back(blocks[i]);
    }
    return res;
}

BigInteger::BigInteger(): blocks({0}), is_positive(true) {}

BigInteger::BigInteger(long long number) {
    is_positive = true;
    if (number < 0) {
        is_positive = false;
        number = -number;
    }
    blocks.clear();
    if (number == 0) {
        blocks.push_back(0);
    }
    while (number > 0) {
        blocks.push_back(number % OSN);
        number /= OSN;
    }
}

BigInteger::operator bool() const {
    if (blocks.size() > 1) return true;
    if (blocks[0]) return true;
    return false;
}

std::string BigInteger::toString() const {
    std::string s = "";
    if (!is_positive) s += '-';
    s += std::to_string(blocks[blocks.size() - 1]);
    for (size_t i = blocks.size() - 1; i > 0; --i) {
        s += Details::block_to_string(blocks[i - 1], OSN);
    }
    return s;
}

BigInteger& BigInteger::operator+=(const BigInteger& number) {
    size_t number_size = number.blocks.size();
    while (blocks.size() < number_size) {
        blocks.push_back(0);
    }
    int delta = (number.is_positive == is_positive)? 1 : -1;
    for (size_t i = 0; i < number_size; ++i) {
        blocks[i] += delta * number.blocks[i];
    }
    normalize();
    normalize_zeros();
    return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& number) {
    return *this += (-number);
}

BigInteger& BigInteger::operator*=(const BigInteger& number) {
    if (blocks.size() < number.blocks.size() + blocks.size()) {
        blocks.resize(number.blocks.size() + blocks.size());
    }

    for (size_t i = blocks.size(); i > 0; --i) {
        for (size_t j = number.blocks.size(); j > 0; --j) {
            if (j != 1) {
                blocks[i + j - 2] += blocks[i - 1] * number.blocks[j - 1];
            }
            else {
                blocks[i + j - 2] = blocks[i - 1] * number.blocks[j - 1];
            }
        }
    }

    is_positive = (is_positive == number.is_positive);
    normalize();
    normalize_zeros();
    return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& number) {
    BigInteger numerator = abs(*this);
    BigInteger denominator = abs(number);
    BigInteger divided;
    BigInteger subtractive;

    std::vector<long long> digits(0);
    blocks.clear();

    long long digit;
    long long bound_digit;
    long long pivot;
    if (numerator.blocks.size() < denominator.blocks.size()) {
        *this = 0;
        return *this;
    }
    for (size_t i = numerator.blocks.size() - denominator.blocks.size() + 1; i > 0; --i) {
        divided = numerator.divide_by_base(i - 1);
        digit = 0;
        bound_digit = OSN;
        while (bound_digit > digit + 1) {
            pivot = (digit + bound_digit) / 2;
            ((divided >= pivot * denominator) ? digit : bound_digit) = pivot;
        }

        subtractive = digit * denominator;
        numerator -= subtractive.multiplicate_by_base(i - 1);
        digits.push_back(digit);
    }

    size_t digit_size = digits.size();
    for (size_t i = 1; i <= digit_size; ++i) {
        blocks.push_back(digits[digit_size - i]);
    }
    if (digit_size == 0) {
        blocks.push_back(0);
    }

    normalize_zeros();
    is_positive = (is_positive == number.is_positive);

    return *this;
}

void BigInteger::halving() {
    long long remainder = 0;
    long long temporary;
    for (size_t i = blocks.size(); i > 0; --i) {
        temporary = blocks[i - 1] + remainder * OSN;
        blocks[i - 1] = (temporary / 2);
        remainder = (temporary % 2);
    }
    normalize_zeros();
}

void BigInteger::doubling() {
    long long remainder = 0;
    long long temporary;
    for (size_t i = 0; i < blocks.size() || remainder; ++i) {
        if (i == blocks.size()) {
            blocks.push_back(0);
        }
        temporary = remainder + blocks[i] * 2;
        blocks[i] = temporary % OSN;
        remainder = temporary / OSN;
    }
    normalize_zeros();
}

BigInteger& BigInteger::operator%=(const BigInteger& number) {
    BigInteger even = *this / number;
    even *= number;
    *this -= even;
    return *this;
}

BigInteger& BigInteger::operator++() {
    return *this += 1;
}

BigInteger BigInteger::operator++(int) {
    BigInteger temporary = *this;
    ++*this;
    return temporary;
}

BigInteger& BigInteger::operator--() {
    return *this -= 1;
}

BigInteger BigInteger::operator--(int) {
    BigInteger temporary = *this;
    --*this;
    return temporary;
}

bool BigInteger::isEven() {
    return blocks[0] % 2 == 0;
}

BigInteger gcd(const BigInteger& first, const BigInteger& second) {
    BigInteger a = abs(first);
    BigInteger b = abs(second);
    BigInteger answer(1);
    if (a == 0 || b == 0) return 0;
    if (a > b) swap(a, b);
    bool evenA, evenB;
    while (a && b) {
        evenA = a.isEven();
        evenB = b.isEven();
        if (evenA && evenB) {
            a.halving();
            b.halving();
            answer.doubling();
            continue;
        }
        if (evenA) {
            a.halving();
            continue;
        }
        if (evenB) {
            b.halving();
            continue;
        }
        if (a > b) swap(a, b);
        b -= a;
    }
    return (a ? a : b) * answer;
}

BigInteger abs(const BigInteger& number) {
    if (number < 0) return -number;
    return number;
}

void swap(BigInteger& first, BigInteger& second) {
    BigInteger temporary = first;
    first = second;
    second = temporary;
}


Rational::Rational(): nominator(0), denominator(1) {}

Rational::Rational(BigInteger number): nominator(number), denominator(1) {}

Rational::Rational(int number): nominator(number), denominator(1) {}

Rational& Rational::operator+=(const Rational& number) {
    nominator = nominator * number.denominator + denominator * number.nominator;
    denominator = denominator * number.denominator;
    reduction();
    return *this;
}

Rational& Rational::operator-=(const Rational& number) {
    return *this += (-number);
}

Rational& Rational::operator*=(const Rational& number) {
    nominator *= number.nominator;
    denominator *= number.denominator;
    reduction();
    return *this;
}

Rational& Rational::operator/=(const Rational& number) {
    nominator *= number.denominator;
    denominator *= number.nominator;
    if (denominator < 0) {
        denominator = -denominator;
        nominator = -nominator;
    }
    reduction();
    return *this;
}

std::string Rational::toString() const {
    std::string result = "";
    result += nominator.toString();
    if (denominator != 1) {
        result += '/';
        result += denominator.toString();
    }
    return result;
}

std::string Rational::asDecimal(size_t precision) const {
    std::string result = "";
    Rational number = *this;
    BigInteger even;
    bool is_dot = false;

    if (number.nominator < 0) {
        result += '-';
        number.nominator = -number.nominator;
    }
    ++precision;

    while (precision > 0) {
        even = number.nominator / number.denominator;
        number -= even;
        number *= 10;
        result += even.toString();
        --precision;

        if (!is_dot && precision > 0) {
            is_dot = true;
            result += '.';
        }
    }
    return result;
}

std::ostream& operator<<(std::ostream& out, const Rational& number) {
    out << number.asDecimal(15);
    return out;
}

std::istream& operator>>(std::istream& in, Rational& number) {
    BigInteger nominator;
    in >> nominator;
    number = Rational(nominator);
    return in;
}

bool operator<(const Rational& first, const Rational& second) {
    return first.nominator * second.denominator < second.nominator * first.denominator;
}

bool operator>(const Rational& first, const Rational& second) {
    return second < first;
}

bool operator==(const Rational& first, const Rational& second) {
    return !(first < second) && !(second < first);
}

bool operator!=(const Rational& first, const Rational& second) {
    return (first < second) || (second < first);
}

bool operator>=(const Rational& first, const Rational& second) {
    return !(first < second);
}

bool operator<=(const Rational& first, const Rational& second) {
    return !(second < first);
}

Rational operator-(const Rational& number) {
    Rational result = number;
    result.nominator = -result.nominator;
    return result;
}

Rational operator+(const Rational& first, const Rational& second) {
    Rational result = first;
    result += second;
    return result;
}

Rational operator-(const Rational& first, const Rational& second) {
    return first + (-second);
}

Rational operator*(const Rational& first, const Rational& second) {
    Rational result = first;
    result *= second;
    return result;
}

Rational operator/(const Rational& first, const Rational& second) {
    Rational result = first;
    result /= second;
    return result;
}

void Rational::reduction() {
    BigInteger coefficient = gcd(nominator, denominator);
    if (coefficient == 0) {
        denominator = 1;
        return;
    }
    nominator /= coefficient;
    denominator /= coefficient;
}

Rational::operator bool() const {
    return nominator != 0;
}

Rational::operator double() const {
    return std::stod(asDecimal(15));
}

template <size_t N>
struct is_prime {
    static const bool value = Details::
            IsPrimeHelper<N, Details::SquareRoot<N>::value>::value;
    static const size_t firstDivisor = Details::
            IsPrimeHelper<N, Details::SquareRoot<N>::value>::firstDivisor;
};

template <>
struct is_prime<1> {
    static const bool value = false;
    static const size_t firstDivisor = 1;
};


template <size_t N>
struct IsPrimeDegree {
    static const bool value = Details::
            IsPrimeDegreeHelper<N, is_prime<N>::firstDivisor>::value;
};

template <>
struct IsPrimeDegree<0> {
    static const bool value = false;
};

template <size_t N>
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

template <size_t N>
const bool is_prime_v = is_prime<N>::value;

template <size_t N>
const bool has_primitive_root_v = has_primitive_root<N>::value;

template <size_t N>
const size_t Residue<N>::_phi = Details::EulerFunction(N);

template <size_t N>
Residue<N>::Residue(int number):
        remainder(number >= 0 ? number % N : (number % int(N)) + N) {}

template <size_t N>
Residue<N>::operator int() const {
    return static_cast<int>(remainder);
}

template <size_t N>
Residue<N>::operator bool() const {
    return remainder != 0;
}


template <size_t N>
Residue<N>& Residue<N>::operator +=(const Residue<N>& another) {
    remainder += another.remainder;
    remainder %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator *=(const Residue<N>& another) {
    remainder *= another.remainder;
    remainder %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator -=(const Residue<N>& another) {
    remainder -= another.remainder;
    remainder += N;
    remainder %= N;
    return *this;
}

template <size_t N>
Residue<N>& Residue<N>::operator /=(const Residue<N>& another) {
    *this *= another.getInverse();
    return *this;
}

template <size_t N>
Residue<N> Residue<N>::operator +(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary += another;
    return temporary;
}

template <size_t N>
Residue<N> Residue<N>::operator -(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary -= another;
    return temporary;
}

template <size_t N>
Residue<N> Residue<N>::operator *(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary *= another;
    return temporary;
}

template <size_t N>
Residue<N> Residue<N>::operator /(const Residue<N>& another) const {
    Residue<N> temporary = *this;
    temporary *= another.getInverse();
    return temporary;
}

template <size_t N>
bool Residue<N>::operator ==(const Residue<N>& another) const {
    return remainder == another.remainder;
}

template <size_t N>
bool Residue<N>::operator !=(const Residue<N>& another) const {
    return !(*this == another);
}

template <size_t N>
Residue<N> Residue<N>::pow(size_t power) const {
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

template <size_t N>
Residue<N> Residue<N>::getInverse() const {
    Details::CheckBoolCompile<is_prime_v<N>>::check();
    int x, y;
    int g = Details::gcdExpanded(remainder, N, x, y);
    if (g == 1) {
        return Residue<N>(x);
    }
    return Residue<N>(0);
}

template <size_t N>
size_t Residue<N>::order() const {
    if (remainder == 1) return 1;
    size_t last_divider = 2;
    for (size_t divider = 2; divider * divider <= _phi; ++divider) {
        if (pow(divider).remainder == 1) {
            return divider;
        }
        last_divider = divider;
    }
    for (size_t divider = last_divider; divider >= 1; --divider) {
        if (pow(_phi / divider).remainder == 1) {
            return _phi / divider;
        }
    }
    return _phi;
}

template <size_t N>
Residue<N> Residue<N>::getPrimitiveRoot() {
    Details::CheckBoolCompile<has_primitive_root_v<N>>::check();
    Residue<N> one(1);
    bool isPrimitiveRoot;
    for (size_t candidateRemainder = 2; candidateRemainder <= N; ++candidateRemainder) {
        isPrimitiveRoot = true;
        Residue<N> candidate(candidateRemainder);
        if (Details::gcd(candidateRemainder, N) == 1 && candidate.pow(_phi / 2) != one) {
            for (size_t power = 2; power <= sqrt(_phi); ++power) {
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

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& number) {
    out << int(number);
    return out;
}

template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& number) {
    int remainder;
    in >> remainder;
    number = Residue<N>(remainder);
    return in;
}

template <typename Field>
std::vector<std::vector<Field>> tableSumAssign(std::vector<std::vector<Field>>& first,
                                               const std::vector<std::vector<Field>>& second) {
    size_t N = first.size();
    size_t M = first[0].size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            first[i][j] += second[i][j];
        }
    }
    return first;
}

template <typename Field>
std::vector<std::vector<Field>> tableSum(const std::vector<std::vector<Field>>& first,
                                         const std::vector<std::vector<Field>>& second) {
    std::vector<std::vector<Field>> answer = first;
    tableSumAssign(answer, second);
    return answer;
}

template <typename Field>
std::vector<std::vector<Field>> tableDifAssign(std::vector<std::vector<Field>>& first,
                                               const std::vector<std::vector<Field>>& second) {
    size_t N = first.size();
    size_t M = first[0].size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            first[i][j] -= second[i][j];
        }
    }
    return first;
}

template <typename Field>
std::vector<std::vector<Field>> tableDif(const std::vector<std::vector<Field>>& first,
                                         const std::vector<std::vector<Field>>& second) {
    std::vector<std::vector<Field>> answer = first;
    tableDifAssign(answer, second);
    return answer;
}

template <typename Field>
std::vector<std::vector<Field>> tableScalarAssign(std::vector<std::vector<Field>>& first,
                                                  const Field& scalar) {
    size_t N = first.size();
    size_t M = first[0].size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            first[i][j] *= scalar;
        }
    }
    return first;
}

template <typename Field>
std::vector<std::vector<Field>> tableProduct(const std::vector<std::vector<Field>>& first,
                                             const std::vector<std::vector<Field>>& second) {
    size_t N = first.size();
    size_t M = second.size();
    size_t K = second[0].size();
    std::vector<std::vector<Field>> Product(N, std::vector<Field>(K, Field(0)));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < K; ++j) {
            for (size_t s = 0; s < M; ++s) {
                Product[i][j] += first[i][s] * second[s][j];
            }
        }
    }
    return Product;
}

template <typename Field>
std::vector<std::vector<Field>> resizeTable(const std::vector<std::vector<Field>>& table, size_t N, size_t M) {
    std::vector<std::vector<Field>> newTable(N, std::vector<Field>(M, Field(0)));
    size_t n = table.size();
    size_t m = table[0].size();
    for (size_t i = 0; i < std::min(n, N); ++i) {
        for (size_t j = 0; j < std::min(m, M); ++j) {
            newTable[i][j] = table[i][j];
        }
    }
    return newTable;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>::Matrix() {
    table.clear();
    for (size_t i = 0; i < N; ++i) {
        table.push_back({});
        for (size_t j = 0; j < M; ++j) {
            table[i].push_back(Field(0));
        }
    }
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>::Matrix(const std::vector<std::vector<Field>>& matrix): table(resizeTable(matrix, N, M)) {}

template <size_t N, size_t M, typename Field>
template<typename AnotherField>
Matrix<N, M, Field>::Matrix(std::initializer_list<std::initializer_list<AnotherField>> matrix) {
    table.clear();
    for (auto i: matrix) {
        table.push_back({});
        for (auto j: i) {
            table.back().push_back(Field(j));
        }
    }
}
template <size_t N, size_t M, typename Field>
const std::vector<std::vector<Field>>& Matrix<N, M, Field>::getTable() const {
    return table;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator+=(const Matrix<N, M, Field>& another) {
    tableSumAssign(table, another.table);
    return *this;
}

template <size_t N, size_t M, typename Field>
Matrix<N, N, Field>& Matrix<N, M, Field>::operator*=(const Matrix<N, M, Field>& another) {
    Details::CheckBoolCompile<N == M>::check();
    *this = *this * another;
    return *this;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator-=(const Matrix<N, M, Field>& another) {
    tableDifAssign(table, another.table);
    return *this;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::operator+(const Matrix<N, M, Field>& second) const{
    Matrix<N, M, Field> answer = *this;
    answer += second;
    return answer;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::operator-(const Matrix<N, M, Field>& second) const {
    Matrix<N, M, Field> answer = *this;
    answer -= second;
    return answer;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=(const Field& scalar) {
    tableScalarAssign(table, scalar);
    return *this;
}

template <size_t N, size_t M, typename Field>
std::vector<Field>& Matrix<N, M, Field>::operator[](size_t index) {
    return table[index];
}

template <size_t N, size_t M, typename Field>
const std::vector<Field>& Matrix<N, M, Field>::operator[](size_t index) const {
    return table[index];
}


template <size_t N, size_t M, size_t K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field>& first,
                              const Matrix<M, K, Field>& second) {
    size_t maxDimension = std::max(std::max(N, M), K);
    if(maxDimension > 64) {
        size_t newSize = Details::upperPowerofTwo(maxDimension);
        std::vector<std::vector<Field>> newFirst = resizeTable(first.getTable(), newSize, newSize);
        std::vector<std::vector<Field>> newSecond = resizeTable(second.getTable(), newSize, newSize);
        return Matrix<N, K, Field>(StrassenProduct(newFirst, newSecond));
    }
    return Matrix<N, K, Field>(tableProduct(first.getTable(), second.getTable()));
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Field& scalar,
                              const Matrix<N, M, Field>& second) {
    Matrix<N, M, Field> Product = second;
    Product *= scalar;
    return Product;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> operator*(const Matrix<N, M, Field>& first,
                              const Field& scalar) {
    Matrix<N, M, Field> Product = first;
    Product *= scalar;
    return Product;
}

template <size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::det() const {
    Details::CheckBoolCompile<N == M>::check();
    Matrix<N, N, Field> matrix = gauss(*this);
    Field determinant(1);
    for (size_t leader = 0; leader < N; ++leader) {
        determinant *= matrix[leader][leader];
    }
    return determinant;
}


template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> gaussExpanded(const Matrix<N, M, Field>& matrix,
                                  Matrix<N, M, Field>& addMatrix,
                                  bool isAddMatrix) {
    Matrix<N, M, Field> temporary = matrix;
    bool isInverted = false;
    for (size_t row = 0, col = 0; row < N && col < M; ++col) {
        for (size_t i = row; i < N; ++i) {
            if (bool(temporary[i][col])) {
                swap(temporary[i], temporary[row]);
                if (isAddMatrix) {
                    swap(addMatrix[i], addMatrix[row]);
                }
                if (i != row) {
                    isInverted = !isInverted;
                }
                break;
            }
        }
        Field leader = temporary[row][col];
        if (!bool(leader)) {
            continue;
        }

        for (size_t i = 0; i < N; ++i) {
            if (i != row && bool(temporary[i][col])) {
                Field coefficient = temporary[i][col];
                coefficient /= leader;
                for (size_t j = 0; j < M; ++j) {
                    temporary[i][j] -= coefficient * temporary[row][j];
                    if (isAddMatrix) {
                        addMatrix[i][j] -= coefficient * addMatrix[row][j];
                    }
                }
            }
        }
        ++row;
    }
    if (isInverted) {
        temporary[0][0] *= Field(-1);
    }
    if (isAddMatrix) {
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                addMatrix[i][j] /= temporary[i][i];
            }
        }
    }
    return temporary;
}

template <size_t N, size_t M, typename Field>
Matrix<N, M, Field> gauss(const Matrix<N, M, Field>& matrix) {
    Matrix<N, M, Field> add;
    return gaussExpanded(matrix, add);
}


template <size_t N, size_t M, typename Field>
std::ostream& operator<<(std::ostream& out,
                         const Matrix<N, M, Field>& matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            out << matrix[i][j] << " ";
        }
        out << "\n";
    }
    return out;
}

template <size_t N, size_t M, typename Field>
std::istream& operator>>(std::istream& in,
                         Matrix<N, M, Field>& matrix) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            in >> matrix[i][j];
        }
    }
    return in;
}

template <size_t N, size_t M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
    Matrix<M, N, Field> answer;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            answer[j][i] = table[i][j];
        }
    }
    return answer;
}

template <size_t N, size_t M, typename Field>
size_t Matrix<N, M, Field>::rank() const {
    Matrix<N, M, Field> matrix = gauss(*this);
    size_t leaders = 0;
    size_t last_leader = 0;
    for (size_t j = 0; j < M; ++j) {
        for (size_t i = last_leader; i < N; ++i) {
            if (bool(matrix[i][j])) {
                ++leaders;
                last_leader = i + 1;
                break;
            }
        }
    }

    return leaders;
}
template <size_t N, size_t M, typename Field>
Field Matrix<N, M, Field>::trace() const {
    Details::CheckBoolCompile<N == M>::check();
    Field traceSum(0);
    for (size_t leader = 0; leader < N; ++leader) {
        traceSum += table[leader][leader];
    }
    return traceSum;
}

template <size_t N, size_t M, typename Field>
void Matrix<N, M, Field>::invert() {
    Details::CheckBoolCompile<N == M>::check();
    Matrix<N, N, Field> Id;
    for (size_t i = 0; i < N; ++i) {
        Id[i][i] = Field(1);
    }
    gaussExpanded(*this, Id, true);
    *this = Id;
}

template <size_t N, size_t M, typename Field>
Matrix<N, N, Field> Matrix<N, M, Field>::inverted() const {
    Details::CheckBoolCompile<N == M>::check();
    Matrix<N, N, Field> matrix = *this;
    matrix.invert();
    return matrix;
}

template <size_t N, size_t M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getRow(size_t index) const {
    return table[index];
}

template <size_t N, size_t M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getColumn(size_t index) const {
    std::vector<Field> column;
    for (size_t i = 0; i < N; ++i) {
        column.push_back(table[i][index]);
    }
    return column;
}

template <size_t N, size_t M, typename Field>
bool Matrix<N, M, Field>::operator==(const Matrix<N, M, Field>& another) const {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (table[i][j] != another[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template <size_t N, size_t M, typename Field>
bool Matrix<N, M, Field>::operator!=(const Matrix<N, M, Field>& another) const {
    return !(*this == another);
}

template <typename Field>
std::vector<std::vector<Field>> StrassenProduct(const std::vector<std::vector<Field>> first,
                                                const std::vector<std::vector<Field>> second) {
    size_t N = first.size();
    if(N < 64) {
        return tableProduct(first, second);
    }
    size_t halfN = N / 2;
    Field zeroElement = Field(0);
    std::vector<std::vector<Field>> A11(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> A12(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> A21(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> A22(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> B11(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> B12(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> B21(halfN, std::vector<Field>(halfN, zeroElement));
    std::vector<std::vector<Field>> B22(halfN, std::vector<Field>(halfN, zeroElement));
    
    for (size_t i = 0; i < halfN; ++i){
        for (size_t j = 0; j < halfN; ++j) {
            A11[i][j] = first[i][j];
            A12[i][j] = first[i][halfN + j];
            A21[i][j] = first[halfN + i][j];
            A22[i][j] = first[halfN + i][halfN + j];
            B11[i][j] = second[i][j];
            B12[i][j] = second[i][halfN + j];
            B21[i][j] = second[halfN + i][j];
            B22[i][j] = second[halfN + i][halfN + j];
        }
    }
    
    std::vector<std::vector<Field>> P1 = StrassenProduct(A11, tableDif(B12, B22));
    std::vector<std::vector<Field>> P2 = StrassenProduct(tableSum(A11, A12), B22);
    std::vector<std::vector<Field>> P3 = StrassenProduct(tableSum(A21, A12), B11);
    std::vector<std::vector<Field>> P4 = StrassenProduct(A22, tableDif(B21, B11));
    std::vector<std::vector<Field>> P5 = StrassenProduct(tableSum(A11, A22), tableSum(B11, B22));
    std::vector<std::vector<Field>> P6 = StrassenProduct(tableDif(A12, A22), tableSum(B21, B22));
    std::vector<std::vector<Field>> P7 = StrassenProduct(tableDif(A11, A21), tableSum(B11, B12));

    std::vector<std::vector<Field>> C11 = tableDif(tableSum(tableSum(P5, P4), P6), P2);
    std::vector<std::vector<Field>> C12 = tableSum(P1, P2);
    std::vector<std::vector<Field>> C21 = tableSum(P3, P4);
    std::vector<std::vector<Field>> C22 = tableDif(tableDif(tableSum(P5, P1), P3), P7);

    std::vector<std::vector<Field>> C(N, std::vector<Field>(N, zeroElement));

    for (size_t i = 0; i < halfN; ++i) {
        for (size_t j = 0; j < halfN; ++j) {
            C[i][j] = C11[i][j];
            C[i][j + halfN] = C12[i][j];
            C[halfN + i][j] = C21[i][j];
            C[halfN + i][halfN + j] = C22[i][j];
        }
    }
    return C;
}

namespace Details {

int gcdExpanded(int first, int second,
                int& coefficientFirst, int& coefficientSecond) {
    if (first == 0) {
        coefficientFirst = 0;
        coefficientSecond = 1;
        return second;
    }
    int coefficientFirstNew;
    int coefficientSecondNew;
    int gcd = gcdExpanded(second % first, first,
                          coefficientFirstNew, coefficientSecondNew);
    coefficientFirst = coefficientSecondNew - (second / first) * coefficientFirstNew;
    coefficientSecond = coefficientFirstNew;
    return gcd;
}

int gcd(int first, int second) {
    int x, y;
    return gcdExpanded(first, second, x, y);
}

size_t EulerFunction(size_t number) {
    size_t _phi = number;
    for (size_t divider = 2; divider * divider <= number; ++divider) {
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

std::string block_to_string(long long number, long long osn) {
    std::string l = std::to_string(number);
    int count_zeros = 0;
    osn /= 10;
    while (number < osn) {
        ++count_zeros;
        osn /= 10;
    }
    if (number == 0) --count_zeros;
    while (count_zeros) {
        l = '0' + l;
        --count_zeros;
    }
    return l;
}

size_t upperPowerofTwo(size_t n) {
    size_t ans = 1;
    size_t one = 1;
    while (ans < n) {
        ans <<= one;
    }
    return ans;
}
} // end of namespace Details
