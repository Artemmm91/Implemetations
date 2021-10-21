#include <iostream>
#include <cstring>

class String{
private:
    size_t size = 0;
    size_t buffer_size = 0;
    char* line = nullptr;

    void change_buffer(size_t new_buffer) {
        buffer_size = new_buffer;
        char* new_line = new char[buffer_size];
        memcpy(new_line, line, size);
        delete[] line;
        line = new_line;
    }

    size_t find(const String& sub_string, size_t start, size_t end) const{
        int delta = start < end ? 1 : -1;
        size_t i = start;
        while(true) {
            if(!memcmp(line + i, sub_string.line, sub_string.length())) return i;
            if(i == end) return size;
            i += delta;
        }
    }

public:
    String(const char* c_string): size(strlen(c_string)), buffer_size(2 * size), line(new char[2 * size]) {
        memcpy(line, c_string, size);
    }

    String(size_t size, char symbol): size(size), buffer_size(2*size), line(new char[2*size]) {
        memset(line, symbol, size);
    }

    String(const char symbol): String(1, symbol) {}

    String(String&& str) {
        std::swap(*this, str);
    }

    explicit String(): size(0), buffer_size(0), line(nullptr){}

    String(const String& string_copy):
            size(string_copy.size), buffer_size(string_copy.buffer_size), line(new char[buffer_size]) {
        memcpy(line, string_copy.line, size);
    }

    ~String() {
        delete[] line;
    }

    bool operator==(const String& string) const{
        if(size != string.length()) return false;
        for(size_t i = 0; i < size; ++i) {
            if(line[i] != string[i]) return false;
        }
        return true;
    }

    String& operator=(const String& string_copy) {
        if(this != &string_copy) {
            delete[] line;
            size = string_copy.size;
            buffer_size = string_copy.buffer_size;
            line = new char[buffer_size];
            memcpy(line, string_copy.line, size);
        }
        return *this;
    }

    String& operator+=(const String& string_add) {
        if(buffer_size < size + string_add.size) {
            change_buffer(buffer_size + string_add.buffer_size);
        }
        memcpy(line + size, string_add.line, string_add.length());
        size += string_add.length();
        return *this;
    }

    String& operator+=(char symbol) {
        push_back(symbol);
        return *this;
    }

    char operator[](size_t index) const{
        return line[index];
    }

    char& operator[](size_t index) {
        return line[index];
    }

    size_t length() const{
        return size;
    }

    void push_back(char new_symbol) {
        if(size == buffer_size) change_buffer(2 * size);
        line[size] = new_symbol;
        ++size;
    }

    void pop_back() {
        if(size != 0) --size;
        if(size * 4 < buffer_size) change_buffer(2 * size);
    }

    char& front() {
        return line[0];
    }

    char front() const{
        return line[0];
    }

    char& back() {
        return line[size - 1];
    }

    char back() const{
        return line[size - 1];
    }

    size_t find(const String& sub_string) const{
        if(sub_string.size > size) return size;
        return find(sub_string, 0, size - sub_string.size);
    }

    // ищет последнее включение подстроки
    size_t rfind(const String& sub_string) const{
        if(sub_string.size > size) return size;
        return find(sub_string, size - sub_string.size + 1, 0);
    }


    String substr(size_t start, size_t count) const{
        String sub_string(count, '0');
        memcpy(sub_string.line, line + start, count);
        return sub_string;
    }

    bool empty() const{
        return size == 0;
    }

    void clear() {
        size = 0;
        delete[] line;
        line = new char[1];
        buffer_size = 1;
    }
};

std::ostream &operator<<(std::ostream& out, const String &string) {
    for(size_t i = 0; i < string.length(); ++i){
        out << string[i];
    }
    return out;
}

std::istream  &operator>>(std::istream& in, String& str) {
    char symbol;
    str.clear();
    in >> std::noskipws;
    while(in >> symbol && symbol != ' ' && symbol != '\n') {
        str.push_back(symbol);
    }
    return in;
}

String operator+(const String& string1, const String& string2) {
    String str = string1;
    str += string2;
    return str;
}

String operator+(const String& string1, char symbol) {
    String str = string1;
    str += symbol;
    return str;
}