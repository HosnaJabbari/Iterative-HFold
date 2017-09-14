#pragma once

#include <string>
#include <vector>

class shape_info
{
private:
    float b_ = -0.6; // intercept
    float m_ = 1.8; // slope
    int seq_len_ = -1;

    bool use_shape_data_ = false;
    std::string shape_file_;

    std::vector<double> data_;

public:
    shape_info() {}

    // based on http://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
    bool is_number(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end()
                && (std::isdigit(*it) || *it == '.' || *it == '-' || *it == '+')
              ) ++it;
        return !s.empty() && it == s.end();
    }

    void set_data(std::string filename);
    void set_b(float b) { b_ = b; }
    void set_m(float m) { m_ = m; }

    double data(int index) { return data_[index]; }
    float b() { return b_; }
    float m() { return m_; }

    void set_sequence_length(int length) {
        seq_len_ = length;
    }

    const int sequence_length() const {
        return seq_len_;
    }

    // prepares filename to be stored in shape_file
    void set_shape_file(std::string filename);

    const std::string shape_file() const{
        return shape_file_;
    }

    const bool use_shape_data() const {
        return use_shape_data_;
    }

};

extern shape_info shape;
