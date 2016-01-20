#ifndef ARCHIVIATORE_CONFIG_H
#define ARCHIVIATORE_CONFIG_H

#include <string>
#include <map>

namespace radarelab {

class Config
{
protected:
    std::map<std::string, std::string> values;

    void set_defaults();

public:
    /// Read config values from the environment
    void read_env();

    /// Read config values from a file
    void read_file(const std::string& fname);

    bool has(const std::string& key) const;

    std::string get(const std::string& key) const;
    std::string get(const std::string& key, const std::string& deflt) const;

    int get_int(const std::string& key) const;
    int get_int(const std::string& key, int deflt) const;

    double get_double(const std::string& key) const;
    double get_double(const std::string& key, double deflt) const;
};

}

#endif
