#include "config.h"
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cctype>

using namespace std;

namespace elaboradar {

namespace {

struct ConfigValue
{
    const char* key; // Config key
    const char* env; // Environment variable name
    const char* def; // Default value

    ConfigValue(const char* key, const char* env, const char* def=0)
        : key(key), env(env), def(def) {}

};

ConfigValue config_layout[] = {
    { "file/first_level", "FIRST_LEVEL_FILE" },
    { "file/vpr_heating", "VPR_HEATING" },
};

/// Trim spaces at both ends of a string
std::string trim(const std::string s)
{
    if (s.empty()) return s;

    size_t start = 0;
    size_t end = s.size() - 1;
    while (start < end && isspace(s[start]))
        ++start;
    while (end > start && isspace(s[end]))
        --end;

    if (start == end) return string();
    if (start == 0 && end == s.size() - 1) return s;
    return s.substr(start, end - start);
}

/// Check if a line is empty or a comment
bool is_empty(const std::string& line)
{
    for (auto i : line)
    {
        if (isspace(i)) continue;
        if (i == '#') return true;
        return false;
    }
    return true;
}

}

void Config::set_defaults()
{
    for (auto v : config_layout)
    {
        if (!v.def) continue;
        values[v.key] = v.def;
    }
}

void Config::read_env()
{
    for (auto v : config_layout)
    {
        if (!v.env) continue;
        const char* val = getenv(v.env);
        if (val == NULL) continue;
        values[v.key] = val;
    }
}

void Config::read_file(const std::string& fname)
{
    ifstream in(fname);
    string line;
    unsigned lineno = 0;

    // Whitelist of valid keys
    unordered_set<string> valid_keys;
    for (auto v : config_layout)
        valid_keys.insert(v.key);

    while (getline(in, line))
    {
        ++lineno;
        if (is_empty(line)) continue;

        size_t pos = line.find('=');
        if (pos == string::npos)
        {
            stringstream s;
            s << fname << ":" << lineno << ": unparsed line (not empty, a comment or key = value)";
            throw runtime_error(s.str());
        }
        string key = trim(line.substr(0, pos));
        string val = trim(line.substr(pos+1));
        if (valid_keys.find(key) == valid_keys.end())
        {
            stringstream s;
            s << fname << ":" << lineno << ": unsupported key '" << key << "'";
            throw runtime_error(s.str());
        }
        values[key] = val;
    }
}

bool Config::has(const std::string& key) const
{
    return values.find(key) != values.end();
}

std::string Config::get(const std::string& key) const
{
    auto res = values.find(key);
    if (res == values.end()) throw std::runtime_error("configuration key " + key + " not found");
    return res->second;
}

std::string Config::get(const std::string& key, const std::string& deflt) const
{
    auto res = values.find(key);
    if (res == values.end()) return deflt;
    return res->second;
}

int Config::get_int(const std::string& key) const
{
    auto res = values.find(key);
    if (res == values.end()) throw std::runtime_error("configuration key " + key + " not found");
    return stoi(res->second);
}

int Config::get_int(const std::string& key, int deflt) const
{
    auto res = values.find(key);
    if (res == values.end()) return deflt;
    return stoi(res->second);
}

double Config::get_double(const std::string& key) const
{
    auto res = values.find(key);
    if (res == values.end()) throw std::runtime_error("configuration key " + key + " not found");
    return stod(res->second);
}

double Config::get_double(const std::string& key, double deflt) const
{
    auto res = values.find(key);
    if (res == values.end()) return deflt;
    return stod(res->second);
}

}
