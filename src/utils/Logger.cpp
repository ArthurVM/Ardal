#include "Logger.hpp"
#include <unordered_map>
#include <string>

namespace ardal {
namespace utils {

LogLevel current_log_level = WARNING; // Default level

void set_log_level(LogLevel level) {
    current_log_level = level;
}

void set_log_level_from_str(const std::string& level_str) {
    static const std::unordered_map<std::string, LogLevel> level_map = {
        {"debug", DEBUG},
        {"info", INFO},
        {"warn", WARNING},
        {"warning", WARNING},
        {"error", ERROR},
        {"critical", CRITICAL},
        {"silent", SILENT}
    };
    auto it = level_map.find(level_str);
    if (it != level_map.end()) {
        current_log_level = it->second;
    }
}

} // namespace utils
} // namespace ardal
