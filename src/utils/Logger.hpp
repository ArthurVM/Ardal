#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

namespace ardal {
namespace utils {

enum LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    CRITICAL,
    SILENT
};

extern LogLevel current_log_level;

void set_log_level(LogLevel level);
void set_log_level_from_str(const std::string& level_str);

} // namespace utils
} // namespace ardal

#endif // LOGGER_HPP
