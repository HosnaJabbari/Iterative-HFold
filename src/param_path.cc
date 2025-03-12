#include "param_path.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>

// Dynamic path (required for conda installation)
char* getParamPath(const char* file_name) {
    static char configPath[512];

    const char* conda_prefix = std::getenv("CONDA_PREFIX");
    if (conda_prefix) {
        std::snprintf(configPath, sizeof(configPath), "%s/share/simfold/params/%s", conda_prefix, file_name);
    } else {
        std::snprintf(configPath, sizeof(configPath), "simfold/params/%s", file_name);
    }

    return configPath;
}
