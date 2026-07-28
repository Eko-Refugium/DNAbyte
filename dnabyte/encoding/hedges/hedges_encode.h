#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include "hedges_config.h"

std::vector<std::string> hedges_encode(
    const std::vector<uint8_t>& data,
    const HedgesConfig& cfg
);