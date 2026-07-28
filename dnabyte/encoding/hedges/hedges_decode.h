#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include "hedges_config.h"


std::vector<uint8_t> hedges_decode(
    const std::vector<std::string>& strands,
    const HedgesConfig& cfg
);
