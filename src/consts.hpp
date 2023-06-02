#pragma once

#include <cstdint>

const unsigned int imageWidth  = 800;
const unsigned int imageHeight = 800;

const double EPSILON = 0.00001;
const uint32_t SOLARIZEDBASE03  = 0x002B36;
const uint32_t SOLARIZEDBASE02  = 0x073642;
const uint32_t SOLARIZEDBASE01  = 0x586e75;
const uint32_t SOLARIZEDBASE00  = 0x657b83;
const uint32_t SOLARIZEDBASE0   = 0x839496;
const uint32_t SOLARIZEDBASE1   = 0x93a1a1;
const uint32_t SOLARIZEDBASE2   = 0xEEE8D5;
const uint32_t SOLARIZEDBASE3   = 0xFDF6E3;
const uint32_t SOLARIZEDYELLOW  = 0xB58900;
const uint32_t SOLARIZEDORANGE  = 0xCB4B16;
const uint32_t SOLARIZEDRED     = 0xDC322F;
const uint32_t SOLARIZEDMAGENTA = 0xD33682;
const uint32_t SOLARIZEDVIOLET  = 0x6C71C4;
const uint32_t SOLARIZEDBLUE    = 0x268BD2;
const uint32_t SOLARIZEDCYAN    = 0x2AA198;
const uint32_t SOLARIZEDGREEN   = 0x859900;


int32_t
sgn(const double x)
{
    return (x > EPSILON) - (x < -EPSILON);
}

int32_t
sgn_strict(const double x)
{
    return (x >= 0.) - (x < 0.);
}
