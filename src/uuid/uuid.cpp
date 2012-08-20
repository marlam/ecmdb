/*
 * Copyright (C) 2011, 2012
 * Martin Lambers <marlam@marlam.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#if HAVE_LIBUUID
# include <uuid/uuid.h>
#else
# include <cstring>
# include <rpc.h>
#endif

#include "uuid.h"


void uuid::generate()
{
#if HAVE_LIBUUID
    uint64_t data[2];
    ::uuid_generate(reinterpret_cast<unsigned char*>(data));
    _d0 = data[0];
    _d1 = data[1];
#else
    UUID u;
    UuidCreate(&u);
    uint64_t data[2];
    std::memcpy(reinterpret_cast<unsigned char *>(data) + 0, &u.Data1, 4);
    std::memcpy(reinterpret_cast<unsigned char *>(data) + 4, &u.Data1, 2);
    std::memcpy(reinterpret_cast<unsigned char *>(data) + 6, &u.Data1, 2);
    std::memcpy(reinterpret_cast<unsigned char *>(data) + 8, &u.Data1, 8);
    _d0 = data[0];
    _d1 = data[1];
#endif
}

std::string uuid::to_string() const
{
#if HAVE_LIBUUID
    uint64_t data[2] = { _d0, _d1 };
    char buf[37];
    ::uuid_unparse(reinterpret_cast<unsigned char*>(data), buf);
    return std::string(buf, 36);
#else
    uint64_t data[2] = { _d0, _d1 };
    UUID u;
    unsigned char *buf;
    std::memcpy(&u.Data1, reinterpret_cast<unsigned char *>(data) + 0, 4);
    std::memcpy(&u.Data1, reinterpret_cast<unsigned char *>(data) + 4, 2);
    std::memcpy(&u.Data1, reinterpret_cast<unsigned char *>(data) + 6, 2);
    std::memcpy(&u.Data1, reinterpret_cast<unsigned char *>(data) + 8, 8);
    UuidToString(&u, &buf);
    std::string s(reinterpret_cast<char *>(buf), 36);
    RpcStringFree(&buf);
    return s;
#endif
}

void uuid::save(std::ostream& os) const
{
    s11n::save(os, _d0);
    s11n::save(os, _d1);
}

void uuid::load(std::istream& is)
{
    s11n::load(is, _d0);
    s11n::load(is, _d1);
}
