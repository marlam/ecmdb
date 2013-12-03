/*
 * Copyright (C) 2011, 2012, 2013
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

#ifndef UUID_H
#define UUID_H

#include <cstdint>
#include <string>

#include "ser.h"


class uuid : public serializable
{
private:
    uint64_t _d0, _d1;

public:
    uuid() : _d0(0), _d1(0) { };

    void generate();
    std::string to_string() const;
    
    bool operator==(const uuid& id) const
    {
        return _d0 == id._d0 && _d1 == id._d1;
    }

    bool operator<(const uuid& id) const
    {
        return _d0 < id._d0 || (_d0 == id._d0 && _d1 < id._d1);
    }

    bool operator>(const uuid& id) const
    {
        return _d0 > id._d0 || (_d0 == id._d0 && _d1 > id._d1);
    }

    // Serialization
    void save(std::ostream& os) const;
    void load(std::istream& is);
};

#endif
