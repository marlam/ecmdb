/*
 * Copyright (C) 2009, 2010, 2011, 2012
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

#ifndef LRU_CACHE_H
#define LRU_CACHE_H

#include <map>

#include "dbg.h"
#include "exc.h"
#include "thread.h"


template<typename ELEMENT_TYPE, typename KEY_TYPE, bool AUTO_SHRINK = true>
class lru_cache
{
private:
    class value                                                 // internal type to store cache elements
    {
    public:
        unsigned long long timestamp;
        const size_t size;
        const ELEMENT_TYPE *element;

        value(unsigned long long t, size_t s, const ELEMENT_TYPE *e) :
            timestamp(t), size(s), element(e)
        {
        }
    };

    // minimum estimated memory requirement for storing one element
    static const size_t _overhead_size = sizeof(value) + sizeof(unsigned long long) + 2 * sizeof(KEY_TYPE);

    size_t _max_size;                                           // max size of cache
    size_t _size;                                               // current size of cache
    unsigned long long _timestamp;                              // current timestamp (assumed to not wrap around)
    std::map<KEY_TYPE, value> _element_map;                     // cache elements
    std::map<unsigned long long, KEY_TYPE> _timestamp_map;      // cache timestamps

    mutex _mutex;                                               // Mutex for locked access

    // remove the least recently used entry
    void remove_lru_element()
    {
        typename std::map<unsigned long long, KEY_TYPE>::iterator timestamp_it = _timestamp_map.begin();
        typename std::map<KEY_TYPE, value>::iterator element_it = _element_map.find(timestamp_it->second);
        assert(element_it != _element_map.end());
        delete element_it->second.element;
        _size -= element_it->second.size;
        _size -= _overhead_size;
        _element_map.erase(element_it);
        _timestamp_map.erase(timestamp_it);
    }

protected:
    // Implement this if you want your cache to be able to fetch unavailable elements itself.
    // In this case, the function should return true if it is able to fetch the element and
    // false otherwise.
    virtual bool fetch_element(const KEY_TYPE& /* key */, ELEMENT_TYPE** /* element */, size_t* /* size */)
    {
        return false;
    }

public:
    lru_cache(size_t max_size) :
        _max_size(max_size), _size(0), _timestamp(0)
    {
    }

    ~lru_cache()
    {
        this->clear();
    }

    const ELEMENT_TYPE* get(const KEY_TYPE& key)
    {
        typename std::map<KEY_TYPE, value>::iterator element_it = _element_map.find(key);
        if (element_it == _element_map.end()) {
            ELEMENT_TYPE* element = NULL;
            size_t element_size;
            if (fetch_element(key, &element, &element_size))
                put(key, element, element_size);
            return element;
        } else {
            unsigned long long old_timestamp = element_it->second.timestamp;
            unsigned long long new_timestamp = atomic::inc_and_fetch(&_timestamp);
            typename std::map<unsigned long long, KEY_TYPE>::iterator timestamp_it = _timestamp_map.find(old_timestamp);
            assert(timestamp_it != _timestamp_map.end());
            _timestamp_map.erase(timestamp_it);
            _timestamp_map.insert(std::pair<unsigned long long, KEY_TYPE>(new_timestamp, key));
            element_it->second.timestamp = new_timestamp;
            return element_it->second.element;
        }
    }

    const ELEMENT_TYPE* locked_get(const KEY_TYPE& key)
    {
        const ELEMENT_TYPE* element;
        _mutex.lock();
        try {
            element = get(key);
        }
        catch (exc& e) {
            _mutex.unlock();
            throw e;
        }
        catch (std::exception& e) {
            _mutex.unlock();
            throw exc(e);
        }
        _mutex.unlock();
        return element;
    }

    void put(const KEY_TYPE& key, const ELEMENT_TYPE* element, size_t size = 0)
    {
        unsigned long long new_timestamp = atomic::inc_and_fetch(&_timestamp);
        typename std::map<KEY_TYPE, value>::iterator element_it = _element_map.find(key);
        if (element_it == _element_map.end()) {
            _element_map.insert(std::pair<KEY_TYPE, value>(key, value(new_timestamp, size, element)));
            _timestamp_map.insert(std::pair<unsigned long long, KEY_TYPE>(new_timestamp, key));
            if (size == 0) {
                // Count number of elements
                _size += 1;
            } else {
                // Count actual size
                _size += size;
                _size += _overhead_size;
            }
            if (AUTO_SHRINK)
                shrink();
        } else {
            unsigned long long old_timestamp = element_it->second.timestamp;
            typename std::map<unsigned long long, KEY_TYPE>::iterator timestamp_it = _timestamp_map.find(old_timestamp);
            assert(timestamp_it != _timestamp_map.end());
            _timestamp_map.erase(timestamp_it);
            _timestamp_map.insert(std::pair<unsigned long long, KEY_TYPE>(new_timestamp, key));
            element_it->second.timestamp = new_timestamp;
            element_it->second.element = element;
        }
    }

    void locked_put(const KEY_TYPE& key, const ELEMENT_TYPE* element, size_t size = 0)
    {
        _mutex.lock();
        try {
            put(key, element, size);
        }
        catch (exc& e) {
            _mutex.unlock();
            throw e;
        }
        catch (std::exception& e) {
            _mutex.unlock();
            throw exc(e);
        }
        _mutex.unlock();
    }

    void clear()
    {
        for (typename std::map<KEY_TYPE, value>::iterator it = _element_map.begin(); it != _element_map.end(); it++) {
            delete it->second.element;
        }
        _element_map.clear();
        _timestamp_map.clear();
        _size = 0;
    }

    void locked_clear()
    {
        _mutex.lock();
        try {
            clear();
        }
        catch (exc& e) {
            _mutex.unlock();
            throw e;
        }
        catch (std::exception& e) {
            _mutex.unlock();
            throw exc(e);
        }
        _mutex.unlock();
    }

    void set_max_size(size_t max_size)
    {
        _max_size = max_size;
        if (AUTO_SHRINK)
            shrink();
    }

    void locked_set_max_size(size_t max_size)
    {
        _mutex.lock();
        try {
            set_max_size(max_size);
        }
        catch (exc& e) {
            _mutex.unlock();
            throw e;
        }
        catch (std::exception& e) {
            _mutex.unlock();
            throw exc(e);
        }
        _mutex.unlock();
    }

    void shrink()
    {
        while (_size > _max_size) {
            remove_lru_element();
        }
    }

    void locked_shrink()
    {
        _mutex.lock();
        try {
            shrink();
        }
        catch (exc& e) {
            _mutex.unlock();
            throw e;
        }
        catch (std::exception& e) {
            _mutex.unlock();
            throw exc(e);
        }
        _mutex.unlock();
    }

#ifndef NDEBUG
    bool check()
    {
        assert(_size <= _max_size);
        assert(_element_map.size() == _timestamp_map.size());
        for (typename std::map<KEY_TYPE, value>::iterator ei = _element_map.begin(); ei != _element_map.end(); ei++) {
            unsigned long long t = ei->second.timestamp;
            typename std::map<unsigned long long, KEY_TYPE>::iterator ti = _timestamp_map.find(t);
            assert(ti != _timestamp_map.end());
            assert(!(ti->second < ei->first) && !(ei->first < ti->second));    // same key
        }
        return true;
    }

    bool locked_check()
    {
        bool r;
        _mutex.lock();
        try {
            r = check();
        }
        catch (exc& e) {
            _mutex.unlock();
            throw e;
        }
        catch (std::exception& e) {
            _mutex.unlock();
            throw exc(e);
        }
        _mutex.unlock();
        return r;
    }
#endif
};

#endif
