// ColumnCache.hpp
/*
Copyright 2025 Arthur V. Morris
*/
#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <functional>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>

namespace py = pybind11;


class ColumnCache {
private:
    size_t _max_cols;
    size_t _packed_rows;
    std::list<size_t> _order;
    std::unordered_map<size_t, std::vector<uint64_t>> _cache;

public:
    ColumnCache(size_t max_bytes, size_t packed_rows)
        : _max_cols(max_bytes / (packed_rows * sizeof(uint64_t))), _packed_rows(packed_rows) {
            // std::cout << "ColumnCache initialized with max_cols: " << _max_cols << ", packed_rows: " << _packed_rows << std::endl;
        }

    const std::vector<uint64_t>& get(size_t col_idx,
                                     const std::function<std::vector<uint64_t>(size_t)>& builder) {
        // std::cout << "ColumnCache get called for column index: " << col_idx << std::endl;
        auto it = _cache.find(col_idx);
        if (it != _cache.end()) {
            // std::cout << "ColumnCache hit for column index: " << col_idx << std::endl;
            _order.remove(col_idx);
            _order.push_back(col_idx);
            return it->second;
        }

        if (_cache.size() >= _max_cols) {
            size_t oldest = _order.front();
            _order.pop_front();
            _cache.erase(oldest);
        }

        auto col = builder(col_idx);
        _order.push_back(col_idx);
        auto& inserted = _cache[col_idx] = std::move(col);
        return inserted;
    }

    void evict(size_t col_idx) {
        _order.remove(col_idx);
        _cache.erase(col_idx);
    }
};