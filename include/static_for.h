#pragma once

template <int First, int Last> struct static_for {
template <typename Fn> void operator()(const Fn& fn) const {
    if constexpr (First < Last) {
        fn(std::integral_constant<int, First>{});
        static_for<First+1, Last>()(fn);
    }
}};

template <int N> struct static_for<N, N> {
template <typename Fn> void operator()(const Fn& fn) const { }
};
