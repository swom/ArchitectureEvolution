// Feed forward neural networks
//
// Hanno 2019


#ifndef ANN_HPP_INCLUDED
#define ANN_HPP_INCLUDED


#include <memory>
#include <type_traits>
#include <algorithm>
#include <array>
#include <tuple>


namespace ann {


  template <size_t size>
  struct state
  {
    static constexpr size_t state_size = size;
  };


  namespace activation {

    //! \brief zero activation function [0] <- u
    struct zero : state<0>
    {
      template <typename T> static constexpr T min = T(0);
      template <typename T> static constexpr T max = T(0);

      template <typename T>
      static T apply(T, const T* __restrict)
      {
        return T(0);
      }
    };


    //! \brief pass-through (no-operation) activation function [u] <- u
    struct identity : state<0>
    {
      template <typename T> static constexpr T min = -std::numeric_limits<T>::max();
      template <typename T> static constexpr T max = +std::numeric_limits<T>::max();

      template <typename T>
      static T apply(T u, const T* __restrict)
      {
        return u;
      }
    };


    //! \brief Hard limit (step) activation
    struct sgn
    {
      //! \brief bipolar hard limit activation [-1,1] <- u
      struct bipolar : state<0>
      {
        template <typename T> static constexpr T min = T(-1);
        template <typename T> static constexpr T max = T(+1);

        template <typename T>
        static T apply(T u, const T * __restrict)
        {
          return u > T(0) ? T(-1) : T(1);
        }
      };


      //! \brief unipolar hard limit activation [0,1] <- u
      struct unipolar : state<0>
      {
        template <typename T> static constexpr T min = T(0);
        template <typename T> static constexpr T max = T(1);

        template <typename T>
        static T apply(T u, const T* __restrict)
        {
          return u > T(0) ? T(0) : T(1);
        }
      };
    };


    //! \brief Rectified linear activation [0,u] <- u
    struct rtlu : state<0>
    {
      template <typename T> static constexpr T min = T(0);
      template <typename T> static constexpr T max = std::numeric_limits<T>::max();

      template <typename T>
      static T apply(T u, const T* __restrict)
      {
        return u = std::max(T(0), u);
      }
    };


    //! \brief hyperbolic activation
    struct tanh
    {
      //! brief bipolar tangent activation [-1,1] <- u
      struct bipolar : state<0>
      {
        template <typename T> static constexpr T min = T(-1);
        template <typename T> static constexpr T max = T(+1);

        template <typename T>
        static T apply(T u, const T * __restrict)
        {
          return std::tanh(u);
        }
      };


      //! brief unipolar tangent activation [0,1] <- u
      struct unipolar : state<0>
      {
        template <typename T> static constexpr T min = T(0);
        template <typename T> static constexpr T max = T(1);

        template <typename T>
        static T apply(T u, const T* __restrict)
        {
          return T(0.5)* (std::tanh(u) + T(1));
        }
      };

    };


    //! \brief arcus tangent activation
    struct arctan
    {
      static constexpr double pi = 3.14159265358979323846;    // still no std::math_constants in C++17 :(

      //! brief bipolar arcus tangent activation [-pi/2,pi/2] <- u
      struct bipolar : state<0>
      {
        template <typename T> static constexpr T min = T(-pi / 2.);
        template <typename T> static constexpr T max = T(+pi / 2.);

        template <typename T>
        static T apply(T u, const T * __restrict)
        {
          return std::atan(u);
        }
      };


      //! brief unipolar arcus tangent activation [0,pi] <- u
      struct unipolar : state<0>
      {
        template <typename T> static constexpr T min = T(0);
        template <typename T> static constexpr T max = T(pi);

        template <typename T>
        static T apply(T u, const T* __restrict)
        {
          return std::atan(u);
        }
      };
    };


    //! \brief sigmoid activation
    struct sig
    {
      //! \brief sigmoid activation [-1, 1] <- u
      //! \tparam n Nominator of the slope parameter
      //! \tparam d Denominator of the slope parameter
      template <int n = 1, int d = 1>
      struct bipolar : state<0>
      {
        template <typename T> static constexpr T min = T(-1);
        template <typename T> static constexpr T max = T(+1);

        template <typename T>
        static T apply(T u, const T * __restrict)
        {
          const T a = static_cast<T>(-n) / d;
          const auto Exp = std::exp(a * u);
          return (T(1) - Exp) / (T(1) + Exp);
        }
      };


      //! \brief sigmoid activation [0, 1] < u
      //! \tparam n Nominator of the slope parameter
      //! \tparam d Denominator of the slope parameter
      template <int n = 1, int d = 1>
      struct unipolar : state<0>
      {
        template <typename T> static constexpr T min = T(0);
        template <typename T> static constexpr T max = T(1);

        template <typename T>
        static T apply(T u, const T* __restrict)
        {
          const T a = static_cast<T>(-n) / d;
          return T(1) / (T(1) + std::exp(a * u));
        }
      };
    };


    //! \brief sigmoid activation with varying shape parameter
    //! activation_state[0]: shape parameter
    struct varsig
    {
      //! \brief sigmoid activation with varying shape parameter [-1,1] <- u
      struct bipolar : state<1>
      {
        template <typename T> static constexpr T min = T(-1);
        template <typename T> static constexpr T max = T(+1);

        template <typename T>
        static T apply(T u, const T * __restrict ps)
        {
          const auto Exp = std::exp(-ps[0] * u);
          return (T(1) - Exp) / (T(1) + Exp);
        }
      };


      //! \brief sigmoid activation with varying shape parameter [0,1] <- u
      struct unipolar : state<1>
      {
        template <typename T> static constexpr T min = T(0);
        template <typename T> static constexpr T max = T(1);

        template <typename T>
        static T apply(T u, const T* __restrict ps)
        {
          return T(1) / (T(1) + std::exp(-ps[0] * u));
        }
      };
    };


  }


  namespace feedback {

    //! \brief no feedback
    struct none : state<0>
    {
      template <typename T>
      static T apply(T u, T* __restrict)
      {
        return u;
      }
    };


    //! \brief direct feedback of \p last_u  * \p f
    //! feedback_state[0]: weight, feedback_state[1]: scratch
    struct direct : state<2>
    {
      template <typename T>
      static T apply(T u, T* __restrict ps)
      {
        return ps[1] = u + ps[1] * ps[0];
      }
    };

  }


  // fully connected neuron
  template <size_t Input,         
            typename Activation,
            typename Feedback = feedback::none,
            bool Biased = true
  >
  struct Neuron
  {
    using neuron_t = Neuron;
    using activation_t = Activation;
    using feedback_t = Feedback;

    // state layout: { [bias weight], #input weights, [activation state], [feedback state] }
    static constexpr size_t input_size = Input;
    static constexpr bool biased = Biased;
    static constexpr size_t ac_begin = (biased ? 1 : 0) + Input;
    static constexpr size_t fb_begin = ac_begin + activation_t::state_size;
    static constexpr size_t state_size = fb_begin + feedback_t::state_size;

    template <typename T>
    static auto feed(const std::array<T, Input>& in, T* __restrict state, size_t /* pos */)
      -> std::enable_if_t<biased, T>
    {
      auto u = state[0];    // bias
      for (size_t i = 1; i <= Input; ++i) {
        u += state[i] * in[i - 1];
      }
      return activation_t::apply(feedback_t::apply(u, state + fb_begin), state + ac_begin);
    }

    template <typename T>
    static auto feed(const std::array<T, Input>& in, T * __restrict state, size_t /* pos */)
      -> std::enable_if_t<!biased, T>
    {
      auto u = T(0);
      for (size_t i = 0; i < Input; ++i) {
        u += state[i] * in[i];
      }
      return activation_t::apply(feedback_t::apply(u, state + fb_begin), state + ac_begin);
    }
  };


  template <size_t Input,
            typename Activation,
            typename Feedback = feedback::none
  >
  using UnbiasedNeuron = Neuron<Input, Activation, Feedback, false>;


  template <typename Neuron, size_t N>
  struct Layer
  {
    using layer_t = Layer;
    using neuron_t = Neuron;

    static constexpr size_t size = N;                                   // number of neurons
    static constexpr size_t input_size = Neuron::input_size;            // inputs
    static constexpr size_t output_size = size;                         // outputs
    static constexpr size_t state_size = size * neuron_t::state_size;   // total state

    template <typename T> static constexpr T min_output = Neuron::template min<T>;
    template <typename T> static constexpr T max_output = Neuron::template max<T>;

    template <typename T>
    static auto feed(const std::array<T, input_size>& in, T* __restrict state)
    {
      std::array<T, output_size> out;
      T* pout = out.data();
      for (size_t i = 0; i < N; ++i, state += neuron_t::state_size) {
        pout[i] = typename neuron_t::feed(in, state, i);
      }
      return out;
    }
  };


  namespace detail {

    template <size_t I, typename L>
    struct accum_state_size
    {
      static constexpr size_t s0 = std::tuple_element_t<I, L>::state_size;
      static constexpr size_t value = s0 + accum_state_size<I - 1, L>::type::value;
    };


    template <typename L>
    struct accum_state_size<0, L> : std::integral_constant<size_t, std::tuple_element_t<0, L>::state_size>
    {};


    template <size_t I, typename L>
    struct accum_state_ofs
    {
      static constexpr size_t value = accum_state_size<I, L>::value - std::tuple_element_t<I, L>::state_size;
    };


    template <typename... L>
    struct network_state_size : accum_state_size<sizeof...(L) - 1, std::tuple<L...>>
    {};


    template <typename L0, typename L1, typename... L>
    struct match_interface_impl
      : std::integral_constant<bool, L0::output_size == L1::input_size && match_interface_impl<L1, L...>::value>
    {
    };


    template <typename L0, typename L1>
    struct match_interface_impl<L0, L1>
      : std::integral_constant<bool, L0::output_size == L1::input_size>
    {
    };


    template <size_t I, typename... L>
    struct match_interface : match_interface_impl<L...>
    {
    };


    template <typename... L>
    struct match_interface<1, L...> : std::integral_constant<bool, true>
    {
    };
  }


  // returned by Network::get_layer, Network::get_neuron
  template <typename T, typename U>
  class Proxy : public U
  {
  public:
    using value_type = std::remove_const_t<T>;
    using type = U;
    static constexpr size_t state_size = type::state_size;

    explicit Proxy(T* state) noexcept : state_(state) {}

    // linear state access functions
    const value_type* cbegin() const { return state_; }
    const value_type* cend() const { return state_ + state_size; }
    const value_type* begin() const { return state_; }
    const value_type* end() const { return state_ + state_size; }
    value_type* begin() { return state_; }
    value_type* end() { return state_ + state_size; }

  private:
    T* state_;
  };


  template <typename T, typename... L>
  class Network
  {
  public:
    static_assert(std::is_arithmetic<T>::value, "Network::value_type shall be arithmetic");
    static_assert(detail::match_interface<sizeof...(L), L...>::value, "Network: layer interfaces don't match");
    static constexpr size_t layers = sizeof...(L);

    using value_type = T;
    using network_t = Network;
    using layer_t = std::tuple<L...>;
    using input_layer_t = std::tuple_element_t<0, layer_t>;
    using output_layer_t = std::tuple_element_t<layers - 1, layer_t>;

    static constexpr size_t output_layer = layers - 1;
    static constexpr size_t input_size = input_layer_t::input_size;
    static constexpr size_t output_size = output_layer_t::output_size;
    static constexpr size_t state_size = detail::network_state_size<L...>::value;

    using input_t = std::array<value_type, input_size>;
    using output_t = std::array<value_type, output_size>;
    
    // creates network with indeterminate state
    Network() : state_{}
    {
    }

    // creates network with state == val
    explicit Network(value_type val)
      : Network()
    {
      for (auto& s : state_) { s = val; };
    }

    // linear state access functions
    const value_type* cbegin() const { return state_.data(); }
    const value_type* cend() const { return state_.data() + state_size; }
    const value_type* begin() const { return state_.data(); }
    const value_type* end() const { return state_.data() + state_size; }
    value_type* begin() { return state_.data(); }
    value_type* end() { return state_.data() + state_size; }

    // Layer access functions
    template <size_t I>
    auto get_layer() 
    {
      static_assert(I < std::tuple_size<layer_t>::value, "NetworkImpl::get_layer index out of range");
      using layerI = std::tuple_element_t<I, layer_t>;
      return Proxy<T,layerI>(state_.data() + detail::accum_state_ofs<I, layer_t>::value);
    }

    // Neuron access functions
    template <size_t I, size_t J>
    auto get_neuron()
    {
      auto layer = get_layer<I>();
      using layer_t  = decltype(layer);
      static_assert(J < layer_t::size, "NetworkImpl::get_neuron neuron index out of range");
      return Proxy<T, layer_t::type::neuron_t>(layer.begin() + J * layer_t::type::neuron_t::state_size);
    }

    // feed forward
    output_t operator()(const input_t & in)
    {
      return do_feed_forward<0>(in);
    }

    // feed forward
    template <typename ... INPUT>
    output_t operator()(INPUT... in)
    {
      static_assert(sizeof...(INPUT) == input_size, "NetworkImpl::operator(): illegal input pack");
      return do_feed_forward<0>({ value_type(in)... });
    }

  private:
    template <size_t I>
    auto do_feed_forward(const std::array<value_type, std::tuple_element_t<I, layer_t>::input_size> & in)
      -> std::enable_if_t<(I == output_layer), output_t>
    {
      return std::tuple_element_t<I, layer_t>::template feed<T>(in, state_.data() + detail::accum_state_ofs<I, layer_t>::value);
    }

    template <size_t I>
    auto do_feed_forward(const std::array<value_type, std::tuple_element_t<I, layer_t>::input_size> & in)
      ->std::enable_if_t<(I < output_layer), output_t>
    {
      return do_feed_forward<I + 1>(std::tuple_element_t<I, layer_t>::template feed<T>(in, state_.data() + detail::accum_state_ofs<I, layer_t>::value));
    }

    std::array<value_type, state_size> state_;
  };


  namespace detail {


    template <size_t I, size_t J, typename network_t, typename Visitor>
    inline auto do_visit_neuron(network_t& network, Visitor&& visitor)
      -> std::enable_if_t<J == std::tuple_element_t<I, typename network_t::layer_t>::size>
    {
    }


    template <size_t I, size_t J, typename network_t, typename Visitor>
    inline auto do_visit_neuron(network_t& network, Visitor&& visitor)
      ->std::enable_if_t<J < std::tuple_element_t<I, typename network_t::layer_t>::size>
    {
      auto n = network. template get_neuron<I, J>();
      using node_t = typename decltype(n)::type;
#if defined(_MSC_VER) && _MSC_VER <= 1900 && !defined(__c2__)
      visitor.operator()(n, I, J);
#else
      visitor.template operator()(n, I, J);
#endif
      do_visit_neuron<I, J + 1>(network, std::forward<Visitor>(visitor));
    }


    template <size_t I, typename network_t, typename Visitor>
    inline auto do_visit_layer(network_t& network, Visitor&& visitor)
      -> std::enable_if_t<I == network_t::output_layer>
    {
      do_visit_neuron<I, 0>(network, std::forward<Visitor>(visitor));
    }


    template <size_t I, typename network_t, typename Visitor>
    inline auto do_visit_layer(network_t& network, Visitor&& visitor)
      ->std::enable_if_t < I < network_t::output_layer>
    {
      do_visit_neuron<I, 0>(network, std::forward<Visitor>(visitor));
      do_visit_layer<I + 1, network_t, Visitor>(network, std::forward<Visitor>(visitor));
    }

  }


  // Apply visitor to all neurons in the network.
  // The visitor shall be a callable c++ object with the following 
  // signature:
  //
  // template <typename Proxy>
  // void operator()(Proxy proxy, size_t layer, size_t node)
  //
  template <typename network_t, typename Visitor>
  void visit_neurons(network_t& network, Visitor&& visitor)
  {
    detail::do_visit_layer<0>(network, std::forward<Visitor>(visitor));
  }


}

#endif
