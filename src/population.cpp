#include "population.h"
#include <cassert>

template< class Ind>
population<Ind>::population(int init_nr_indiv,
                            double mut_rate,
                            double mut_step,
                            std::vector<int> net_arch
                            ):
    m_vec_indiv(static_cast<unsigned int>(init_nr_indiv),
                individual{ind_param{net_param{net_arch, linear, net_arch}}}),
    m_vec_new_indiv(static_cast<unsigned int>(init_nr_indiv),
                    individual{ind_param{net_param{net_arch, linear, net_arch}}}),
    m_mut_rate_weight{mut_rate},
    m_mut_step{mut_step}
{}

bool operator==(const pop_param& lhs, const pop_param& rhs)
{
    bool n_inds = lhs.number_of_inds == rhs.number_of_inds;
    bool mut_rate_weights = lhs.mut_rate_weight == rhs.mut_rate_weight;
    bool mut_steps = lhs.mut_step == rhs.mut_step;
    bool mut_rate_activations = lhs.mut_rate_activation == rhs.mut_rate_activation;
    bool mut_rate_duplication = lhs.mut_rate_duplication == rhs.mut_rate_duplication;
    bool n_trials = lhs.n_trials == rhs.n_trials;

    return n_inds && mut_rate_weights && mut_steps && mut_rate_activations && mut_rate_duplication && n_trials;
}

namespace pop {
std::vector<double> adjust_distances(std::vector<double> distances)
{
    for(double& dist : distances)
    {
        dist += 0.0000000000000001;
    }
    return distances;
}

template< class Ind>
bool all_nets_equals_to(const population<Ind>& p, const typename Ind::net_t& n)
{
    return std::all_of(p.get_inds().begin(), p.get_inds().end(),
                       [n](const Ind& i)
    {return i.get_net() == n;});
}


std::vector<double> create_rescaled_fitness_vec(const std::vector<double>& distance_from_target,
                                                double selection_strength)
{
    std::vector<double> fitness_inds(distance_from_target.size());

#pragma omp parallel for
    for(int i = 0; i < distance_from_target.size(); i++)
    {
        auto ind_fit = std::exp(-selection_strength * distance_from_target[i]);

        fitness_inds[i] = ind_fit;
    }

    return fitness_inds;
}

void check_and_correct_dist(std::vector<double>& distance_from_target, double& min_distance)
{
    if(min_distance == 0)
    {
        distance_from_target =  adjust_distances(distance_from_target);

        min_distance = *std::min_element(distance_from_target.begin(),
                                         distance_from_target.end());
    }
}

template< class Ind>
const Ind &get_nth_ind(const population<Ind>& p, size_t ind_index)
{
    return p.get_inds()[ind_index];
}

template< class Ind>
Ind &get_nth_ind(population<Ind>& p, size_t ind_index)
{
    return p.get_inds()[ind_index];
}

template< class Ind>
const typename Ind::net_t& get_nth_ind_net(const population<Ind>& p, size_t ind_index)
{
    return get_nth_ind(p, ind_index).get_net();
}

std::vector<double> rescale_dist_to_fit(std::vector<double> distance_from_target,
                                        double selection_strength)
{
    auto fitness_inds = create_rescaled_fitness_vec(distance_from_target, selection_strength);

    return fitness_inds;
}

}

population<> produce_simple_pop()
{

    pop_param p_p;
    p_p.number_of_inds = 2;

    net_param n_p;
    n_p.max_arc = {1,1};
    n_p.net_arc = {1,1};

    ind_param i_p;
    i_p.net_par = n_p;

    population p{p_p,{}};

    std::mt19937_64 rng;
    p.get_inds_nonconst().at(0).mutate(1, 1, rng, 0, 0);

    return p;
}

#ifndef NDEBUG
void test_population() noexcept
{
    {
        int nelement = 10;
        population pop{10};
        assert (static_cast<int>(pop.get_inds().size()) == nelement );
    }
    //A population has a member variable called m_mut_step
    //And m_mut_rate
    //By default initialized to 0.01 (mut_rate)
    // and 0.1 (mut_step)
    {
        population p;
        assert(are_equal_with_tolerance(p.get_mut_rate_weight(), 0.01));
        assert(are_equal_with_tolerance(p.get_mut_step(), 0.1));

        auto mut_rate = 5.0;
        population p2{0, mut_rate};
        assert(are_equal_with_tolerance(p2.get_mut_rate_weight(), mut_rate));

        auto mut_step = 5.0;
        population p3{0, 0, mut_step};
        assert(are_equal_with_tolerance(p3.get_mut_step(), mut_step));
    }

    ///Population can be initialized with network architecture for inds
    {
        std::vector<int> net_arch{1,33,3,1};
        net_param n_p{net_arch, linear, net_arch};
        population p{{1, 0, 0, 0, 0}, {{n_p}}};
        assert(pop::get_nth_ind_net(p, 0) == network{n_p});
    }

    //Population has a buffer_vector for the new_population, with size equal to number of inds
    {
        population p;
        assert(p.get_new_inds().size() == p.get_inds().size());
    }

#define FIX_ISSUE_32
#ifdef FIX_ISSUE_32
    ///Individuals with higher fitness are preferentially selected for the next generation
    {
        int n_inds = 2;
        size_t first_ind = 0;
        size_t second_ind = 1;
        population p{n_inds};
        std::mt19937_64 rng;

        //make first ind net recognizable
        auto new_net =  change_all_weights_values_and_activations(pop::get_nth_ind_net(p,first_ind), 123456);
        pop::change_nth_ind_net(p, first_ind, new_net);

        pop::set_nth_ind_fitness(p, first_ind, 1);
        pop::set_nth_ind_fitness(p, second_ind, 0);

        pop::reproduce(p, rng);

        assert(pop::all_nets_equals_to(p, new_net));
    }
#endif

    //#define FIX_ISSUE_37

    {
        net_param net_par;

        ind_param i_p{net_par};

        int number_of_inds = 132;
        double mut_rate = 0.314;
        double mut_step = 0.1414;

        pop_param p_p{number_of_inds, mut_rate, mut_step, mut_rate, mut_rate};

        population p{p_p, i_p};

        for(const auto& ind : p.get_inds())
        {
            assert(ind.get_net() == network{net_par});
        }

        assert(are_equal_with_tolerance(double(p.get_inds().size()), number_of_inds) &
               are_equal_with_tolerance(p.get_mut_rate_weight(), mut_rate) &
               are_equal_with_tolerance(p.get_mut_step(), mut_step));
    }

    ///It is possible to obtain a vector of robustness measures from
    /// a vector of individuals
    {
        double robust_weight = 10;
        double frail_weight = 0;
        int n_mutations = 10;
        std::mt19937_64 rng;
        auto rng2 = rng;

        population robust_p;
        robust_p.change_all_inds_weights(robust_weight);

        population frail_p;
        assert(pop::all_inds_weights_have_value(frail_p, frail_weight));

        auto robust_inds_robustness = pop::calc_mutation_sensibility_all_inds(robust_p, n_mutations, rng);
        auto frail_inds_robustness = pop::calc_mutation_sensibility_all_inds(frail_p, n_mutations, rng2);

        assert(pairwise_comparison_for_majority(robust_inds_robustness,
                                                frail_inds_robustness));
    }

    ///Individuals are ranked based on fitness
    {
       auto p = produce_simple_pop();

        std::vector<double> input{1};
        double optimal_value = 1;
        double sel_str = 1;

        pop::calc_fitness(p, optimal_value, sel_str, input);

        assert(pop::all_fitnesses_are_not_equal(p.get_inds()));
        assert(pop::is_sorted_by_fitness(p.get_inds()) &&
               pop::is_sorted_by_rank(p.get_inds()));
    }

    ///The first rank is 0, the last is equal to the number of inds - 1
    {

    }
}
#endif
