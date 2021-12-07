#include "population.h"

#include <cassert>

template< mutation_type M>
population<M>::population(int init_nr_indiv,
                       double mut_rate,
                       double mut_step,
                       std::vector<int> net_arch
                       ):
    m_vec_indiv(static_cast<unsigned int>(init_nr_indiv),
                individual{ind_param{net_param{net_arch}}}),
    m_vec_new_indiv(static_cast<unsigned int>(init_nr_indiv),
                    individual{ind_param{net_param{net_arch}}}),
    m_mut_rate{mut_rate},
    m_mut_step{mut_step}
{}


template< mutation_type M>
population<M>::population(const pop_param &p_p,const ind_param& i_p):
    m_vec_indiv(static_cast<unsigned int>(p_p.number_of_inds)),
    m_vec_new_indiv(static_cast<unsigned int>(p_p.number_of_inds)),
    m_mut_rate{p_p.mut_rate},
    m_mut_step{p_p.mut_step}
{
    for(auto& ind : m_vec_indiv){ind = individual{i_p};}
    for(auto& ind : m_vec_new_indiv){ind = individual{i_p};}
}

std::vector<double> adjust_distances(std::vector<double> distances)
{
    for(double& dist : distances)
    {
        dist += 0.0000000000000001;
    }
    return distances;
}

template<mutation_type M>
bool all_nets_equals_to(const population<M>& p, const network<M>& n)
{
    return std::all_of(p.get_inds().begin(), p.get_inds().end(),
                       [n](const individual<M>& i)
    {return i.get_net() == n;});
}

template<mutation_type M>
rndutils::mutable_discrete_distribution<>  create_mut_dist_fit(population<M>& p)
{
    rndutils::mutable_discrete_distribution<> mut_dist;

    mut_dist.mutate_transform(p.get_inds().begin(),
                              p.get_inds().end(),
                              [](const individual<M>& i){return i.get_fitness();});
    return  mut_dist;
}

std::vector<double> create_rescaled_fitness_vec(std::vector<double> distance_from_target,
                                                double selection_strength)
{
    std::vector<double> fitness_inds;
    for(size_t i = 0; i != distance_from_target.size(); i++)
    {
        auto ind_fit = std::exp(-selection_strength * distance_from_target[i]);

        fitness_inds.push_back(ind_fit);
    }

    return fitness_inds;
}

template<mutation_type M>
void change_nth_ind_net(population<M>& p, size_t ind_index, network<M> n)
{
    p.change_nth_ind_net(ind_index, n);
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

template<mutation_type M>
void select_new_pop(population<M>& p,
                    const rndutils::mutable_discrete_distribution<>& mut_dist,
                    std::mt19937_64& rng)
{
    for( size_t i = 0; i != p.get_inds().size(); i++)
    {
        auto selected_ind_index = mut_dist(rng);
        auto selected_ind = p.get_inds()[selected_ind_index];
        p.get_new_inds()[i] = selected_ind;
        p.get_new_inds()[i].mutate(p.get_mut_rate(),
                                   p.get_mut_step(),
                                   rng);
    }
}

template< mutation_type M>
void swap_new_with_old_pop(population<M>& p)
{
    p.get_inds().swap(p.get_new_inds());
}

template< mutation_type M>
const individual<M> &get_nth_ind(const population<M>& p, size_t ind_index)
{
    return p.get_inds()[ind_index];
}

template< mutation_type M>
individual<M> &get_nth_ind(population<M>& p, size_t ind_index)
{
    return p.get_inds()[ind_index];
}

template<mutation_type M>
const network<M> &get_nth_ind_net(const population<M>& p, size_t ind_index)
{
    return get_nth_ind(p, ind_index).get_net();
}

std::vector<double> rescale_dist_to_fit(std::vector<double> distance_from_target,
                                        double selection_strength)
{
    auto fitness_inds = create_rescaled_fitness_vec(distance_from_target, selection_strength);

    return fitness_inds;
}

template<mutation_type M>
void reproduce(population<M>& p, std::mt19937_64& rng)
{
    auto mut_dist = create_mut_dist_fit(p);

    select_new_pop(p, mut_dist, rng);

    swap_new_with_old_pop(p);
}

template<mutation_type M>
void set_nth_ind_fitness (population<M>& p, size_t ind_index, double fitness)
{
    auto& ind = p.get_inds()[ind_index];
    ind.set_fitness(fitness);
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
        assert(are_equal_with_tolerance(p.get_mut_rate(), 0.01));
        assert(are_equal_with_tolerance(p.get_mut_step(), 0.1));

        auto mut_rate = 5.0;
        population p2{0, mut_rate};
        assert(are_equal_with_tolerance(p2.get_mut_rate(), mut_rate));

        auto mut_step = 5.0;
        population p3{0 ,0, mut_step};
        assert(are_equal_with_tolerance(p3.get_mut_step(), mut_step));
    }

    ///Population can be initialized with network architecture for inds
    {
        std::vector<int> net_arch{1,33,3,1};
        population p{1, 0, 0, net_arch};
        assert(get_nth_ind_net(p, 0) == network{net_arch});
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
        auto new_net =  change_all_weights_values_and_activations(get_nth_ind_net(p,first_ind), 123456);
        change_nth_ind_net(p, first_ind, new_net);

        set_nth_ind_fitness(p, first_ind, 1);
        set_nth_ind_fitness(p, second_ind, 0);

        reproduce(p, rng);

        assert(all_nets_equals_to(p, new_net));
    }
#endif

    //#define FIX_ISSUE_37

    {
        net_param net_par;

        ind_param i_p{net_par};

        int number_of_inds = 132;
        double mut_rate = 0.314;
        double mut_step = 0.1414;

        pop_param p_p{number_of_inds, mut_rate, mut_step};

        population p{p_p, i_p};

        for(const auto& ind : p.get_inds())
        {
            assert(ind.get_net() == network{net_par});
        }

        assert(are_equal_with_tolerance(p.get_inds().size(), number_of_inds) &
               are_equal_with_tolerance(p.get_mut_rate(), mut_rate) &
               are_equal_with_tolerance(p.get_mut_step(), mut_step));
    }

}
#endif
