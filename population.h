#ifndef POPULATION_H
#define POPULATION_H

#include "individual.h"
#include "rndutils.hpp"
#include <vector>


struct pop_param
{
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(pop_param,
                                   number_of_inds,
                                   mut_rate,
                                   mut_step)
    int number_of_inds;
    double mut_rate;
    double mut_step;
};

template <class Ind = individual<>>
class population
{
public:

    using ind_t = Ind;

    population(int init_nr_indiv = 1,
               double mut_rate = 0.01,
               double mut_step = 0.1,
               std::vector<int> net_arch = {1,2,1});

    population(const pop_param &p_p,const ind_param& i_p):
        m_vec_indiv(p_p.number_of_inds, Ind{i_p}),
        m_vec_new_indiv(p_p.number_of_inds, Ind{i_p}),
        m_mut_rate{p_p.mut_rate},
        m_mut_step{p_p.mut_step}
    {}

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(population,
                                   m_vec_indiv,
                                   m_mut_rate,
                                   m_mut_step);

    ///Changes the network of the nth individual to a given network
    template<class Net>
    void change_nth_ind_net(size_t ind_index, const Net& n){
        m_vec_indiv[ind_index].change_net(n);
    }


    ///Get const ref to vector of individuals
    const std::vector<Ind>& get_inds() const noexcept{return m_vec_indiv;}

    ///Get ref to vector of individuals
    std::vector<Ind>& get_inds() noexcept{return m_vec_indiv;}

    ///Returns the ref tot the mutable fitness distribution
    rndutils::mutable_discrete_distribution<>& get_fitness_dist() noexcept{return m_fitness_dist;}
    ///Get const ref to vector of individuals
    const std::vector<Ind>& get_new_inds() const noexcept{return m_vec_new_indiv;}

    ///Get ref to vector of individuals
    std::vector<Ind>& get_new_inds() noexcept{return m_vec_new_indiv;}

    ///Return mutation rate
    double get_mut_rate() const noexcept {return m_mut_rate;}

    ///Return mutation step
    double get_mut_step() const noexcept {return m_mut_step;}

private:

    std::vector<Ind> m_vec_indiv;
    std::vector<Ind> m_vec_new_indiv;
    double m_mut_rate;
    double m_mut_step;
    rndutils::mutable_discrete_distribution<> m_fitness_dist;

};

///Checks that 2 populations are equal
template< class Ind>
bool operator== (const population<Ind>& lhs, const population<Ind>& rhs)
{
    bool inds = lhs.get_inds() == rhs.get_inds();
    bool mut_rate = are_equal_with_tolerance(lhs.get_mut_rate(), rhs.get_mut_rate());
    bool mut_step = are_equal_with_tolerance(lhs.get_mut_step(), rhs.get_mut_step());

    return inds && mut_rate && mut_step;
}
namespace pop {



///Calculates the avg_fitness of the population
template< class Ind>
double avg_fitness(const population<Ind>& p){
    auto fitnesses = extract_fitnesses(p.get_inds());
    return calc_mean(fitnesses);
}

///Rescales the distance fro the target of an ind
///to a fitness value between 0  and 1
std::vector<double> rescale_dist_to_fit(std::vector<double> distance_from_target,
                                        double selection_strength);

///Checks that all individuals in the pop have the same input
template< class Ind>
bool all_individuals_have_same_input(const population<Ind> &p)
{
    std::vector<double> input_first_individual = p.get_inds()[0].get_input_values();

    for(auto& ind : p.get_inds()){
        if(input_first_individual != ind.get_input_values()){
            return false;
        }
    }
    return true;
}

///Assign inputs to a population
template<class Pop>
void assign_new_inputs_to_inds(Pop &p, const std::vector<double> &inputs)
{
    for(auto& ind : p.get_inds()){
        ind.assign_input(inputs);
    }
}

///Calculates the distance from the output of one individual's network
/// to the optimal output given a series of inputs
template<class Ind>
std::vector<double> calc_dist_from_target(const std::vector<Ind>& inds, double env_value)
{
    std::vector<double> distance_from_target;

    for(const auto& ind : inds)
    {
        auto sqr_distance = calc_sqr_distance(ind, env_value);
        distance_from_target.push_back(sqr_distance);
    }

    return distance_from_target;
}

///Sets the fitness of the nth ind in the population
template< class Ind>
void set_nth_ind_fitness (population<Ind>& p, size_t ind_index, double fitness)
{
    auto& ind = p.get_inds()[ind_index];
    ind.set_fitness(fitness);
}

///Sets the fitness of the individuals to the one contained in the fitness vector
template< class Ind>
void set_fitness_inds(population<Ind>& p, const std::vector<double>& fitness_vector)
{
    assert(p.get_inds().size() == fitness_vector.size());

    for(size_t i = 0; i != fitness_vector.size(); i++)
    {
        set_nth_ind_fitness(p, i, fitness_vector[i]);
    }
}

///Calculates the fitness of inds in pop given a target env_value
template< class Ind>
population<Ind>& calc_fitness(population<Ind>& p, const double& optimal_value,const double &sel_str)
{

    std::vector<double> distance_from_target = calc_dist_from_target(p.get_inds(), optimal_value);

    auto fitness_vector = rescale_dist_to_fit(distance_from_target, sel_str);

    set_fitness_inds(p, fitness_vector);

    return p;
}
///changes the net of the nth individual to a given net
template< class Ind>
void change_nth_ind_net(population<Ind>& p, size_t ind_index, typename Ind::net_t n);

///Creates a mutable distribution from whihc to draw inds based on fitness
template< class Ind>
rndutils::mutable_discrete_distribution<> create_mut_dist_fit(population<Ind>& p)
{
    rndutils::mutable_discrete_distribution<> mut_dist;

    mut_dist.mutate_transform(p.get_inds().begin(),
                              p.get_inds().end(),
                              [](const Ind& i){return i.get_fitness();});
    return  mut_dist;
}

///Extracts a vector of the fitnesses of individuals into a double vectors
template<class Ind>
std::vector<double> extract_fitnesses(const std::vector<Ind>& inds)
{
    std::vector<double> fitnesses;
    for(const auto& ind : inds)
    {
        fitnesses.push_back(ind.get_fitness());
    }
    return fitnesses;
}

///Gets the best n individuals in a pop
template< class Ind>
std::vector<Ind> get_best_n_inds(const population<Ind>& p, int nth)
{
    auto inds = p.get_inds();
    std::nth_element(inds.begin(), inds.begin() + nth, inds.end(),
                     [](const Ind& lhs, const Ind& rhs)
    {return lhs.get_fitness() > rhs.get_fitness();});

    return std::vector<Ind>(inds.begin(), inds.begin() + nth);
}


template< class Ind>
const Ind& get_nth_ind(const population<Ind>& p, size_t ind_index);

template< class Ind>
Ind& get_nth_ind(population<Ind>& p, size_t ind_index);

///Returns the fitness of the nth individual
template< class Ind>
double get_nth_ind_fitness(const population<Ind>& p, const size_t& ind_index)
{
    return p.get_inds()[ind_index].get_fitness();
}

template< class Ind>
const typename Ind::net_t& get_nth_ind_net(const population<Ind>& p, size_t ind_index);

///Select inds for new pop from old pop based on mutable dist
/// and mutates them
template< class Ind>
void select_new_pop(population<Ind>& p,
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

///Swaps a vector of new_inds with the vector of old inds
template< class Ind>
void swap_new_with_old_pop(population<Ind>& p)
{
    p.get_inds().swap(p.get_new_inds());
}

///Reproduces inds with a probability proportional to their fitness
template< class Ind>
void reproduce(population<Ind>& p, std::mt19937_64& rng)
{
    auto mut_dist = create_mut_dist_fit(p);

    select_new_pop(p, mut_dist, rng);

    swap_new_with_old_pop(p);
}

///Calculates the standard deviation
template< class Ind>
double var_fitness(const population<Ind> &p){
    auto inds = p.get_inds();
    auto fitnesses = extract_fitnesses(inds);
    return calc_stdev(fitnesses);
}

///Returns the input of the nth individual
template< class Ind>
const std::vector<double> &get_nth_individual_input(const population<Ind> &p, const int n)
{
  return get_nth_ind(p, n).get_input_values();
}


void test_population() noexcept;

}
#endif // POPULATION_H
