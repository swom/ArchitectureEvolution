#include "population.h"

#include <cassert>

population::population(int init_nr_indiv,
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


population::population(const pop_param &p_p,const ind_param& i_p):
    m_vec_indiv(static_cast<unsigned int>(p_p.number_of_inds)),
    m_vec_new_indiv(static_cast<unsigned int>(p_p.number_of_inds)),
    m_mut_rate{p_p.mut_rate},
    m_mut_step{p_p.mut_step}
{
    for(auto& ind : m_vec_indiv){ind = individual{i_p};}
    for(auto& ind : m_vec_new_indiv){ind = individual{i_p};}
}



bool operator== (const population& lhs, const population& rhs)
{
    bool inds = lhs.get_inds() == rhs.get_inds();
    bool mut_rate = are_equal_with_tolerance(lhs.get_mut_rate(), rhs.get_mut_rate());
    bool mut_step = are_equal_with_tolerance(lhs.get_mut_step(), rhs.get_mut_step());

    return inds && mut_rate && mut_step;
}

double avg_fitness(const population& p)
{
    auto fitnesses = extract_fitnesses(p.get_inds());
    return calc_mean(fitnesses);
}

std::vector<double> adjust_distances(std::vector<double> distances)
{
    for(double& dist : distances)
    {
        dist += 0.0000000000000001;
    }
    return distances;
}

bool all_nets_equals_to(const population& p, const network& n)
{
    return std::all_of(p.get_inds().begin(), p.get_inds().end(),
                       [n](const individual& i)
    {return i.get_net() == n;});
}

std::vector<double> calc_dist_from_target(const std::vector<individual>& inds, double env_value)
{
    std::vector<double> distance_from_target;

    for(const auto& ind : inds)
    {
        auto sqr_distance = calc_sqr_distance(ind, env_value);
        distance_from_target.push_back(sqr_distance);
    }

    return distance_from_target;
}

void population::calc_fitness(const double& optimal_value,const double &sel_str)
{

    std::vector<double> distance_from_target = calc_dist_from_target(get_inds(), optimal_value);

    auto fitness_vector = rescale_dist_to_fit(distance_from_target, sel_str);

    set_fitness_inds(*this, fitness_vector);
}

void population::assign_inputs_to_inds(const std::vector<double> &inputs)
{
  for(auto& ind : m_vec_indiv){
      ind.assign_input(inputs);
  }
}

rndutils::mutable_discrete_distribution<>  create_mut_dist_fit(population& p)
{
    rndutils::mutable_discrete_distribution<> mut_dist;

    mut_dist.mutate_transform(p.get_inds().begin(),
                              p.get_inds().end(),
                              [](const individual& i){return i.get_fitness();});
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


void check_and_correct_dist(std::vector<double>& distance_from_target, double& min_distance)
{
    if(min_distance == 0)
    {
        distance_from_target =  adjust_distances(distance_from_target);

        min_distance = *std::min_element(distance_from_target.begin(),
                                         distance_from_target.end());
    }
}

std::vector<double> extract_fitnesses(const std::vector<individual>& inds)
{
    std::vector<double> fitnesses;
    for(const auto& ind : inds)
    {
        fitnesses.push_back(ind.get_fitness());
    }
    return fitnesses;
}

void select_new_pop(population& p,
                    const rndutils::mutable_discrete_distribution<>& mut_dist,
                    std::mt19937_64& rng)
{
    for( size_t i = 0; i != p.get_inds().size(); i++)
    {
        auto selected_ind_index = mut_dist(rng);
        auto selected_ind = p.get_inds()[selected_ind_index];
        p.change_nth_new_ind(i, selected_ind);
        individual mutated_ind = p.get_new_inds()[i] ;
        mutated_ind.mutate(p.get_mut_rate(),
                            p.get_mut_step(),
                            rng);

        p.change_nth_new_ind(i, mutated_ind);

    }
}

void swap_new_with_old_pop(population& p)
{
    std::vector<individual> inds = p.get_inds();
    std::vector<individual> new_inds = p.get_new_inds();
    inds.swap(new_inds);
    p.change_vector_individuals(inds);
    p.change_vector_new_inds(new_inds);
}

std::vector<individual> get_best_n_inds(const population& p, int nth)
{
    auto inds = p.get_inds();
    std::nth_element(inds.begin(), inds.begin() + nth, inds.end(),
                     [](const individual& lhs, const individual& rhs)
    {return lhs.get_fitness() > rhs.get_fitness();});

    return std::vector<individual>(inds.begin(), inds.begin() + nth);
}

const individual& get_nth_ind(const population& p, size_t ind_index)
{
    return p.get_inds()[ind_index];
}

double get_nth_ind_fitness(const population& p, const size_t& ind_index)
{
    return p.get_inds()[ind_index].get_fitness();
}

const network& get_nth_ind_net(const population& p, size_t ind_index)
{
    return get_nth_ind(p, ind_index).get_net();
}

void population::change_nth_ind_net(size_t ind_index, const network& n)
{
    m_vec_indiv[ind_index].change_net(n);
}

std::vector<double> rescale_dist_to_fit(std::vector<double> distance_from_target,
                                        double selection_strength)
{
    auto fitness_inds = create_rescaled_fitness_vec(distance_from_target, selection_strength);

    return fitness_inds;
}

void population::reproduce(std::mt19937_64& rng)
{
    auto mut_dist = create_mut_dist_fit(*this);

    select_new_pop(*this, mut_dist, rng);

    swap_new_with_old_pop(*this);
}

void set_fitness_inds(population& p, const std::vector<double>& fitness_vector)
{
    assert(p.get_inds().size() == fitness_vector.size());

    for(size_t i = 0; i != fitness_vector.size(); i++)
    {
        set_nth_ind_fitness(p, i, fitness_vector[i]);
    }
}

void set_nth_ind_fitness (population& p, size_t ind_index, double fitness)
{
    auto ind = p.get_inds()[ind_index];
    ind.set_fitness(fitness);
    p.change_nth_individual(ind_index, ind);
}

double var_fitness(const population& p)
{
    auto inds = p.get_inds();
    auto fitnesses = extract_fitnesses(inds);
    return calc_stdev(fitnesses);
}

bool all_individuals_have_same_input(const population &p)
{
  std::vector<double> input_first_individual = p.get_inds()[0].get_input_values();

  for(auto& ind : p.get_inds()){
      if(input_first_individual != ind.get_input_values()){
          return false;
        }
    }
  return true;
}

const std::vector<double> &get_nth_individual_input(const population &p, const int n)
{
  return get_nth_ind(p, n).get_input_values();
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
        auto new_net =  change_all_weights(get_nth_ind_net(p,first_ind), 123456);
        p.change_nth_ind_net(first_ind, new_net);

        set_nth_ind_fitness(p, first_ind, 1);
        set_nth_ind_fitness(p, second_ind, 0);

        p.reproduce(rng);

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
