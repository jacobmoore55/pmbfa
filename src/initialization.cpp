#include "pmbfa_types.h"

namespace pmbfa_init {


struct Candidate {
  Rcpp::List params;
  double score;
  Candidate(const Rcpp::List& p, double s) : params(p), score(s) {}
};

inline int option_int(const Rcpp::List& options, const char* name,
                      int fallback, int minimum) {
  if (!options.containsElementNamed(name) || Rf_isNull(options[name]))
    return fallback;
  double value = Rcpp::as<double>(options[name]);
  if (!std::isfinite(value) || value < minimum ||
      value > std::numeric_limits<int>::max() || std::floor(value) != value)
    Rcpp::stop("%s must be an integer >= %d.", name, minimum);
  return static_cast<int>(value);
}

inline double option_positive(const Rcpp::List& options, const char* name,
                               double fallback) {
  if (!options.containsElementNamed(name) || Rf_isNull(options[name]))
    return fallback;
  double value = Rcpp::as<double>(options[name]);
  if (!(value > 0.0) || !std::isfinite(value))
    Rcpp::stop("%s must be positive and finite.", name);
  return value;
}

inline int random_index(int size) {
  if (size <= 1) return 0;
  return std::min(size - 1, static_cast<int>(R::runif(0.0, size)));
}

// Calculate centers without constructing an observation-sized temporary array.
Rcpp::NumericMatrix centers(const Rcpp::NumericVector& x, int n, int p,
                            int N, int G, const std::vector<int>& assignment,
                            std::vector<int>& counts) {
  const int cells = n * p;
  Rcpp::NumericMatrix mean(cells, G);
  counts.assign(G, 0);
  std::vector<double> overall(cells, 0.0);
  for (int obs = 0; obs < N; ++obs) {
    const double* observation = x.begin() + static_cast<R_xlen_t>(obs) * cells;
    int group = assignment[obs];
    if (group >= 0) {
      ++counts[group];
      double* destination = &mean(0, group);
      for (int cell = 0; cell < cells; ++cell) destination[cell] += observation[cell];
    }
    for (int cell = 0; cell < cells; ++cell) overall[cell] += observation[cell] / N;
  }
  for (int group = 0; group < G; ++group) {
    double* destination = &mean(0, group);
    for (int cell = 0; cell < cells; ++cell)
      destination[cell] = counts[group] ? destination[cell] / counts[group] : overall[cell];
  }
  return mean;
}

// Balanced nearest-neighbor starts, retaining supplied labels.  With partial
// labels the allocation is proportional to known memberships, with at least
// one unknown assigned to a missing group whenever this is feasible.
std::vector<int> nearest_assignment(const Rcpp::NumericVector& x,
                                     int n, int p, int N, int G,
                                     const std::vector<int>& known) {
  const int cells = n * p;
  std::vector<int> assignment = known;
  std::vector<int> remaining, counts;
  remaining.reserve(N);
  for (int obs = 0; obs < N; ++obs)
    if (assignment[obs] < 0) remaining.push_back(obs);
  Rcpp::NumericMatrix mean = centers(x, n, p, N, G, assignment, counts);
  const int unknown_count = remaining.size();
  const int known_count = N - unknown_count;
  std::vector<int> allocation(G, 0);
  if (known_count == 0) {
    for (int group = 0; group < G - 1; ++group) allocation[group] = N / G;
    allocation[G - 1] = N - (G - 1) * (N / G);
  } else {
    int allocated = 0;
    std::vector<std::pair<double, int>> remainders;
    for (int group = 0; group < G; ++group) {
      double target = static_cast<double>(unknown_count) * counts[group] / known_count;
      allocation[group] = static_cast<int>(std::floor(target));
      allocated += allocation[group];
      remainders.emplace_back(target - allocation[group], group);
    }
    std::stable_sort(remainders.begin(), remainders.end(),
      [](const std::pair<double, int>& left, const std::pair<double, int>& right) {
        return left.first > right.first;
      });
    for (int extra = allocated; extra < unknown_count; ++extra)
      ++allocation[remainders[extra - allocated].second];
    for (int group = 0; group < G; ++group) {
      if (counts[group] > 0 || allocation[group] > 0) continue;
      int donor = -1;
      for (int other = 0; other < G; ++other) {
        int spare = allocation[other] - (counts[other] == 0 ? 1 : 0);
        if (spare > 0 && (donor < 0 || allocation[other] > allocation[donor]))
          donor = other;
      }
      if (donor >= 0) {
        --allocation[donor];
        ++allocation[group];
      }
    }
  }
  std::vector<std::pair<double, int>> distances;
  distances.reserve(unknown_count);
  for (int group = 0; group < G; ++group) {
    if (allocation[group] == 0 || remaining.empty()) continue;
    const double* reference;
    if (counts[group] > 0) reference = &mean(0, group);
    else reference = x.begin() + static_cast<R_xlen_t>(remaining[random_index(remaining.size())]) * cells;
    distances.clear();
    for (int obs : remaining) {
      const double* observation = x.begin() + static_cast<R_xlen_t>(obs) * cells;
      double distance = 0.0;
      for (int cell = 0; cell < cells; ++cell) {
        double difference = observation[cell] - reference[cell];
        distance += difference * difference;
      }
      distances.emplace_back(distance, obs);
    }
    int take = std::min(allocation[group], static_cast<int>(remaining.size()));
    std::partial_sort(distances.begin(), distances.begin() + take, distances.end(),
      [](const std::pair<double, int>& left, const std::pair<double, int>& right) {
        return left.first < right.first ||
          (left.first == right.first && left.second < right.second);
      });
    for (int index = 0; index < take; ++index) assignment[distances[index].second] = group;
    remaining.erase(std::remove_if(remaining.begin(), remaining.end(),
      [&assignment](int obs) { return assignment[obs] >= 0; }), remaining.end());
  }
  return assignment;
}

// Repair empty k-means groups using only unlabeled observations.  If labels
// make a nonempty group impossible, its center subsequently uses the mean.
void repair_empty(std::vector<int>& assignment, const std::vector<int>& known,
                   int G) {
  std::vector<int> counts(G, 0);
  for (int group : assignment) if (group >= 0) ++counts[group];
  for (int group = 0; group < G; ++group) {
    if (counts[group] > 0) continue;
    std::vector<int> donors;
    for (int obs = 0; obs < static_cast<int>(assignment.size()); ++obs)
      if (known[obs] < 0 && assignment[obs] >= 0 && counts[assignment[obs]] > 1)
        donors.push_back(obs);
    if (donors.empty()) continue;
    int obs = donors[random_index(donors.size())];
    --counts[assignment[obs]];
    assignment[obs] = group;
    ++counts[group];
  }
}

Rcpp::List make_params(const Rcpp::NumericVector& x, int n, int p, int N,
                        int G, int q, int r, const std::string& family,
                        const Rcpp::NumericVector& labels,
                        const std::vector<int>& assignment,
                        double variance_floor, double nu_min, double shape_min) {
  const int cells = n * p;
  std::vector<int> counts;
  Rcpp::NumericMatrix mean = centers(x, n, p, N, G, assignment, counts);
  Rcpp::NumericMatrix sigma(n, G), psi(p, G);
  std::fill(sigma.begin(), sigma.end(), std::max(1.0, variance_floor));
  std::fill(psi.begin(), psi.end(), std::max(1.0, variance_floor));
  Rcpp::NumericVector lambda(static_cast<R_xlen_t>(G) * n * q);
  Rcpp::NumericVector delta(static_cast<R_xlen_t>(G) * p * r);
  lambda.attr("dim") = Rcpp::IntegerVector::create(G, n, q);
  delta.attr("dim") = Rcpp::IntegerVector::create(G, p, r);
  // The screening model has common loading matrices, just as in the original.
  for (int factor = 0; factor < q; ++factor)
    for (int row = 0; row < n; ++row) {
      double value = R::runif(-1.0, 1.0);
      for (int group = 0; group < G; ++group)
        lambda[group + static_cast<R_xlen_t>(G) * (row + n * factor)] = value;
    }
  for (int factor = 0; factor < r; ++factor)
    for (int col = 0; col < p; ++col) {
      double value = R::runif(-1.0, 1.0);
      for (int group = 0; group < G; ++group)
        delta[group + static_cast<R_xlen_t>(G) * (col + p * factor)] = value;
    }
  Rcpp::NumericVector pi(G);
  double weight_sum = 0.0;
  for (int group = 0; group < G; ++group) {
    pi[group] = std::max(1, counts[group]);
    weight_sum += pi[group];
  }
  pi = pi / weight_sum;
  Rcpp::List args(G);
  for (int group = 0; group < G; ++group) {
    Rcpp::NumericMatrix location(n, p);
    std::copy(&mean(0, group), &mean(0, group) + cells, location.begin());
    Rcpp::List component = Rcpp::List::create(Rcpp::_["M"] = location);
    if (family != "Normal") {
      Rcpp::NumericMatrix skew(n, p);
      std::fill(skew.begin(), skew.end(), 0.1);
      component["A"] = skew;
    }
    if (family == "ST") component["nu"] = std::max(10.0, nu_min);
    if (family == "GH") {
      component["omega"] = 50.0;
      component["lambda"] = 2.0;
    }
    if (family == "NIG") component["gam"] = std::max(1.0, shape_min);
    if (family == "VG") component["gamma"] = std::max(1.0, shape_min);
    args[group] = component;
  }
  Rcpp::List moddat = Rcpp::List::create(
    Rcpp::_["N"] = N, Rcpp::_["G"] = G, Rcpp::_["n"] = n,
    Rcpp::_["p"] = p, Rcpp::_["q"] = q, Rcpp::_["r"] = r,
    Rcpp::_["name"] = family == "Normal" ? "CCC-CCC" : "CCC-CCC-U-C",
    Rcpp::_["family"] = family, Rcpp::_["class"] = Rcpp::clone(labels));
  Rcpp::List gmp = Rcpp::List::create(
    Rcpp::_["Sigmak"] = sigma, Rcpp::_["Psik"] = psi,
    Rcpp::_["Lambdak"] = lambda, Rcpp::_["Deltak"] = delta,
    Rcpp::_["pik"] = pi);
  return Rcpp::List::create(Rcpp::_["moddat"] = moddat,
                            Rcpp::_["gmp"] = gmp, Rcpp::_["args"] = args);
}

double fit_score(const Rcpp::List& fit,double minimum_component) {
  if (!fit.containsElementNamed("loglik") || Rf_isNull(fit["loglik"]))
    return -std::numeric_limits<double>::infinity();
  Rcpp::NumericVector history = fit["loglik"];
  if (history.size() == 0 || !std::isfinite(history[history.size() - 1]))
    return -std::numeric_limits<double>::infinity();
  if(fit.containsElementNamed("diagnostics")){
    Rcpp::List diagnostics=fit["diagnostics"];
    if(diagnostics.containsElementNamed("effective_sizes")){
      Rcpp::NumericVector sizes=diagnostics["effective_sizes"];
      for(double size:sizes)if(!std::isfinite(size)||size<=minimum_component)
        return -std::numeric_limits<double>::infinity();
    }
  }
  return history[history.size() - 1];
}

void keep_candidate(std::vector<Candidate>& bank, const Rcpp::List& params,
                      double score, int capacity) {
  if (!std::isfinite(score)) return;
  if (static_cast<int>(bank.size()) < capacity) bank.emplace_back(params, score);
  else {
    auto worst = std::min_element(bank.begin(), bank.end(),
      [](const Candidate& left, const Candidate& right) { return left.score < right.score; });
    if (score > worst->score) *worst = Candidate(params, score);
  }
}


Rcpp::List initialize(
    const Rcpp::NumericVector& x, const int G, const int q, const int r,
    const std::string family,
    Rcpp::List aargs,
    const Rcpp::NumericVector labels) {
  using namespace pmbfa_init;
  if (Rf_isNull(aargs)) aargs = Rcpp::List::create();
  Rcpp::IntegerVector dimensions = x.attr("dim");
  if (dimensions.size() != 3 || dimensions[0] < 1 || dimensions[1] < 1 || dimensions[2] < 1)
    Rcpp::stop("x must be a nonempty numeric array with dimensions (n, p, N).");
  int n = dimensions[0], p = dimensions[1], N = dimensions[2];
  if (G < 1 || G > N || q < 1 || q >= n || r < 1 || r >= p)
    Rcpp::stop("Require 1 <= G <= N, 1 <= q < n, and 1 <= r < p.");
  if (family != "Normal" && family != "ST" && family != "GH" &&
      family != "NIG" && family != "VG") Rcpp::stop("Unknown distribution family.");
  for (double value : x) if (!std::isfinite(value)) Rcpp::stop("x must contain finite values.");
  if (labels.size() != 0 && labels.size() != N)
    Rcpp::stop("labels must be empty or have length N.");
  Rcpp::NumericVector fixed_labels(N, 0.0);
  std::vector<int> known(N, -1);
  if (labels.size()) {
    for (int obs = 0; obs < N; ++obs) {
      double label = labels[obs];
      if (!std::isfinite(label) || label < 0.0 || label > G || label != std::floor(label))
        Rcpp::stop("labels must be integers in 0:G, with 0 denoting unknown membership.");
      fixed_labels[obs] = label;
      if (label > 0.0) known[obs] = static_cast<int>(label) - 1;
    }
  }
  int ncandid = option_int(aargs, "ncandid", 200, 1);
  int nstart = option_int(aargs, "nstart", 20, 1);
  int nrep = option_int(aargs, "nrep", 1, 1);
  int init_iter = option_int(aargs, "init.iter", 5, 0);
  if (nrep > nstart || nstart > ncandid)
    Rcpp::stop("Require nrep <= nstart <= ncandid.");
  double variance_floor = option_positive(aargs, "variance.floor", 0.1);
  double nu_min = option_positive(aargs, "nu.min", 4.1);
  double shape_min = option_positive(aargs, "shape.min", 0.5);
  double minimum_component=pm::option(aargs,"min.component",6.0);
  if(!std::isfinite(minimum_component)||minimum_component<0)
    Rcpp::stop("min.component must be finite and nonnegative.");
  std::string method = "nearest";
  if (aargs.containsElementNamed("init_method") && !Rf_isNull(aargs["init_method"]))
    method = Rcpp::as<std::string>(aargs["init_method"]);
  if (method != "nearest" && method != "kmeans")
    Rcpp::stop("init_method must be 'nearest' or 'kmeans'.");
  Rcpp::List scoring_options = Rcpp::clone(aargs);
  scoring_options["max.iter"] = 0;
  scoring_options["verbose"] = false;
  std::vector<Candidate> screened, refined;
  screened.reserve(nstart);
  refined.reserve(nrep);
  Rcpp::List fallback;
  int failed_screening = 0, failed_refinement = 0, kmeans_failures = 0;
  Rcpp::NumericMatrix flat;
  if (method == "kmeans") {
    flat = Rcpp::NumericMatrix(N, n * p);
    for (int obs = 0; obs < N; ++obs)
      for (int cell = 0; cell < n * p; ++cell)
        flat(obs, cell) = x[cell + static_cast<R_xlen_t>(n) * p * obs];
  }
  for (int candidate = 0; candidate < ncandid; ++candidate) {
    Rcpp::checkUserInterrupt();
    std::vector<int> assignment = nearest_assignment(x, n, p, N, G, known);
    if (method == "kmeans") {
      std::vector<int> counts;
      Rcpp::NumericMatrix means = centers(x, n, p, N, G, assignment, counts);
      Rcpp::NumericMatrix starting_centers(G, n * p);
      for (int group = 0; group < G; ++group)
        for (int cell = 0; cell < n * p; ++cell)
          starting_centers(group, cell) = means(cell, group);
      try {
        Rcpp::Function kmeans = Rcpp::Environment::namespace_env("stats")["kmeans"];
        Rcpp::Language km_call(kmeans, Rcpp::_["x"] = flat,
          Rcpp::_["centers"] = starting_centers, Rcpp::_["iter.max"] = 20,
          Rcpp::_["algorithm"] = "Lloyd");
        Rcpp::Function suppress = Rcpp::Environment::base_env()["suppressWarnings"];
        Rcpp::Language quiet_call(suppress, km_call);
        Rcpp::List result = quiet_call.eval();
        Rcpp::IntegerVector cluster = result["cluster"];
        for (int obs = 0; obs < N; ++obs)
          assignment[obs] = known[obs] >= 0 ? known[obs] : cluster[obs] - 1;
        repair_empty(assignment, known, G);
      } catch (Rcpp::internal::InterruptedException&) { throw; }
        catch (std::exception&) { ++kmeans_failures; }
    }
    Rcpp::List params = make_params(x, n, p, N, G, q, r, family, fixed_labels,
                                   assignment, variance_floor, nu_min, shape_min);
    if (fallback.size() == 0) fallback = params;
    try {
      Rcpp::List fit = family == "Normal" ? pm::fit(x, params, scoring_options, true)
                                           : pm::fit(x, params, scoring_options, false);
      double score = fit_score(fit,minimum_component);
      if (std::isfinite(score)) {
        Rcpp::List accepted = fit.containsElementNamed("params") ?
          Rcpp::as<Rcpp::List>(fit["params"]) : params;
        keep_candidate(screened, accepted, score, nstart);
      } else ++failed_screening;
    } catch (Rcpp::internal::InterruptedException&) { throw; }
      catch (std::exception&) { ++failed_screening; }
  }
  std::stable_sort(screened.begin(), screened.end(),
    [](const Candidate& left, const Candidate& right) { return left.score > right.score; });
  scoring_options["max.iter"] = init_iter;
  for (const Candidate& candidate : screened) {
    Rcpp::checkUserInterrupt();
    try {
      Rcpp::List fit = family == "Normal" ? pm::fit(x, candidate.params, scoring_options, true)
                                           : pm::fit(x, candidate.params, scoring_options, false);
      double score = fit_score(fit,minimum_component);
      if (std::isfinite(score)) {
        Rcpp::List accepted = fit.containsElementNamed("params") ?
          Rcpp::as<Rcpp::List>(fit["params"]) : candidate.params;
        keep_candidate(refined, accepted, score, nrep);
      } else ++failed_refinement;
    } catch (Rcpp::internal::InterruptedException&) { throw; }
      catch (std::exception&) { ++failed_refinement; }
  }
  // A failed short fit must not erase a successfully scored raw initialization.
  bool used_raw_fallback = refined.empty();
  if (refined.empty()) {
    for (const Candidate& candidate : screened)
      keep_candidate(refined, candidate.params, candidate.score, nrep);
  }
  std::stable_sort(refined.begin(), refined.end(),
    [](const Candidate& left, const Candidate& right) { return left.score > right.score; });
  Rcpp::List output(nrep);
  Rcpp::NumericVector scores(nrep, NA_REAL);
  for (int start = 0; start < nrep; ++start) {
    if (refined.empty()) output[start] = Rcpp::clone(fallback);
    else {
      int index = start % refined.size();
      output[start] = Rcpp::clone(refined[index].params);
      scores[start] = refined[index].score;
    }
  }
  output.attr("initialization") = Rcpp::List::create(
    Rcpp::_["method"] = method, Rcpp::_["ncandid"] = ncandid,
    Rcpp::_["screened"] = static_cast<int>(screened.size()),
    Rcpp::_["failed_screening"] = failed_screening,
    Rcpp::_["failed_refinement"] = failed_refinement,
    Rcpp::_["kmeans_failures"] = kmeans_failures,
    Rcpp::_["unscored_fallback"] = refined.empty(),
    Rcpp::_["raw_fallback"] = used_raw_fallback,
    Rcpp::_["scores"] = scores);
  return output;
}

} // namespace pmbfa_init
