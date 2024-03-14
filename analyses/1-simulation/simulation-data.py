#!/usr/bin/env python
# AUTHOR: Stephen Blaskowski
# CREATE DATE: 13 December 2022

# simulated replicate resampling experiment

# imports
from concurrent.futures import ProcessPoolExecutor
import datetime
import json
import numpy as np
from pathlib import Path
import pandas as pd
import scipy
from sklearn.model_selection import ParameterGrid
import tensorly as tl
from tensorly import check_random_state
from tensorly.cp_tensor import CPTensor
from barnacle import (
    SparseCP, 
    simulated_sparse_tensor, 
    visualize_3d_tensor, 
    plot_factors_heatmap, 
    recovery_relevance, 
    pairs_precision_recall
)
from barnacle.tensors import SparseCPTensor
from tlab.cp_tensor import store_cp_tensor, load_cp_tensor
from tlviz.visualisation import optimisation_diagnostic_plots
from tlviz.model_evaluation import relative_sse, core_consistency
from tlviz.factor_tools import factor_match_score, cosine_similarity, degeneracy_score
from tlviz.multimodel_evaluation import similarity_evaluation
import xarray as xr


# function to return random replicate labelings
def generate_replicate_labels(sample_names, random_state=None, replicate_map=None):
    '''Generates random replicate labels to align with an input vector of sample names.
    
    Parameters
    ----------
    sample_names : np.ndarray
        Input array of sample names. Must be sorted in ascending order.
    random_state : {None, int, numpy.random.RandomState}
        Default is None.
    replicate_map : {None, dict}
        Map of integer labels to preferred replicate labels.
        Example:
            {0:'A', 1:'B', 2:'C'}
        Default is None.
            
    Returns
    -------
    replicate_labels : np.ndarray
        Array of randomly generated replicate lables, to be aligned with input `sample_names`.
    
    '''
    # check that input is a numpy array
    if not type(sample_names) is np.ndarray:
        raise AssertionError('`sample_names` must be a numpy.ndarray')
    # check sample_names are sorted
    if not np.all(sample_names[:-1] <= sample_names[1:]):
        raise AssertionError('`sample_names` must be sorted in ascending order')
    # get counts of each sample name
    names, counts = np.unique(sample_names, return_counts=True)
    # get random state
    rns = check_random_state(random_state)
    # generate replicate labels
    replicate_labels = [rns.choice(np.arange(counts.max()), size=c, replace=False) for c in counts]
    replicate_labels = np.concatenate(replicate_labels)
    # map preferred replicate labels
    if replicate_map is not None:
        mapped_replicate_labels = [replicate_map[i] for i in replicate_labels]
        replicate_labels = np.array(mapped_replicate_labels)
    # return result
    return replicate_labels


# function to separate out subtensors of each replicate
def separate_replicates(dataset, coordinates, data_variable, replicate_label='replicate'):
    '''Separates data from each replicate set into its own independent DataArray.
    
    Parameters
    ----------
    dataset : xarray.Dataset
        Dataset with replicates
    coordinates : list of str
        Coordinates to be preserved in each replicate set.
    data_variable : str
        Name of data variable to be selected in each replicate set.
    replicate_label : str, default is 'replicate'
        Label of replicate field in `dataset`. 
 
    Returns
    -------
    replicate_sets : dict of xarray.DataArrays
        Set of replicate DataArrays, each keyed on its replicate label.
    
    '''
    # get list of replicate labels
    replicates = np.unique(dataset[replicate_label].to_numpy())
    # pull out each replicate subset
    subsets = list()
    for rep in replicates:
        # pull out only data with of the particular replicate
        df = dataset.where(dataset[replicate_label] == rep, drop=True).to_dataframe().reset_index()
        # reform DataArray with specified coordinates and data_variable
        rep_da = xr.DataArray.from_series(df.set_index(coordinates)[data_variable])
        # add to dict of replicate subsets
        subsets.append(rep_da)
    # arrange in dictionary and return
    return dict(zip(replicates, subsets))


# function to select common indices between two datasets
def select_common_indices(dataset_1, dataset_2, coordinates):
    '''Finds common indices between two datasets.
    
    Parameters
    ----------
    dataset_1 : xarray.Dataset
        Dataset with common coordinates to be compared.
    dataset_2 : xarray.Dataset
        Dataset with common coordinates to be compared.
    coordinates : list of str
        Coordinates to be compared between `dataset_1` and `dataset_2`.
        
    Returns
    -------
    common_index_labels : list of numpy.Arrays
        Common indices' labels, one per coordinate provided.
    indices_1 : list of numpy.Arrays of ints
        Common integer indices, one per coordinate provided.
    indices_2 : list of numpy.Arrays of ints
        Common integer indices, one per coordinate provided.
    '''
    # initialize outputs
    if len(coordinates) > 1:
        common_index_labels = {}
        indices_1 = {}
        indices_2 = {}
    # loop through coordinates
    for coord in coordinates:
        # get shared coordinate labels
        shared_labels = np.intersect1d(
            dataset_1.indexes[coord], 
            dataset_2.indexes[coord], 
            assume_unique=True, 
            return_indices=False
        )
        # get dataset 1 index
        _, index_1, _ = np.intersect1d(
            dataset_1.indexes[coord], 
            shared_labels, 
            assume_unique=True, 
            return_indices=True
        )
        # get dataset 2 index
        _, index_2, _ = np.intersect1d(
            dataset_2.indexes[coord], 
            shared_labels, 
            assume_unique=True, 
            return_indices=True
        )
        # store labels and indices
        if len(coordinates) > 1: 
            common_index_labels[coord] = shared_labels
            indices_1[coord] = index_1
            indices_2[coord] = index_2
    # return results
    if len(coordinates) > 1: 
        return common_index_labels, indices_1, indices_2
    else:
        return shared_labels, index_1, index_2
    
    
# function to select subset of indices in cp tensor
def subset_cp_tensor(cp_tensor, subset_indices):
    '''Selects subset of cp_tensor based on provided indices
    
    Parameters
    ----------
    cp_tensor : tensorly.CPTensor
        CPTensor object with (weights, factors).
    subset_indices : dict(int: index-like)
        Dictionary with mode as key and value an integer index of 
        the positions to be downselected from `cp_tensor`.
        Example: {1: [0, 1, 3, 4, 5, 8]}
        
    Returns
    -------
    subset_cp : tensorly.CPTensor
        Subset CPTensor.
    '''
    weights, factors = cp_tensor
    new_factors = factors.copy()
    for mode, index in subset_indices.items():
        new_factors[mode] = factors[mode][index]
    return(CPTensor((weights, new_factors)))


def fit_save_model(model, data, path, fit_params):
    '''Helper function that takes an instantiated model and data as input,
    fits the model to the data, and returns the fit model. Optionally, the model
    and its settings can be saved to an input file path.
    
    Parameters
    ----------
    model : barnacle.SparseCP
        Instantiated and parameterized SparseCP model.
    data : numpy.ndarray
        Input data tensor.
    path : pathlib.Path
        Path directory where output will be saved. If path=None, no data will be saved.
        If a legitimate filepath is provided, the fit model, in addition to parameters 
        will be saved at the provided filepath.
    fit_params : dict
        Keyword arguments to be passed to the SparseCP.fit_transform() method. 
        Pass empty dictionary if no kwargs are to be passed. 
            
    Returns
    -------
    model : barnacle.SparseCP
        Fit model.
    '''
    if path is not None:
        # make path directories if they don't exist yet
        if not path.exists():
            path.mkdir(parents=True)
        # save parameters
        if model._best_cp_index is not None:
            raise AssertionError('The `model` passed has already been fit')
        else:
            with open(path / 'model_parameters.txt', 'w') as f:
                f.write(json.dumps(model.__dict__, indent=4))
    _ = model.fit_transform(data, return_losses=False, **fit_params)
    # save best fit model
    if path is not None:
        store_cp_tensor(model.decomposition_, path / 'fitted_model.h5')
    # return model
    return model


# function to count number of nonzero components in a cp tensor
def nonzero_components(cp, return_trimmed_cp=False):
    accumulator = np.ones_like(cp.weights)
    for f in cp.factors:
        accumulator *= f.sum(axis=0)
    if return_trimmed_cp:
        raise NotImplementedError
    else:
        return (accumulator != 0.0).sum()
    
    
# function to calculate f1 score from composite precision & recall scores
def composite_f1(precision, recall):
    '''
    Calculates F1 score from precision and recall.'''
    numerator = precision + recall
    if numerator == 0:
        return 0
    else:
        return (2 * precision * recall) / numerator


# main experiment script
def main():
    
    # output directory and experiment parameters
    base_dir = Path('../../data/1-simulation')
    filepath_fit_data = base_dir / 'fitting_data.csv'
    filepath_cv_data = base_dir / 'cv_data.csv'
    n_simulations = 100
    replicates = ['A', 'B', 'C']
    n_replicates = len(replicates) 
    fitting_results = []
    cv_results = []
    
    # simulation parameters
    sim_shape_pools = [
        np.arange(10, 1001), 
        np.arange(10, 101), 
        np.arange(10, 101)
    ]
    sim_rank_pool = np.arange(2, 11)
    sim_density_dist = scipy.stats.uniform()
    sim_noise_level_dist = scipy.stats.uniform(-1, 2)    # log10 of noise_level drawn from this
    sim_dist_list = [
        scipy.stats.uniform(loc=-1, scale=2), 
        scipy.stats.uniform(), 
        scipy.stats.uniform()
    ]
    sim_params = []
    
    # set random states
    seed = 1987
    rns = check_random_state(seed)
    sim_density_dist.random_state = rns
    sim_noise_level_dist.random_state = rns
    
    # model parameters
    model_params = {
        'rank': [int(i) for i in np.linspace(1, 12, 12)], 
        'lambdas': [[i, 0.0, 0.0] for i in [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]], 
        'nonneg_modes': [[1, 2]],
        'tol': [1e-5], 
        'n_iter_max': [2000], 
        'n_initializations': [5]
    }
    param_grid = list(ParameterGrid(model_params))
    
    # begin experiment
    for sim_id in range(n_simulations):
        output_dir = base_dir / 'sim{}'.format(sim_id)
        # make path directories if they don't exist yet
        if not output_dir.exists():
            output_dir.mkdir(parents=True)
        
        # set simulation parameters
        sim_p = {}
        sim_p['id'] = sim_id
        sim_p['shape'] = [int(rns.choice(pool)) for pool in sim_shape_pools]
        sim_p['rank'] = int(rns.choice(sim_rank_pool))
        min_densities = [1/s for s in sim_p['shape']]
        # draw densities until they are above the minimum necessary so that
        # all factor matricies are full rank
        densities_above_mins = False
        while not densities_above_mins:
            sim_p['densities'] = [d for d in sim_density_dist.rvs(len(sim_p['shape']))]
            densities_above_mins = np.all([d > min_densities[i] for i, d in enumerate(sim_p['densities'])]) 
        sim_p['noise_level'] = 10 ** sim_noise_level_dist.rvs()
        
        # re-seed simulated data until all factor matrices are full rank
        full_rank = False
        while not full_rank:
            # generate new seed from rns
            sim_p['seed'] = rns.randint(2**32)
            # generate SparseCP simulated data tensor
            sim_tensor = simulated_sparse_tensor(
                shape=sim_p['shape'],                
                rank=sim_p['rank'],                         
                densities=sim_p['densities'], 
                factor_dist_list=sim_dist_list, 
                random_state=sim_p['seed']
            )
            # check that all factors are full rank
            full_rank = np.all([np.linalg.matrix_rank(f) == sim_p['rank'] for f in sim_tensor.factors])
        # save parameters of successful simulation
        sim_params.append(sim_p)
        with open(output_dir / 'simulation_parameters.txt', 'w') as f:
            f.write(json.dumps(sim_p, indent=4))
        print('Beginning simulation {}: {}'.format(sim_id, sim_tensor), flush=True)
        
        # store ground truth simulated dataset
        store_cp_tensor(sim_tensor, output_dir / 'simulation_ground_truth.h5')
        
        # save plot of simulation factors
        heatmap_params = {'vmin':-1, 'vmax':1, 'cmap':'coolwarm', 'center':0}
        fig, ax = plot_factors_heatmap(
            tl.cp_normalize(sim_tensor).factors, 
            mask_thold=[0, 0], 
            ratios=False, 
            heatmap_kwargs=heatmap_params)
        fig.savefig(output_dir / 'simulation_factors.png')     
        
        # generate simulated replicates
        for rep in replicates:
            print('\nFitting replicate {} models\n'.format(rep), flush=True)
            
            # make a new directory
            path = output_dir / 'replicate{}'.format(rep)
            if not path.exists():
                path.mkdir(parents=True)
                
            # generate unique noisy tensor
            tensor = sim_tensor.to_tensor(
                noise_level=sim_p['noise_level'], 
                sparse_noise=True, 
                random_state=rns
            )
            
            # save simulated replicate data
            da = xr.DataArray(tensor)
            da.to_netcdf(path / 'simulation_data_replicate_{}.nc'.format(rep))
    
            # instantiate models, with reproducible random seed for each
            models = [SparseCP(**p, random_state=rns.randint(2**32)) for p in param_grid]

            # assemble job parameters
            job_params = (
                models, 
                [tensor for m in models], 
                [path / 'rank{}'.format(m.rank) / 'lambda{}'.format(m.lambdas[0]) for m in models], 
                [{'threads': 1, 'verbose': 0} for m in models]
            )
                
            # run jobs 
            executor = ProcessPoolExecutor()
            fit_models = executor.map(fit_save_model, *job_params)
                
            # iterate through models in order
            for model in fit_models:
                
                # calculate metrics
                rank = model.rank
                lamb = model.lambdas[0]
                best_init = model._best_cp_index
                loss = model.loss_[-1]
                cvg_iter = len(model.loss_)
                data_sse = relative_sse(model.decomposition_, tensor)
                model_sse = relative_sse(model.decomposition_, sim_tensor.to_tensor())
                fms = factor_match_score(model.decomposition_, sim_tensor, consider_weights=False, allow_smaller_rank=True)
                degeneracy = degeneracy_score(model.decomposition_)
                cc = core_consistency(model.decomposition_, tensor)
                can_fms = [factor_match_score(model.decomposition_, c, consider_weights=False, allow_smaller_rank=True) for c in model.candidates_]
                can_sse = [relative_sse(c, tensor) for c in model.candidates_]
                
                # record metrics
                fitting_results.append(
                    {
                        'datetime': datetime.datetime.now(), 
                        'simulation_id': sim_id, 
                        'replicate': rep, 
                        'rank': rank, 
                        'lambda': lamb, 
                        'best_init': best_init, 
                        'loss': loss, 
                        'convergence_iterations': cvg_iter, 
                        'data_sse': data_sse, 
                        'model_sse': model_sse, 
                        'fms': fms, 
                        'degeneracy': degeneracy, 
                        'core_consistency': cc, 
                        'candidate_fms': can_fms, 
                        'candidate_sse': can_sse
                    }
                )
                
                # print some metrics
                print('rank:{}, lambda:{}, sse:{}, fms:{}'.format(
                    rank, 
                    lamb, 
                    data_sse, 
                    fms
                ), flush=True)
            
                # save data
                fitting_df = pd.DataFrame(fitting_results)
                fitting_df.to_csv(filepath_fit_data, index=False)
            
            # shut down executor
            executor.shutdown()

    # collect cross validation results
    for sim_id in range(n_simulations):
        # set path of simulated data
        sim_path = base_dir / 'sim{}'.format(sim_id)
        # read in ground truth simulation data
        true_cp = load_cp_tensor(sim_path / 'simulation_ground_truth.h5')
        print('\nBeginning cross validataion of simulation {}: {}'.format(sim_id, true_cp), flush=True)
        # pull noise level from parameter file
        with open(sim_path / 'simulation_parameters.txt') as f:
            sim_params = json.load(f)
            noise_level = sim_params['noise_level']
            sim_densities = sim_params['densities']
            gene_sparsity = sim_params['densities'][0] * sim_params['shape'][0]
        # read in simulated data
        rep_data = {}
        for rep in replicates:
            rep_path = sim_path / 'replicate{}'.format(rep) 
            rep_data[rep] = xr.open_dataarray(rep_path / 'simulation_data_replicate_{}.nc'.format(rep))
        # iterate through all parameter combos
        for params in param_grid:
            # get all the models
            cps = {}
            expt_path = 'rank{}/lambda{}'.format(params['rank'], params['lambdas'][0])
            for rep in replicates:
                cp_path = sim_path / 'replicate{}'.format(rep) / expt_path
                cps[rep] = load_cp_tensor(cp_path / 'fitted_model.h5')
            for modeled_rep in replicates:
                for comparison_rep in replicates:
                    # calculate fms & cosine similiary scores against other fit models
                    if modeled_rep < comparison_rep:
                        model_fms = factor_match_score(
                            cps[modeled_rep], 
                            cps[comparison_rep], 
                            consider_weights=False, 
                            allow_smaller_rank=True
                        )
                        model_cossim = cosine_similarity(
                            cps[modeled_rep].factors[0], 
                            cps[comparison_rep].factors[0]
                        )
                        model_sup_cossim = cosine_similarity(
                            (cps[modeled_rep].factors[0] != 0), 
                            (cps[comparison_rep].factors[0] != 0)
                        )
                    else:
                        model_fms = model_cossim = model_sup_cossim = np.nan
                    # calculate fms, cosine, & cluster scores against true factors
                    if modeled_rep == comparison_rep:
                        true_fms = factor_match_score(
                            true_cp, 
                            cps[modeled_rep], 
                            consider_weights=False, 
                            allow_smaller_rank=True
                        )
                        true_cossim = cosine_similarity(
                            true_cp.factors[0], 
                            cps[modeled_rep].factors[0]
                        )
                        true_sup_cossim = cosine_similarity(
                            (true_cp.factors[0] != 0), 
                            (cps[modeled_rep].factors[0] != 0)
                        )
                        clusters_true = SparseCPTensor(true_cp).get_clusters(0, boolean=True)
                        clusters_fit = SparseCPTensor(cps[modeled_rep]).get_clusters(0, boolean=True)
                        recovery, relevance = recovery_relevance(clusters_true, clusters_fit)
                        precision, recall = pairs_precision_recall(clusters_true, clusters_fit)
                    else:
                        true_fms = true_cossim = true_sup_cossim = np.nan
                        recovery = relevance = precision = recall = np.nan
                    # calculate relative sse
                    rel_sse = relative_sse(cps[modeled_rep], rep_data[comparison_rep].data)
                    # keep results
                    cv_results.append(
                        {
                            'simulation_id': sim_id, 
                            'simulation_rank': true_cp.rank, 
                            'simulation_shape': true_cp.shape, 
                            'simulation_densities': sim_densities, 
                            'simulation_mean_gene_sparsity': gene_sparsity, 
                            'noise_level': noise_level,
                            'rank': params['rank'], 
                            'lambda': params['lambdas'][0], 
                            'modeled_replicate': modeled_rep, 
                            'comparison_replicate': comparison_rep, 
                            'replicate_pair': '{}, {}'.format(modeled_rep, comparison_rep), 
                            'n_components': nonzero_components(cps[modeled_rep]), 
                            'mean_gene_sparsity': 
                                (cps[modeled_rep].factors[0] != 0.0).sum(axis=0).mean(), 
                            'relative_sse': rel_sse, 
                            'model_fms': model_fms, 
                            'true_fms': true_fms, 
                            'model_factor0_cosine_similarity': model_cossim, 
                            'true_factor0_cosine_similarity': true_cossim, 
                            'model_factor0_support_cosine_similarity': model_sup_cossim, 
                            'true_factor0_support_cosine_similarity': true_sup_cossim, 
                            'recovery': recovery, 
                            'relevance': relevance, 
                            'precision': precision, 
                            'recall': recall, 
                            'f1': composite_f1(precision, recall)
                        }
                    ) 
                    # store results in dataframe
                    cv_df = pd.DataFrame(cv_results)
        # save data for each simulation separately
        cv_df[cv_df['simulation_id'] == sim_id].to_csv(sim_path / 'cv_data.csv', index=False)
    # save copy of all results combined
    cv_df.to_csv(filepath_cv_data, index=False)

if __name__ == "__main__":
  main()
