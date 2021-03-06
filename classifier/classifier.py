#!/usr/bin/env python
import argparse
import sys
import uuid

import calour as ca
import calour_utils as cu

import numpy as np
import pandas as pd
import calour


def classifier_performance_matrix(exp: ca.Experiment, use_subset_features=True, shuffle=False, shuffle_source=False):
    '''test the cross cohort classifier performance

    Parameters
    ----------
    exp: ca.Experiment with the different cohort subsets
    use_subset_features: bool, optional
        True to filter features for classifier before each cross-cohort classification (feature intersection of the 2 cohorts)
        False to use all the exp features for classifiers
    shuffle: bool, optional
        if True, randomize the HC/disease labels of validation cohort before predicting using the model (for null hypothesis values)
    shuffle_source: bool, optional
        if True, randomize the HC/disease labels of training cohort before predicting (for null hypothesis values)

    Returns
    -------
    ROC pandas.dataframe
        2d numpy.array. row is the training experiment, column is the testing experiment, value is the ROC.
        indexes and columns are the cohort ids (exp field from exp.sample_metadata)
    balanced_accuracy pandas.dataframe
        2d numpy.array. row is the training experiment, column is the testing experiment, value is the balanced accuracy score.
        indexes and columns are the cohort ids (exp field from exp.sample_metadata)
    '''
    ca.set_log_level('ERROR')
    num_exp = len(exp.sample_metadata['exp'].unique())
    print('processing %d experiments' % num_exp)
    roc_mat = np.zeros([num_exp, num_exp])
    accuracy_mat = np.zeros([num_exp, num_exp])

    ids1 = []
    for idx1, (id1, exp1) in enumerate(exp.iterate('exp')):
        if shuffle_source:
            exp1.sample_metadata['type'] = exp1.sample_metadata['type'].sample(frac=1).values
        exp1 = exp1.filter_sum_abundance(0, strict=True)
        ids1.append(id1.split('/')[-1])
        # if not subset of features, train once on the exp1 dataset
        if not use_subset_features:
            model = cu.classify_fit(exp1, 'type')
        for idx2, (id2, exp2) in enumerate(exp.iterate('exp')):
            exp2 = exp2.filter_sum_abundance(0, strict=True)
            cexp1 = exp1.filter_ids(exp2.feature_metadata.index)
            cexp2 = exp2.filter_ids(cexp1.feature_metadata.index)
            # if shuffle, mix the HC/disease of exp2
            if shuffle:
                cexp2.sample_metadata['type'] = cexp2.sample_metadata['type'].sample(frac=1).values

            # if same experiment, so do training/validation
            if idx1 == idx2:
                # keep 2/3 of samples
                cexp1 = cexp1.downsample('exp', keep=2 * int(len(cexp1.sample_metadata) / 3))
                cexp2 = cexp2.filter_ids(cexp1.sample_metadata.index, axis='s', negate=True)

            # if using subset features, train the model after subsetting
            if use_subset_features:
                model = cu.classify_fit(cexp1, 'type')

            # now predict on exp2 and measure performance
            res = cu.classify_predict(cexp2, 'type', model, plot_it=False)
            try:
                # we wrap this in a try since set may include only one class
                roc_auc = cu.classify_get_roc(res)
            except:
                roc_auc = None
            roc_mat[idx1, idx2] = roc_auc

            try:
                accuracy = cu.classify_get_accuracy(res)
            except:
                accuracy = None
            accuracy_mat[idx1, idx2] = accuracy

    ca.set_log_level('INFO')
    resdf = pd.DataFrame(roc_mat, index=ids1, columns=ids1)
    resdf_accuracy = pd.DataFrame(accuracy_mat, index=ids1, columns=ids1)
    return resdf, resdf_accuracy


def main(argv):
    parser = argparse.ArgumentParser(description='metaanalysis cross-classifier', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help='name of input biom table')
    parser.add_argument('-m', '--map', help='name of input mapping file')
    parser.add_argument('-o', '--output', help='name of output file')

    parser.add_argument('--use-subset', help="Use only subset of features present in both cohorts for classifier training", action='store_true', default=True)
    parser.add_argument('--shuffle', help="Shuffle testing labels", action='store_true', default=False)
    parser.add_argument('--shuffle-source', help="Shuffle training labels", action='store_true', default=False)
    parser.add_argument('--uname', help="add unique signature to name", action='store_true', default=False)

    args = parser.parse_args(argv)

    cname = args.output
    if args.uname:
        cname += '_'
        cname = cname + str(uuid.uuid1())
        # cname += str(int(time.time() * 1000000))
    print('started processing file %s' % cname)

    # load the experiment
    print('loading experiment %s' % args.input)
    exp = ca.read_amplicon(args.input, args.map, min_reads=1000, normalize=10000)
    print(exp)

    # run the classifier
    print('running the classifier')
    resdf_roc, resdf_accuracy = classifier_performance_matrix(exp=exp, use_subset_features=args.use_subset, shuffle=args.shuffle, shuffle_source=args.shuffle_source)

    # save the results
    print('saving to %s' % cname)
    resdf_roc.to_csv(cname + '_roc' + '.csv')
    resdf_accuracy.to_csv(cname + '_accuracy' + '.csv')


if __name__ == "__main__":
    main(sys.argv[1:])
