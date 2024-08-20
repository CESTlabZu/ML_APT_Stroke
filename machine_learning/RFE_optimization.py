#!/usr/bin/env python
# coding: utf-8
__author__ = "Malvika Viswanathan, Leqi Yin, Yashwant Kurmi, Zhongliang Zu"
__maintainer__ = "Malvika Viswanathan, Zhongliang Zu"
__email__ = "zhongliang.zu@vumc.org"


import numpy as np
from numpy import mean
from numpy import std
from sklearn.model_selection import cross_val_score
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestRegressor
from sklearn.pipeline import Pipeline
from sklearn.metrics import get_scorer_names
from matplotlib import pyplot
import scipy.io


# load matrix with training and target data
# mat = scipy.io.loadmat('...')

# define X_train
# X_train = mat[...]
# y = mat[...]

def spectrum(a,w,n):
    step=50
    offset= [i for i in range(-1500,1550, step)]
    rffreq=np.array([-4000, -3500, -3000, -2500]+ offset+ [2500, 3000,3500,4000])
    y = a/(1+(((rffreq-n)**2)/((0.5*(w))**2)))
    return y



def get_models():
    models = dict()
    for i in range(1,68):
        rfe = RFE(estimator=RandomForestRegressor(), n_features_to_select=i)
        model = RandomForestRegressor()
        models[str(i)] = Pipeline(steps=[('s',rfe),('m',model)])
    return models

# evaluate a give model using cross-validation
def evaluate_model(model, X, y):
     scores = cross_val_score(model, X, y, scoring='neg_mean_absolute_error', n_jobs=-1, error_score='raise')
     return scores
 

# define dataset
X, y = X_train, y
# get the models to evaluate
models = get_models()
# evaluate the models and store results
results, names = list(), list()
for name, model in models.items():  
    print(name)
    scores = evaluate_model(model, X, y)
    results.append(scores)
    names.append(name)
    print('>%s %.3f (%.3f)' % (name, mean(scores), std(scores)))
# plot model performance for comparison
pyplot.boxplot(results, labels=names, showmeans=True)
pyplot.show()

min_score=np.argmin(scores, axis=0)

# define the method
rfe = RFE(estimator=RandomForestRegressor(), n_features_to_select=min_score)
# fit the model
rfe.fit(X_train, y)
# transform the data


for i in range(X.shape[1]):
    if( rfe.ranking_[i] == 1):
        print('Column: %d, Selected %s, Rank: %.3f' % (i, rfe.support_[i], rfe.ranking_[i]))
        
